################################################################################
#                                                                              #
#                   SYNGAP1 VARIANT VIEWER - Shiny Application                 #
#                                                                              #
#  Purpose: Interactive genome browser for visualizing SYNGAP1 genetic         #
#           variants with patient data and research asset tracking             #
#                                                                              #
#  Architecture:                                                               #
#    1. Package Dependencies (~line 33)                                        #
#                                                                              #
#    1B. ClinVar Data Loading (~line 72)                                       #
#       - Fetch SYNGAP1 variants from NCBI ClinVar via E-utilities API         #
#       - Cache results locally with 7-day expiry                              #
#                                                                              #
#    2-6. Data Structures, Helpers & API Integration (~lines 441-828)          #
#       - Amino acid code lookup table                                         #
#       - Protein nomenclature formatting functions                            #
#       - cDNA coordinate extraction                                           #
#       - Ensembl REST API queries with permanent disk cache                   #
#                                                                              #
#    7. Data Loading & Preprocessing (~line 830)                               #
#       - Load patient variant data from CSV                                   #
#       - Convert cDNA coordinates to genome positions via Ensembl API         #
#       - Merge positions into variant data frame                              #
#                                                                              #
#    8. UI Definition (~line 1046)                                             #
#       - Sidebar with dynamic variant type buttons                            #
#       - Static ClinVar track button + log-spaced size filter slider (1–50,000 KB) #
#       - Checkboxes for filtering by research assets                          #
#       - Main panel with IGV genome browser                                   #
#                                                                              #
#    9-9B. GFF3 Track Formatting (~lines 1195-1428)                            #
#       - Format internal variant data as GFF3 for IGV                        #
#       - Format ClinVar data as GFF3 for IGV                                 #
#                                                                              #
#    10. Server Logic (~line 1430)                                             #
#       - Reactive filtering based on user selections                          #
#       - Dynamic track loading/updating for internal variants                 #
#       - ClinVar track loading on button click                                #
#       - ClinVar track reloading on size slider change                        #
#       - IGV browser rendering and interaction                                #
#                                                                              #
#  Design Philosophy:                                                          #
#    - Separate data processing from UI to improve performance                 #
#    - Cache external API calls to work offline and avoid rate limits          #
#    - Use reactive programming for smooth filter interactions                 #
#    - Generate tracks dynamically based on available variant types            #
#                                                                              #
################################################################################

################################################################################
# SECTION 1: PACKAGE DEPENDENCIES
################################################################################

# Core Shiny framework for web application
library(shiny)

# igvShiny: Wrapper for Integrative Genomics Viewer (IGV) in R
# Why: Provides interactive genome browser functionality directly in Shiny
library(igvShiny)

# Data manipulation with tidyverse syntax
# Why: Clean, readable syntax for filtering and transforming variant data
library(dplyr)

# String pattern matching and manipulation
# Why: Parse variant nomenclature (e.g., "c.333del", "p.R135X")
library(stringr)

# HTTP requests to REST APIs
# Why: Query Ensembl API for coordinate conversions
library(httr)

# JSON parsing
# Why: Parse Ensembl API responses
library(jsonlite)

# Function memoization (caching)
# Why: Cache expensive API calls to disk for instant retrieval on subsequent runs
library(memoise)

# Extended Shiny widgets
# Why: sliderTextInput() provides a logarithmic-style slider by accepting an
#      explicit vector of stop values, giving fine-grained control at the low
#      end of the ClinVar size filter without sacrificing upper-range coverage
library(shinyWidgets)

################################################################################
# CITATION
################################################################################
# Shannon P, Gladki A, Scigocka K (2024). igvShiny: igvShiny: a wrapper of 
# Integrative Genomics Viewer (IGV - an interactive tool for visualization 
# and exploration integrated genomic data). R package version 1.2.0, 
# https://gladkia.github.io/igvShiny/, https://github.com/gladkia/igvShiny.

################################################################################
# SECTION 1B: CLINVAR DATA LOADING — LIVE NCBI E-UTILITIES FETCH WITH CACHE
################################################################################

# ClinVar data is fetched live from NCBI and cached locally for 7 days.
#
# Why live fetch instead of a static file?
#   - ClinVar is updated continuously; new SYNGAP1 submissions appear regularly.
#   - A static file requires manual re-download and redeployment to stay current.
#   - The E-utilities API is free, does not require an account, and is stable.
#
# Cache strategy (time-based, not permanent):
#   Unlike the Ensembl coordinate cache (which is permanent because hg38
#   positions never change), ClinVar classifications and submissions DO change.
#   We use a 7-day expiry: on first run or after 7 days, the app re-fetches;
#   otherwise it loads the cached .rds file instantly.
#
#   Cache files:
#     cache/clinvar/clinvar_data.rds   — parsed data frame
#     cache/clinvar/last_updated.txt   — ISO timestamp of last successful fetch
#
# NCBI E-utilities used (no API key required; rate limit = 3 req/sec):
#   1. esearch  — retrieve all ClinVar variation IDs for SYNGAP1[gene]
#   2. esummary — retrieve full variant metadata for each ID (batched, 200/call)
#
# Note on API key:
#   If your organization obtains an NCBI API key in the future, add it as:
#     api_key = "YOUR_KEY_HERE"
#   to the GET() calls below, which raises the rate limit to 10 req/sec.

# Ensure cache directory exists
dir.create("cache/clinvar", recursive = TRUE, showWarnings = FALSE)

CLINVAR_CACHE_RDS       <- "cache/clinvar/clinvar_data.rds"
CLINVAR_CACHE_TIMESTAMP <- "cache/clinvar/last_updated.txt"
CLINVAR_CACHE_DAYS      <- 7          # Re-fetch after this many days
CLINVAR_BATCH_SIZE      <- 200        # IDs per esummary call (NCBI recommended max)
CLINVAR_RATE_DELAY      <- 0.4        # Seconds between requests (stays under 3/sec)
NCBI_BASE               <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"

#' Check Whether the ClinVar Cache is Still Fresh
#'
#' @description
#' Returns TRUE if a valid cache file exists and was written within
#' CLINVAR_CACHE_DAYS days. Returns FALSE if the cache is absent,
#' unreadable, or expired — triggering a fresh API fetch.
clinvar_cache_is_fresh <- function() {
  if (!file.exists(CLINVAR_CACHE_RDS) ||
      !file.exists(CLINVAR_CACHE_TIMESTAMP)) return(FALSE)
  
  timestamp_str <- tryCatch(
    readLines(CLINVAR_CACHE_TIMESTAMP, n = 1, warn = FALSE),
    error = function(e) return("")
  )
  if (nchar(trimws(timestamp_str)) == 0) return(FALSE)
  
  last_updated <- tryCatch(
    as.POSIXct(timestamp_str, format = "%Y-%m-%dT%H:%M:%S", tz = "UTC"),
    error = function(e) return(NA)
  )
  if (is.na(last_updated)) return(FALSE)
  
  age_days <- as.numeric(difftime(Sys.time(), last_updated, units = "days"))
  return(age_days < CLINVAR_CACHE_DAYS)
}

#' Fetch All SYNGAP1 ClinVar Variation IDs via esearch
#'
#' @description
#' Queries NCBI esearch for all ClinVar entries associated with SYNGAP1.
#' Returns a character vector of variation IDs (e.g. c("12345", "67890", ...)).
#' Returns NULL on network failure so the caller can fall back to cache.
fetch_clinvar_ids <- function() {
  message("ClinVar: querying NCBI esearch for SYNGAP1 variation IDs...")
  
  url <- paste0(
    NCBI_BASE, "esearch.fcgi",
    "?db=clinvar",
    "&term=SYNGAP1%5Bgene%5D",   # SYNGAP1[gene] — URL-encoded
    "&retmax=10000",              # Upper bound; SYNGAP1 currently has ~1900
    "&retmode=json"
  )
  
  resp <- tryCatch(
    GET(url, httr::timeout(30)),   # Qualified to avoid namespace collision with other packages
    error = function(e) {
      warning("ClinVar esearch network error: ", conditionMessage(e))
      NULL
    }
  )
  
  if (is.null(resp)) return(NULL)
  
  if (http_error(resp)) {
    warning("ClinVar esearch HTTP error: status ", status_code(resp),
            " — ", content(resp, "text", encoding = "UTF-8"))
    return(NULL)
  }
  
  parsed <- tryCatch(
    fromJSON(content(resp, "text", encoding = "UTF-8")),
    error = function(e) {
      warning("ClinVar esearch JSON parse error: ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(parsed)) return(NULL)
  
  ids <- parsed$esearchresult$idlist
  message("ClinVar: found ", length(ids), " variation IDs.")
  return(ids)
}

#' Parse a Single esummary Record into a One-Row Data Frame
#'
#' @description
#' Extracts the fields we need from one ClinVar esummary JSON result object.
#' Any field that is absent or malformed in the JSON is set to NA gracefully.
#'
#' Field mapping verified against live NCBI esummary API response (2026-02):
#'   - Assembly field is "assembly_name" (not "Assembly")
#'   - Somatic field is "clinical_impact_classification" (not "somatic_clinical_impact")
#'   - molecular_consequence_list entries are plain strings (not named sub-objects)
#'   - trait_set lives inside germline_classification (not at record top level)
#'   - Canonical SPDI is at variation_set[[1]]$canonical_spdi
#'   - Variant name is at variation_set[[1]]$variation_name (also in top-level title)
#'   - Gene symbol is at genes[[1]]$symbol
#'
#' @param rec One list element from the esummary "result" object
#' @return Single-row data frame with all ClinVar metadata columns
parse_esummary_record <- function(rec) {
  
  # Helper: safely extract a scalar value, return NA if absent/empty
  safe <- function(x, default = NA_character_) {
    if (is.null(x) || length(x) == 0) return(default)
    v <- trimws(as.character(x[[1]]))
    if (length(v) == 0 || is.na(v) || v == "") return(default)
    v
  }
  
  # --- Variant name, gene, basic fields ---
  vset         <- rec$variation_set[[1]]
  name         <- safe(vset$variation_name)
  gene         <- safe(rec$genes[[1]]$symbol)
  protein_chg  <- safe(rec$protein_change)
  accession    <- safe(rec$accession)
  variation_id <- safe(rec$uid)
  variant_type <- safe(rec$obj_type)
  canonical_spdi <- safe(vset$canonical_spdi)
  
  # molecular_consequence_list entries are plain character strings
  mol_conseq <- if (!is.null(rec$molecular_consequence_list) &&
                    length(rec$molecular_consequence_list) > 0) {
    safe(rec$molecular_consequence_list[[1]])
  } else { NA_character_ }
  
  # allele_id and dbsnp_id are not present in esummary — set to NA
  # (available in efetch full XML if needed in future)
  allele_id <- NA_character_
  dbsnp_id  <- NA_character_
  
  # --- Condition: trait_set lives inside germline_classification ---
  germ      <- rec$germline_classification
  germ_sig  <- safe(germ$description)
  germ_date <- safe(germ$last_evaluated)
  germ_rev  <- safe(germ$review_status)
  
  condition <- if (!is.null(germ$trait_set) && length(germ$trait_set) > 0) {
    trait_names <- sapply(germ$trait_set, function(t) {
      # trait_name can be either a direct string or a list with a name field
      if (is.character(t)) safe(t)
      else safe(t$trait_name)
    })
    trait_names <- trait_names[!is.na(trait_names)]
    if (length(trait_names) > 0) paste(trait_names, collapse = "; ") else NA_character_
  } else { NA_character_ }
  
  # --- Somatic clinical impact (field is "clinical_impact_classification") ---
  som       <- rec$clinical_impact_classification
  som_sig   <- safe(som$description)
  som_date  <- safe(som$last_evaluated)
  som_rev   <- safe(som$review_status)
  
  # --- Oncogenicity ---
  onco      <- rec$oncogenicity_classification
  onco_sig  <- safe(onco$description)
  onco_date <- safe(onco$last_evaluated)
  onco_rev  <- safe(onco$review_status)
  
  # --- GRCh38 genomic coordinates ---
  # variation_loc is a list of location objects, one per assembly.
  # Match on assembly_name == "GRCh38" (field verified from live response).
  chr38   <- NA_character_
  start38 <- NA_real_
  end38   <- NA_real_
  
  locs <- vset$variation_loc
  if (!is.null(locs) && length(locs) > 0) {
    for (loc in locs) {
      if (!is.null(loc$assembly_name) && loc$assembly_name == "GRCh38") {
        chr38   <- safe(loc$chr)
        start38 <- suppressWarnings(as.numeric(safe(loc$start)))
        end38   <- suppressWarnings(as.numeric(safe(loc$stop)))
        # SNVs have start == stop; give 1 bp width so IGV renders them visibly
        if (!is.na(start38) && !is.na(end38) && start38 == end38) {
          end38 <- end38 + 1L
        }
        break
      }
    }
  }
  
  # Return one-row data frame with column names matching create_clinvar_gff3_data()
  data.frame(
    Name                                          = name,
    Gene.s.                                       = gene,
    Protein.change                                = protein_chg,
    Condition.s.                                  = condition,
    Accession                                     = accession,
    VariationID                                   = variation_id,
    AlleleID.s.                                   = allele_id,
    dbSNP.ID                                      = dbsnp_id,
    Canonical.SPDI                                = canonical_spdi,
    Variant.type                                  = variant_type,
    Molecular.consequence                         = mol_conseq,
    Germline.classification                       = germ_sig,
    Germline.date.last.evaluated                  = germ_date,
    Germline.review.status                        = germ_rev,
    Somatic.clinical.impact                       = som_sig,
    Somatic.clinical.impact.date.last.evaluated   = som_date,
    Somatic.clinical.impact.review.status         = som_rev,
    Oncogenicity.classification                   = onco_sig,
    Oncogenicity.date.last.evaluated              = onco_date,
    Oncogenicity.review.status                    = onco_rev,
    clinvar_chr                                   = if (!is.na(chr38)) paste0("chr", chr38) else NA_character_,
    clinvar_start                                 = start38,
    clinvar_end                                   = end38,
    stringsAsFactors = FALSE
  )
}

#' Fetch Full ClinVar Metadata for a Vector of Variation IDs via esummary
#'
#' @description
#' Calls NCBI esummary in batches of CLINVAR_BATCH_SIZE IDs, parses each
#' result record, and returns a single combined data frame.
#'
#' @param ids Character vector of ClinVar variation IDs from fetch_clinvar_ids()
#' @return Data frame of all variants, or NULL on complete failure
fetch_clinvar_summaries <- function(ids) {
  batches    <- split(ids, ceiling(seq_along(ids) / CLINVAR_BATCH_SIZE))
  n_batches  <- length(batches)
  all_rows   <- vector("list", n_batches)
  
  message("ClinVar: fetching summaries in ", n_batches,
          " batches of up to ", CLINVAR_BATCH_SIZE, " IDs each...")
  
  for (i in seq_along(batches)) {
    message("  Batch ", i, " / ", n_batches, "...")
    
    id_string <- paste(batches[[i]], collapse = ",")
    url <- paste0(
      NCBI_BASE, "esummary.fcgi",
      "?db=clinvar",
      "&id=", id_string,
      "&retmode=json"
    )
    
    resp <- tryCatch(GET(url), error = function(e) NULL)
    if (is.null(resp) || http_error(resp)) {
      warning("  esummary batch ", i, " failed — skipping.")
      Sys.sleep(CLINVAR_RATE_DELAY)
      next
    }
    
    parsed <- tryCatch(
      fromJSON(content(resp, "text", encoding = "UTF-8"), simplifyVector = FALSE),
      error = function(e) NULL
    )
    if (is.null(parsed) || is.null(parsed$result)) {
      warning("  esummary batch ", i, " returned unparseable JSON — skipping.")
      Sys.sleep(CLINVAR_RATE_DELAY)
      next
    }
    
    # parsed$result is a named list; "uids" is the list of IDs, the rest are records
    uids    <- parsed$result$uids
    records <- parsed$result[names(parsed$result) != "uids"]
    
    batch_rows <- lapply(records, function(rec) {
      tryCatch(parse_esummary_record(rec), error = function(e) NULL)
    })
    batch_rows <- Filter(Negate(is.null), batch_rows)
    
    if (length(batch_rows) > 0) {
      all_rows[[i]] <- do.call(rbind, batch_rows)
    }
    
    Sys.sleep(CLINVAR_RATE_DELAY)   # Respect NCBI rate limit (3 req/sec)
  }
  
  combined <- do.call(rbind, Filter(Negate(is.null), all_rows))
  if (is.null(combined) || nrow(combined) == 0) return(NULL)
  
  message("ClinVar: successfully parsed ", nrow(combined), " variant records.")
  return(combined)
}

#' Load ClinVar Data — from Cache if Fresh, Otherwise Fetch from NCBI
#'
#' @description
#' Top-level entry point called once at app startup.
#'
#' Decision logic:
#'   1. If cache exists and is < 7 days old → load .rds (instant, no network)
#'   2. Otherwise → fetch from NCBI E-utilities, save to cache
#'   3. If fetch fails AND stale cache exists → use stale cache with a warning
#'   4. If fetch fails AND no cache → return empty data frame (app still loads)
#'
#' The 7-day window balances freshness against startup time.
#' ClinVar is updated monthly for most genes; 7-day polling is more than adequate.
#'
#' @return Data frame ready for use by create_clinvar_gff3_data()
load_clinvar_data <- function() {
  
  # ── Path 1: Fresh cache ──────────────────────────────────────────────────
  if (clinvar_cache_is_fresh()) {
    message("ClinVar: loading from cache (< ", CLINVAR_CACHE_DAYS, " days old).")
    return(readRDS(CLINVAR_CACHE_RDS))
  }
  
  # ── Path 2: Fetch from NCBI ──────────────────────────────────────────────
  message("ClinVar: cache absent or expired — fetching from NCBI E-utilities...")
  
  ids  <- fetch_clinvar_ids()
  data <- if (!is.null(ids) && length(ids) > 0) fetch_clinvar_summaries(ids) else NULL
  
  if (!is.null(data) && nrow(data) > 0) {
    # Save fresh data and update timestamp
    saveRDS(data, CLINVAR_CACHE_RDS)
    writeLines(
      format(Sys.time(), format = "%Y-%m-%dT%H:%M:%S", tz = "UTC"),
      CLINVAR_CACHE_TIMESTAMP
    )
    message("ClinVar: cache updated successfully.")
    return(data)
  }
  
  # ── Path 3: Fetch failed — fall back to stale cache if available ─────────
  if (file.exists(CLINVAR_CACHE_RDS)) {
    warning(
      "ClinVar: NCBI fetch failed. Falling back to stale cache. ",
      "Data may be outdated."
    )
    return(readRDS(CLINVAR_CACHE_RDS))
  }
  
  # ── Path 4: No data at all ───────────────────────────────────────────────
  warning(
    "ClinVar: NCBI fetch failed and no cache found. ",
    "ClinVar track will be empty."
  )
  return(data.frame())
}

# Load (or fetch) ClinVar data at startup.
# On first run this takes ~30-60 seconds depending on network speed.
# On subsequent runs within 7 days it loads from cache in < 1 second.
clinvar_data <- load_clinvar_data()

################################################################################
# SECTION 2: DATA STRUCTURES & CONSTANTS
################################################################################

# Amino Acid Code Mapping Table
# Purpose: Convert single-letter amino acid codes to three-letter codes
# Why: Three-letter codes are more readable for biologists
# Example: "R135X" → "Arg135ter" (Arginine at position 135 to stop codon)
aa_code_map <- c(
  A = "Ala",  # Alanine
  C = "Cys",  # Cysteine
  D = "Asp",  # Aspartate
  E = "Glu",  # Glutamate
  F = "Phe",  # Phenylalanine
  G = "Gly",  # Glycine
  H = "His",  # Histidine
  I = "Ile",  # Isoleucine
  K = "Lys",  # Lysine
  L = "Leu",  # Leucine
  M = "Met",  # Methionine
  N = "Asn",  # Asparagine
  P = "Pro",  # Proline
  Q = "Gln",  # Glutamine
  R = "Arg",  # Arginine
  S = "Ser",  # Serine
  T = "Thr",  # Threonine
  V = "Val",  # Valine
  W = "Trp",  # Tryptophan
  Y = "Tyr",  # Tyrosine
  X = "ter"   # Termination (stop codon)
)

################################################################################
# SECTION 3: HELPER FUNCTIONS - PROTEIN NOMENCLATURE
################################################################################

#' Convert Single-Letter to Three-Letter Amino Acid Codes
#'
#' @description
#' Converts protein change notation from compact single-letter format to 
#' more readable three-letter format following HGVS nomenclature standards.
#'
#' @param protein_change_str String containing protein change (e.g., "p.R135X")
#'
#' @return String with three-letter amino acid codes (e.g., "Arg135ter")
#'
#' @details
#' Handles multiple protein change formats:
#'   - Simple substitutions: R135X → Arg135ter
#'   - Frameshifts: K138Hfs*11 → Lys138Hisfs*11
#'   - Synonymous: T120= → Thr120=
#'   - Deletions/insertions: maintained as-is
#'
#' Design Decision: Why do this conversion?
#'   - More readable for clinical researchers
#'   - Standard in clinical genetics publications
#'   - Reduces ambiguity (e.g., "Ser" vs "Cys" when using S vs C)
#'
#' @examples
#' swap_one_letter_to_three_letter("p.R135X")     # → "Arg135ter"
#' swap_one_letter_to_three_letter("R135X")       # → "Arg135ter" (p. prefix optional)
#' swap_one_letter_to_three_letter("K138Hfs*11")  # → "Lys138Hisfs*11"
swap_one_letter_to_three_letter <- function(protein_change_str) {
  # Remove the "p." prefix if present (HGVS standard prefix for protein changes)
  # Why: We want to work with just "R135X" not "p.R135X" for easier parsing
  protein_change_str <- gsub("^p\\.", "", protein_change_str)
  
  # Define regex pattern to match protein change components
  # Pattern breakdown:
  #   ^([A-Z])           - Start amino acid (single letter)
  #   ([0-9]+)           - Position number (e.g., 135)
  #   ([A-Z](?:fs\\*\\d+)|[A-Z]|=|del|ins|dup|[a-z]+\\*\\d+)$ - End change
  #
  # Why this pattern? Covers main variant types:
  #   - Simple changes: R135X (Arg→stop)
  #   - Frameshifts: K138Hfs*11
  #   - Synonymous: T120=
  #   - Special cases: del, ins, dup
  match <- str_match(protein_change_str, "^([A-Z])([0-9]+)([A-Z](?:fs\\*\\d+)|[A-Z]|=|del|ins|dup|[a-z]+\\*\\d+)$")
  
  # If pattern matches, proceed with conversion
  if (!is.na(match[1,1])) {
    # Extract components
    start_aa <- match[1,2]  # First amino acid (e.g., "R")
    position <- match[1,3]  # Position (e.g., "135")
    change <- match[1,4]    # Change (e.g., "X" for stop)
    
    # Convert start amino acid to three-letter code
    # Design: Use lookup table for fast, reliable conversion
    if (start_aa %in% names(aa_code_map)) {
      start_aa_three <- aa_code_map[start_aa]
    } else {
      # If not in map (shouldn't happen), keep original
      start_aa_three <- start_aa
    }
    
    # Handle the change component
    # Check if change is amino acid + frameshift/stop notation
    # Pattern: ([A-Z]|\*) for amino acid or stop, (fs\*\d+|\*\d+)? for frameshift
    change_match <- str_match(change, "^([A-Z]|\\*)(fs\\*\\d+|\\*\\d+)?$")
    
    if (!is.na(change_match[1,1])) {
      # Extract the amino acid and frameshift components
      change_aa <- change_match[1,2]    # Amino acid (e.g., "X" or "H")
      fs_part <- change_match[1,3]      # Frameshift notation (e.g., "fs*11")
      
      # Convert change amino acid to three-letter code
      if (change_aa %in% names(aa_code_map)) {
        change_aa_three <- aa_code_map[change_aa]
      } else {
        change_aa_three <- change_aa
      }
      
      # Reconstruct with frameshift notation if present
      # Example: "H" + "fs*11" → "Hisfs*11"
      change_three <- paste0(change_aa_three, ifelse(is.na(fs_part), "", fs_part))
    } else {
      # If change is 'fs*number' without preceding amino acid,
      # or other special notation like "=", "del", "ins"
      # Keep as-is since no amino acid to convert
      change_three <- change
    }
    
    # Reconstruct the full protein change string
    # Example: "Arg" + "135" + "ter" = "Arg135ter"
    protein_change_three <- paste0(start_aa_three, position, change_three)
    return(protein_change_three)
  } else {
    # If pattern doesn't match (e.g., complex variants, structural changes),
    # return original string unchanged
    # Design Decision: Fail gracefully rather than error out
    return(protein_change_str)
  }
}

#' Extract and Format Protein Change
#'
#' @description
#' Wrapper function that handles NA values and calls the conversion function.
#'
#' @param protein_change Raw protein change string from CSV
#'
#' @return Formatted protein change with three-letter codes, or NA if input is NA
#'
#' @details
#' This is a simple wrapper to handle edge cases before conversion.
#' Separates data validation from conversion logic (single responsibility principle).
extract_protein_change <- function(protein_change) {
  # Handle missing/empty values gracefully
  # Why: Not all variants have protein predictions (e.g., intronic variants)
  if (is.na(protein_change) || protein_change == "") return(protein_change)
  
  # Swap one-letter codes for three-letter codes
  protein_change_three <- swap_one_letter_to_three_letter(protein_change)
  return(protein_change_three)
}

################################################################################
# SECTION 4: HELPER FUNCTIONS - COORDINATE EXTRACTION
################################################################################

#' Extract cDNA Coordinate from Variant String
#'
#' @description
#' Parses HGVS cDNA notation to extract the primary coordinate number.
#' Handles both single positions and ranges.
#'
#' @param variant HGVS cDNA notation (e.g., "c.333del", "c.190_200del")
#'
#' @return Numeric cDNA coordinate, or NA if unparseable
#'
#' @details
#' Handles multiple formats:
#'   - Simple: c.333del → 333
#'   - Range: c.190_200del → 190 (uses start position)
#'   - Intronic: c.190-2A>G → 190 (extracts exonic reference)
#'   - Invalid: c.GAIN → NA (no numeric coordinate)
#'
#' Design Decision: Why extract coordinates?
#'   - Needed to query Ensembl API for genome positions
#'   - cDNA coordinates are relative to transcript
#'   - IGV browser needs absolute genome coordinates
#'   - This bridges the gap between clinical nomenclature and visualization
#'
#' @examples
#' extract_cDNA_coordinate("c.333del")      # → 333
#' extract_cDNA_coordinate("c.190_200del")  # → 190 (start of range)
#' extract_cDNA_coordinate("c.190-2A>G")    # → 190 (intronic, uses exon ref)
#' extract_cDNA_coordinate("c.GAIN")        # → NA (no coordinate)
extract_cDNA_coordinate <- function(variant) {
  # Try to match simple format: c.NUMBER (e.g., c.333, c.490)
  # Pattern: c\. matches literal "c.", ([0-9]+) captures the number
  match <- str_match(variant, "c\\.([0-9]+)")
  
  if (!is.na(match[1, 2])) {
    # Successfully matched - return the coordinate as numeric
    # Why numeric? Needed for API query and sorting
    return(as.numeric(match[1, 2]))
  } else {
    # Didn't match simple format - try range format
    # Handle ranges like c.190_200 by taking the first part of the range
    # Why first position? 
    #   - Represents start of affected region
    #   - Ensembl API can map ranges, but we simplify to single point
    #   - IGV will show variant at this position (good enough for visualization)
    range_match <- str_match(variant, "c\\.([0-9]+)_([0-9]+)")
    
    if (!is.na(range_match[1, 2])) {
      # Return the first coordinate in the range as numeric
      return(as.numeric(range_match[1, 2]))
    } else {
      # Couldn't parse coordinate - return NA
      # This happens for:
      #   - Structural variants (e.g., "GAIN", "Entire coding sequence")
      #   - Complex rearrangements
      #   - Malformed strings
      # Design: Return NA rather than error - let downstream handle gracefully
      return(NA)
    }
  }
}

################################################################################
# SECTION 5: API INTEGRATION - ENSEMBL REST API
################################################################################

#' Query Ensembl REST API for Genome Positions
#'
#' @description
#' Converts cDNA coordinates (relative to transcript) to genome coordinates
#' (absolute positions on chromosome) using the Ensembl REST API.
#'
#' @param cDNA_coordinate Numeric cDNA position (e.g., 333)
#' @param transcript_id Ensembl transcript ID (default: ENST00000418600 for SYNGAP1)
#'
#' @return List with seqid (chromosome), start, and end positions
#'
#' @details
#' API Endpoint: https://rest.ensembl.org/map/cdna/{transcript}/{coordinate}..{coordinate}
#' 
#' Example Query:
#'   Input:  cDNA_coordinate = 333, transcript = ENST00000418600
#'   URL:    https://rest.ensembl.org/map/cdna/ENST00000418600/333..333
#'   Output: chr6:33,425,796-33,425,796 (hg38)
#'
#' Why this design?
#'   - Separates API logic from caching (see memoise wrapper below)
#'   - Returns consistent structure even on errors (list with NAs)
#'   - Logs errors but doesn't crash app (defensive programming)
#'
#' Rate Limits:
#'   - 15 requests/second per IP
#'   - 55,000 requests/hour per IP
#'   - Mitigated by memoization (cache results permanently)
#'
#' Error Handling:
#'   - Invalid coordinates → returns NA positions
#'   - Network errors → returns NA positions
#'   - Malformed responses → returns NA positions
#'   - Always returns valid structure for downstream processing
#'
#' @examples
#' get_genome_positions_function(333)     # → list(seqid="chr6", start=33425796, end=33425796)
#' get_genome_positions_function(NA)      # → list(seqid="chr6", start=NA, end=NA)
#' get_genome_positions_function(99999)   # → list(seqid="chr6", start=NA, end=NA) [out of bounds]
get_genome_positions_function <- function(cDNA_coordinate, transcript_id = "ENST00000418600") {
  # Handle NA coordinates immediately
  # Why check first? Avoids unnecessary API call for invalid data
  if (is.na(cDNA_coordinate)) {
    return(list(seqid = "chr6", start = NA, end = NA))
  }
  
  # Construct API request
  # Ensembl server base URL
  server <- "https://rest.ensembl.org"
  
  # Build endpoint path
  # Format: /map/cdna/{transcript_id}/{start}..{end}
  # Using same number for start and end gets a single position
  # Could use range (e.g., "333..335") for deletions, but single point is simpler
  ext <- paste0("/map/cdna/", transcript_id, "/", cDNA_coordinate, "..", cDNA_coordinate)
  
  # Make HTTP GET request
  # content_type("application/json") tells Ensembl we want JSON response
  # Why GET? Read-only operation, no data modification
  response <- GET(paste0(server, ext), content_type("application/json"))
  
  # Check if request succeeded
  # http_error() returns FALSE for status codes 200-299
  if (!http_error(response)) {
    # Parse JSON response
    # Ensembl returns complex nested structure - we need the "mappings" field
    result <- fromJSON(content(response, as = "text", encoding = "UTF-8"))
    
    # Verify we got mapping data
    # Structure check: Does result have "mappings" and is it non-empty?
    # Why check? Some coordinates may be valid but unmappable (intronic regions)
    if ("mappings" %in% names(result) && nrow(result$mappings) > 0) {
      # Extract genome coordinates from first mapping
      # [1] takes first row - usually only one mapping, but API can return multiple
      # for transcripts with alternative splicing
      start <- result$mappings$start[1]
      end <- result$mappings$end[1]
      
      # Return as named list
      # seqid = "chr6" hardcoded because SYNGAP1 is on chromosome 6
      # Could extract from API response, but hardcoding is simpler and faster
      return(list(seqid = "chr6", start = start, end = end))
    } else {
      # Mappings field missing or empty - coordinate might be intronic or invalid
      # Log warning for debugging but don't crash
      message("No mapping data found for cDNA coordinate ", cDNA_coordinate)
      return(list(seqid = "chr6", start = NA, end = NA))
    }
  } else {
    # HTTP error occurred (4xx or 5xx status code)
    # Could be: rate limit (429), not found (404), server error (500), etc.
    # Log the error with status code for debugging
    message("Error fetching genome position for coordinate ", cDNA_coordinate, 
            ". Status code: ", status_code(response))
    
    # Return NA positions - app will skip this variant in visualization
    # Design: Graceful degradation - show what we can, skip what we can't
    return(list(seqid = "chr6", start = NA, end = NA))
  }
}

################################################################################
# SECTION 6: PERFORMANCE OPTIMIZATION - CACHING WITH MEMOIZATION
################################################################################

#' Setup Cache Directory
#'
#' @description
#' Creates a disk-based cache for storing API responses.
#' 
#' Design Decision: Why disk cache instead of memory?
#'   - Persists across R sessions (app restarts don't lose cache)
#'   - Can handle large amounts of data without RAM limits
#'   - Shared across different runs of the app
#'   - Enables offline operation after first run
#'
#' Cache Structure:
#'   cache/genome_positions/
#'     ├── [hash1].rds  ← Response for coordinate 28
#'     ├── [hash2].rds  ← Response for coordinate 333
#'     └── ...          (one file per unique coordinate)
#'
#' Performance Impact:
#'   - First run:  ~143 API calls × 300ms = 43 seconds
#'   - Cached run: ~143 file reads × 2ms = 0.3 seconds
#'   - Speedup: 150x faster!

# Initialize disk cache
# cachem::cache_disk() creates cache directory if it doesn't exist
genome_position_cache <- cachem::cache_disk("cache/genome_positions")

#' Memoised Version of API Function
#'
#' @description
#' Wraps get_genome_positions_function with memoization to cache results.
#'
#' How Memoization Works:
#'   1. Function called with arguments (e.g., cDNA_coordinate=333)
#'   2. memoise() creates hash key from function name + arguments
#'   3. Checks if cache file exists for this key
#'   4a. Cache hit:  Read .rds file, return value (2ms)
#'   4b. Cache miss: Call original function, save result, return value (300ms)
#'
#' Design Decision: Why wrap in memoise?
#'   - Transparent caching - no changes to calling code
#'   - Automatic cache invalidation if function changes (hash includes function body)
#'   - Works with any function (generic solution)
#'   - Handles concurrent access safely
#'
#' Cache Persistence:
#'   - Files never expire (genomic coordinates don't change)
#'   - Survives R restarts
#'   - Can manually clear by deleting cache directory
#'
#' Usage:
#'   # First call (no cache) - queries API
#'   pos1 <- get_genome_positions(333)  # Takes 300ms
#'   
#'   # Second call (cached) - reads from disk
#'   pos2 <- get_genome_positions(333)  # Takes 2ms (150x faster!)
get_genome_positions <- memoise(get_genome_positions_function, 
                                cache = genome_position_cache)

################################################################################
# SECTION 7: DATA LOADING & PREPROCESSING
################################################################################

#' Load Variant Data
#'
#' @description
#' Loads patient variant data from CSV file.
#'
#' Design Decision: Why load at startup rather than in server?
#'   - Data doesn't change during app session
#'   - Loading once is more efficient than per-user
#'   - Makes debugging easier (data loaded before UI interactions)
#'   - Preprocessing happens once, not per user connection
#'
#' File Format: updatedCitizen191.csv
#'   - 153 rows (variants)
#'   - 33 columns (metadata, patient info, research assets)
#'   - Key columns:
#'     - SYNGAP1.variant: cDNA notation (e.g., "c.333del")
#'     - variant.type: Variant classification
#'     - predicted.protein: Protein change
#'     - biorepository, iPSC.line, mouse.line: Research assets
#'     - X..patients.in.Citizen.Health: Number of patients

# Load CSV data
# stringsAsFactors = FALSE prevents automatic factor conversion
# Why? We want to manipulate strings freely without factor level constraints
variant_data <- read.csv("updatedCitizen191.csv", stringsAsFactors = FALSE)

#' Standardize Variant Type Names
#'
#' @description
#' Normalizes variant type nomenclature for consistency.
#'
#' Problem: CSV has inconsistent naming
#'   - "missense VUS" vs "VUS - missense"
#'   - "frameshift" vs "frameshift deletion" vs "frameshift insertion"
#'   - Mixed case
#'
#' Solution: Map all variations to canonical names
#'
#' Why important?
#'   - Consistent button labels in UI
#'   - Reliable filtering
#'   - Correct color assignments
#'   - Clean track names in IGV

variant_data <- variant_data %>%
  mutate(
    # First, convert all to lowercase for case-insensitive matching
    variant.type = tolower(variant.type),
    
    # Then map to canonical names using case_when (like switch statement)
    # Order matters: more specific patterns first
    variant.type = case_when(
      # Missense variants (amino acid substitutions)
      variant.type %in% c('missense') ~ 'missense',
      
      # Missense VUS (Variants of Uncertain Significance)
      # Why separate? Different clinical interpretation, shown in orange
      variant.type %in% c('missense vus', 'vus - missense') ~ 'missense-VUS',
      
      # Nonsense variants (premature stop codons)
      variant.type %in% c('nonsense', 'nonsense vus') ~ 'nonsense',
      
      # Frameshift variants (reading frame disruptions)
      # Groups all frameshift subtypes together
      variant.type %in% c('frameshift', 'frameshift deletion', 'frameshift insertion') ~ 'frameshift',
      
      # Indels (small insertions/deletions maintaining frame)
      variant.type %in% c('indel', 'insertion deletion') ~ 'indel',
      
      # Copy number gains
      variant.type %in% c('gain', 'gain exon 3 lp') ~ 'gain',
      
      # VUS (general uncertain significance)
      variant.type %in% c('vus', 'vus - gain entire seq') ~ 'vus',
      
      # Intronic variants (outside coding sequence)
      variant.type %in% c('intronic') ~ 'intronic',
      
      # Copy number losses
      variant.type %in% c('loss') ~ 'loss',
      
      # Default: keep original if no match
      # Handles new types without code changes
      TRUE ~ variant.type
    )
  )

# Convert to factor for efficient storage and categorical operations
# Why factor? 
#   - Memory efficient (stores as integers with labels)
#   - Enforces valid values (can't accidentally add typos)
#   - Useful for grouping and counting
variant_data$variant.type <- factor(variant_data$variant.type)

#' Precompute Derived Fields
#'
#' @description
#' Extract and format data needed for visualization before app starts.
#'
#' Design Philosophy: Do expensive work once at startup
#'   - Parsing happens once, not every time a user filters
#'   - Results stored in data frame for fast access
#'   - Users get instant UI response
#'
#' What's being precomputed:
#'   1. cDNA coordinates (for API queries)
#'   2. Formatted protein changes (for display)

variant_data <- variant_data %>%
  mutate(
    # Extract numeric cDNA coordinate from variant string
    # sapply applies function to each row (vectorized operation)
    # Example: "c.333del" → 333
    cDNA_coordinate = sapply(SYNGAP1.variant, extract_cDNA_coordinate),
    
    # Format protein changes with three-letter amino acids
    # Example: "p.R135X" → "Arg135ter"
    protein_change = sapply(predicted.protein, extract_protein_change)
  )

#' Fetch Genome Positions via API
#'
#' @description
#' Convert all cDNA coordinates to genome positions using Ensembl API.
#'
#' Critical Performance Section!
#' This is where memoization saves 40+ seconds on subsequent runs.
#'
#' Strategy: Batch prefetch at startup
#'   - Get unique coordinates (143 from 153 variants - some share coordinates)
#'   - Query API once per unique coordinate
#'   - Cache results for future runs
#'   - Merge back into main data frame
#'
#' Why prefetch vs on-demand?
#'   ✓ All data ready before user interaction
#'   ✓ No delays when clicking buttons
#'   ✓ Progress visible in console during first load
#'   ✓ Easier to handle errors (all at once vs scattered)
#'   ✗ Slight delay on first startup (acceptable tradeoff)

# Get list of unique cDNA coordinates
# Why unique()? Same coordinate may appear in multiple variants
# Example: c.490C>T appears in 8 patients, but only query API once
unique_coords <- unique(variant_data$cDNA_coordinate)

# Fetch positions for all unique coordinates
# lapply returns a list of results, one per coordinate
# get_genome_positions is memoised, so:
#   - First run: Makes 143 API calls (cached for future)
#   - Later runs: Reads 143 cached files (instant!)
#
# Progress: Watch console during first run to see API queries
# You'll see messages like "No mapping data found for..." for invalid coordinates
positions <- lapply(unique_coords, function(coord) {
  get_genome_positions(coord, transcript_id = "ENST00000418600")
})

# Convert list of positions to data frame for merging
# sapply extracts each field (start, end) from the list
# Creates a data frame with columns: cDNA_coordinate, seqid, start, end
positions_df <- data.frame(
  cDNA_coordinate = unique_coords,
  seqid = "chr6",  # SYNGAP1 is on chromosome 6 (hardcoded for simplicity)
  start = sapply(positions, function(pos) pos$start),
  end = sapply(positions, function(pos) pos$end),
  stringsAsFactors = FALSE
)

# Merge genome positions back into main variant data
# Left join: Keep all variants, add positions where available
# by = "cDNA_coordinate" matches rows between tables
# all.x = TRUE keeps variants even if position lookup failed (they get NA)
#
# Result: variant_data now has columns: ..., cDNA_coordinate, seqid, start, end
# This is what IGV needs to display variants on the genome!
variant_data <- merge(variant_data, positions_df, by = "cDNA_coordinate", all.x = TRUE)

#' Define Color Scheme for Variant Types
#'
#' @description
#' Color palette for different mutation types in IGV browser.
#'
#' Design Considerations:
#'   - Colorblind-friendly (avoid red/green for critical distinctions)
#'   - Sufficient contrast for visibility
#'   - Semantic meaning where possible
#'   - Consistent with scientific conventions
#'
#' Color Choices:
#'   - Blue (#3E6CB5): Missense - common, moderate impact
#'   - Orange (#FFA500): Missense-VUS - uncertain, needs attention
#'   - Green (#19AE69): Nonsense - clearly pathogenic
#'   - Purple (#6F3592): Frameshift - severe disruption
#'   - Light Blue (#D0DDEE): Indel - variable impact
#'   - Light Green (#CBE7D4): Gain - structural variant
#'   - Lavender (#D4CAE1): VUS - uncertain
#'   - Gray (#4E4E4E): Intronic/Loss - likely benign or uncertain

color_table <- list(
  missense  = "#3E6CB5",      # Blue
  "missense-VUS" = "#FFA500",  # Orange (attention-grabbing for VUS)
  nonsense  = "#19AE69",      # Green
  frameshift = "#6F3592",     # Purple
  indel     = "#D0DDEE",      # Light Blue
  gain      = "#CBE7D4",      # Light Green
  vus       = "#D4CAE1",      # Lavender
  intronic  = "#4E4E4E",      # Gray
  loss      = "#4E4E4E",      # Gray
  clinvar   = "#00BCD4"       # Teal - visually distinct from all internal track colors
)

################################################################################
# SECTION 8: USER INTERFACE DEFINITION
################################################################################

#' Shiny UI Layout
#'
#' @description
#' Defines the visual layout and interactive elements of the web application.
#'
#' Architecture: Two-column layout
#'   - Left sidebar (30% width): Controls and filters
#'   - Right main panel (70% width): IGV genome browser
#'
#' Design Philosophy:
#'   - Put controls on left (common web convention)
#'   - Maximize space for visualization
#'   - Group related controls together
#'   - Clear visual hierarchy

ui <- shinyUI(fluidPage(
  # Application title
  # Appears at top of page
  titlePanel("SVV for SYNGAP1 Variant Viewer"),
  
  # Two-column layout with sidebar and main content
  sidebarLayout(
    
    # ============================================================
    # LEFT SIDEBAR: Controls
    # ============================================================
    sidebarPanel(
      # Section header for variant type buttons
      h3("Add Tracks"),
      
      # Dynamic UI element: Buttons generated based on available variant types
      # Why dynamic? 
      #   - Number of variant types varies by dataset
      #   - Filters change which types are visible
      #   - Automatically adapts to data
      # Rendered by server (see output$trackButtons below)
      uiOutput("trackButtons"),
      
      # ----------------------------------------
      # ClinVar Track Button (static)
      # ----------------------------------------
      # This button is static (not part of the dynamic variant-type loop)
      # because ClinVar is an external reference dataset, not a filtered
      # subset of the internal Citizen Health variant data.
      # It always appears regardless of which filters are active.
      hr(style = "border-top: 1px solid #ccc; margin: 8px 0;"),
      actionButton(
        inputId = "addTrack_clinvar",
        label   = "ClinVar"
      ),
      
      # ClinVar size filter slider
      # Filters ClinVar variants shown on the track by their genomic span.
      #
      # Why sliderTextInput() instead of sliderInput()?
      #   ClinVar spans an enormous size range: ~1 bp SNVs up to ~46,632 KB
      #   whole-chromosome CNVs. A linear slider across that range would make
      #   the low end (where most interesting variants live) nearly impossible
      #   to control. sliderTextInput() accepts an explicit vector of stop
      #   values, so we can space them logarithmically: fine steps at the bottom
      #   (1, 2, 3 … 10 KB) and progressively coarser steps toward the top
      #   (5,000 KB steps above 10,000 KB). Every stop is a "round" KB value
      #   that makes clinical sense.
      #
      # Stop sequence (45 positions total):
      #   1–9 KB      : step 1 KB   (fine control for point variants / small indels)
      #   10–90 KB    : step 10 KB  (single-exon to multi-exon deletions)
      #   100–900 KB  : step 100 KB (sub-megabase SVs)
      #   1,000–9,000 KB : step 1,000 KB (megabase-scale CNVs)
      #   10,000–50,000 KB : step 5,000 KB (cytogenetic-band CNVs)
      #
      # Default: 100 KB — excludes whole-gene CNVs while retaining variants
      #   within the SYNGAP1 locus (~37 kb). Users drag right to reveal
      #   progressively larger structural variants.
      #
      # Behavior: If the ClinVar track is already displayed, moving the slider
      # immediately reloads it with only variants at or below the new threshold.
      sliderTextInput(
        inputId  = "clinvar_size_kb",
        label    = "Max ClinVar variant size (KB)",
        choices  = c(
          seq(1,     9,     by = 1),      #  1–9 KB    : step 1 KB
          seq(10,    90,    by = 10),     # 10–90 KB   : step 10 KB
          seq(100,   900,   by = 100),    # 100–900 KB : step 100 KB
          seq(1000,  9000,  by = 1000),   # 1–9 MB     : step 1,000 KB
          seq(10000, 50000, by = 5000)    # 10–50 MB   : step 5,000 KB
        ),
        selected = 100,
        grid     = FALSE,
        width    = "100%"
      ),
      
      # Spacer for visual separation
      # Improves readability by grouping buttons vs checkboxes
      
      # ========================================
      # Research Asset Filters
      # ========================================
      # Why these filters?
      #   - Help researchers find variants with available resources
      #   - Enables focused analysis on well-characterized variants
      #   - Supports research planning and collaboration
      
      # Filter 1: Biorepository samples
      # Shows only variants where biological samples are stored
      # Use case: "Which variants can we get DNA/RNA for?"
      div(
        class = "form-group",
        checkboxInput("biorepository_filter", 
                      "Has biorepository samples", 
                      value = FALSE)  # Default unchecked (show all)
      ),
      
      # Filter 2: iPSC cell lines
      # Shows only variants with induced pluripotent stem cell lines
      # Use case: "Which variants can we study in cell culture?"
      div(
        class = "form-group",
        checkboxInput("iPSC_line_filter", 
                      "Cell line available", 
                      value = FALSE)
      ),
      
      # Filter 3: Mouse models
      # Shows only variants with corresponding mouse lines
      # Use case: "Which variants have animal models for testing therapies?"
      div(
        class = "form-group",
        checkboxInput("Mouse_line_filter", 
                      "Mouse line available", 
                      value = FALSE)
      ),
      
      # Set sidebar width to 30% of page (3 out of 12 columns)
      # Why 30%? 
      #   - Enough space for controls without crowding
      #   - Leaves 70% for genome browser (where the action is)
      width = 3
    ),
    
    # ============================================================
    # RIGHT MAIN PANEL: Visualization
    # ============================================================
    mainPanel(
      # IGV genome browser widget
      # This is where variants are visualized on chromosome 6
      #
      # Parameters:
      #   - id: 'igvShiny_0' - unique identifier for this instance
      #   - height: "700px" - tall enough to see multiple tracks
      #
      # Rendered by server (see output$igvShiny_0 below)
      igvShinyOutput('igvShiny_0', height = "700px"),
      
      # Set main panel width to 70% of page (9 out of 12 columns)
      width = 9
    )
  ),
  
  # ============================================================
  # FOOTER: Citation
  # ============================================================
  # Horizontal rule for visual separation
  hr(),
  
  # Citation div with centered text and small font
  # Why include citation?
  #   - Give credit to igvShiny developers
  #   - Help users find original tool for their own projects
  #   - Professional appearance
  tags$div(
    style = "text-align: center; font-size: 12px;",
    tags$p("Citation:"),
    tags$p(
      "Shannon P, Gladki A, Scigocka K (2024). ",
      tags$em("igvShiny: igvShiny: a wrapper of Integrative Genomics Viewer "),
      "(IGV - an interactive tool for visualization and exploration integrated genomic data). ",
      "R package version 1.2.0, ",
      tags$a(href = "https://gladkia.github.io/igvShiny/", 
             "https://gladkia.github.io/igvShiny/"), ", ",
      tags$a(href = "https://github.com/gladkia/igvShiny", 
             "https://github.com/gladkia/igvShiny.")
    )
  )
))

################################################################################
# SECTION 9: DATA FORMATTING FOR IGV
################################################################################

#' Create GFF3-formatted Data for IGV Tracks
#'
#' @description
#' Converts variant data to GFF3 format required by IGV browser.
#'
#' @param data Filtered variant data frame
#' @param variant_type Type of variant for this track (e.g., "missense")
#'
#' @return Data frame in GFF3 format for igvShiny
#'
#' @details
#' GFF3 Format Requirements:
#'   1. seqid    - Chromosome (chr6)
#'   2. source   - Data source (not used, set to ".")
#'   3. type     - Feature type (variant type)
#'   4. start    - Genomic start position
#'   5. end      - Genomic end position
#'   6. score    - Confidence score (not used, set to ".")
#'   7. strand   - DNA strand (+/-)
#'   8. phase    - Coding phase (not used for variants, set to ".")
#'   9. attributes - Semicolon-separated key=value pairs
#'
#' Attributes Field:
#'   Contains metadata displayed when user clicks variant:
#'   - ID: Unique identifier
#'   - Name: Display name (protein change if available, else cDNA)
#'   - Description: Full cDNA notation
#'   - Number of patients: Clinical relevance indicator
#'
#' Design Decision: Why GFF3?
#'   - Standard genomics format (widely understood)
#'   - Supported natively by IGV
#'   - Flexible attributes system
#'   - Human-readable text format

create_gff3_data <- function(data, variant_type) {
  data <- data %>%
    # Remove variants without genome positions
    # Why filter? IGV can't display variants without coordinates
    # These are typically structural variants or failed API queries
    filter(!is.na(start)) %>%
    
    mutate(
      # Column 1: Chromosome
      seqid = "chr6",
      
      # Column 2: Source
      # "." indicates no specific source
      # Could use "patient_database" or "clinical_variants" if needed
      source = ".",
      
      # Column 3: Feature type
      # Uses the variant type (missense, nonsense, etc.)
      # This is what gets colored in IGV based on color_table
      type = variant_type,
      
      # Column 6: Score
      # "." indicates no score
      # Could add CADD scores, conservation scores, etc. if available
      score = ".",
      
      # Column 7: Strand
      # "+" for forward strand (SYNGAP1 is on forward strand of chr6)
      # If gene were on reverse strand, this would be "-"
      strand = "+",
      
      # Column 8: Phase
      # "." for not applicable (phase is for CDS features, not variants)
      phase = ".",
      
      # Column 9: Attributes (most important for display!)
      # Build different attribute strings based on data availability
      attributes = ifelse(
        # Check if protein change is available
        is.na(protein_change) | protein_change == "",
        
        # NO protein change - use cDNA notation only
        # This happens for intronic variants, synonymous variants, etc.
        paste0(
          "ID=", SYNGAP1.variant,           # Unique ID
          ";Name=", SYNGAP1.variant,        # Display name
          ";Description=cDNA: ", SYNGAP1.variant,  # Tooltip text
          ";Number of patients=", X..patients.in.Citizen.Health  # Patient count
        ),
        
        # YES protein change - show protein change as name
        # This is what users typically want to see
        paste0(
          "ID=", protein_change,            # Unique ID (protein level)
          ";Name=", protein_change,         # Display name (e.g., "Arg135ter")
          ";Description=cDNA: ", SYNGAP1.variant,  # Still show cDNA in tooltip
          ";Number of patients=", X..patients.in.Citizen.Health  # Patient count
        )
      )
    ) %>%
    # Select and order columns to match GFF3 specification
    # Order is critical - IGV expects columns in this exact order
    select(seqid, source, type, start, end, score, strand, phase, attributes)
  
  return(data)
}

################################################################################
# SECTION 9B: CLINVAR DATA FORMATTING FOR IGV
################################################################################

#' Create GFF3-formatted Data from ClinVar API Data for IGV Tracks
#'
#' @description
#' Converts the parsed ClinVar data frame into GFF3 format required by IGV.
#' All available ClinVar metadata fields are packed into the GFF3 attributes
#' column so they appear in IGV's feature detail pop-up when a user clicks
#' a variant.
#'
#' @param data The clinvar_data data frame (fetched at startup via NCBI E-utilities)
#' @param max_size_kb Numeric. Maximum variant span in kilobases to include in
#'   the track. Variants whose (end - start) exceeds this threshold are excluded.
#'   Defaults to 100 KB. Driven by the "Max ClinVar variant size" slider in the
#'   UI; pass as.numeric(input$clinvar_size_kb) from the server.
#'
#' @return Data frame in GFF3 format for use with loadGFF3TrackFromLocalData()
#'
#' @details
#' GFF3 columns produced:
#'   seqid      - "chr6" (all SYNGAP1 variants are on chromosome 6)
#'   source     - "ClinVar"
#'   type       - "clinvar" (drives teal color via color_table)
#'   start      - GRCh38 start position (parsed from GRCh38Location)
#'   end        - GRCh38 end position
#'   score      - "."
#'   strand     - "+"
#'   phase      - "."
#'   attributes - All ClinVar metadata fields as key=value pairs
#'
#' Attributes included (all available ClinVar fields):
#'   Name, ProteinChange, Condition, Accession, VariationID, AlleleID,
#'   dbSNP_ID, CanonicalSPDI, VariantType, MolecularConsequence,
#'   GermlineClassification, GermlineLastEvaluated, GermlineReviewStatus,
#'   SomaticClinicalImpact, SomaticLastEvaluated, SomaticReviewStatus,
#'   OncogenicityClassification, OncogenicityLastEvaluated,
#'   OncogenicityReviewStatus
#'
#' Design Decision: Why include all fields?
#'   - ClinVar metadata is rich and clinically meaningful
#'   - Users can click a variant to see full clinical context
#'   - Avoids information loss from the source dataset
#'   - Mirrors IGV's standard ClinVar track behavior

create_clinvar_gff3_data <- function(data, max_size_kb = 100) {
  
  # Guard: return empty GFF3 frame if ClinVar data failed to load.
  # This happens when both the NCBI fetch and the cache are unavailable,
  # leaving clinvar_data as a zero-column data.frame().
  # Without this guard, dplyr::filter() crashes trying to reference
  # clinvar_start on a frame that has no columns at all.
  required_cols <- c("clinvar_start", "clinvar_end", "clinvar_chr")
  if (nrow(data) == 0 || !all(required_cols %in% names(data))) {
    message("ClinVar: no data available — track will be empty.")
    return(data.frame())
  }
  
  # Helper: sanitize a value for GFF3 attribute string
  # Replaces semicolons and equals signs which are GFF3 delimiters,
  # and trims whitespace. Returns "." for NA/empty values.
  clean_attr <- function(x) {
    x <- as.character(x)
    x <- trimws(x)
    x[is.na(x) | x == ""] <- "."
    # Escape GFF3 reserved characters inside values
    x <- gsub(";", "%3B", x, fixed = TRUE)
    x <- gsub("=", "%3D", x, fixed = TRUE)
    x
  }
  
  # Convert KB threshold to base pairs for coordinate arithmetic
  max_size_bp <- max_size_kb * 1000
  
  data %>%
    # Remove variants without valid GRCh38 coordinates
    # These are typically older submissions that predate hg38 mapping
    filter(!is.na(clinvar_start) & !is.na(clinvar_end)) %>%
    
    # Remove variants whose span exceeds the user-specified size threshold.
    #
    # Why: Large copy number variants (e.g. chr6:156,974-46,789,291) span
    # millions of base pairs and physically cover the entire SYNGAP1 locus.
    # At every zoom level they sit on top of all smaller point variants,
    # intercepting every click and making the underlying variants unreachable.
    #
    # The threshold is controlled by the "Max variant size" slider in the UI
    # (default 100 kb). Lowering it progressively hides larger SVs and CNVs,
    # letting the user focus on point variants and small indels.
    #   - SYNGAP1 itself is ~37 kb, so 100 kb already excludes whole-gene CNVs.
    #   - Large CNVs spanning cytogenetic bands are not interpretable at this
    #     zoom level and are better reviewed directly in ClinVar's own browser.
    filter((clinvar_end - clinvar_start) <= max_size_bp) %>%
    
    mutate(
      # --- Required GFF3 columns ---
      seqid  = clinvar_chr,          # Chromosome (e.g. "chr6")
      source = "ClinVar",            # Data provenance label
      type   = "clinvar",            # Drives teal color via color_table
      start  = clinvar_start,
      end    = clinvar_end,
      score  = ".",
      strand = "+",                  # SYNGAP1 is on forward strand
      phase  = ".",
      
      # --- Attributes column: all ClinVar metadata ---
      # IGV displays these as a table when the user clicks a feature.
      # Fields are semicolon-separated key=value pairs per GFF3 spec.
      attributes = paste0(
        "ID=",                              clean_attr(Accession),
        ";Name=",                           clean_attr(Name),
        ";ProteinChange=",                  clean_attr(Protein.change),
        ";Condition=",                      clean_attr(Condition.s.),
        ";Accession=",                      clean_attr(Accession),
        ";VariationID=",                    clean_attr(VariationID),
        ";AlleleID=",                       clean_attr(AlleleID.s.),
        ";dbSNP_ID=",                       clean_attr(dbSNP.ID),
        ";CanonicalSPDI=",                  clean_attr(Canonical.SPDI),
        ";VariantType=",                    clean_attr(Variant.type),
        ";MolecularConsequence=",           clean_attr(Molecular.consequence),
        ";GermlineClassification=",         clean_attr(Germline.classification),
        ";GermlineLastEvaluated=",          clean_attr(Germline.date.last.evaluated),
        ";GermlineReviewStatus=",           clean_attr(Germline.review.status),
        ";SomaticClinicalImpact=",          clean_attr(Somatic.clinical.impact),
        ";SomaticLastEvaluated=",           clean_attr(Somatic.clinical.impact.date.last.evaluated),
        ";SomaticReviewStatus=",            clean_attr(Somatic.clinical.impact.review.status),
        ";OncogenicityClassification=",     clean_attr(Oncogenicity.classification),
        ";OncogenicityLastEvaluated=",      clean_attr(Oncogenicity.date.last.evaluated),
        ";OncogenicityReviewStatus=",       clean_attr(Oncogenicity.review.status)
      )
    ) %>%
    # Select GFF3 columns in the required order
    select(seqid, source, type, start, end, score, strand, phase, attributes)
}

################################################################################
# SECTION 10: SERVER LOGIC (REACTIVE PROGRAMMING)
################################################################################

#' Shiny Server Function
#'
#' @description
#' Defines the reactive logic that powers the application.
#'
#' Reactive Programming Concepts:
#'   - Reactive values: Change over time (e.g., checkbox state)
#'   - Reactive expressions: Auto-recalculate when inputs change
#'   - Observers: Execute side effects when inputs change
#'
#' Server Architecture:
#'   1. State management (track which buttons were clicked)
#'   2. Data filtering (respond to checkbox changes)
#'   3. UI generation (create buttons based on filtered data)
#'   4. Event handling (respond to button clicks)
#'   5. Track updates (refresh IGV when filters change)
#'   6. IGV initialization (set up genome browser)

server <- function(input, output, session) {
  
  # ============================================================
  # STATE MANAGEMENT
  # ============================================================
  
  #' Track State Tracker
  #'
  # @description
  # Keeps track of which variant type tracks have been added to IGV.
  #
  # Why needed?
  #   - Remember user's track selections across filter changes
  #   - Avoid re-adding tracks that are already displayed
  #   - Enable smart updates (only refresh visible tracks)
  #
  # Data structure: Named list
  #   Example: list(missense = TRUE, nonsense = TRUE, frameshift = FALSE)
  #   TRUE = track is currently displayed
  #   FALSE or missing = track not displayed
  #
  # reactiveVal() creates a reactive variable that can be read and set
  # Like a reactive version of a regular variable
  added_tracks <- reactiveVal(list())
  
  # ClinVar Track State Tracker
  #
  # Separate from added_tracks because ClinVar is not part of the dynamic
  # variant-type loop. This flag lets the slider observer know whether to
  # reload the ClinVar track when the size threshold changes.
  #   TRUE  = ClinVar track is currently displayed in IGV
  #   FALSE = ClinVar track has not been added yet (slider changes are no-ops)
  clinvar_track_added <- reactiveVal(FALSE)
  
  # ============================================================
  # DATA FILTERING (REACTIVE EXPRESSION)
  # ============================================================
  
  #' Filtered Variant Data
  #'
  # @description
  # Reactive expression that filters variant data based on checkbox selections.
  #
  # Reactive Expression Behavior:
  #   - Automatically recalculates when inputs change
  #   - Caches result until inputs change (efficient)
  #   - Can be used by multiple observers/outputs
  #
  # Why reactive()? 
  #   - Avoid duplicating filter logic
  #   - Automatic dependency tracking
  #   - Efficient re-computation (only when needed)
  #
  # Filter Logic:
  #   - Start with all variants
  #   - Apply each filter if checkbox is checked
  #   - Filters are cumulative (AND logic)
  #   - Result: Variants matching all selected criteria
  
  filtered_data <- reactive({
    # Start with full dataset
    data <- variant_data
    
    # Filter 1: Biorepository samples
    # Only applies if checkbox is checked
    if (input$biorepository_filter) {
      # Keep only rows where biorepository column has data
      # !is.na() checks for non-missing values
      # != "" checks for non-empty strings
      # Why both? Some missing data is NA, some is empty string
      data <- data %>% filter(!is.na(biorepository) & biorepository != "")
    }
    
    # Filter 2: iPSC cell lines
    if (input$iPSC_line_filter) {
      # Keep only rows where iPSC.line column has data
      # Note: Column name has dot (iPSC.line) not underscore
      data <- data %>% filter(!is.na(iPSC.line) & iPSC.line != "")
    }
    
    # Filter 3: Mouse models
    if (input$Mouse_line_filter) {
      # Keep only rows where mouse.line column has data
      data <- data %>% filter(!is.na(mouse.line) & mouse.line != "")
    }
    
    # Return filtered data
    # Result depends on which checkboxes are checked:
    #   - None checked: All 153 variants
    #   - One checked: Subset with that resource
    #   - Multiple checked: Intersection (variants with ALL resources)
    data
  })
  
  # ============================================================
  # DYNAMIC UI GENERATION
  # ============================================================
  
  #' Generate Variant Type Buttons
  #'
  # @description
  # Creates action buttons dynamically based on filtered variant types.
  #
  # Why dynamic?
  #   - Number and types of variants change with filters
  #   - Button list should reflect what's actually displayable
  #   - Automatic adaptation to data changes
  #
  # renderUI() creates UI elements reactively
  #   - Re-runs when filtered_data() changes
  #   - Updates button list in real-time
  
  output$trackButtons <- renderUI({
    
    # Get unique variant types from filtered data
    # As filters change, this list changes too
    types <- unique(filtered_data()$variant.type)
    
    # Exclude problematic types
    # Why exclude?
    #   - "gain exon 3" might be malformed or edge case
    #   - "extra copy" is structural variant without coordinates
    # Better to skip than show broken tracks
    types <- types[types != "gain exon 3"]
    types <- types[types != "extra copy"]
    
    # Create one button per variant type
    # lapply returns a list of button widgets
    buttons <- lapply(types, function(type) {
      # Create action button
      # ID: Unique identifier using variant type
      #     make.names() ensures valid R names (no spaces/special chars)
      #     Example: "missense-VUS" → "missense.VUS"
      # Label: Display text (what user sees)
      #        Shows the variant type name
      actionButton(inputId = paste0("addTrack_", make.names(type)), 
                   label = type)
    })
    
    # Convert list of buttons to tagList
    # Why tagList? Shiny's way of combining multiple UI elements
    # Renders as a vertical stack of buttons
    do.call(tagList, buttons)
  })
  
  # ============================================================
  # EVENT HANDLING: Button Clicks
  # ============================================================
  
  #' Handle Track Addition Button Clicks
  #'
  # @description
  # Responds to user clicking variant type buttons by adding IGV tracks.
  #
  # Challenge: How to observe dynamically created buttons?
  # Solution: observe() with loop over all variant types
  #
  # Design Pattern: Closure with local()
  #   - Loop variable 'type' would be shared across iterations
  #   - local() creates isolated scope
  #   - Each observeEvent gets its own copy of 't'
  #   - Prevents bugs from variable sharing
  
  observe({
    # Get all variant types (not just filtered ones)
    # Why all? Button IDs exist even when hidden by filters
    types <- unique(variant_data$variant.type)
    
    # Create observer for each variant type
    for (type in types) {
      # local() creates isolated scope to avoid closure issues
      local({
        # Copy loop variable to local scope
        # Each iteration gets its own 't' variable
        # Without this, all observers would share the same 'type'
        t <- type
        
        # Create observer for this specific button
        # Fires when button is clicked
        observeEvent(input[[paste0("addTrack_", make.names(t))]], {
          
          # Update state: Mark this track as added
          # Get current state
          added_tracks_list <- added_tracks()
          # Set this type to TRUE
          added_tracks_list[[t]] <- TRUE
          # Save updated state
          added_tracks(added_tracks_list)
          
          # Get data for this variant type
          # Uses filtered_data() so respects active filters
          track_data <- filtered_data() %>% filter(variant.type == t)
          
          # Check if we have any data to display
          if (nrow(track_data) > 0) {
            # YES - variants to show
            
            # Convert to GFF3 format
            gff3_data <- create_gff3_data(track_data, variant_type = t)
            
            # Load track into IGV
            # Why loadGFF3TrackFromLocalData?
            #   - Loads from R data frame (not file)
            #   - Fast (no disk I/O)
            #   - Dynamic (updates in real-time)
            loadGFF3TrackFromLocalData(
              session,                          # Shiny session object
              id = "igvShiny_0",               # IGV widget ID
              trackName = t,                   # Track name (e.g., "missense")
              tbl = gff3_data,                 # Data in GFF3 format
              colorTable = color_table,        # Color scheme
              colorByAttribute = "type",       # Color by variant type
              trackHeight = 40,                # Height in pixels
              displayMode = "EXPANDED",        # Show all variants
              visibilityWindow = 1e8           # Show at all zoom levels (100Mb)
            )
          } else {
            # NO - filtered out all variants of this type
            # Load empty track to show "no data" rather than error
            # Why load empty track?
            #   - User feedback (shows track exists but filtered out)
            #   - Prevents confusion (track name visible, just no data)
            #   - Consistent behavior (all clicks produce visible result)
            loadGFF3TrackFromLocalData(
              session,
              id = "igvShiny_0",
              trackName = t,
              tbl = data.frame(),              # Empty data frame
              colorTable = color_table,
              colorByAttribute = "type",
              trackHeight = 40,
              displayMode = "EXPANDED",
              visibilityWindow = 1e8
            )
          }
          
          # Re-center view on SYNGAP1 gene after adding track
          # Why? Ensures user sees the gene region, not random chromosome region
          # "SYNGAP1" is a gene symbol - IGV knows where it is
          showGenomicRegion(session, id = "igvShiny_0", region = "SYNGAP1")
          
        }, ignoreInit = TRUE)  # Don't fire on initial app load, only on clicks
      })
    }
  })
  
  # ============================================================
  # CLINVAR TRACK BUTTON HANDLER
  # ============================================================
  
  #' Handle ClinVar Track Button Click
  #'
  # @description
  # Loads the ClinVar variant track into IGV when the user clicks "ClinVar".
  #
  # Design Notes:
  #   - Handled separately from the dynamic variant-type loop because ClinVar
  #     is an external reference dataset, not a filtered subset of internal data.
  #   - Not affected by the biorepository / cell-line / mouse-line checkboxes.
  #   - Uses pre-parsed clinvar_data (built at startup) for instant response.
  #   - Passes the current slider value to create_clinvar_gff3_data() so the
  #     initial load already respects whatever size threshold is set.
  #   - Sets clinvar_track_added(TRUE) so the slider observer knows to react
  #     to future threshold changes.
  #   - Track name "ClinVar" is used so it appears as a labeled track in IGV.
  
  observeEvent(input$addTrack_clinvar, {
    
    # Mark ClinVar track as active so the slider observer can reload it
    clinvar_track_added(TRUE)
    
    # Build GFF3 from parsed ClinVar data, applying the current size threshold.
    # as.numeric() is required because sliderTextInput() returns a character.
    clinvar_gff3 <- create_clinvar_gff3_data(clinvar_data,
                                             max_size_kb = as.numeric(input$clinvar_size_kb))
    
    if (nrow(clinvar_gff3) > 0) {
      # Load ClinVar track into IGV
      # trackName = "ClinVar" sets the label visible in the IGV track panel
      loadGFF3TrackFromLocalData(
        session,
        id                = "igvShiny_0",
        trackName         = "ClinVar",
        tbl               = clinvar_gff3,
        colorTable        = color_table,
        colorByAttribute  = "type",       # "clinvar" type → teal (#00BCD4)
        trackHeight       = 50,           # Slightly taller than internal tracks
        displayMode       = "EXPANDED",   # Show all variants (not collapsed)
        visibilityWindow  = 1e8           # Visible at all zoom levels
      )
    } else {
      # Fallback: load empty track if all rows lacked GRCh38 coordinates
      # Unlikely with a fresh ClinVar download, but handles gracefully
      loadGFF3TrackFromLocalData(
        session,
        id                = "igvShiny_0",
        trackName         = "ClinVar",
        tbl               = data.frame(),
        colorTable        = color_table,
        colorByAttribute  = "type",
        trackHeight       = 50,
        displayMode       = "EXPANDED",
        visibilityWindow  = 1e8
      )
    }
    
    # Re-center view on SYNGAP1 after adding track
    showGenomicRegion(session, id = "igvShiny_0", region = "SYNGAP1")
    
  }, ignoreInit = TRUE)
  
  # ============================================================
  # CLINVAR SIZE SLIDER HANDLER
  # ============================================================
  
  #' Reload ClinVar Track When Size Threshold Changes
  #'
  # @description
  # Responds to changes in the "Max ClinVar variant size" slider by reloading
  # the ClinVar track with only variants at or below the new KB threshold.
  #
  # Only fires if ClinVar track has already been added (clinvar_track_added()
  # is TRUE) — avoids a no-op reload on app startup when the slider is
  # initialized but no track exists yet.
  #
  # Design: Same reload pattern as the filter-change handler used for internal
  # variant tracks — replace the existing track with a freshly filtered one.
  # IGV handles the swap without flicker.
  
  observeEvent(input$clinvar_size_kb, {
    
    # Only reload if the ClinVar track is actually displayed
    if (!clinvar_track_added()) return()
    
    # Re-build GFF3 with the updated size threshold.
    # as.numeric() is required because sliderTextInput() returns a character.
    clinvar_gff3 <- create_clinvar_gff3_data(clinvar_data,
                                             max_size_kb = as.numeric(input$clinvar_size_kb))
    
    if (nrow(clinvar_gff3) > 0) {
      loadGFF3TrackFromLocalData(
        session,
        id                = "igvShiny_0",
        trackName         = "ClinVar",
        tbl               = clinvar_gff3,
        colorTable        = color_table,
        colorByAttribute  = "type",
        trackHeight       = 50,
        displayMode       = "EXPANDED",
        visibilityWindow  = 1e8
      )
    } else {
      # All ClinVar variants exceed the threshold — show empty track
      loadGFF3TrackFromLocalData(
        session,
        id                = "igvShiny_0",
        trackName         = "ClinVar",
        tbl               = data.frame(),
        colorTable        = color_table,
        colorByAttribute  = "type",
        trackHeight       = 50,
        displayMode       = "EXPANDED",
        visibilityWindow  = 1e8
      )
    }
    
  }, ignoreInit = TRUE)  # Don't fire on startup before a track exists
  
  # ============================================================
  # FILTER CHANGE HANDLING
  # ============================================================
  
  #' Update Tracks When Filters Change
  #'
  # @description
  # Refreshes visible tracks when user changes filter checkboxes.
  #
  # Challenge: Filters change data, but tracks are already displayed
  # Solution: Observe filtered_data() and reload all visible tracks
  #
  # Behavior:
  #   - Check biorepository filter → Tracks update to show only those variants
  #   - Uncheck filter → Tracks update to show all variants again
  #
  # Design Decision: Reload vs Hide
  #   - Could hide filtered-out variants (faster)
  #   - Instead, reload entire track (cleaner)
  #   - Trade: Performance for correctness
  #   - Dataset is small (153 variants) so reload is fast enough
  
  observeEvent(filtered_data(), {
    # Get list of currently displayed tracks
    added_tracks_list <- added_tracks()
    
    # Update each displayed track
    # Only updates tracks that user has added (not all possible tracks)
    for (t in names(added_tracks_list)) {
      # Check if this track is supposed to be visible
      if (added_tracks_list[[t]]) {
        # Get filtered data for this variant type
        track_data <- filtered_data() %>% filter(variant.type == t)
        
        # Check if filtered data contains any variants
        if (nrow(track_data) > 0) {
          # YES - show filtered variants
          gff3_data <- create_gff3_data(track_data, variant_type = t)
          
          # Reload track with new data
          # This replaces the existing track
          # IGV handles the update smoothly (no flicker)
          loadGFF3TrackFromLocalData(
            session,
            id = "igvShiny_0",
            trackName = t,
            tbl = gff3_data,
            colorTable = color_table,
            colorByAttribute = "type",
            trackHeight = 40,
            displayMode = "EXPANDED",
            visibilityWindow = 1e8
          )
        } else {
          # NO - all variants filtered out
          # Show empty track
          loadGFF3TrackFromLocalData(
            session,
            id = "igvShiny_0",
            trackName = t,
            tbl = data.frame(),  # Empty
            colorTable = color_table,
            colorByAttribute = "type",
            trackHeight = 40,
            displayMode = "EXPANDED",
            visibilityWindow = 1e8
          )
        }
      }
    }
  }, ignoreInit = TRUE)  # Don't fire on startup, only when filters actually change
  
  # ============================================================
  # IGV INITIALIZATION
  # ============================================================
  
  #' Render IGV Genome Browser
  #'
  # @description
  # Initializes the IGV browser widget when app starts.
  #
  # renderIgvShiny() is called once when app loads
  # Creates the IGV instance that will display tracks
  #
  # Configuration:
  #   - Genome: hg38 (GRCh38, current human reference)
  #   - Initial view: SYNGAP1 gene locus
  #   - Display mode: SQUISHED (compact track view)
  #   - Tracks: Empty initially (added by user button clicks)
  
  output$igvShiny_0 <- renderIgvShiny({
    # Parse and validate genome specification
    # genomeName: "hg38" specifies human genome build 38
    # initialLocus: "SYNGAP1" centers view on SYNGAP1 gene
    #   - IGV recognizes gene symbols
    #   - Automatically zooms to gene boundaries
    #   - Alternative: Could use coordinates "chr6:33,420,000-33,448,000"
    genomeOptions <- parseAndValidateGenomeSpec(
      genomeName = "hg38",      # Human genome build 38 (GRCh38)
      initialLocus = "SYNGAP1"  # Start viewing at SYNGAP1 gene
    )
    
    # Create IGV instance
    # displayMode: "SQUISHED" makes tracks compact vertically
    #   - Alternatives: "EXPANDED" (tall), "COLLAPSED" (minimal)
    #   - SQUISHED is good default: Readable but space-efficient
    # tracks: list() starts with no tracks
    #   - User adds tracks by clicking buttons
    #   - Clean initial view (not overwhelming)
    igvShiny(
      genomeOptions,             # Genome configuration
      displayMode = "SQUISHED",  # Compact track display
      tracks = list()            # Start with no tracks
    )
  })
}

################################################################################
# SECTION 11: APPLICATION LAUNCH
################################################################################

#' Launch Shiny Application
#'
# @description
# Starts the web server and opens the application in browser.
#
# shinyApp() combines UI and server into runnable app
# When you run this script, it:
#   1. Loads all the data (CSV reading, API queries, etc.)
#   2. Starts a web server (usually on http://localhost:XXXX)
#   3. Opens your default browser to the app
#   4. Waits for user interaction
#   5. Responds to clicks, filters, etc.
#   6. Stops when you close the browser or press Escape in R console
#
# Run this file with: source("app.R") or click "Run App" in RStudio

shinyApp(ui = ui, server = server)

################################################################################
# END OF APPLICATION
################################################################################
#
# Summary of Architecture:
#
# 1. DATA PIPELINE
#    CSV → Parse → API Query → Cache → Merge → Ready for Display
#
# 2. USER INTERACTIONS
#    Click Button → Filter Data → Format GFF3 → Load Track → Update IGV
#    Check Filter → Filter Data → Reload Tracks → Update IGV
#    Move Size Slider → Re-filter ClinVar by span → Reload ClinVar Track → Update IGV
#
# 3. PERFORMANCE OPTIMIZATIONS
#    - Memoization: Ensembl coordinate lookups cached permanently to disk (150x speedup)
#    - Time-limited cache: ClinVar data re-fetched from NCBI every 7 days
#    - Preprocessing: Parse data once at startup, not per interaction
#    - Reactive expressions: Auto-recompute only when inputs change
#    - Batch operations: Prefetch all coordinates before UI loads
#
# 4. KEY DESIGN DECISIONS
#    - Separate data processing from UI rendering
#    - Cache Ensembl coordinate lookups permanently; ClinVar data with 7-day expiry
#    - Use reactive programming for smooth interactions
#    - Generate UI dynamically based on data
#    - Graceful degradation when data is missing
#    - Clear separation of concerns (functions, sections)
#
# 5. FUTURE IMPROVEMENTS
#    - Add genome positions directly to CSV (eliminate API dependency)
#    - Add download button for filtered variants
#    - Add variant search functionality
#    - Add zoom to variant feature
#    - Add track reordering
#    - Add patient detail modal on variant click
#
################################################################################