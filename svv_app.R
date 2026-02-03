################################################################################
#                                                                              #
#                   SYNGAP1 VARIANT VIEWER - Shiny Application                 #
#                                                                              #
#  Purpose: Interactive genome browser for visualizing SYNGAP1 genetic         #
#           variants with patient data and research asset tracking             #
#                                                                              #
#  Architecture:                                                               #
#    1. Data Loading & Preprocessing (lines 1-172)                             #
#       - Load patient variant data from CSV                                   #
#       - Convert cDNA coordinates to genome positions via Ensembl API         #
#       - Cache results for fast subsequent loads                              #
#                                                                              #
#    2. UI Definition (lines 187-229)                                          #
#       - Sidebar with dynamic variant type buttons                            #
#       - Checkboxes for filtering by research assets                          #
#       - Main panel with IGV genome browser                                   #
#                                                                              #
#    3. Server Logic (lines 266-398)                                           #
#       - Reactive filtering based on user selections                          #
#       - Dynamic track loading/updating                                       #
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

################################################################################
# CITATION
################################################################################
# Shannon P, Gladki A, Scigocka K (2024). igvShiny: igvShiny: a wrapper of 
# Integrative Genomics Viewer (IGV - an interactive tool for visualization 
# and exploration integrated genomic data). R package version 1.2.0, 
# https://gladkia.github.io/igvShiny/, https://github.com/gladkia/igvShiny.

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
  loss      = "#4E4E4E"       # Gray
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
#
# 3. PERFORMANCE OPTIMIZATIONS
#    - Memoization: API calls cached to disk (150x speedup)
#    - Preprocessing: Parse data once at startup, not per interaction
#    - Reactive expressions: Auto-recompute only when inputs change
#    - Batch operations: Prefetch all coordinates before UI loads
#
# 4. KEY DESIGN DECISIONS
#    - Separate data processing from UI rendering
#    - Cache external API calls permanently
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