# SYNGAP1 Variant Viewer

An interactive R Shiny application for visualizing SYNGAP1 genetic variants using the Integrative Genomics Viewer (IGV). This tool enables researchers to explore variant data on a genome browser with dynamic filtering by available research resources.

Developed for **[CURE SYNGAP1](https://curesyngap1.org/)** - a global group of families committed to accelerating the science to cure SYNGAP1 & to supporting each other.

![License](https://img.shields.io/badge/license-MIT-blue.svg)
![R Version](https://img.shields.io/badge/R-%E2%89%A54.0.0-blue)
![Status](https://img.shields.io/badge/status-active-success.svg)

---

## Table of Contents

- [Overview](#overview)
- [Features](#features)
- [Installation](#installation)
- [Usage](#usage)
- [Architecture](#architecture)
- [File Structure](#file-structure)
- [Performance & Caching](#performance--caching)
- [Maintenance](#maintenance)
- [Troubleshooting](#troubleshooting)
- [Citation](#citation)
- [Acknowledgments](#acknowledgments)
- [License](#license)

---

## Overview

The SYNGAP1 Variant Viewer is a web-based application designed to help researchers visualize and analyze genetic variants in the SYNGAP1 gene. SYNGAP1 variants are associated with neurodevelopmental disorders, and this tool facilitates:

- **Variant visualization** on chromosome 6 using an embedded genome browser
- **ClinVar integration**: live NCBI ClinVar data fetched automatically and displayed as a toggleable reference track
- **Filtering by research assets** (biorepository samples, iPSC cell lines, mouse models)
- **Interactive exploration** with color-coded variant types
- **Coordinate conversion** from cDNA (transcript-relative) to genomic positions (genome-relative)

The application currently displays **153 variants** tracked through Citizen Health and partner research organizations.

### Key Technologies

- **R Shiny**: Web application framework
- **igvShiny**: R wrapper for Integrative Genomics Viewer
- **Ensembl REST API**: Genomic coordinate conversion
- **NCBI E-utilities**: Live ClinVar variant data retrieval
- **Reactive programming**: Real-time UI updates

---

## Features

### Core Functionality

- ✅ **Interactive Genome Browser**: Embedded IGV displaying variants on chromosome 6
- ✅ **ClinVar Reference Track**: Toggleable track showing all SYNGAP1 ClinVar submissions, fetched live from NCBI and cached for 7 days; variant size filterable via a log-spaced slider (1–50,000 KB, default 100 KB)
- ✅ **Multiple Variant Types**: Missense, nonsense, frameshift, indel, intronic, and structural variants
- ✅ **Dynamic Track Loading**: Add/remove variant tracks by type
- ✅ **Color-Coded Display**: Each variant type has a distinct color for easy identification
- ✅ **Smart Filtering**: Filter variants by available research resources
- ✅ **Offline Capability**: Works without internet after initial setup (via caching)
- ✅ **Variant Metadata**: Display variant-associated variant counts and research asset availability

### Filtering Options

1. **Biorepository samples**: Show only variants with stored biological samples
2. **Cell line availability**: Show only variants with iPSC cell lines
3. **Mouse line availability**: Show only variants with mouse models

Filters are combinatorial (AND logic) - checking multiple filters shows variants that meet ALL criteria.

### Variant Display

Each variant shows:
- Protein change (e.g., "Arg135ter") or cDNA notation (e.g., "c.333del")
- Number of patients with this variant
- Genomic position on chromosome 6
- Variant type classification

---

## Installation

### Prerequisites

- **R version**: ≥ 4.0.0
- **Operating System**: Windows, macOS, or Linux
- **Internet connection**: Required for initial setup only

### Step 1: Clone the Repository

```bash
git clone https://github.com/cmcneil-02/syngap1-variant-viewer.git
cd syngap1-variant-viewer
```

### Step 2: Install R Package Dependencies

The app requires several R packages. Install them using the provided script:

```r
# Run the automated installation script
source("svv_packages.R")
```

This will install:
- `shiny` - Web application framework
- `igvShiny` - IGV genome browser wrapper (from Bioconductor)
- `dplyr` - Data manipulation
- `stringr` - String processing
- `httr` - HTTP requests for API
- `jsonlite` - JSON parsing
- `memoise` - Function caching
- `cachem` - Cache management
- `shinyWidgets` - Extended UI widgets (log-spaced ClinVar size slider)

### Step 3: Create Cache Directory

The cache directories will be created automatically on first run, but you can create them manually:

```r
dir.create("cache/genome_positions", recursive = TRUE)
dir.create("cache/clinvar", recursive = TRUE)
```

### Step 4: Run the Application

```r
# In R console
shiny::runApp("svv_app.R")

# Or in RStudio
# Open svv_app.R and click "Run App" button
```

The app will open in your default web browser.

---

## Usage

### Starting the Application

1. Ensure all dependencies are installed (run `svv_packages.R` if needed)
2. Run: `shiny::runApp("svv_app.R")`
3. The app opens in your web browser at `http://localhost:XXXX`

### Basic Workflow

**Step 1: View the genome browser**
- The IGV browser loads centered on the SYNGAP1 gene
- Reference genome: hg38 (GRCh38)
- Chromosome: chr6

**Step 2: Add variant tracks**
- Click buttons in the left sidebar to add variant type tracks
- Example: Click "missense" to display all missense variants
- Each track appears as a horizontal layer in IGV
- Click **"ClinVar"** (below a separator) to overlay all SYNGAP1 ClinVar submissions as a teal reference track; this track is independent of the research-asset filters
- Use the **"Max ClinVar variant size"** slider beneath the ClinVar button to control which variants appear by their genomic span. The slider is log-spaced (1–50,000 KB) so you get fine-grained control at the low end (1 KB steps for point variants) and coarser steps at the high end (5,000 KB steps for large CNVs). Default is 100 KB. Moving the slider immediately reloads the track if it is already displayed.

**Step 3: Apply filters (optional)**
- Check "Has biorepository samples" to show only variants with stored samples
- Check "Cell line available" to show only variants with iPSC lines
- Check "Mouse line available" to show only variants with mouse models
- Tracks update automatically when filters change

**Step 4: Explore variants**
- Zoom in/out using IGV controls
- Click on variants to see details
- Tracks show variant names (protein changes or cDNA notation)

### First Run vs. Subsequent Runs

**First Run** (with internet):
- Takes ~70-120 seconds total
- Queries NCBI E-utilities to fetch all SYNGAP1 ClinVar entries (~30-60 seconds); result cached to `cache/clinvar/` for 7 days
- Queries Ensembl API to convert cDNA coordinates to genome positions (~40-60 seconds); result cached permanently to `cache/genome_positions/`

**Subsequent Runs within 7 days** (can be offline):
- Takes <2 seconds
- Reads both the ClinVar data and genome positions from disk cache
- No internet required

**After 7 days:**
- ClinVar cache is considered stale; the app re-fetches from NCBI on startup (~30-60 seconds)
- Genome position cache is permanent and never needs refreshing

---

## Architecture

### Application Structure

The app follows a three-layer architecture:

```
┌─────────────────────────────────────────┐
│ DATA LAYER                              │
│ • Fetch ClinVar data (NCBI, 7-day cache)│
│ • Load updatedCitizen191.csv            │
│ • Parse variant nomenclature            │
│ • Convert coordinates (API + cache)     │
│ • Format protein changes                │
└─────────────────┬───────────────────────┘
                  ↓
┌─────────────────────────────────────────┐
│ PRESENTATION LAYER (UI)                 │
│ • Sidebar controls                      │
│ • IGV browser widget                    │
│ • Dynamic variant-type buttons          │
│ • Static ClinVar track button           │
└─────────────────┬───────────────────────┘
                  ↓
┌─────────────────────────────────────────┐
│ LOGIC LAYER (Server)                    │
│ • Reactive filtering                    │
│ • Track management                      │
│ • Event handling (incl. ClinVar button) │
└─────────────────────────────────────────┘
```

### Key Components

**1. ClinVar Data Loading**
- Source: NCBI E-utilities API (esearch + esummary)
- Fetch: All SYNGAP1 ClinVar variation IDs, then full metadata in batches of 200
- Cache: `cache/clinvar/clinvar_data.rds` — expires after 7 days, then re-fetched automatically
- Fallback: If NCBI is unreachable, app uses stale cache (or loads empty track as last resort)
- Display: GFF3 track colored teal (#00BCD4); all ClinVar metadata visible on click in IGV

**2. Coordinate Conversion**
- Input: cDNA notation (e.g., "c.333del")
- Extraction: Numeric coordinate (333)
- API Query: Ensembl REST API
- Output: Genome position (chr6:33,425,796)
- Caching: Permanent disk storage via `memoise`

**2. Coordinate Conversion**
- Input: cDNA notation (e.g., "c.333del")
- Extraction: Numeric coordinate (333)
- API Query: Ensembl REST API
- Output: Genome position (chr6:33,425,796)
- Caching: Permanent disk storage via `memoise`

**3. Protein Nomenclature Formatting**
- Converts single-letter amino acid codes to three-letter codes
- Example: "p.R135X" → "Arg135ter"
- Improves readability for biologists

**4. Reactive Filtering**
- Filter changes trigger automatic data recalculation
- Tracks reload with filtered data
- No manual refresh needed
- Note: ClinVar track is unaffected by research-asset filters (it is a reference dataset, not internal data)

**5. GFF3 Track Generation**
- Converts both internal variant data and ClinVar data to GFF3 format
- IGV-compatible genomic feature format
- Internal tracks include: ID, Name, Description, variant count
- ClinVar track includes: all ClinVar metadata fields (germline classification, review status, conditions, etc.)

### Data Flow

```
NCBI E-utilities (esearch + esummary)
  ↓
Cache to cache/clinvar/ (7-day expiry)
  ↓
[On button click: format as GFF3 → ClinVar track in IGV]

updatedCitizen191.csv
  ↓
Extract cDNA coordinates
  ↓
Query Ensembl API (permanently cached)
  ↓
Merge genome positions
  ↓
User selects filters
  ↓
Filter variant data
  ↓
Generate GFF3 tracks
  ↓
Display in IGV browser
```

For detailed architecture documentation, see **`svv_architecture.Rmd`**.

---

## File Structure

```
syngap1-variant-viewer/
├── svv_app.R                       # Fully annotated application code
├── svv_architecture.Rmd            # Detailed architecture documentation
├── svv_packages.R                  # Automated dependency installation
├── README.md                       # This file
├── updatedCitizen191.csv          # Variant data (153 variants, 33 columns)
├── cache/
│   ├── genome_positions/          # Ensembl coordinate cache (permanent, auto-generated)
│   │   ├── [hash1].rds
│   │   ├── [hash2].rds
│   │   └── ... (~143 files)
│   └── clinvar/                   # ClinVar data cache (7-day expiry, auto-generated)
│       ├── clinvar_data.rds       # Parsed ClinVar data frame
│       └── last_updated.txt       # ISO timestamp of last successful fetch
└── LICENSE                        # License file
```

### File Descriptions

**svv_app.R**
- Main application with extensive inline comments (500+ lines of documentation)
- Production-ready code with explanations of every function and design decision
- Use this file to run the application
- Excellent resource for understanding the codebase

**svv_architecture.Rmd**
- R Markdown document explaining the overall architecture
- Design philosophy and development timeline
- Key design decisions and tradeoffs
- Common programming patterns used
- Read this to understand the big picture

**svv_packages.R**
- Automated script to install all required R packages
- Handles both CRAN and Bioconductor dependencies
- Run once during initial setup

**updatedCitizen191.csv**
- Variant data from Citizen Health
- 153 variants across 33 data columns
- Includes: variant notation, type, patient counts, research assets
- **Required** for application to function

**cache/genome_positions/**
- Auto-generated directory storing API response cache
- Contains ~143 .rds files (one per unique cDNA coordinate)
- Enables offline operation after first run
- Can be regenerated by deleting directory and rerunning app

---

## Performance & Caching

The app maintains two separate caches to avoid redundant network calls.

### Ensembl Coordinate Cache (Permanent)

Genomic positions for cDNA coordinates are fetched once from the Ensembl REST API and stored permanently, since hg38 positions never change.

```r
# Permanent memoized cache
get_genome_positions <- memoise(
  get_genome_positions_function, 
  cache = cachem::cache_disk("cache/genome_positions")
)
```

**How it works:**
1. First call to `get_genome_positions(333)` queries the Ensembl API (300ms)
2. Result saved to `cache/genome_positions/[hash].rds`
3. Subsequent calls read from disk (2ms) — **150× faster**
4. Cache persists indefinitely across R sessions

### ClinVar Cache (7-Day Expiry)

ClinVar submissions change regularly, so the app re-fetches from NCBI after 7 days.

**How it works:**
1. On startup, the app checks `cache/clinvar/last_updated.txt`
2. If cache is < 7 days old, `clinvar_data.rds` is loaded instantly (<1 second)
3. If cache is absent or expired, the app fetches from NCBI (esearch + esummary, ~30-60 seconds), then saves to `cache/clinvar/`
4. If the fetch fails but a stale cache exists, the stale cache is used with a warning

**To force a re-fetch before 7 days:**
```r
unlink("cache/clinvar/last_updated.txt")
# Restart the app — it will treat the cache as expired
```

### Performance Metrics

| Operation | First Run | Cached Run |
|-----------|-----------|------------|
| Single Ensembl coordinate | 300ms | 2ms |
| 143 Ensembl coordinates | ~43 seconds | ~0.3 seconds |
| Ensembl cache speedup | — | **150×** |
| ClinVar fetch (~1,900 variants) | ~30-60 seconds | <1 second |

### Cache Management

**When to clear the Ensembl cache:**
- Switched genome builds (hg38 → hg19)
- Changed transcript IDs
- Cache corruption

```r
# Delete Ensembl coordinate cache
unlink("cache/genome_positions/*.rds")
```

**When to clear the ClinVar cache:**
- To force an immediate re-fetch before the 7-day window expires

```r
# Force ClinVar re-fetch on next startup
unlink("cache/clinvar/last_updated.txt")
```

**Cache sizes:**
- Ensembl: ~2KB per coordinate × 143 coordinates ≈ 300KB total
- ClinVar: ~1-3MB for the full SYNGAP1 entry set (single .rds file)

---

## Maintenance

### Updating Variant Data

To update the variant dataset:

1. Export new data from your master database (Excel)
2. Ensure column names match the current structure
3. Replace `updatedCitizen191.csv`
4. Restart the application

The app will automatically:
- Parse new variants
- Query Ensembl API for any new coordinates (or use cache for existing ones)
- Update the visualization

### Updating ClinVar Data

ClinVar data is refreshed automatically every 7 days on startup — no manual action is needed. To force an immediate refresh:

```r
unlink("cache/clinvar/last_updated.txt")
# Then restart the app
```

### Adding New Variant Types

If your data includes new variant classifications:

1. Add to variant type normalization (~line 878 in `svv_app.R`)
2. Add color to color table (~line 1033)
3. No other changes needed — app generates buttons dynamically

Example:
```r
# In variant type normalization
variant.type = case_when(
  variant.type %in% c('splice site') ~ 'splice',
  # ... existing mappings ...
)

# In color table
color_table <- list(
  splice = "#FF5733",  # New color for splice variants
  # ... existing colors ...
)
```

### Changing Genome Build

Currently configured for **hg38**. To switch to hg19:

1. Update `genomeName = "hg19"` in the `renderIgvShiny` block (~line 1843 in `svv_app.R`)
2. Clear Ensembl cache: `unlink("cache/genome_positions/*.rds")`
3. Restart app (will rebuild cache for hg19)

Note: The ClinVar GRCh38 coordinate filter in `create_clinvar_gff3_data` would also need updating if switching builds.

### Changing Transcript

Currently uses **ENST00000418600** (SYNGAP1 canonical transcript).

To change:
1. Update `transcript_id = "ENST00000XXXXXX"` in the coordinate fetch loop (~line 989 in `svv_app.R`)
2. Clear Ensembl cache: `unlink("cache/genome_positions/*.rds")`
3. Restart app

---

## Troubleshooting

### Common Issues

**Problem: `Error in library(igvShiny) : there is no package called 'igvShiny'`**

**Solution:**
```r
source("svv_packages.R")
# Or manually:
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("igvShiny")
```

---

**Problem: ClinVar track is empty or missing**

**Solution:** This can happen if the NCBI fetch failed on startup and no cache existed. Check the R console for ClinVar-related warnings. Try:
```r
unlink("cache/clinvar/last_updated.txt")
# Then restart the app
```
If the problem persists, verify your internet connection and that `https://eutils.ncbi.nlm.nih.gov` is reachable.

---

**Problem: ClinVar track shows outdated data**

**Solution:** The cache refreshes automatically every 7 days. To force an immediate re-fetch:
```r
unlink("cache/clinvar/last_updated.txt")
# Restart the app
```

---

**Problem: App loads slowly on first run**

**Solution:** This is expected! First run takes ~70-120 seconds to:
- Fetch all SYNGAP1 ClinVar entries from NCBI (~30-60 seconds)
- Query Ensembl API for cDNA→genome coordinate conversion (~40-60 seconds)

Both caches are then stored to disk; subsequent runs start in <2 seconds.

---

**Problem: No variants displayed on tracks**

**Solution:** Check if filters are too restrictive. Uncheck all filter boxes to see all variants, then reapply filters incrementally.

---

**Problem: Variants missing genome positions**

**Solution:** Some variants (structural, intronic) may not have valid cDNA coordinates. These show as NA in positions and are skipped in visualization. This is expected behavior.

---

**Problem: `Error: API rate limit exceeded`**

**Solution:** You've made >15 requests/second to Ensembl. Wait 60 seconds and restart. After first successful run, cache prevents this issue.

---

**Problem: Cache directory not found**

**Solution:**
```r
dir.create("cache/genome_positions", recursive = TRUE)
```

---

### Logging & Debugging

Enable console logging to track API calls:

During first run, watch for messages like:
- `"ClinVar: querying NCBI esearch for SYNGAP1 variation IDs..."` — ClinVar fetch started
- `"ClinVar: found X variation IDs."` — esearch succeeded
- `"ClinVar: fetching summaries in N batches..."` — esummary batches in progress
- `"ClinVar: successfully parsed X variant records."` — ClinVar data ready
- `"ClinVar: loading from cache (< 7 days old)."` — cache hit, no fetch needed
- `"No mapping data found for cDNA coordinate X"` — Ensembl coordinate couldn't be mapped (expected for some variants)
- `"Error fetching genome position for coordinate Y"` — Ensembl API error (will skip)

ClinVar warnings to watch for:
- `"ClinVar: NCBI fetch failed. Falling back to stale cache."` — network issue; app continues with old data
- `"ClinVar: NCBI fetch failed and no cache found."` — ClinVar track will be empty this session

These messages are informational and don't indicate problems unless all variants fail.

---

## Citation

### Software Citations

If you use this tool in your research, please cite:

**This application:**
```
SYNGAP1 Variant Viewer
Developed for Cure SYNGAP1
https://github.com/cmcneil-02/syngap1-variant-viewer
```

**igvShiny (required dependency):**
```
Shannon P, Gladki A, Scigocka K (2024). igvShiny: igvShiny: a wrapper of 
Integrative Genomics Viewer (IGV - an interactive tool for visualization 
and exploration integrated genomic data). R package version 1.2.0, 
https://gladkia.github.io/igvShiny/, https://github.com/gladkia/igvShiny
```

**Ensembl REST API:**
```
Yates AD, et al. (2020). Ensembl 2020. Nucleic Acids Research, 
48(D1), D682-D688. doi: 10.1093/nar/gkz966
```

**NCBI E-utilities (ClinVar data):**
```
Sayers EW, et al. (2022). Database resources of the National Center for 
Biotechnology Information. Nucleic Acids Research, 50(D1), D20-D26. 
doi: 10.1093/nar/gkab1112
```

### BibTeX Entries

```bibtex
@software{syngap1_viewer,
  title = {SYNGAP1 Variant Viewer},
  author = {CURE SYNGAP1 Team},
  organization = {Cure SYNGAP1},
  year = {2025},
  url = {https://github.com/cmcneil-02/syngap1-variant-viewer}
}

@Manual{igvShiny,
  title = {igvShiny: igvShiny wrapper of Integrative Genomics Viewer},
  author = {Paul Shannon and Arkadiusz Gladki and Karolina Scigocka},
  year = {2024},
  note = {R package version 1.2.0},
  url = {https://github.com/gladkia/igvShiny}
}
```

---

## Acknowledgments

### Organizations

**Cure SYNGAP1**
- Website: https://curesyngap1.org/
- Mission: Finding a cure for SYNGAP1-related disorders
- This application was developed to support SYNGAP1 research and patient advocacy efforts

**Citizen Health**
- Patient variant data source and registry partner

### Software & Tools

- **igvShiny developers**: Paul Shannon, Arkadiusz Gladki, Karolina Scigocka
- **Ensembl**: Genome annotation and REST API services
- **NCBI**: ClinVar database and E-utilities API
- **R Shiny team**: Application framework
- **SYNGAP1 research community**: Collaborative data sharing and research efforts

### Contributors

- Collin McNeil

---

## License

- MIT License - See LICENSE file for details

This software is provided for research and academic use.

---

## Contact & Support

**Issues & Questions:**
- Open an issue on GitHub: https://github.com/cmcneil-02/SYNGAP1-Variant-Viewer/issues

**CURE SYNGAP1:**
- Website: https://curesyngap1.org/

**Collaboration Inquiries:**
- [Contact method to be determined]

---

## Data Source

Variant data is sourced from **Citizen Health** and partner research organizations. The dataset includes:
- 153 SYNGAP1 variants
- Patient counts per variant
- Research asset availability (biorepository samples, cell lines, mouse models)
- Clinical and functional annotations

Data is publicly available through Citizen Health and Cure SYNGAP1.

---

## Version History

### v2.1.0 (2026)
- Added live ClinVar integration via NCBI E-utilities (esearch + esummary)
- ClinVar data fetched automatically on startup and cached locally for 7 days
- ClinVar variants displayed as a toggleable teal reference track in IGV
- All ClinVar metadata (germline classification, review status, conditions, etc.) accessible on variant click
- Log-spaced size filter slider (1–50,000 KB, default 100 KB) controls which ClinVar variants are shown by genomic span; reloads track live on change

### v1.0.0 (2025)
- Initial release
- Support for 153 variants
- Three filter types (biorepository, iPSC, mouse models)
- IGV browser integration
- Ensembl API with permanent caching
- Color-coded variant types
- Reactive filtering

---

**Community Contributions Welcome!**

---

**Last Updated:** 02-17-2026  
**Repository:** cmcneil-02  
**Documentation Version:** 2.1  
**Maintained by:** CURE SYNGAP1 Team