# SYNGAP1 Variant Viewer

An interactive R Shiny application for visualizing SYNGAP1 genetic variants using the Integrative Genomics Viewer (IGV). This tool enables researchers to explore patient variant data on a genome browser with dynamic filtering by available research resources.

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
- **Filtering by research assets** (biorepository samples, iPSC cell lines, mouse models)
- **Interactive exploration** with color-coded variant types
- **Coordinate conversion** from cDNA (transcript-relative) to genomic positions (genome-relative)

The application currently displays **153 patient variants** tracked through Citizen Health and partner research organizations.

### Key Technologies

- **R Shiny**: Web application framework
- **igvShiny**: R wrapper for Integrative Genomics Viewer
- **Ensembl REST API**: Genomic coordinate conversion
- **Reactive programming**: Real-time UI updates

---

## Features

### Core Functionality

- ✅ **Interactive Genome Browser**: Embedded IGV displaying variants on chromosome 6
- ✅ **Multiple Variant Types**: Missense, nonsense, frameshift, indel, intronic, and structural variants
- ✅ **Dynamic Track Loading**: Add/remove variant tracks by type
- ✅ **Color-Coded Display**: Each variant type has a distinct color for easy identification
- ✅ **Smart Filtering**: Filter variants by available research resources
- ✅ **Offline Capability**: Works without internet after initial setup (via caching)
- ✅ **Patient Metadata**: Display variant-associated patient counts and research asset availability

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

### Step 3: Create Cache Directory

The cache directory will be created automatically on first run, but you can create it manually:

```r
dir.create("cache/genome_positions", recursive = TRUE)
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
- Takes ~40-60 seconds
- Queries Ensembl API to convert cDNA coordinates to genome positions
- Caches all results to disk in `cache/genome_positions/`
- Creates ~143 cache files

**Subsequent Runs** (can be offline):
- Takes <1 second
- Reads positions from cache
- No internet required

---

## Architecture

### Application Structure

The app follows a three-layer architecture:

```
┌─────────────────────────────────────────┐
│ DATA LAYER                              │
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
│ • Dynamic button generation             │
└─────────────────┬───────────────────────┘
                  ↓
┌─────────────────────────────────────────┐
│ LOGIC LAYER (Server)                    │
│ • Reactive filtering                    │
│ • Track management                      │
│ • Event handling                        │
└─────────────────────────────────────────┘
```

### Key Components

**1. Coordinate Conversion**
- Input: cDNA notation (e.g., "c.333del")
- Extraction: Numeric coordinate (333)
- API Query: Ensembl REST API
- Output: Genome position (chr6:33,425,796)
- Caching: Permanent disk storage via `memoise`

**2. Protein Nomenclature Formatting**
- Converts single-letter amino acid codes to three-letter codes
- Example: "p.R135X" → "Arg135ter"
- Improves readability for biologists

**3. Reactive Filtering**
- Filter changes trigger automatic data recalculation
- Tracks reload with filtered data
- No manual refresh needed

**4. GFF3 Track Generation**
- Converts variant data to GFF3 format
- IGV-compatible genomic feature format
- Includes attributes: ID, Name, Description, Patient count

### Data Flow

```
updatedCitizen191.csv
  ↓
Extract cDNA coordinates
  ↓
Query Ensembl API (cached)
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
│   └── genome_positions/          # API response cache (auto-generated)
│       ├── [hash1].rds
│       ├── [hash2].rds
│       └── ... (~143 files)
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
- Patient variant data from Citizen Health
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

### Memoization Strategy

The app uses **memoization** to cache expensive Ensembl API calls:

```r
# Cached function
get_genome_positions <- memoise(
  get_genome_positions_function, 
  cache = cachem::cache_disk("cache/genome_positions")
)
```

**How it works:**
1. First call to `get_genome_positions(333)` queries the API (300ms)
2. Result saved to `cache/genome_positions/[hash].rds`
3. Subsequent calls read from disk (2ms) - **150× faster**
4. Cache persists across R sessions

### Performance Metrics

| Operation | First Run | Cached Run |
|-----------|-----------|------------|
| Single coordinate | 300ms | 2ms |
| 143 coordinates | ~43 seconds | ~0.3 seconds |
| Speedup | — | **150×** |

### Cache Management

**When to clear cache:**
- Switched genome builds (hg38 → hg19)
- Changed transcript IDs
- Cache corruption

**How to clear:**
```r
# Delete all cache files
unlink("cache/genome_positions/*.rds")
```

**Cache size:**
- ~2KB per variant coordinate
- 143 variants = ~300KB total (negligible)

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
- Query API for any new coordinates (or use cache for existing ones)
- Update the visualization

### Adding New Variant Types

If your data includes new variant classifications:

1. Add to variant type normalization (lines 131-147 in `svv_app.R`)
2. Add color to color table (lines 175-185)
3. No other changes needed - app generates buttons dynamically

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

1. Update line 395 in `svv_app.R`: `genomeName = "hg19"`
2. Clear cache: `unlink("cache/genome_positions/*.rds")`
3. Restart app (will rebuild cache for hg19)

### Changing Transcript

Currently uses **ENST00000418600** (SYNGAP1 canonical transcript).

To change:
1. Update line 96 in `svv_app.R`: `transcript_id = "ENST00000XXXXXX"`
2. Clear cache
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

**Problem: App loads slowly on first run**

**Solution:** This is expected! First run takes ~40-60 seconds to query Ensembl API and build cache. Subsequent runs are <1 second.

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
- `"No mapping data found for cDNA coordinate X"` - Coordinate couldn't be mapped (expected for some variants)
- `"Error fetching genome position for coordinate Y"` - API error (will retry or skip)

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

### v1.0.0 (2025)
- Initial release
- Support for 153 variants
- Three filter types (biorepository, iPSC, mouse models)
- IGV browser integration
- Ensembl API with caching
- Color-coded variant types
- Reactive filtering

---

**Community Contributions Welcome!**

---

**Last Updated:** 02-03-2026  
**Repository:** cmcneil-02  
**Documentation Version:** 1.0  
**Maintained by:** CURE SYNGAP1 Team