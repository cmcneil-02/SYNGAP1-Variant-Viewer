# SYNGAP1 Variant Viewer - Installation Script
# This script installs all required packages for the R Shiny app

cat("========================================\n")
cat("SYNGAP1 Variant Viewer - Setup Script\n")
cat("========================================\n\n")

# Function to check if package is installed
check_package <- function(pkg) {
  return(requireNamespace(pkg, quietly = TRUE))
}

# 1. Install BiocManager if needed
cat("Step 1/4: Checking BiocManager...\n")
if (!check_package("BiocManager")) {
  cat("  Installing BiocManager...\n")
  install.packages("BiocManager", repos = "https://cran.r-project.org")
} else {
  cat("  BiocManager already installed ✓\n")
}

# 2. Install CRAN packages
cat("\nStep 2/4: Installing CRAN packages...\n")
cran_packages <- c("shiny", "dplyr", "stringr", "httr", "jsonlite", "memoise", "cachem", "shinyWidgets")

for (pkg in cran_packages) {
  if (!check_package(pkg)) {
    cat(paste0("  Installing ", pkg, "...\n"))
    install.packages(pkg, repos = "https://cran.r-project.org")
  } else {
    cat(paste0("  ", pkg, " already installed ✓\n"))
  }
}

# 3. Install igvShiny from Bioconductor
cat("\nStep 3/4: Installing igvShiny from Bioconductor...\n")
if (!check_package("igvShiny")) {
  cat("  Installing igvShiny...\n")
  BiocManager::install("igvShiny", update = FALSE, ask = FALSE)
} else {
  cat("  igvShiny already installed ✓\n")
}

# 4. Create cache directories
cat("\nStep 4/4: Creating cache directories...\n")
if (!dir.exists("cache/genome_positions")) {
  dir.create("cache/genome_positions", recursive = TRUE)
  cat("  Created cache/genome_positions/ ✓\n")
} else {
  cat("  cache/genome_positions/ already exists ✓\n")
}
if (!dir.exists("cache/clinvar")) {
  dir.create("cache/clinvar", recursive = TRUE)
  cat("  Created cache/clinvar/ ✓\n")
} else {
  cat("  cache/clinvar/ already exists ✓\n")
}

# Verify installation
cat("\n========================================\n")
cat("Installation Complete! Verifying...\n")
cat("========================================\n\n")

required_packages <- c("shiny", "igvShiny", "dplyr", "stringr", "httr", "jsonlite", "memoise", "shinyWidgets")
all_installed <- TRUE

for (pkg in required_packages) {
  if (check_package(pkg)) {
    version <- tryCatch({
      as.character(packageVersion(pkg))
    }, error = function(e) {
      "unknown"
    })
    cat(sprintf("✓ %s (version %s)\n", pkg, version))
  } else {
    cat(sprintf("✗ %s - FAILED TO INSTALL\n", pkg))
    all_installed <- FALSE
  }
}

cat("\n========================================\n")
if (all_installed) {
  cat("SUCCESS! All packages installed.\n")
  cat("You can now run the app with:\n")
  cat("  shiny::runApp('app.R')\n")
} else {
  cat("WARNING! Some packages failed to install.\n")
  cat("Please review the error messages above.\n")
}
cat("========================================\n")