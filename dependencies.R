#!/usr/bin/env Rscript

# This file lists all package dependencies for the app
# Run this script to install all required packages

# List of required packages
required_packages <- c(
  "shiny",
  "shinydashboard",
  "DT",
  "plotly",
  "ggplot2",
  "pheatmap",
  "reshape2",
  "dplyr",
  "tidyr",
  "ggpubr",
  "ggthemes",
  "ggrepel",
  "stringr",
  "data.table",
  "corrplot",
  "circlize",
  "RColorBrewer",
  "pROC",
  "caret",
  "heatmaply",
  "shinyBS",
  "shinyWidgets",
  "tibble",
  "ggiraph",
  "wesanderson"
)

# Function to check and install missing packages
install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
    cat(paste0("Installed package: ", pkg, "\n"))
  } else {
    cat(paste0("Package already installed: ", pkg, "\n"))
  }
}

# Install missing packages
for (pkg in required_packages) {
  install_if_missing(pkg)
}

cat("All required packages are installed.\n")
