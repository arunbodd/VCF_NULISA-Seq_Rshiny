# init.R
#
# This script is executed during the application startup.
# It ensures all required packages are installed on the server.

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
  "wesanderson",
  "markdown"
)

# Function to check and install missing packages
install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cran.rstudio.com/")
  }
}

# Install missing packages
for (pkg in required_packages) {
  install_if_missing(pkg)
}
