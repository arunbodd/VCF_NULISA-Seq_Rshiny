# Deployment script for VCF vs Control Proteomics Analysis app

# Load required packages
required_packages <- c("rsconnect", "markdown", "shiny", "shinydashboard", "dplyr", "ggplot2", 
                      "DT", "pROC", "reshape2", "plotly", "ggrepel", "stringr", "tidyr", 
                      "pheatmap", "RColorBrewer", "ComplexHeatmap", "circlize", "STRINGdb", 
                      "igraph", "ggpubr", "ggthemes", "data.table", "corrplot")

# Check and install missing packages
missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]
if (length(missing_packages) > 0) {
  cat("Installing missing packages:", paste(missing_packages, collapse = ", "), "\n")
  install.packages(missing_packages)
}

# Load rsconnect
library(rsconnect)

# IMPORTANT: Never store credentials in your code!
# Instead, use one of these secure approaches:

# OPTION 1: Use the rsconnect::setAccountInfo() function interactively before running this script
# Run this in the console, NOT in this script:
# rsconnect::setAccountInfo(name='YOUR_ACCOUNT_NAME', token='YOUR_TOKEN', secret='YOUR_SECRET')

# OPTION 2: Store credentials in environment variables
# You can set these in your .Renviron file or using Sys.setenv()
if (Sys.getenv("SHINYAPPS_NAME") != "" && 
    Sys.getenv("SHINYAPPS_TOKEN") != "" && 
    Sys.getenv("SHINYAPPS_SECRET") != "") {
  
  rsconnect::setAccountInfo(
    name = Sys.getenv("SHINYAPPS_NAME"),
    token = Sys.getenv("SHINYAPPS_TOKEN"),
    secret = Sys.getenv("SHINYAPPS_SECRET")
  )
} else {
  message("Please set your shinyapps.io credentials before deploying.")
  message("You can do this by running rsconnect::setAccountInfo() or by setting environment variables.")
  message("See the README for more information.")
  # Uncomment to stop execution if credentials are missing
  # stop("Missing shinyapps.io credentials")
}

# Check if required files exist
required_files <- c("app.R", "config.R", "about.md", "README.md")
missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  stop("Missing required files: ", paste(missing_files, collapse = ", "))
}

# Check if data directory exists and contains required files
if (!dir.exists("data")) {
  stop("Data directory not found. Please create a 'data' directory with the required data files.")
}

data_files <- c(
  "data/Controlonly_samples.csv",
  "data/VCP_samples.csv",
  "data/Control_VCP_metadata.csv"
)

missing_data_files <- data_files[!file.exists(data_files)]

if (length(missing_data_files) > 0) {
  stop("Missing data files: ", paste(missing_data_files, collapse = ", "))
}

# Deploy the app with data directory included
tryCatch({
  deployApp(
    appName = "VCF-Proteomics-Analysis",
    appTitle = "VCF vs Control Proteomics Analysis",
    appFiles = c(
      "app.R",
      "config.R",
      "about.md",
      "README.md",
      "data"  # Include the entire data directory
    ),
    launch.browser = TRUE
  )
  
  # Print success message
  cat("\nDeployment completed!\n")
  cat("Your app should be available at: https://arunbodd.shinyapps.io/VCF-Proteomics-Analysis/\n")
}, error = function(e) {
  cat("\nError during deployment:\n")
  print(e)
  cat("\nTroubleshooting tips:\n")
  cat("1. Check that all required files exist\n")
  cat("2. Ensure your data files are in the 'data' directory\n")
  cat("3. Verify your shinyapps.io credentials are correct\n")
  cat("4. Check your internet connection\n")
  cat("5. Try running the app locally first: shiny::runApp()\n")
})
