#!/usr/bin/env Rscript

# Script to deploy the Proteomics Analysis App to RConnect/shinyapps.io

# Load the rsconnect package
library(rsconnect)

# Ensure required packages are installed locally
if (!requireNamespace("markdown", quietly = TRUE)) {
  install.packages("markdown")
}

# Set the application name
app_name <- "VCP_Proteomics_Analysis"

# Deploy the application
deployApp(
  appDir = ".",                   # Current directory
  appName = app_name,             # Application name
  appTitle = "Proteomics Analysis App",
  appFiles = c(                   # Files to include in the deployment
    "app.R",
    "about.md",
    "init.R",
    "data/Control_VCP_metadata.csv",
    "data/Controlonly_samples.csv", 
    "data/VCP_samples.csv"
  ),
  appPrimaryDoc = "app.R",        # Main application file
  launch.browser = TRUE,          # Open browser after deployment
  forceUpdate = TRUE              # Force update if the app already exists
)

# Print success message
cat(paste0("Application '", app_name, "' has been deployed successfully.\n"))
