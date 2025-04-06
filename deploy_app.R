# Deployment script for VCF vs Control Proteomics Analysis app

# Load required packages
if (!require("rsconnect")) {
  install.packages("rsconnect")
}
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

# Deploy the app with data directory included
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
