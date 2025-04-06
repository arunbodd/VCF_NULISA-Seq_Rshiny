# Script to check data file structure
# Load necessary libraries
library(dplyr)

# File paths
data_dir <- "data"
control_file <- file.path(data_dir, "Controlonly_samples.csv")
vcp_file <- file.path(data_dir, "VCP_samples.csv")
metadata_file <- file.path(data_dir, "Control_VCP_metadata.csv")

# Check if files exist
cat("Checking if files exist:\n")
cat("Control file exists:", file.exists(control_file), "\n")
cat("VCP file exists:", file.exists(vcp_file), "\n")
cat("Metadata file exists:", file.exists(metadata_file), "\n")

# Read files and check structure
if (file.exists(control_file)) {
  control_data <- read.csv(control_file, header = TRUE, nrows = 5)
  cat("\nControl file structure:\n")
  str(control_data)
  cat("\nControl file column names:\n")
  print(colnames(control_data))
}

if (file.exists(vcp_file)) {
  vcp_data <- read.csv(vcp_file, header = TRUE, nrows = 5)
  cat("\nVCP file structure:\n")
  str(vcp_data)
  cat("\nVCP file column names:\n")
  print(colnames(vcp_data))
}

if (file.exists(metadata_file)) {
  metadata <- read.csv(metadata_file, header = TRUE, nrows = 5)
  cat("\nMetadata file structure:\n")
  str(metadata)
  cat("\nMetadata file column names:\n")
  print(colnames(metadata))
}
