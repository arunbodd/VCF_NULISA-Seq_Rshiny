# Configuration file for VCF vs Control Proteomics Analysis app

# Path to data directory (use full path)
DATA_DIR <- "data"  # Data directory within the app folder

# Data file names
CONTROL_SAMPLES_FILE <- "Controlonly_samples.csv"
VCF_SAMPLES_FILE <- "VCP_samples.csv"
METADATA_FILE <- "Control_VCP_metadata.csv"

# Full paths to data files
CONTROL_SAMPLES_PATH <- file.path(DATA_DIR, CONTROL_SAMPLES_FILE)
VCF_SAMPLES_PATH <- file.path(DATA_DIR, VCF_SAMPLES_FILE)
METADATA_PATH <- file.path(DATA_DIR, METADATA_FILE)

# Function to check if data files exist
check_data_files <- function() {
  missing_files <- c()
  
  if (!file.exists(CONTROL_SAMPLES_PATH)) {
    missing_files <- c(missing_files, CONTROL_SAMPLES_PATH)
  }
  
  if (!file.exists(VCF_SAMPLES_PATH)) {
    missing_files <- c(missing_files, VCF_SAMPLES_PATH)
  }
  
  if (!file.exists(METADATA_PATH)) {
    missing_files <- c(missing_files, METADATA_PATH)
  }
  
  return(missing_files)
}
