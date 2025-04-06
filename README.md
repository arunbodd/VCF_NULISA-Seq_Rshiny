# VCF vs Control Proteomics Analysis

A Shiny application for analyzing proteomics data comparing VCF and Control samples.

## Features

- **Dataset Overview**: View sample information and basic statistics
- **PCA Analysis**: Interactive PCA plot with customizable options
- **Differential Expression Analysis**: Volcano plot and table of differentially expressed proteins
- **Heatmap Visualization**: Customizable heatmap of top differentially expressed proteins
- **ROC Analysis**: Single and multi-protein ROC curve analysis with LOOCV
- **Protein Visualization**: Boxplots, violin plots, and line plots for individual proteins

## Requirements

The application requires the following R packages:

```r
library(shiny)
library(shinydashboard)
library(dplyr)
library(ggplot2)
library(DT)
library(pROC)
library(reshape2)
library(plotly)
library(ggrepel)
library(stringr)
library(tidyr)
library(pheatmap)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
library(STRINGdb)
library(igraph)
library(ggpubr)
library(ggthemes)
library(data.table)
library(corrplot)
```

## Running the Application

To run the application locally:

1. Clone this repository:
```bash
git clone https://github.com/arunbodd/VCF_NULISA-Seq_Rshiny.git
cd VCF_NULISA-Seq_Rshiny
```

2. Ensure your data files are in the `data/` directory:
   - `data/Controlonly_samples.csv`: Control sample data
   - `data/VCP_samples.csv`: VCF sample data
   - `data/Control_VCP_metadata.csv`: Sample metadata

3. Run the app:
```r
shiny::runApp()
```

## Live Demo

A live demo of the application is available at:

**[VCF vs Control Proteomics Analysis App](https://arunbodd.shinyapps.io/VCF-Proteomics-Analysis/)**

## Data Format Requirements

1. **Control Samples CSV**:
   - Must include columns: `Control.Plasma`, `Target`, and `NPQ`

2. **VCF Samples CSV**:
   - Must include columns: `VCP.Plasma`, `Target`, and `NPQ`
   - Note: The file is still named VCP_samples.csv but the app refers to this as VCF

3. **Metadata CSV**:
   - Must include a `Sample` column that matches the sample names in the data files

## Deployment

If you wish to deploy your own version of this app:

1. Install the rsconnect package:
```r
install.packages("rsconnect")
```

2. Set up your shinyapps.io account credentials using one of these secure methods:

   **Method A: Interactive setup (recommended for personal use)**
   ```r
   # Run this in the R console, NOT in a script that will be committed to version control
   rsconnect::setAccountInfo(
     name='YOUR_ACCOUNT_NAME',
     token='YOUR_TOKEN',
     secret='YOUR_SECRET'
   )
   ```

   **Method B: Environment variables (recommended for teams/CI)**
   ```r
   # Add to your .Renviron file (never commit this file to git)
   SHINYAPPS_NAME="YOUR_ACCOUNT_NAME"
   SHINYAPPS_TOKEN="YOUR_TOKEN"
   SHINYAPPS_SECRET="YOUR_SECRET"
   ```

3. Use the provided deployment script:
```r
source("deploy_app.R")
```

This will deploy the app with the data from your data/ directory.

### Security Best Practices

- **NEVER commit API keys or secrets to version control**
- Use environment variables or interactive setup for credentials
- Consider using `.gitignore` to exclude any files containing sensitive information
- For team deployments, use CI/CD secrets management

## Download Options

The application provides download options for all plots in various formats:
- PNG
- PDF
- TIFF
- JPEG

## Future Development

There are plans to convert this application to Python using Streamlit for better user experience and visualization capabilities.
