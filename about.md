# VCF vs Control Proteomics Analysis

This application provides interactive visualization and analysis tools for proteomics data comparing VCF (Valosin-Containing Protein) and Control samples.

## Features

- **Dataset Overview**: View sample information and basic statistics about the dataset
- **PCA Analysis**: Interactive Principal Component Analysis to visualize sample clustering and variance explained
- **Differential Expression**: Identify proteins with significant differences between VCF and Control samples using linear models
- **Heatmap**: Hierarchical clustering of samples and proteins with customizable visualization options
- **Gene Visualization**: Detailed visualization of individual protein expression with boxplots, violin plots, and line plots
- **ROC Analysis**: Evaluate the diagnostic potential of proteins using Leave-One-Out Cross-Validation

## Data

The application uses proteomics data from plasma samples, comparing VCF samples with Control samples. The data is preprocessed and normalized to allow for meaningful comparisons between groups.

## Analysis Methods

- **Principal Component Analysis (PCA)**: Reduces the dimensionality of the data to visualize patterns and clustering
- **Differential Expression**: Linear models are used to identify statistically significant differences between groups
- **Hierarchical Clustering**: Groups samples and proteins based on expression similarity
- **Statistical Tests**: Wilcoxon rank-sum tests are used to compare expression levels between groups
- **ROC Analysis**: Leave-One-Out Cross-Validation with logistic regression to evaluate diagnostic potential

## How to Use

1. Navigate through the tabs on the left sidebar to explore different analyses
2. Adjust parameters using the controls on the left side of each analysis tab
3. Interact with plots by hovering, zooming, and clicking on elements
4. Download data and plots using the download buttons where available

## Data

The application uses proteomics data from the NULISA-Seq platform, comparing VCF samples with Control samples. The data files are stored in the `data/` directory:

- `Controlonly_samples.csv`: Control sample data
- `VCP_samples.csv`: VCF sample data (Note: The file is still named VCP_samples.csv but the app refers to this as VCF)
- `Control_VCP_metadata.csv`: Sample metadata

## Contact

For questions or issues, please contact the application administrator.

## Citation

If you use this application in your research, please cite:

[Citation information will be added here]
