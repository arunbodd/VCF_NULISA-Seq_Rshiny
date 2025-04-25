# Release Notes

## v1.0.0 (April 25, 2025)

### Initial Release

This is the initial release of the VCP NULISA-Seq Rshiny application, a tool for analyzing proteomics data comparing VCP (Valosin-Containing Protein) and Control samples.

### Features

- **Dataset Overview**: View sample information and basic statistics
- **PCA Analysis**: Interactive PCA plot with customizable options
- **Differential Expression Analysis**: Volcano plot and table of differentially expressed proteins
- **Heatmap Visualization**: Customizable heatmap of top differentially expressed proteins
- **ROC Analysis**: Single and multi-protein ROC curve analysis with LOOCV
- **Protein Visualization**: Boxplots, violin plots, and line plots for individual proteins

### Technical Improvements

- Standardized naming convention across the application (VCP)
- Fixed deployment link in documentation
- Improved error handling for data loading
- Enhanced visualization options with downloadable plots

### Known Issues

- Large datasets may experience performance issues
- Some visualizations may require adjustments for optimal display on different screen sizes

### Future Plans

- Convert application to Python using Streamlit
- Add additional statistical analysis methods
- Implement batch correction features
- Enhance protein network visualization
