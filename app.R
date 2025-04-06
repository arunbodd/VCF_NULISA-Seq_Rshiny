#!/usr/bin/env Rscript

library(shiny)
library(shinydashboard)
library(DT)
library(plotly)
library(ggplot2)
library(pheatmap)
library(reshape2)
library(dplyr)
library(tidyr)
library(ggpubr)
library(ggthemes)
library(ggrepel)
library(stringr)
library(data.table)
library(corrplot)
library(circlize)
library(RColorBrewer)
library(pROC)
library(caret)
library(heatmaply)
library(shinyBS)
library(shinyWidgets)
library(tibble)
library(ggiraph)
library(wesanderson)

# Global variables declaration
utils::globalVariables(c(
  "Protein", "Sample", "Target", "NPQ", "Expression", "Pvalue", "FDR", 
  "regression_coefficient", "std_error", "t_value", "significance_category",
  "PC1", "PC2", "Specificity", "Sensitivity", "specificity", "sensitivity",
  "binary_group", "1 - specificity", "1 - Specificity"
))

# Helper function to ensure the measurement column is named "NPQ"
fix_names <- function(df, default_measure_col = "NPQ") {
  # Trim any extra whitespace in column names
  colnames(df) <- trimws(colnames(df))
  if(ncol(df) < 3) stop("Input file must have at least three columns: Sample, Target, and a measurement column.")
  
  # Print column names for debugging
  cat("Column names before fixing:", paste(colnames(df), collapse=", "), "\n")
  
  # If the measurement column is not present, try to detect it
  if(!(default_measure_col %in% colnames(df))) {
    # Try case-insensitive match first
    candidate_cols <- grep(paste0("^", default_measure_col, "$"), colnames(df), ignore.case = TRUE, value = TRUE)
    if(length(candidate_cols) >= 1) {
      idx <- which(tolower(colnames(df)) == tolower(candidate_cols[1]))
      colnames(df)[idx] <- default_measure_col
      cat("Found case-insensitive match for", default_measure_col, "at index", idx, "\n")
    } else {
      # Try pattern matching
      candidate <- grep("npq|value|expression", tolower(colnames(df)), value = TRUE)
      if(length(candidate) >= 1) {
        idx <- which(tolower(colnames(df)) == candidate[1])
        colnames(df)[idx] <- default_measure_col
        cat("Found pattern match for", default_measure_col, "at index", idx, "\n")
      } else {
        # Fallback: assume the third column is the measurement column
        colnames(df)[3] <- default_measure_col
        cat("Using fallback: renamed column 3 to", default_measure_col, "\n")
      }
    }
  }
  
  # Print column names after fixing
  cat("Column names after fixing:", paste(colnames(df), collapse=", "), "\n")
  
  return(df)
}

# Use a relative path for the data directory
data_dir <- file.path(getwd(), "data")
required_files <- c(
  "Controlonly_samples.csv",
  "VCP_samples.csv",
  "Control_VCP_metadata.csv"
)

if (!dir.exists(data_dir)) {
  stop(paste("Data directory", data_dir, "does not exist. Please create it and add the required data files."))
}

missing_files <- required_files[!file.exists(file.path(data_dir, required_files))]
if (length(missing_files) > 0) {
  stop(paste("The following required data files are missing:", 
             paste(file.path(data_dir, missing_files), collapse = ", ")))
}

# Define UI
ui <- dashboardPage(
  dashboardHeader(title = "Proteomics Analysis App"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Data Overview", tabName = "dataoverview", icon = icon("table")),
      menuItem("PCA", tabName = "pca", icon = icon("chart-area")),
      menuItem("Differential Expression", tabName = "diffexpr", icon = icon("chart-bar")),
      menuItem("Heatmap", tabName = "heatmap", icon = icon("th")),
      menuItem("ROC Analysis", tabName = "roc", icon = icon("chart-line")),
      menuItem("LOOCV Analysis", tabName = "loocv", icon = icon("chart-line")),
      menuItem("Distribution Plots", tabName = "distplots", icon = icon("chart-bar")),
      menuItem("About", tabName = "about", icon = icon("info-circle"))
    ),
    hr(),
    radioButtons("downloadFormat", "Download Format:",
                 choices = c("PNG", "JPEG", "PDF"),
                selected = "Default")
  ),
  dashboardBody(
    tabItems(
      tabItem(tabName = "dataoverview",
              fluidRow(
                box(width = 12, title = "Data Overview", status = "primary",
                    helpText("This tab displays the differential protein expression results."),
                    fluidRow(
                      column(3, 
                             sliderInput("significance", "Significance Threshold:",
                                         min = 0.001, max = 0.1, value = 0.05, step = 0.001)),
                      column(3, 
                             numericInput("foldchange", "Fold Change Threshold:",
                                          value = 1.5, min = 1, max = 10, step = 0.1)),
                      column(3, 
                             actionButton("applyFilters", "Apply Filters", 
                                          icon = icon("filter"), 
                                          class = "btn-primary")),
                      column(3,
                             downloadButton("downloadData", "Download Data"))
                    ),
                    br(),
                    DTOutput("diffExprTable")
                )
              )
      ),
      tabItem(tabName = "pca",
              fluidRow(
                div(id = "pcaTab",
                    box(width = 8, title = "PCA Plot", status = "primary",
                        plotlyOutput("pcaPlot", height = "500px")
                    )
                ),
                box(width = 4, title = "PCA Options", status = "info",
                    sliderInput("pcaPointSize", "Point Size:", min = 1, max = 10, value = 3, step = 0.5),
                    checkboxInput("showLabels", "Show Sample Labels", FALSE),
                    selectInput("pcaColorScheme", "Color Scheme:",
                               choices = c("Default", "Zissou1", "Darjeeling1", "FantasticFox1", "Moonrise1", "Royal1", "BottleRocket1"),
                               selected = "Default"),
                    helpText("Principal Component Analysis (PCA) visualizes the main sources of variation in the data.")
                )
              )
      ),
      tabItem(tabName = "diffexpr",
              fluidRow(
                div(id = "volcanoTab",
                    box(width = 8, title = "Volcano Plot", status = "primary",
                        plotlyOutput("volcanoPlot", height = "500px")
                    )
                ),
                box(width = 4, title = "Volcano Plot Options", status = "info",
                    sliderInput("volcanoPointSize", "Point Size:", min = 1, max = 10, value = 3, step = 0.5),
                    checkboxInput("showVolcanoLabels", "Show Protein Labels for Significant Proteins", FALSE),
                    selectInput("volcanoColorScheme", "Color Scheme:",
                               choices = c("Default", "Zissou1", "Darjeeling1", "FantasticFox1", "Moonrise1", "Royal1", "BottleRocket1"),
                               selected = "Default"),
                    helpText("Volcano plot shows statistical significance vs. magnitude of change.")
                )
              )
      ),
      tabItem(tabName = "heatmap",
              fluidRow(
                div(id = "heatmapTab",
                    box(width = 8, title = "Heatmap", status = "primary",
                        plotlyOutput("heatmapPlot", height = "600px"),
                        textOutput("heatmapInfo")
                    )
                ),
                box(width = 4, title = "Heatmap Options", status = "info",
                    checkboxInput("clusterRows", "Cluster Rows", TRUE),
                    checkboxInput("clusterColumns", "Cluster Columns", TRUE),
                    sliderInput("heatmapFDR", "FDR Threshold:",
                               min = 0.001, max = 0.1, value = 0.05, step = 0.001),
                    sliderInput("heatmapFC", "Fold Change Threshold:",
                               min = 1, max = 5, value = 1.5, step = 0.1),
                    selectInput("heatmapColorScheme", "Color Scheme:",
                                choices = c("RdBu", "RdYlBu", "Spectral", "YlGnBu", "YlOrRd", "viridis", "magma", "plasma"),
                                selected = "RdBu"),
                    helpText("Heatmap visualizes expression patterns across samples for top differentially expressed proteins.")
                )
              )
      ),
      tabItem(tabName = "distplots",
              fluidRow(
                div(id = "distplotTab",
                    box(width = 8, title = "Distribution Plots", status = "primary",
                        plotOutput("distPlot", height = "600px"),
                        downloadButton("downloadDistPlot", "Download Plot")
                    ),
                    box(width = 4, title = "Plot Options", status = "info",
                        fluidRow(
                          column(4, 
                                 selectInput("distplotType", "Plot Type:",
                                            choices = c("Boxplot", "Violin Plot", "Expression Plot"),
                                            selected = "Boxplot")
                          ),
                          column(4,
                                 selectInput("distplotProtein", "Select Protein:", choices = NULL)
                          ),
                          column(4,
                                 conditionalPanel(
                                   condition = "input.distplotType == 'Boxplot' || input.distplotType == 'Violin Plot'",
                                   checkboxInput("showDistPoints", "Show Individual Points", TRUE)
                                 )
                          )
                        ),
                        selectInput("colorScheme", "Color Scheme:",
                                   choices = c("Default", "Zissou1", "Darjeeling1","FantasticFox1", "Moonrise1", "Royal1", "BottleRocket1"),
                                   selected = "Default"),
                        helpText("Visualize protein expression distributions across samples.")
                    )
                )
              )
      ),
      tabItem(tabName = "roc",
              fluidRow(
                div(id = "rocTab",
                    box(width = 8, title = "ROC Curves", status = "primary",
                        plotlyOutput("rocCurve", height = "600px")
                    )
                ),
                box(width = 4, title = "ROC Options", status = "info",
                    selectInput("rocprotein", "Select Protein:", choices = NULL),
                    helpText("ROC curves evaluate the classification performance of a protein.")
                )
              )
      ),
      tabItem(tabName = "loocv",
              fluidRow(
                div(id = "loocvTab",
                    box(width = 12, title = "LOOCV Analysis", status = "primary",
                        DTOutput("loocvTable"),
                        downloadButton("downloadLOOCV", "Download Results")
                    )
                )
              )
      ),
      tabItem(tabName = "about",
              fluidRow(
                box(width = 12, title = "About", status = "primary",
                    includeMarkdown("about.md")
                )
              )
      )
    )
  )
)

# Define server logic
server <- function(input, output, session) {
  
  # Reactive expression to load and process data
  data <- reactive({
    cat("Loading data files...\n")
    
    ## --- Load Control Data ---
    exprs_control <- read.csv(file.path(data_dir, "Controlonly_samples.csv"), header = TRUE, sep = ",", stringsAsFactors = FALSE)
    cat("Control samples loaded. Dimensions:", nrow(exprs_control), "x", ncol(exprs_control), "\n")
    
    # Print original column names for debugging
    cat("Original Control column names:", paste(colnames(exprs_control), collapse=", "), "\n")
    
    # Trim column names and rename first two columns
    colnames(exprs_control) <- trimws(colnames(exprs_control))
    if (!("Sample" %in% colnames(exprs_control))) {
      cat("Renaming first column from", colnames(exprs_control)[1], "to Sample\n")
      colnames(exprs_control)[1] <- "Sample"
    }
    if (!("Target" %in% colnames(exprs_control))) {
      cat("Renaming second column from", colnames(exprs_control)[2], "to Target\n")
      colnames(exprs_control)[2] <- "Target"
    }
    
    # Ensure NPQ column exists
    exprs_control <- fix_names(exprs_control, "NPQ")
    
    # Check if NPQ column exists after fixing
    if (!("NPQ" %in% colnames(exprs_control))) {
      stop("NPQ column not found in Control data after fixing column names")
    }
    
    # Print column names before melting
    cat("Control column names before melting:", paste(colnames(exprs_control), collapse=", "), "\n")
    
    tryCatch({
      # Ensure NPQ column is numeric
      exprs_control$NPQ <- as.numeric(exprs_control$NPQ)
      
      # Create a wide-format data frame directly
      Cntrl_reshaped_data <- reshape2::dcast(exprs_control, Target ~ Sample, value.var = "NPQ")
      cat("Control reshaped data dimensions:", nrow(Cntrl_reshaped_data), "x", ncol(Cntrl_reshaped_data), "\n")
    }, error = function(e) {
      cat("Error in Control data reshaping:", e$message, "\n")
      # Alternative approach if the above fails
      cat("Trying alternative reshaping approach for Control data...\n")
      control_matrix <- matrix(exprs_control$NPQ, 
                               nrow = length(unique(exprs_control$Target)),
                               dimnames = list(unique(exprs_control$Target), unique(exprs_control$Sample)))
      Cntrl_reshaped_data <- as.data.frame(control_matrix)
      Cntrl_reshaped_data$Target <- rownames(Cntrl_reshaped_data)
      Cntrl_reshaped_data <- Cntrl_reshaped_data[, c("Target", setdiff(colnames(Cntrl_reshaped_data), "Target"))]
      cat("Control reshaped data dimensions (alternative):", nrow(Cntrl_reshaped_data), "x", ncol(Cntrl_reshaped_data), "\n")
    })
    
    ## --- Load VCP Data ---
    exprs_vcp <- read.csv(file.path(data_dir, "VCP_samples.csv"), header = TRUE, sep = ",", stringsAsFactors = FALSE)
    cat("VCP samples loaded. Dimensions:", nrow(exprs_vcp), "x", ncol(exprs_vcp), "\n")
    
    # Print original column names for debugging
    cat("Original VCP column names:", paste(colnames(exprs_vcp), collapse=", "), "\n")
    
    colnames(exprs_vcp) <- trimws(colnames(exprs_vcp))
    if (!("Sample" %in% colnames(exprs_vcp))) {
      cat("Renaming first column from", colnames(exprs_vcp)[1], "to Sample\n")
      colnames(exprs_vcp)[1] <- "Sample"
    }
    if (!("Target" %in% colnames(exprs_vcp))) {
      cat("Renaming second column from", colnames(exprs_vcp)[2], "to Target\n")
      colnames(exprs_vcp)[2] <- "Target"
    }
    
    # Ensure NPQ column exists
    exprs_vcp <- fix_names(exprs_vcp, "NPQ")
    
    # Check if NPQ column exists after fixing
    if (!("NPQ" %in% colnames(exprs_vcp))) {
      stop("NPQ column not found in VCP data after fixing column names")
    }
    
    # Print column names before melting
    cat("VCP column names before melting:", paste(colnames(exprs_vcp), collapse=", "), "\n")
    
    tryCatch({
      exprs_vcp$NPQ <- as.numeric(exprs_vcp$NPQ)
      VCP_reshaped_data <- reshape2::dcast(exprs_vcp, Target ~ Sample, value.var = "NPQ")
      cat("VCP reshaped data dimensions:", nrow(VCP_reshaped_data), "x", ncol(VCP_reshaped_data), "\n")
    }, error = function(e) {
      cat("Error in VCP data reshaping:", e$message, "\n")
      cat("Trying alternative reshaping approach for VCP data...\n")
      vcp_matrix <- matrix(exprs_vcp$NPQ, 
                           nrow = length(unique(exprs_vcp$Target)),
                           dimnames = list(unique(exprs_vcp$Target), unique(exprs_vcp$Sample)))
      VCP_reshaped_data <- as.data.frame(vcp_matrix)
      VCP_reshaped_data$Target <- rownames(VCP_reshaped_data)
      VCP_reshaped_data <- VCP_reshaped_data[, c("Target", setdiff(colnames(VCP_reshaped_data), "Target"))]
      cat("VCP reshaped data dimensions (alternative):", nrow(VCP_reshaped_data), "x", ncol(VCP_reshaped_data), "\n")
    })
    
    ## --- Load Metadata ---
    metadata <- read.csv(file.path(data_dir, "Control_VCP_metadata.csv"), header = TRUE, sep = ",")
    cat("Metadata loaded. Dimensions:", nrow(metadata), "x", ncol(metadata), "\n")
    colnames(metadata) <- trimws(colnames(metadata))
    # Assume metadata has columns 'Sample' and 'Treatment'
    metadata$Sample <- toupper(gsub("[^[:alnum:]]", "", metadata$Sample))
    
    # Clean sample names in expression data but preserve original format
    # Store original sample names for display
    original_cntrl_samples <- colnames(Cntrl_reshaped_data)[-1]
    original_vcp_samples <- colnames(VCP_reshaped_data)[-1]
    
    # Create a mapping of original to processed sample names
    sample_name_map <- data.frame(
      Original = c(original_cntrl_samples, original_vcp_samples),
      Processed = toupper(gsub("[^[:alnum:]]", "", c(original_cntrl_samples, original_vcp_samples))),
      Group = c(rep("Control", length(original_cntrl_samples)), rep("VCP", length(original_vcp_samples)))
    )
    
    # Now proceed with the standard processing for merging
    colnames(Cntrl_reshaped_data)[-1] <- toupper(gsub("[^[:alnum:]]", "", colnames(Cntrl_reshaped_data)[-1]))
    colnames(VCP_reshaped_data)[-1] <- toupper(gsub("[^[:alnum:]]", "", colnames(VCP_reshaped_data)[-1]))
    
    ## --- Merge Expression Data ---
    cat("Merging control and VCP data...\n")
    merged_exprn <- merge(Cntrl_reshaped_data, VCP_reshaped_data, by = "Target", all = TRUE)
    cat("Merged data dimensions:", nrow(merged_exprn), "x", ncol(merged_exprn), "\n")
    
    na_count <- sum(is.na(merged_exprn))
    if (na_count > 0) {
      cat("WARNING: Found", na_count, "NA values in merged data\n")
    }
    
    # Automatically remove APOE protein (case-insensitive)
    merged_exprn <- merged_exprn %>% filter(!grepl("^APOE$", Target, ignore.case = TRUE))
    cat("APOE protein removed. New dimensions:", nrow(merged_exprn), "x", ncol(merged_exprn), "\n")
    
    # Preserve protein names by creating a 'Protein' column,
    # then set rownames and remove the original Target column.
    rownames(merged_exprn) <- merged_exprn$Target
    merged_exprn$Protein <- merged_exprn$Target
    merged_exprn <- merged_exprn[, c("Protein", setdiff(colnames(merged_exprn), c("Target", "Protein")))]
    # Convert all sample columns to numeric (skip the Protein column)
    merged_exprn[, setdiff(names(merged_exprn), "Protein")] <- lapply(merged_exprn[, setdiff(names(merged_exprn), "Protein")], as.numeric)
    
    ## --- Check Metadata vs. Expression Sample Names ---
    expr_sample_names <- colnames(merged_exprn)[-1]  # exclude the Protein column
    if (!all(metadata$Sample %in% expr_sample_names)) {
      warning("Some sample names in metadata do not match the expression data. Only matching samples will be used.")
      metadata <- metadata %>% filter(Sample %in% expr_sample_names)
    }
    
    # Create Group column based on sample names.
    control_samples <- setdiff(colnames(Cntrl_reshaped_data), "Target")
    vcp_samples <- setdiff(colnames(VCP_reshaped_data), "Target")
    
    metadata <- metadata %>%
      mutate(Group = case_when(
        Sample %in% control_samples ~ "Control",
        Sample %in% vcp_samples ~ "VCP",
        TRUE ~ NA_character_
      )) %>%
      filter(!is.na(Group))
    
    cat("Number of samples in metadata after grouping:", nrow(metadata), "\n")
    print(table(metadata$Group))
    
    ordered_metadata <- metadata %>% arrange(Group)
    
    cat("Final merged expression data dimensions:", nrow(merged_exprn), "x", ncol(merged_exprn), "\n")
    cat("Metadata dimensions:", nrow(metadata), "x", ncol(metadata), "\n")
    
    # Create long format version for distribution plots
    merged_exprn_long <- merged_exprn %>% 
      pivot_longer(cols = -Protein, names_to = "Sample", values_to = "NPQ") %>%
      left_join(ordered_metadata, by = "Sample")
    
    list(
      merged_exprn = merged_exprn,
      merged_exprn_long = merged_exprn_long,
      metadata = metadata,
      ordered_metadata = ordered_metadata,
      sample_name_map = sample_name_map
    )
  })
  
  # Perform PCA
  pca_data <- reactive({
    req(data())
    cat("Performing PCA...\n")
    pca_result <- prcomp(t(data()$merged_exprn[,-1]), scale. = TRUE)  # exclude Protein column
    cat("PCA complete\n")
    pca_result
  })
  
  # Differential expression analysis
  diff_expr <- reactive({
    req(data())
    cat("Starting differential expression analysis\n")
    
    merged_exprn <- data()$merged_exprn
    ordered_metadata <- data()$ordered_metadata
    
    if (is.null(merged_exprn) || nrow(merged_exprn) == 0) {
      cat("ERROR: No expression data available\n")
      return(NULL)
    }
    if (is.null(ordered_metadata) || nrow(ordered_metadata) == 0) {
      cat("ERROR: No metadata available\n")
      return(NULL)
    }
    
    # Use rownames (which are protein names)
    allProteins <- rownames(merged_exprn)
    if (length(allProteins) == 0) {
      cat("ERROR: No proteins found in expression data\n")
      return(NULL)
    }
    
    cat("Number of proteins to analyze:", length(allProteins), "\n")
    
    allProteins_modelstats <- vector(mode = 'list', length = length(allProteins))
    names(allProteins_modelstats) <- allProteins
    proteins_processed <- 0
    
    for (idx in seq_along(allProteins)) {
      p <- allProteins[idx]
      protein_expr <- merged_exprn[p, setdiff(colnames(merged_exprn), "Protein")]
      model_data <- data.frame(
        Expression = as.numeric(protein_expr),
        Group = ordered_metadata$Group,
        Sample = ordered_metadata$Sample
      )
      model_data <- model_data[complete.cases(model_data), ]
      
      if (nrow(model_data) < 3) {
        cat("WARNING: Not enough data points for protein", p, "\n")
        next
      }
      if (length(unique(model_data$Group)) < 2) {
        cat("WARNING: Not enough groups for protein", p, "\n")
        next
      }
      
      tryCatch({
        model <- lm(Expression ~ Group, data = model_data)
        model_summary <- summary(model)
        allProteins_modelstats[[idx]] <- list(
          protein = p,
          regression_coefficient = model_summary$coefficients[2, 1],
          std_error = model_summary$coefficients[2, 2],
          t_value = model_summary$coefficients[2, 3],
          p_value = model_summary$coefficients[2, 4]
        )
        proteins_processed <- proteins_processed + 1
        if (proteins_processed %% 100 == 0) {
          cat("Processed", proteins_processed, "proteins\n")
        }
      }, error = function(e) {
        cat("ERROR processing protein", p, ":", e$message, "\n")
      })
    }
    
    allProteins_modelstats <- allProteins_modelstats[!sapply(allProteins_modelstats, is.null)]
    if (length(allProteins_modelstats) == 0) {
      cat("ERROR: No valid statistical results\n")
      return(NULL)
    }
    
    allProteins_modelstats_df <- do.call(rbind, lapply(allProteins_modelstats, function(x) {
      data.frame(
        Protein = x$protein,
        regression_coefficient = x$regression_coefficient,
        std_error = x$std_error,
        t_value = x$t_value,
        Pvalue = x$p_value
      )
    }))
    
    allProteins_modelstats_df <- allProteins_modelstats_df %>%
      mutate(FDR = p.adjust(Pvalue, method = "fdr"))
    
    cat("Differential expression analysis complete. Proteins analyzed:", nrow(allProteins_modelstats_df), "\n")
    allProteins_modelstats_df
  })
  
  top_proteins <- reactive({
    req(diff_expr())
    diff_expr() %>% 
      filter(FDR < 0.05 & abs(regression_coefficient) > log2(1.5)) %>% 
      arrange(regression_coefficient)
  })
  
  # Update protein selection dropdowns when data is loaded
  observe({
    req(diff_expr())
    protein_choices <- diff_expr()$Protein
    updateSelectInput(session, "distplotProtein", choices = protein_choices)
    updateSelectInput(session, "rocprotein", choices = protein_choices)
    updateSelectInput(session, "loocvprotein", choices = protein_choices)
  })
  
  filtered_diff_expr <- reactive({
    req(diff_expr(), input$significance, input$foldchange)
    diff_expr() %>%
      filter(FDR < input$significance & abs(regression_coefficient) > log2(input$foldchange))
  })
  
  output$diffExprTable <- renderDT({
    req(filtered_diff_expr())
    datatable(filtered_diff_expr(), 
              options = list(pageLength = 10, scrollX = TRUE),
              filter = 'top',
              rownames = FALSE) %>%
      formatRound(columns = c('regression_coefficient', 'std_error', 't_value', 'Pvalue', 'FDR'), digits = 4)
  })
  
  output$downloadData <- downloadHandler(
    filename = function() {
      paste("differential_expression_", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(filtered_diff_expr(), file, row.names = FALSE)
    }
  )
  
  output$pcaPlot <- renderPlotly({
    req(pca_data())
    pca_df <- as.data.frame(pca_data()$x)
    pca_df$Sample <- rownames(pca_df)
    meta <- data()$ordered_metadata
    pca_df <- merge(pca_df, meta, by = "Sample")
    
    # Get original sample names and create a mapping for display
    sample_map <- data()$sample_name_map
    
    # Create display names for hover text
    pca_df$DisplayName <- sapply(pca_df$Sample, function(samp) {
      idx <- which(sample_map$Processed == samp)
      if(length(idx) > 0) {
        return(sub("_Emory.*", "", sample_map$Original[idx[1]]))
      } else {
        return(samp)
      }
    })
    
    if(input$pcaColorScheme != "Default") {
      palette <- wesanderson::wes_palette(name = input$pcaColorScheme, type = "discrete")
      p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Group, text = DisplayName)) +
        geom_point(size = input$pcaPointSize) +
        labs(
          title = "PCA Plot",
          x = paste0("PC1 (", round(summary(pca_data())$importance[2,1] * 100, 2), "% Variance)"),
          y = paste0("PC2 (", round(summary(pca_data())$importance[2,2] * 100, 2), "% Variance)")
        ) +
        theme_minimal() +
        scale_color_manual(values = palette)
    } else {
      p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Group, text = DisplayName)) +
        geom_point(size = input$pcaPointSize) +
        labs(
          title = "PCA Plot",
          x = paste0("PC1 (", round(summary(pca_data())$importance[2,1] * 100, 2), "% Variance)"),
          y = paste0("PC2 (", round(summary(pca_data())$importance[2,2] * 100, 2), "% Variance)")
        ) +
        theme_minimal() +
        scale_color_manual(values = c("Control" = "#66c2a5", "VCP" = "#fc8d62"))
    }
    
    if (input$showLabels) {
      p <- p + geom_text_repel(aes(label = Sample), size = 3)
    }
    
    ggplotly(p, tooltip = c("text", "color")) %>% 
      layout(legend = list(orientation = "v", y = 0.9, x = 1.05, xanchor = "left"))
  })
  
  output$volcanoPlot <- renderPlotly({
    req(diff_expr())
    
    volcanoPlotDat <- diff_expr() %>%
      filter(!is.na(Pvalue) & !is.infinite(Pvalue) & !is.na(regression_coefficient) & !is.infinite(regression_coefficient)) %>%
      mutate(significance_category = case_when(
        Pvalue < input$significance & regression_coefficient > log2(input$foldchange) ~ "Significant & Upregulated",
        Pvalue < input$significance & regression_coefficient < -log2(input$foldchange) ~ "Significant & Downregulated",
        TRUE ~ "Not Significant"
      ))
    
    p <- ggplot(volcanoPlotDat, aes(x = regression_coefficient, y = -log10(Pvalue), color = significance_category,
                                    text = paste("Protein:", Protein, "<br>log2 Fold Change:", round(regression_coefficient, 3),
                                                 "<br>-log10 p-value:", round(-log10(Pvalue), 3)))) +
      geom_point(size = input$volcanoPointSize, alpha = 0.6) +
      geom_hline(yintercept = -log10(input$significance), linetype = "dashed") +
      geom_vline(xintercept = c(-log2(input$foldchange), log2(input$foldchange)), linetype = "dashed") +
      labs(
        title = "Volcano Plot",
        x = "log2 Fold Change",
        y = "-log10 P-value"
      ) +
      theme_minimal() +
      scale_color_manual(values = c(
        "Significant & Upregulated" = "#fc8d62",
        "Significant & Downregulated" = "#66c2a5",
        "Not Significant" = "#cccccc"
      ))
    
    if(input$volcanoColorScheme != "Default") {
      palette <- wesanderson::wes_palette(name = input$volcanoColorScheme, type = "discrete")
      p <- p + scale_color_manual(values = palette)
    }
    
    if (input$showVolcanoLabels) {
      sig_proteins <- volcanoPlotDat %>% filter(Pvalue < input$significance & abs(regression_coefficient) > log2(input$foldchange))
      p <- p + geom_text_repel(data = sig_proteins, aes(label = Protein), size = 3, box.padding = 0.5, point.padding = 0.2)
    }
    
    ggplotly(p, tooltip = "text") %>% layout(legend = list(orientation = "v", y = 0.9, x = 1.05, xanchor = "left"))
  })
  
  output$heatmapInfo <- renderText({
    req(diff_expr())
    top_proteins <- diff_expr() %>% 
      filter(FDR < input$heatmapFDR & abs(regression_coefficient) > log2(input$heatmapFC)) %>% 
      pull(Protein)
    paste("Displaying", length(top_proteins), "proteins that meet the threshold criteria.")
  })
  
  output$heatmapPlot <- renderPlotly({
    req(diff_expr(), data())
    
    # Diagnostic output to check filtering
    cat("FDR threshold:", input$heatmapFDR, "\n")
    cat("Fold change threshold:", input$heatmapFC, "\n")
    
    # Count proteins meeting criteria before filtering
    total_proteins <- nrow(diff_expr())
    cat("Total proteins before filtering:", total_proteins, "\n")
    
    # Count proteins meeting just FDR criteria
    fdr_filtered <- diff_expr() %>% filter(FDR < input$heatmapFDR)
    cat("Proteins meeting FDR <", input$heatmapFDR, ":", nrow(fdr_filtered), "\n")
    
    # Count proteins meeting just fold change criteria
    fc_filtered <- diff_expr() %>% filter(abs(regression_coefficient) > log2(input$heatmapFC))
    cat("Proteins meeting fold change >", input$heatmapFC, ":", nrow(fc_filtered), "\n")
    
    # Get proteins meeting both criteria
    top_proteins <- diff_expr() %>% 
      filter(FDR < input$heatmapFDR & abs(regression_coefficient) > log2(input$heatmapFC)) %>% 
      arrange(regression_coefficient) %>% pull(Protein)
    
    cat("Proteins meeting both criteria:", length(top_proteins), "\n")
    
    validate(
      need(length(top_proteins) > 0, "No proteins meet the current threshold criteria. Please adjust the FDR or fold change thresholds.")
    )
    
    top_proteins_expression <- data()$merged_exprn[top_proteins, setdiff(colnames(data()$merged_exprn), "Protein")]
    top_proteins_expression_scaled <- t(scale(t(top_proteins_expression)))
    
    # Get original sample names and create a mapping for display
    sample_map <- data()$sample_name_map
    
    # Create a named vector for column renaming
    col_names <- colnames(top_proteins_expression_scaled)
    display_names <- sapply(col_names, function(col) {
      idx <- which(sample_map$Processed == col)
      if(length(idx) > 0) {
        return(sub("_Emory.*", "", sample_map$Original[idx[1]]))
      } else {
        return(col)
      }
    })
    
    # Rename columns for display
    colnames(top_proteins_expression_scaled) <- display_names
    
    annotation_df <- data()$metadata %>%
      select(Sample, Group) %>% column_to_rownames("Sample")
      
    # Update row names in annotation to match display names
    new_row_names <- sapply(rownames(annotation_df), function(row) {
      idx <- which(sample_map$Processed == row)
      if(length(idx) > 0) {
        return(sub("_Emory.*", "", sample_map$Original[idx[1]]))
      } else {
        return(row)
      }
    })
    
    # Create new annotation dataframe with updated row names
    annotation_df_new <- annotation_df
    rownames(annotation_df_new) <- new_row_names
    
    heatmaply(
      top_proteins_expression_scaled,
      Rowv = input$clusterRows,
      Colv = input$clusterColumns,
      col_side_colors = annotation_df_new,
      hclust_method = "complete",
      dist_method = "euclidean",
      main = "Heatmap of Differentially Expressed Proteins",
      xlab = "Samples",
      ylab = "Proteins",
      colors = brewer.pal(9, input$heatmapColorScheme)
    )
  })
  
  output$distPlot <- renderPlot({
    req(input$distplotProtein)
    
    # Debug data loading
    if(is.null(data())) {
      message("Data object is NULL")
      return(ggplot() + 
               geom_text(aes(0,0, label="Data object is NULL - check data loading")) + 
               theme_void())
    }
    
    if(!exists("merged_exprn_long", where=data())) {
      message("merged_exprn_long not found in data object. Available objects: ", 
              paste(names(data()), collapse=", "))
      return(ggplot() + 
               geom_text(aes(0,0, label="Data structure issue - check data processing")) + 
               theme_void())
    }
    
    plot_data <- tryCatch({
      # Get the raw data first
      raw_data <- data()$merged_exprn_long %>%
        filter(Protein == input$distplotProtein)
      
      # Join with sample name map to get original names
      sample_map <- data()$sample_name_map
      
      processed_data <- raw_data %>%
        left_join(sample_map, by = c("Sample" = "Processed")) %>%
        mutate(DisplayName = ifelse(!is.na(Original), 
                                    sub("_Emory.*", "", Original), 
                                    Sample)) %>%
        # Use DisplayName for plotting but keep Group
        mutate(Sample = DisplayName) %>%
        select(-Original, -DisplayName, -Group.y) %>%
        rename(Group = Group.x)
      
      message("Processed sample names: ", paste(unique(processed_data$Sample), collapse=", "))
      
      processed_data
    }, error = function(e) {
      message("Error processing sample names: ", e$message)
      # Fallback to original data if mapping fails
      data()$merged_exprn_long %>%
        filter(Protein == input$distplotProtein)
    })
    
    if(is.null(plot_data)) {
      return(ggplot() + 
               geom_text(aes(0,0, label="Error filtering data")) + 
               theme_void())
    }
    
    if(nrow(plot_data) == 0) {
      return(ggplot() + 
               geom_text(aes(0,0, label="No data for selected protein")) + 
               theme_void())
    }
    
    stat_test <- t.test(NPQ ~ Group, data = plot_data)
    
    if(input$distplotType == "Boxplot") {
      p <- ggplot(plot_data, aes(x = Group, y = NPQ, fill = Group)) +
        geom_boxplot(outlier.shape = NA) +
        {if(input$showDistPoints) geom_jitter(width = 0.1, alpha = 0.5)} +
        labs(title = paste("Boxplot for", input$distplotProtein),
             subtitle = paste0("p-value: ", signif(stat_test$p.value, 3), " (t-test)")) +
        theme_minimal() +
        theme(legend.position = "none")
    } 
    else if(input$distplotType == "Violin Plot") {
      p <- ggplot(plot_data, aes(x = Group, y = NPQ, fill = Group)) +
        geom_violin(trim = FALSE) +
        geom_boxplot(width = 0.1, fill = "white", alpha = 0.5) +
        {if(input$showDistPoints) geom_jitter(width = 0.1, alpha = 0.5)} +
        labs(title = paste("Violin Plot for", input$distplotProtein),
             subtitle = paste0("p-value: ", signif(stat_test$p.value, 3), " (t-test)")) +
        theme_minimal() +
        theme(legend.position = "none")
    }
    else { # Expression Plot
      p <- ggplot(plot_data, aes(x = Group, y = NPQ, color = Group)) +
        geom_boxplot(width = 0.5, alpha = 0.2) +
        geom_jitter(width = 0.2, size = 3, aes(text = Sample)) +
        stat_summary(fun = mean, geom = "point", size = 4, shape = 23, fill = "white") +
        labs(title = paste("Expression Plot for", input$distplotProtein),
             subtitle = paste0("p-value: ", signif(stat_test$p.value, 3), " (t-test)"),
             y = "Expression (NPQ)") +
        theme_minimal() +
        theme(legend.position = "none")
    }
    
    # Apply color scheme
    if(input$colorScheme != "Default") {
      palette <- wesanderson::wes_palette(name = input$colorScheme, type = "discrete")
      if(input$distplotType == "Expression Plot") {
        p <- p + scale_color_manual(values = palette)
      } else {
        p <- p + scale_fill_manual(values = palette)
      }
    }
    
    p
  })
  
  output$downloadDistPlot <- downloadHandler(
    filename = function() {
      paste(tolower(gsub(" ", "_", input$distplotType)), "_", 
            input$distplotProtein, "_", Sys.Date(), ".", 
            tolower(input$downloadFormat), sep = "")
    },
    content = function(file) {
      plot_data <- tryCatch({
        # Get the raw data first
        raw_data <- data()$merged_exprn_long %>%
          filter(Protein == input$distplotProtein)
        
        # Join with sample name map to get original names
        sample_map <- data()$sample_name_map
        
        processed_data <- raw_data %>%
          left_join(sample_map, by = c("Sample" = "Processed")) %>%
          mutate(DisplayName = ifelse(!is.na(Original), 
                                      sub("_Emory.*", "", Original), 
                                      Sample)) %>%
          # Use DisplayName for plotting but keep Group
          mutate(Sample = DisplayName) %>%
          select(-Original, -DisplayName, -Group.y) %>%
          rename(Group = Group.x)
        
        message("Processed sample names: ", paste(unique(processed_data$Sample), collapse=", "))
        
        processed_data
      }, error = function(e) {
        message("Error processing sample names: ", e$message)
        # Fallback to original data if mapping fails
        data()$merged_exprn_long %>%
          filter(Protein == input$distplotProtein)
      })
      
      if(is.null(plot_data) || nrow(plot_data) == 0) {
        return(NULL)
      }
      
      stat_test <- t.test(NPQ ~ Group, data = plot_data)
      
      if(input$distplotType == "Boxplot") {
        p <- ggplot(plot_data, aes(x = Group, y = NPQ, fill = Group)) +
          geom_boxplot(outlier.shape = NA) +
          {if(input$showDistPoints) geom_jitter(width = 0.1, alpha = 0.5)} +
          labs(title = paste("Boxplot for", input$distplotProtein),
               subtitle = paste0("p-value: ", signif(stat_test$p.value, 3), " (t-test)")) +
          theme_minimal() +
          theme(legend.position = "none")
      } 
      else if(input$distplotType == "Violin Plot") {
        p <- ggplot(plot_data, aes(x = Group, y = NPQ, fill = Group)) +
          geom_violin(trim = FALSE) +
          geom_boxplot(width = 0.1, fill = "white", alpha = 0.5) +
          {if(input$showDistPoints) geom_jitter(width = 0.1, alpha = 0.5)} +
          labs(title = paste("Violin Plot for", input$distplotProtein),
               subtitle = paste0("p-value: ", signif(stat_test$p.value, 3), " (t-test)")) +
          theme_minimal() +
          theme(legend.position = "none")
      }
      else { # Expression Plot
        p <- ggplot(plot_data, aes(x = Group, y = NPQ, color = Group)) +
          geom_boxplot(width = 0.5, alpha = 0.2) +
          geom_jitter(width = 0.2, size = 3, aes(text = Sample)) +
          stat_summary(fun = mean, geom = "point", size = 4, shape = 23, fill = "white") +
          labs(title = paste("Expression Plot for", input$distplotProtein),
               subtitle = paste0("p-value: ", signif(stat_test$p.value, 3), " (t-test)"),
               y = "Expression (NPQ)") +
          theme_minimal() +
          theme(legend.position = "none")
      }
      
      # Apply color scheme
      if(input$colorScheme != "Default") {
        palette <- wesanderson::wes_palette(name = input$colorScheme, type = "discrete")
        if(input$distplotType == "Expression Plot") {
          p <- p + scale_color_manual(values = palette)
        } else {
          p <- p + scale_fill_manual(values = palette)
        }
      }
      
      if (tolower(input$downloadFormat) == "png") {
        ggsave(file, plot = p, device = "png", width = 10, height = 8, dpi = 300)
      } else if (tolower(input$downloadFormat) == "jpeg") {
        ggsave(file, plot = p, device = "jpeg", width = 10, height = 8, dpi = 300)
      } else {
        ggsave(file, plot = p, device = "pdf", width = 10, height = 8)
      }
    }
  )
  
  output$rocCurve <- renderPlotly({
    req(diff_expr(), data(), input$rocprotein)
    
    protein_expr <- as.numeric(data()$merged_exprn[input$rocprotein, setdiff(colnames(data()$merged_exprn), "Protein")])
    sample_names <- colnames(data()$merged_exprn)[-1]
    group_info <- data()$ordered_metadata$Group
    
    protein_data <- data.frame(
      Sample = sample_names,
      Group = group_info,
      Expression = protein_expr,
      Protein = input$rocprotein,
      stringsAsFactors = FALSE
    )
    
    validate(
      need(nrow(protein_data) > 0, "No data available for this protein"),
      need(!all(is.na(protein_data$Expression)), "Expression values are all NA for this protein"),
      need(n_distinct(protein_data$Group) >= 2, "Need at least two groups for ROC analysis")
    )
    
    protein_data <- protein_data[!is.na(protein_data$Expression), ]
    
    validate(
      need(sum(protein_data$Group == "VCP") > 0, "No positive cases for ROC analysis"),
      need(sum(protein_data$Group == "Control") > 0, "No negative cases for ROC analysis")
    )
    
    tryCatch({
      protein_data$binary_group <- ifelse(protein_data$Group == "VCP", 1, 0)
      roc_obj <- roc(protein_data$binary_group, protein_data$Expression)
      auc_value <- auc(roc_obj)
      
      roc_data <- data.frame(
        specificity = roc_obj$sensitivities,
        sensitivity = roc_obj$specificities
      )
      
      # Create a static ggplot version for downloads
      static_p <- ggplot(roc_data, aes(x = 1 - specificity, y = sensitivity)) +
        geom_line(color = "red", size = 1) +
        geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
        labs(
          title = paste("ROC Curve for", input$rocprotein),
          x = "False Positive Rate (1 - Specificity)",
          y = "True Positive Rate (Sensitivity)"
        ) +
        annotate("text", x = 0.7, y = 0.1, 
                 label = paste("AUC =", round(auc_value, 3)),
                 size = 5, color = "red", fontface = "bold") +
        theme_minimal() +
        coord_equal()
      
      # Store the static plot for download handler to use
      output$static_roc_plot <- renderPlot({
        static_p
      })
      
      # Create an interactive plotly version for display
      p <- plot_ly(roc_data, x = ~(1 - specificity), y = ~sensitivity, type = 'scatter', mode = 'lines',
                   line = list(color = 'red', width = 2), name = input$rocprotein) %>%
        add_trace(x = c(0, 1), y = c(0, 1), mode = 'lines', line = list(color = 'gray', width = 1, dash = 'dash'),
                  name = 'Random', showlegend = FALSE) %>%
        layout(
          title = paste("ROC Curve for", input$rocprotein),
          xaxis = list(title = "False Positive Rate (1 - Specificity)"),
          yaxis = list(title = "True Positive Rate (Sensitivity)"),
          annotations = list(
            list(
              x = 0.7,
              y = 0.1,
              text = paste("AUC =", round(auc_value, 3)),
              showarrow = FALSE,
              font = list(size = 14, color = "red"),
              bgcolor = "rgba(255, 255, 255, 0.8)",
              bordercolor = "rgba(0, 0, 0, 0.2)",
              borderwidth = 1
            )
          ),
          shapes = list(
            list(
              type = "rect",
              x0 = 0, x1 = 1, y0 = 0, y1 = 1,
              line = list(color = "rgba(0,0,0,0)")
            )
          )
        )
      
      p
    }, error = function(e) {
      cat("Error in ROC analysis:", e$message, "\n")
      plotly_empty(type = "scatter", mode = "markers") %>%
        layout(title = paste("Error in ROC analysis for", input$rocprotein),
               annotations = list(
                 x = 0.5,
                 y = 0.5,
                 text = e$message,
                 showarrow = FALSE
               ))
    })
  })
 
  
  # LOOCV function with skipping of problematic folds
  calc_loocv_auc <- function(data, protein) {
    protein_data <- data.frame(
      Expression = as.numeric(data$merged_exprn[protein, setdiff(colnames(data$merged_exprn), "Protein")]),
      Group = data$ordered_metadata$Group
    )
    protein_data$binary_group <- ifelse(protein_data$Group == "VCP", 1, 0)
    loocv_predicted_prob <- rep(NA, nrow(protein_data))
    loocv_true_class <- protein_data$binary_group
    
    for (i in 1:nrow(protein_data)) {
      train_data <- protein_data[-i, ]
      test_data <- protein_data[i, ]
      if (length(unique(train_data$binary_group)) < 2) {
        next  # Skip fold if training set lacks class diversity
      }
      model <- glm(binary_group ~ Expression, data = train_data, family = "binomial")
      loocv_predicted_prob[i] <- predict(model, newdata = test_data, type = "response")
    }
    
    valid_idx <- !is.na(loocv_predicted_prob)
    if(sum(valid_idx) < 2) return(NULL)
    
    roc_loocv <- roc(loocv_true_class[valid_idx], loocv_predicted_prob[valid_idx])
    auc_loocv <- auc(roc_loocv)
    
    list(
      auc = auc_loocv,
      roc_obj = roc_loocv,
      predictions = data.frame(
        true_class = loocv_true_class[valid_idx],
        predicted_prob = loocv_predicted_prob[valid_idx]
      )
    )
  }
  
  output$loocvTable <- renderDT({
    req(diff_expr(), data())
    
    # Get significant proteins
    sig_proteins <- filtered_diff_expr() %>%
      arrange(regression_coefficient) %>% pull(Protein)
    
    # Calculate LOOCV for each protein
    loocv_results_list <- lapply(sig_proteins, function(protein) {
      result <- tryCatch({
        loocv_result <- calc_loocv_auc(data(), protein)
        if(is.null(loocv_result)) return(NULL)
        
        data.frame(
          Protein = protein,
          AUC = round(loocv_result$auc, 3),
          Sensitivity = round(mean(loocv_result$roc_obj$sensitivities), 3),
          Specificity = round(mean(loocv_result$roc_obj$specificities), 3)
        )
      }, error = function(e) {
        message("Error calculating LOOCV for protein ", protein, ": ", e$message)
        NULL
      })
      return(result)
    })
    
    # Combine results and remove NULL entries
    loocv_results_df <- do.call(rbind, loocv_results_list[!sapply(loocv_results_list, is.null)])
    
    # Return results table
    datatable(loocv_results_df, 
              options = list(
                pageLength = 10,
                lengthMenu = c(10, 25, 50, 100),
                scrollX = TRUE
              ),
              rownames = FALSE)
  })

  output$downloadLOOCV <- downloadHandler(
    filename = function() {
      paste("loocv_results_all_", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      # Get significant proteins
      sig_proteins <- filtered_diff_expr() %>%
        arrange(regression_coefficient) %>% pull(Protein)
      
      # Calculate LOOCV for each protein
      loocv_results_list <- lapply(sig_proteins, function(protein) {
        result <- tryCatch({
          loocv_result <- calc_loocv_auc(data(), protein)
          if(is.null(loocv_result)) return(NULL)
          
          data.frame(
            Protein = protein,
            AUC = round(loocv_result$auc, 3),
            Sensitivity = round(mean(loocv_result$roc_obj$sensitivities), 3),
            Specificity = round(mean(loocv_result$roc_obj$specificities), 3)
          )
        }, error = function(e) {
          NULL
        })
        return(result)
      })
      
      # Combine results and remove NULL entries
      loocv_results_df <- do.call(rbind, loocv_results_list[!sapply(loocv_results_list, is.null)])
      
      write.csv(loocv_results_df, file, row.names = FALSE)
    }
)
  
  # Tooltip observers using shinyBS
  observe({
    bsTooltip("pcaTab", "PCA shows overall variance in the dataset and how samples cluster.", placement = "top", trigger = "hover")
    bsTooltip("volcanoTab", "Volcano plot shows the relationship between significance and fold change.", placement = "top", trigger = "hover")
    bsTooltip("heatmapTab", "Heatmap visualizes expression patterns across samples for top differentially expressed proteins.", placement = "top", trigger = "hover")
    bsTooltip("distplotTab", "Distribution plots show expression patterns for individual proteins.", placement = "top", trigger = "hover")
    bsTooltip("rocTab", "ROC curves evaluate classification performance of a protein.", placement = "top", trigger = "hover")
    bsTooltip("loocvTab", "LOOCV assesses classification robustness by leaving one sample out.", placement = "top", trigger = "hover")
  })
}

shinyApp(ui = ui, server = server)