# Install required packages if not already installed
required_packages <- c("dplyr", "ggplot2", "car", "effsize")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

# Load required libraries
library(car)      # for Levene's test
library(effsize)  # for effect size calculations
library(dplyr)
library(ggplot2)

statdata <- read.csv("C:/Users/samaw/Downloads/LEP_Data.csv",header=T)
summary(statdata)

# Statistical Analysis: Greek Speakers in Greece vs Abroad
# Comparing Control group vs Diaspora groups collectively
# Create binary grouping variable
statdata$Group_Binary <- ifelse(statdata$Group == "Ctrl", "Control", "Diaspora")

# Convert Group_Binary to factor
statdata$Group_Binary <- factor(statdata$Group_Binary, levels = c("Control", "Diaspora"))

# Function to perform appropriate statistical tests
perform_stat_test <- function(data, variable, group_var = "Group_Binary") {
  ##cat("\n", rep("=",60), "\n")
  #cat("ANALYSIS FOR:", variable, "\n")
  #cat(rep("=",60), "\n")
  
  # Remove missing values for this analysis
  clean_data <- data[!is.na(data[[variable]]), ]
  
  if(variable %in% c("Gender", "Birth")) {
    # categorical variables - Chi-square test
    #cat("Variable type: categorical\n")
    #cat("Test used: Chi-square test (or Fisher's exact test if expected frequencies < 5)\n\n")
    
    # Create contingency table
    cont_table <- table(clean_data[[variable]], clean_data[[group_var]])
    print("Contingency Table:")
    print(cont_table)
    print(addmargins(cont_table))
    
    # Check expected frequencies
    expected <- chisq.test(cont_table)$expected
    #cat("\nExpected frequencies:\n")
    print(expected)
    
    # Perform appropriate test
    if(any(expected < 5)) {
      #cat("\nUsing Fisher's Exact Test (some expected frequencies < 5)\n")
      test_result <- fisher.test(cont_table)
      print(test_result)
    } else {
      #cat("\nUsing Chi-square Test\n")
      test_result <- chisq.test(cont_table)
      print(test_result)
    }
    
  } else {
    # Continuous/ordinal variables
    #cat("Variable type: Continuous/Ordinal\n")
    
    # Descriptive statistics by group
    desc_stats <- clean_data %>%
      group_by(!!sym(group_var)) %>%
      summarise(
        n = n(),
        mean = mean(!!sym(variable), na.rm = TRUE),
        median = median(!!sym(variable), na.rm = TRUE),
        sd = sd(!!sym(variable), na.rm = TRUE),
        min = min(!!sym(variable), na.rm = TRUE),
        max = max(!!sym(variable), na.rm = TRUE),
        .groups = 'drop'
      )
    
    #cat("\nDescriptive Statistics:\n")
    print(desc_stats)
    
    # Shapiro-Wilk test for normality (if n <= 50 per group)
    control_data <- clean_data[clean_data[[group_var]] == "Control", variable]
    expat_data <- clean_data[clean_data[[group_var]] == "Diaspora", variable]
    
    #cat("\nNormality Tests (Shapiro-Wilk):\n")
    if(length(control_data) <= 50) {
      shapiro_control <- shapiro.test(control_data)
      #cat("Control group: W =", round(shapiro_control$statistic, 4), 
         # ", p =", round(shapiro_control$p.value, 4), "\n")
    }
    
    if(length(expat_data) <= 50) {
      shapiro_expat <- shapiro.test(expat_data)
      #cat("Diaspora group: W =", round(shapiro_expat$statistic, 4), 
         # ", p =", round(shapiro_expat$p.value, 4), "\n")
    }
    
    # Levene's test for equal variances
    if(require(car, quietly = TRUE)) {
      levene_test <- car::leveneTest(clean_data[[variable]], clean_data[[group_var]])
      #cat("\nLevene's Test for Equal Variances:\n")
      print(levene_test)
    }
    
    # Perform Mann-Whitney U test (recommended for non-normal data)
    #cat("\nMann-Whitney U Test (Wilcoxon rank-sum test):\n")
    #cat("(Recommended for non-normal distributions)\n")
    wilcox_result <- wilcox.test(clean_data[[variable]] ~ clean_data[[group_var]], 
                                 exact = FALSE, correct = TRUE)
    print(wilcox_result)
    
    # Also perform t-test for comparison
    #cat("\nWelch's t-test (for comparison, assumes normality):\n")
    t_result <- t.test(clean_data[[variable]] ~ clean_data[[group_var]], 
                       var.equal = FALSE)
    print(t_result)
    
    # Effect size (Cohen's d equivalent for Mann-Whitney)
    if(require(effsize, quietly = TRUE)) {
      cohen_d <- effsize::cohen.d(clean_data[[variable]], clean_data[[group_var]])
      #cat("\nEffect size (Cohen's d):\n")
      print(cohen_d)
    }
  }
}


# Variables to analyze
variables_to_test <- c("Age", "Gender", "Birth", "GenN", "YrsAbroad", "ExpGr", "ExpFr")

# Perform analysis for each variable
for(var in variables_to_test) {
  perform_stat_test(statdata, var)
}

# Summary table of p-values
#cat("\n", rep("=",80), "\n")
#cat("SUMMARY OF P-VALUES\n")
#cat(rep("=",80), "\n")

summary_results <- data.frame(
  Variable = character(),
  Test_Used = character(),
  P_Value = numeric(),
  Significant = character(),
  stringsAsFactors = FALSE
)

# You'll need to manually fill this based on the results above
#cat("After running the analysis, create a summary table with:\n")
#cat("- Variable name\n")
#cat("- Test used (Chi-square/Fisher's/Mann-Whitney/t-test)\n")
#cat("- P-value\n")
#cat("- Significance (p < 0.05)\n")

# Visualization function
create_visualizations <- function(data, variables) {
  for(var in variables) {
    if(var %in% c("Gender", "Birth")) {
      # Bar plot for categorical variables
      p <- ggplot(data, aes(x = !!sym(var), fill = Group_Binary)) +
        geom_bar(position = "dodge") +
        labs(title = paste("Distribution of", var, "by Group"),
             x = var, y = "Count", fill = "Group") +
        theme_minimal()
      print(p)
    } else {
      # Box plot for continuous variables
      p <- ggplot(data, aes(x = Group_Binary, y = !!sym(var), fill = Group_Binary)) +
        geom_boxplot(alpha = 0.7) +
        geom_jitter(width = 0.2, alpha = 0.5) +
        labs(title = paste("Distribution of", var, "by Group"),
             x = "Group", y = var) +
        theme_minimal() +
        theme(legend.position = "none")
      print(p)
    }
  }
}

# Create visualizations
#cat("\nCreating visualizations...\n")
create_visualizations(statdata, variables_to_test)

# Education Level Analysis
analyze_education <- function(data) {
  #cat(rep("=", 80), "\n")
  #cat("EDUcatION LEVEL ANALYSIS\n")
  #cat(rep("=", 80), "\n\n")
  
  # Remove missing values
  clean_data <- data[!is.na(data$Edu), ]
  #cat("Sample size after removing missing values:", nrow(clean_data), "\n")
  
  # Create contingency table
  #cat("\n1. CONTINGENCY TABLE\n")
  #cat(rep("-", 40), "\n")
  cont_table <- table(clean_data$Edu, clean_data$Group_Binary)
  print("Raw Contingency Table:")
  print(cont_table)
  
  #cat("\nContingency Table with Margins:\n")
  print(addmargins(cont_table))
  
  # Proportions within each group
  #cat("\n2. PROPORTIONS WITHIN GROUPS\n")
  #cat(rep("-", 40), "\n")
  prop_table <- prop.table(cont_table, margin = 2) # proportions by column (group)
  #cat("Proportions within each group:\n")
  print(round(prop_table, 3))
  
  # Overall proportions
  #cat("\nOverall proportions:\n")
  overall_prop <- prop.table(cont_table)
  print(round(overall_prop, 3))
  
  # Expected frequencies for chi-square test
  #cat("\n3. EXPECTED FREQUENCIES\n")
  #cat(rep("-", 40), "\n")
  expected <- chisq.test(cont_table)$expected
  print(round(expected, 2))
  
  # Check if chi-square test assumptions are met
  min_expected <- min(expected)
  cells_less_than_5 <- sum(expected < 5)
  total_cells <- length(expected)
  
  #cat("\nChi-square test assumptions:\n")
  #cat("Minimum expected frequency:", round(min_expected, 2), "\n")
  #cat("Cells with expected frequency < 5:", cells_less_than_5, "out of", total_cells, "\n")
  #cat("Percentage of cells < 5:", round((cells_less_than_5/total_cells)*100, 1), "%\n")
  
  # Perform appropriate statistical test
  #cat("\n4. STATISTICAL TESTS\n")
  #cat(rep("-", 40), "\n")
  
  if(cells_less_than_5/total_cells > 0.2) {
    cat("WARNING: >20% of cells have expected frequency < 5\n")
    cat("Using Fisher's Exact Test (Monte Carlo simulation)\n\n")
    
    # Fisher's exact test with simulation for larger tables
    fisher_result <- fisher.test(cont_table, simulate.p.value = TRUE, B = 10000)
    print(fisher_result)
    
  } else {
    cat("Using Pearson's Chi-square Test\n\n")
    chi_result <- chisq.test(cont_table)
    print(chi_result)
    
    # Also provide Fisher's test for comparison
    cat("\nFisher's Exact Test (for comparison):\n")
    fisher_result <- fisher.test(cont_table, simulate.p.value = TRUE, B = 10000)
    print(fisher_result)
  }
  
  # Effect size (Cramér's V)
  #cat("\n5. EFFECT SIZE\n")
  #cat(rep("-", 40), "\n")
  chi_stat <- chisq.test(cont_table)$statistic
  n <- sum(cont_table)
  df_min <- min(nrow(cont_table) - 1, ncol(cont_table) - 1)
  cramers_v <- sqrt(chi_stat / (n * df_min))
  
  #cat("Cramér's V:", round(cramers_v, 3), "\n")
  #cat("Effect size interpretation:\n")
  #cat("  Small: 0.1, Medium: 0.3, Large: 0.5\n")
  
  if(cramers_v < 0.1) {
    effect_interp <- "negligible"
  } else if(cramers_v < 0.3) {
    effect_interp <- "small"
  } else if(cramers_v < 0.5) {
    effect_interp <- "medium"
  } else {
    effect_interp <- "large"
  }
  #cat("  Current effect size:", effect_interp, "\n")
  
  # Detailed breakdown by education level
  #cat("\n6. DETAILED BREAKDOWN BY EDUcatION LEVEL\n")
  #cat(rep("-", 50), "\n")
  
  detailed_stats <- clean_data %>%
    group_by(Edu, Group_Binary) %>%
    summarise(count = n(), .groups = 'keep') %>%
    group_by(Edu) %>%
    mutate(
      total_in_edu_level = sum(count),
      percentage_in_group = round((count / total_in_edu_level) * 100, 1)
    ) %>%
    ungroup()
  
  print(detailed_stats)
  
  # Post-hoc analysis: pairwise comparisons if overall test is significant
  if(exists("chi_result") && chi_result$p.value < 0.05) {
    #cat("\n7. POST-HOC ANALYSIS\n")
    #cat(rep("-", 40), "\n")
    #cat("Overall test is significant. Consider examining standardized residuals:\n\n")
    
    # Standardized residuals
    std_residuals <- chisq.test(cont_table)$stdres
    #cat("Standardized residuals (|z| > 2 suggests significant contribution):\n")
    print(round(std_residuals, 2))
    
    # Highlight significant cells
    significant_cells <- which(abs(std_residuals) > 2, arr.ind = TRUE)
    if(nrow(significant_cells) > 0) {
      #cat("\nSignificant deviations from expected (|z| > 2):\n")
      for(i in 1:nrow(significant_cells)) {
        row_idx <- significant_cells[i, 1]
        col_idx <- significant_cells[i, 2]
        edu_level <- rownames(std_residuals)[row_idx]
        group <- colnames(std_residuals)[col_idx]
        residual <- std_residuals[row_idx, col_idx]
        direction <- ifelse(residual > 0, "more than expected", "fewer than expected")
        #cat(paste0("  ", edu_level, " in ", group, " group: ", direction, 
                 #  " (z = ", round(residual, 2), ")\n"))
      }
    }
  }
  
  return(list(
    contingency_table = cont_table,
    proportions = prop_table,
    test_result = if(exists("chi_result")) chi_result else fisher_result,
    cramers_v = cramers_v,
    effect_size = effect_interp
  ))
}

# Create visualizations for education
create_education_plots <- function(data) {
  # Remove missing values
  clean_data <- data[!is.na(data$Edu), ]
  
  # Define education level order (assuming this is the intended hierarchy)
  edu_order <- c("LYC", "TEC", "STU", "BAC", "MAS", "PHD")
  clean_data$Edu <- factor(clean_data$Edu, levels = edu_order)
  
  #cat("\n8. VISUALIZATIONS\n")
  #cat(rep("-", 40), "\n")
  
  
  # 1. Grouped bar chart
  p2 <- ggplot(clean_data, aes(x = Edu, fill = Group_Binary)) +
    geom_bar(position = "dodge") +
    labs(title = "Education Levels by Group (Grouped)",
         subtitle = paste("Total sample size:", nrow(clean_data)),
         x = "Education Level", 
         y = "Count",
         fill = "Group") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_manual(values = c("Control" = "#2E86AB", "Expatriate" = "#A23B72"))
  
  print(p2)
  
  
  # 2. Group comparison plot
  p4 <- ggplot(clean_data, aes(x = Group_Binary, fill = Edu)) +
    geom_bar(position = "fill") +
    labs(title = "Education Distribution in the Control and Diaspora Groups",
         subtitle = "Shows the education composition of each group",
         x = "Group", 
         y = "Proportion",
         fill = "Education Level") +
    theme_minimal() +
    scale_y_continuous(labels = scales::percent_format()) +
    scale_fill_brewer(palette = "Set3")
  
  print(p4)
  
  #cat("\nEducation level codes (assumed hierarchy):\n")
  #cat("LYC = Lyceum (High School)\n")
  #cat("TEC = Technical Education\n") 
  #cat("STU = Student (Current)\n")
  #cat("BAC = Bachelor's Degree\n")
  #cat("MAS = Master's Degree\n")
  #cat("PHD = Doctoral Degree\n")
}

# Run the analysis
#cat("Starting Education Level Analysis...\n\n")
edu_results <- analyze_education(statdata)

# Create visualizations
create_education_plots(statdata)

# Summary
#cat("\n", rep("=", 80), "\n")
#cat("SUMMARY\n")
#cat(rep("=", 80), "\n")
#cat("Education level analysis completed.\n")
#cat("Check the p-value from the statistical test above to determine significance.\n")
#cat("Cramér's V effect size:", round(edu_results$cramers_v, 3), "(", edu_results$effect_size, "effect )\n")
#cat("Examine standardized residuals and visualizations for interpretation.\n")






##########################





SpontData <- read.csv("C:/Users/samaw/Downloads/Spont_Data.csv")
summary(SpontData)

# SpontData Analysis: Correlations and Group Differences
# Analyzing relationships between independent and dependent variables
# with multiple comparison corrections

library(corrplot)
library(Hmisc)
library(psych)

# Install required packages if not available
required_packages <- c("dplyr", "ggplot2", "corrplot", "Hmisc", "psych", "car", "effsize")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)


# Define variable groups
independent_vars <- c("Group", "Age", "Gender", "Edu", "Yrs_Abroad", "GenN", "Exp_GR", "Exp_FR", "InDeg")
dependent_vars <- c("Subject_Drop", "Opp_Drop_Rate", "Length_Drop_Rate", 
                    "Overt_Subject_Pro", "Opp_S_Prn_Ratio", "Length_SPrnR")
continuous_vars <- c("Age", "Yrs_Abroad", "GenN", "Exp_GR", "Exp_FR", "InDeg", dependent_vars)
categorical_vars <- c("Group", "Gender", "Edu")

#cat(rep("=", 80), "\n")
#cat("SPONTDATA COMPREHENSIVE ANALYSIS\n")
#cat(rep("=", 80), "\n\n")

#cat("Sample size:", nrow(SpontData), "\n")
#cat("Independent variables:", length(independent_vars), "\n")
#cat("Dependent variables:", length(dependent_vars), "\n\n")

# 1. CORRELATION MATRIX ANALYSIS
#cat("1. CORRELATION MATRIX ANALYSIS\n")
#cat(rep("-", 50), "\n")

# Focus only on independent vs dependent variable correlations
# Remove dependent-dependent correlations from analysis
independent_continuous <- c("Age", "Yrs_Abroad", "GenN", "Exp_GR", "Exp_FR", "InDeg")
all_analysis_vars <- c(independent_continuous, dependent_vars)

# Remove rows with missing values for correlation analysis
analysis_data <- SpontData[, all_analysis_vars]
complete_analysis <- analysis_data[complete.cases(analysis_data), ]

#cat("Complete cases for correlation analysis:", nrow(complete_analysis), "\n")
#cat("Analyzing correlations between independent and dependent variables only\n\n")

# Spearman correlation (non-parametric)
#cat("Using Spearman correlation (robust to non-normal distributions)\n")
spear_corr <- cor(complete_analysis, method = "spearman", use = "complete.obs")
spear_cor_test <- rcorr(as.matrix(complete_analysis), type = "spearman")

# Extract only independent-dependent correlations
ind_dep_correlations <- spear_corr[independent_continuous, dependent_vars]
ind_dep_p_values <- spear_cor_test$P[independent_continuous, dependent_vars]

#cat("Independent vs Dependent Variable Correlations:\n")
print(round(ind_dep_correlations, 3))

#cat("\nP-values for Independent vs Dependent correlations:\n")
print(round(ind_dep_p_values, 3))

# Apply Benjamini-Hochberg correction only to independent-dependent correlations
p_values_vector <- as.vector(ind_dep_p_values)
p_values_vector <- p_values_vector[!is.na(p_values_vector)]
n_comparisons <- length(p_values_vector)
adjusted_p <- p.adjust(p_values_vector, method = "BH")

#cat("\nNumber of correlation tests performed (independent vs dependent only):", n_comparisons, "\n")
#cat("Multiple comparison correction: Benjamini-Hochberg (FDR)\n")

# Create adjusted p-value matrix
adjusted_p_matrix <- matrix(adjusted_p, nrow = length(independent_continuous), ncol = length(dependent_vars))
rownames(adjusted_p_matrix) <- independent_continuous
colnames(adjusted_p_matrix) <- dependent_vars

#cat("\nAdjusted p-values (Benjamini-Hochberg) - Independent vs Dependent only:\n")
print(round(adjusted_p_matrix, 3))

# Identify significant correlations after correction
significant_cors <- which(adjusted_p_matrix < 0.05, arr.ind = TRUE)
if(nrow(significant_cors) > 0) {
  cat("\nSignificant correlations after multiple comparison correction (p < 0.05):\n")
  for(i in 1:nrow(significant_cors)) {
    row_idx <- significant_cors[i, 1]
    col_idx <- significant_cors[i, 2]
    var1 <- rownames(adjusted_p_matrix)[row_idx]
    var2 <- colnames(adjusted_p_matrix)[col_idx]
    cor_val <- ind_dep_correlations[row_idx, col_idx]
    p_val <- adjusted_p_matrix[row_idx, col_idx]
    cat(sprintf("  %s - %s: r = %.3f, adjusted p = %.3f\n", var1, var2, cor_val, p_val))
  }
} else {
  cat("\nNo significant correlations after multiple comparison correction.\n")
}

# Visualization of correlation matrix
#cat("\nCreating correlation matrix visualization...\n")

# Create a focused correlation matrix showing only independent vs dependent relationships
# This will be a rectangular matrix (not square)
corrplot(ind_dep_correlations, 
         method = "color",
         tl.cex = 0.8,
         tl.col = "black",
         tl.srt = 45,
         title = "Independent vs Dependent Variable Correlations",
         mar = c(0,0,2,0),
         is.corr = FALSE)  # Important: tells corrplot this isn't a square correlation matrix

# Add significance indicators
sig_indicator <- ifelse(adjusted_p_matrix < 0.05, "*", 
                        ifelse(adjusted_p_matrix < 0.1, ".", ""))

# Create a plot with correlation values and significance
corrplot(ind_dep_correlations, 
         method = "color",
         tl.cex = 0.8,
         tl.col = "black",
         tl.srt = 45,
         title = "Independent vs Dependent Correlations with Significance",
         mar = c(0,0,2,0),
         addCoef.col = "black",
         number.cex = 0.7,
         is.corr = FALSE)

cat("Significance levels: * p < 0.05 (after correction), . p < 0.1 (after correction)\n")

# 2. GROUP DIFFERENCES ANALYSIS
cat("\n\n2. GROUP DIFFERENCES ANALYSIS\n")
cat(rep("-", 50), "\n")

analyze_group_differences <- function(data, dependent_variables) {
  results_list <- list()
  p_values_for_correction <- c()
  test_names <- c()
  
  for(dv in dependent_variables) {
    cat("\n", rep("-", 30), "\n")
    cat("DEPENDENT VARIABLE:", dv, "\n")
    cat(rep("-", 30), "\n")
    
    # Check if variable exists and is numeric
    if(!dv %in% names(data)) {
      cat("Variable", dv, "not found in dataset. Skipping...\n")
      next
    }
    
    if(!is.numeric(data[[dv]])) {
      cat("Variable", dv, "is not numeric. Converting to numeric...\n")
      data[[dv]] <- as.numeric(as.character(data[[dv]]))
    }
    
    # Remove missing values and ensure Group variable exists
    if(!"Group" %in% names(data)) {
      cat("Group variable not found. Skipping analysis.\n")
      return(NULL)
    }
    
    clean_data <- data[!is.na(data[[dv]]) & !is.na(data$Group), ]
    
    # Check if we have enough data
    if(nrow(clean_data) < 4) {
      cat("Insufficient data for", dv, ". Skipping...\n")
      next
    }
    
    # Check if we have at least 2 groups
    if(length(unique(clean_data$Group)) < 2) {
      cat("Less than 2 groups found for", dv, ". Skipping...\n")
      next
    }
    
    # Descriptive statistics by group
    desc_stats <- clean_data %>%
      group_by(Group) %>%
      summarise(
        n = n(),
        mean = mean(!!sym(dv), na.rm = TRUE),
        median = median(!!sym(dv), na.rm = TRUE),
        sd = sd(!!sym(dv), na.rm = TRUE),
        min = min(!!sym(dv), na.rm = TRUE),
        max = max(!!sym(dv), na.rm = TRUE),
        IQR = IQR(!!sym(dv), na.rm = TRUE),
        .groups = 'drop'
      )
    
    cat("Descriptive Statistics by Group:\n")
    print(desc_stats)
    
    # Shapiro-Wilk tests for normality by group
    cat("\nNormality Tests (Shapiro-Wilk):\n")
    groups <- unique(clean_data$Group)
    for(group in groups) {
      group_data <- clean_data[clean_data$Group == group, dv]
      group_data <- group_data[!is.na(group_data)]  # Remove NAs
      if(length(group_data) >= 3 && length(group_data) <= 5000) {
        tryCatch({
          shapiro_result <- shapiro.test(group_data)
          cat(sprintf("%s group: W = %.4f, p = %.4f\n", 
                      group, shapiro_result$statistic, shapiro_result$p.value))
        }, error = function(e) {
          cat(sprintf("%s group: Could not perform Shapiro test\n", group))
        })
      } else {
        cat(sprintf("%s group: Sample size (%d) not suitable for Shapiro test\n", 
                    group, length(group_data)))
      }
    }
    
    # Mann-Whitney U test (recommended for non-normal data)
    cat("\nMann-Whitney U Test:\n")
    tryCatch({
      wilcox_result <- wilcox.test(as.numeric(clean_data[[dv]]) ~ factor(clean_data$Group), 
                                   exact = FALSE, correct = TRUE)
      print(wilcox_result)
      
      # Store p-value for multiple comparison correction
      p_values_for_correction <- c(p_values_for_correction, wilcox_result$p.value)
      test_names <- c(test_names, paste("Group difference:", dv))
      
    }, error = function(e) {
      cat("Error in Mann-Whitney U test for", dv, ":", e$message, "\n")
      next
    })
    
    # Effect size calculation - FIXED VERSION
    cat("\nEffect size calculation:\n")
    tryCatch({
      # Manual Cohen's d calculation (most reliable)
      group_names <- unique(clean_data$Group)
      if(length(group_names) == 2) {
        group1_data <- as.numeric(clean_data[clean_data$Group == group_names[1], dv])
        group2_data <- as.numeric(clean_data[clean_data$Group == group_names[2], dv])
        
        # Remove any remaining NAs
        group1_data <- group1_data[!is.na(group1_data)]
        group2_data <- group2_data[!is.na(group2_data)]
        
        if(length(group1_data) > 0 && length(group2_data) > 0) {
          mean1 <- mean(group1_data)
          mean2 <- mean(group2_data)
          sd1 <- sd(group1_data)
          sd2 <- sd(group2_data)
          n1 <- length(group1_data)
          n2 <- length(group2_data)
          
          # Pooled standard deviation
          pooled_sd <- sqrt(((n1-1)*sd1^2 + (n2-1)*sd2^2)/(n1+n2-2))
          
          # Cohen's d
          cohens_d <- (mean1 - mean2) / pooled_sd
          
          cat(sprintf("Cohen's d = %.3f\n", cohens_d))
          cat("Effect size interpretation: ")
          if(abs(cohens_d) < 0.2) cat("Negligible effect\n")
          else if(abs(cohens_d) < 0.5) cat("Small effect\n")
          else if(abs(cohens_d) < 0.8) cat("Medium effect\n")
          else cat("Large effect\n")
          
          # Store results
          results_list[[dv]] <- list(
            descriptives = desc_stats,
            test_result = wilcox_result,
            effect_size = list(estimate = cohens_d, 
                               interpretation = ifelse(abs(cohens_d) < 0.2, "negligible",
                                                       ifelse(abs(cohens_d) < 0.5, "small",
                                                              ifelse(abs(cohens_d) < 0.8, "medium", "large"))))
          )
        }
      }
    }, error = function(e) {
      cat("Could not calculate effect size for", dv, "- Error:", e$message, "\n")
      # Store results without effect size
      results_list[[dv]] <- list(
        descriptives = desc_stats,
        test_result = wilcox_result,
        effect_size = NULL
      )
    })
    
    # Also perform t-test for comparison
    cat("\nWelch's t-test (for comparison):\n")
    tryCatch({
      t_result <- t.test(as.numeric(clean_data[[dv]]) ~ factor(clean_data$Group), 
                         var.equal = FALSE)
      print(t_result)
    }, error = function(e) {
      cat("Error in t-test for", dv, ":", e$message, "\n")
    })
  }
  
  # Apply multiple comparison correction
  if(length(p_values_for_correction) > 0) {
    cat("\n\n", rep("=", 60), "\n")
    cat("MULTIPLE COMPARISON CORRECTION\n")
    cat(rep("=", 60), "\n")
    
    cat("Number of tests performed:", length(p_values_for_correction), "\n")
    cat("Method: Benjamini-Hochberg (FDR control)\n\n")
    
    adjusted_p_values <- p.adjust(p_values_for_correction, method = "BH")
    
    # Create summary table
    summary_table <- data.frame(
      Test = test_names,
      Raw_P_Value = round(p_values_for_correction, 4),
      Adjusted_P_Value = round(adjusted_p_values, 4),
      Significant = ifelse(adjusted_p_values < 0.05, "Yes", "No"),
      stringsAsFactors = FALSE
    )
    
    cat("Summary of Group Difference Tests:\n")
    print(summary_table)
    
    # Highlight significant results
    significant_tests <- summary_table[summary_table$Significant == "Yes", ]
    if(nrow(significant_tests) > 0) {
      cat("\nSignificant group differences after correction:\n")
      for(i in 1:nrow(significant_tests)) {
        cat(sprintf("  %s (adjusted p = %.4f)\n", 
                    significant_tests$Test[i], significant_tests$Adjusted_P_Value[i]))
      }
    } else {
      cat("\nNo significant group differences after multiple comparison correction.\n")
    }
    
    return(list(results = results_list, summary = summary_table))
  } else {
    cat("No valid tests were performed.\n")
    return(list(results = results_list, summary = data.frame()))
  }
}

# Now run the analysis with the corrected function
cat("Running group differences analysis with corrected function...\n")
group_analysis <- analyze_group_differences(SpontData, dependent_vars)

# 3. catEGORICAL VARIABLE ANALYSIS
#cat("\n\n3. catEGORICAL VARIABLE ANALYSIS\n")
#cat(rep("-", 50), "\n")

analyze_categorical_effects <- function(data, categorical_var, dependent_variables) {
  #cat("Analyzing effects of", categorical_var, "on dependent variables\n")
  #cat(rep("-", 40), "\n")
  
  p_values <- c()
  test_names <- c()
  
  for(dv in dependent_variables) {
    #cat("\nDV:", dv, "vs", categorical_var, "\n")
    
    # Remove missing values
    clean_data <- data[!is.na(data[[dv]]) & !is.na(data[[categorical_var]]), ]
    
    if(categorical_var == "Gender" || categorical_var == "Group") {
      # Two groups - Mann-Whitney U test
      if(length(unique(clean_data[[categorical_var]])) == 2) {
        test_result <- wilcox.test(clean_data[[dv]] ~ clean_data[[categorical_var]])
        #cat(sprintf("  Mann-Whitney U: p = %.4f\n", test_result$p.value))
        p_values <- c(p_values, test_result$p.value)
        test_names <- c(test_names, paste(categorical_var, "effect on", dv))
      }
    } else if(categorical_var == "Edu") {
      # Multiple groups - Kruskal-Wallis test
      test_result <- kruskal.test(clean_data[[dv]] ~ factor(clean_data[[categorical_var]]))
      #cat(sprintf("  Kruskal-Wallis: p = %.4f\n", test_result$p.value))
      p_values <- c(p_values, test_result$p.value)
      test_names <- c(test_names, paste(categorical_var, "effect on", dv))
    }
  }
  
  # Apply multiple comparison correction
  if(length(p_values) > 0) {
    adjusted_p <- p.adjust(p_values, method = "BH")
    
    #cat("\nMultiple comparison correction results:\n")
    correction_table <- data.frame(
      Test = test_names,
      Raw_P = round(p_values, 4),
      Adjusted_P = round(adjusted_p, 4),
      Significant = ifelse(adjusted_p < 0.05, "Yes", "No")
    )
    print(correction_table)
    
    return(correction_table)
  }
}

# Analyze categorical variables
if("Gender" %in% names(SpontData)) {
  gender_analysis <- analyze_categorical_effects(SpontData, "Gender", dependent_vars)
}

if("Edu" %in% names(SpontData)) {
  edu_analysis <- analyze_categorical_effects(SpontData, "Edu", dependent_vars)
}

# 4. VISUALIZATION OF KEY RELATIONSHIPS
#cat("\n\n4. VISUALIZATION OF KEY RELATIONSHIPS\n")
#cat(rep("-", 50), "\n")

# Box plots for group differences
create_group_boxplots <- function(data, dependent_variables) {
  for(dv in dependent_variables) {
    p <- ggplot(data, aes(x = Group, y = !!sym(dv), fill = Group)) +
      geom_boxplot(alpha = 0.7) +
      geom_jitter(width = 0.2, alpha = 0.6) +
      labs(title = paste("Group Differences in", dv),
           x = "Group", y = dv) +
      theme_minimal() +
      theme(legend.position = "none") +
      scale_fill_manual(values = c("Control" = "#2E86AB", "Francophone" = "#A23B72"))
    
    print(p)
  }
}

#cat("Creating boxplots for group comparisons...\n")
create_group_boxplots(SpontData, dependent_vars)

# Correlation heatmap focusing on significant relationships
#cat("\nCreating focused independent-dependent correlation heatmap...\n")

# Create visualization specifically for independent vs dependent relationships
if(nrow(complete_analysis) > 0) {
  # Create a comprehensive plot showing all independent-dependent correlations
  corrplot(ind_dep_correlations, 
           method = "color",
           tl.cex = 0.7,
           tl.col = "black",
           tl.srt = 45,
           title = "All Independent vs Dependent Variable Correlations",
           mar = c(0,0,2,0),
           addCoef.col = "black",
           number.cex = 0.6,
           is.corr = FALSE,
           col = colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))(200))
} else {
  cat("Insufficient complete cases for correlation matrix.\n")
}

# 5. SUMMARY AND RECOMMENDATIONS
cat("\n\n", rep("=", 80), "\n")
cat("ANALYSIS SUMMARY AND RECOMMENDATIONS\n")
cat(rep("=", 80), "\n")

cat("1. All correlation and group difference tests used non-parametric methods\n")
cat("   (Spearman correlations, Mann-Whitney U, Kruskal-Wallis)\n\n")

cat("2. Multiple comparison corrections applied using Benjamini-Hochberg method\n")
cat("   to control False Discovery Rate (FDR)\n\n")

cat("3. Sample size: n =", nrow(SpontData), "(10 per group)\n")
cat("   Note: Small sample size may limit power to detect effects\n\n")

cat("4. Check the results above for:\n")
cat("   - Significant correlations after correction\n")
cat("   - Significant group differences after correction\n")
cat("   - Effect sizes for practical significance\n\n")

cat("5. Visualizations provided for:\n")
cat("   - Correlation matrices\n")
cat("   - Group comparison boxplots\n")
cat("   - Key variable relationships\n\n")

cat("Analysis completed successfully!\n")