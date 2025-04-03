# Setup -------------------------------------------------------------------
setwd("/Users/gastaldi/Documents/ownCloud/GitHub/MMQT-longDTA/longDTA_all_MMQT_txt")

### Load necessary libraries
set.seed(42)
library(dplyr)
library(openxlsx)
library(MASS)
library(factoextra)
library(readxl)
library(tidyverse)
library(writexl)
library(lme4)
library(lmerTest)
library(broom)
library(ggstatsplot)
library(palmerpenguins)
library(psych)
library(reshape2)
library(effsize)
library(emmeans)

# Merge individual txt files and first quality control --------------------
# Get a list of all .txt files in the working directory
files <- list.files(pattern = "*.txt")
# Initialize an empty data frame for the merged data
merged_data <- data.frame()
# Loop through the files
for (file in files) {
  # Read the file
  data <- read.table(file, header = TRUE, sep = "\t")
  # Add the ImageID column
  data$ImageID <- substr(file, 1, 11)
  # Add the data to the merged_data data frame
  merged_data <- bind_rows(merged_data, data)
}
### Add Info from manual quality control 
setwd("/Users/gastaldi/Documents/ownCloud/GitHub/MMQT-longDTA/longDTA_MMQT_Analysis")
longDTA_MMQT_QC <- read.xlsx("longDTA_MMQT_QC.xlsx")
# Merge the two data frames on the "ImageID" column
merged_data <- merge(merged_data, longDTA_MMQT_QC, by = "ImageID")
# Write the merged data to a new Excel file
write.xlsx(merged_data, "merged_data.xlsx")

### filter cells by Quality=ok, mergedSoma =0,nucleiN = 1,thrDistCells 15
# stackDistM_1, 2, 4, 5 are XY; stackDistM_3 and 6 are Z 
# groupShortestDist > 15, branchN > 0

filtered_data <- merged_data %>%
  filter(Quality == "ok", nucleiN == 1, mergedSoma == 0, stackDistM_1 > 15, stackDistM_2 > 15, stackDistM_4 > 15, stackDistM_5 > 15, stackDistM_3 > 5, stackDistM_6 > 5, groupShortestDist > 15, branchN > 0)
# Write the merged data to a new Excel file
write.xlsx(filtered_data, "filtered_data.xlsx")

# Parameter selection -----------------------------------------------------
setwd("longDTA_MMQT_Analysis")
selected_data <- filtered_data %>%
  dplyr::select(sphericity,
                circularity,
                branchN,
                solidity,
                nodesEndingRatio,
                branchNodesN,
                branchEndNodesPerBranch,
                branchNodesPerBranch,
                branchLengthSkel_4,
                branchSegmentsPerBranch,
                betweenness_4,
                closeness_4,
                closenessSeed,
                branchLengthRatio_4,
                Group,
                Timepoint,
                AnimalID,
                ImageID,
                soma)
write.xlsx(selected_data, "selected_data.xlsx")

# Correlation Matrix ------------------------------------------------------
setwd("longDTA_MMQT_Analysis\\")


numeric_df <- filtered_data %>% select_if(is.numeric)
# Identify constant columns
constant_columns <- sapply(numeric_df, function(col) length(unique(col)) == 1)
# Remove constant columns
numeric_df <- numeric_df[, !constant_columns]
# Remove NA
numeric_df <- na.omit(numeric_df)
# Compute correlation matrix
cor_matrix <- cor(numeric_df)
# Write the merged data to a new Excel file
cor_df <- as.data.frame(cor_matrix)
write.xlsx(cor_df, "cor_matrix.xlsx")

# PCA Analysis ------------------------------------------------------------
setwd("longDTA_MMQT_Analysis")

#Load filered dataset
filtered_data <- read_excel("filtered_data.xlsx")
filtered_data <- filtered_data %>% mutate(CellID = paste(AnimalID, soma, sep="_"))

# Branching PCA
branching_data <- filtered_data %>%
  dplyr::select(nodesEndingRatio, 
                branchN, 
                branchNodesN, 
                branchNodesPerBranch, 
                branchEndNodesPerBranch, 
                branchSegmentsPerBranch, 
                branchLengthSkel_4, 
                branchLengthRatio_4, 
                closenessSeed, 
                closeness_4, 
                betweenness_4,
                soma, 
                CellID,
                AnimalID,
                Group,
                Timepoint)

# Select only the numerical columns for PCA
branching_PCA <- branching_data %>%
  select_if(is.numeric) %>% prcomp(scale. = TRUE) 

# Write the data to a new Excel file
branching_PCA_xlsx <- cbind(branching_data, branching_PCA$x)
write.xlsx(branching_PCA_xlsx, "branching_PCA.xlsx") 

# Select only the numerical columns for PCA
numerical_data <- selected_data %>%
  select_if(is.numeric)
# Run PCA
pca_result <- prcomp(numerical_data, scale. = TRUE)
# Get the principal component values
principal_components <- pca_result$x
# Add the principal component values to the selected_data
selected_data_with_PCA <- cbind(selected_data, principal_components)

# Write the data to a new Excel file
write.xlsx(selected_data_with_PCA, "selected_data_with_PCs.xlsx")
#Load PCA dataset
selected_data_with_PCA <- read_excel("selected_data_with_PCs.xlsx")
# Plot PCA
fviz_eig(pca_result, addlabels = T)
fviz_pca_biplot(branching_PCA, 
                label = "var", 
                habillage = filtered_data$Group,
                col.var = "black",
                )

# Composite Score ---------------------------------------------------------
setwd("longDTA_MMQT_Analysis")

selected_data <- read_excel("selected_data.xlsx")
selected_data_with_morphology_score <- selected_data
selected_data_with_morphology_score$Timepoint <- as.factor(selected_data_with_morphology_score$Timepoint)
selected_data_with_morphology_score$soma <- as.factor(selected_data_with_morphology_score$soma)
selected_data_with_morphology_score <- selected_data_with_morphology_score %>% mutate(CellID = paste(AnimalID, soma, sep="_"))

# Cronbachs alpha Morphology Score
morphology_data <- selected_data_with_morphology_score %>%
  dplyr::select(sphericity, 
                circularity,
                solidity,
                nodesEndingRatio, 
                branchN, 
                branchNodesN, 
                branchNodesPerBranch, 
                branchEndNodesPerBranch, 
                branchSegmentsPerBranch, 
                branchLengthSkel_4, 
                branchLengthRatio_4, 
                closenessSeed, 
                closeness_4, 
                betweenness_4)
morphology_transformed <- morphology_data %>% select_if(is.numeric)  %>% scale() 
positive_correl_morphology <- c("sphericity", "circularity", "solidity", "nodesEndingRatio", "closenessSeed", "closeness_4")
negative_correl_morphology <- c("branchN", "betweenness_4", "branchNodesN", "branchNodesPerBranch", "branchEndNodesPerBranch", "branchSegmentsPerBranch", "branchLengthSkel_4", "branchLengthRatio_4")
morphology_Cronbach <- alpha(morphology_transformed, keys = negative_correl_morphology)
write_xlsx(morphology_Cronbach$total, "morphology_Cronbach.xlsx")
# Great Cronbachs alpha = 0.91, create Morphology composite score
morphology_transformed <- morphology_data %>% select_if(is.numeric)  %>% scale() %>% as.data.frame() 
morphology_transformed[negative_correl_morphology] <- -morphology_transformed[negative_correl_morphology]
selected_data_with_morphology_score$Morphology_score <- rowMeans(morphology_transformed, na.rm = TRUE)
write_xlsx(selected_data_with_morphology_score, "selected_data_with_morphology_score.xlsx")



# Remove Outliers ---------------------------------------------------------
#Input = selected_data_with_morphology_score; Output = stats_df_clean
setwd("longDTA_MMQT_Analysis")
selected_data_with_morphology_score <- read_excel("selected_data_with_morphology_score.xlsx")
stats_df <- as.data.frame(selected_data_with_morphology_score)
stats_df$GroupByTime <- paste(stats_df$Group, stats_df$Timepoint, sep = "_")
stats_df$GroupByTime <- as.factor(stats_df$GroupByTime)

stats_df[] <- lapply(stats_df, function(x) {
  if(is.numeric(x)) return(x)
  factor(x)
})
find_outliers <- function(x) {
  qnt <- quantile(x, probs=c(.25, .75), na.rm = TRUE)
  H <- 1.5 * IQR(x, na.rm = TRUE)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  return(y)
}
stats_df_split <- split(stats_df, stats_df$GroupByTime)
stats_df_clean <- lapply(stats_df_split, function(df) {
  numeric_columns <- sapply(df, is.numeric)
  df[numeric_columns] <- lapply(df[numeric_columns], find_outliers)
  return(df)
})
stats_df_clean <- do.call(rbind, stats_df_clean)
stats_df_clean <- na.omit(stats_df_clean)
write_xlsx(stats_df_clean, "stats_df_clean.xlsx")

# Export Figures individual -----------------------------------------------
setwd("longDTA_MMQT_Analysis")
stats_df_clean <- read_excel("stats_df_clean.xlsx")

plotting_df <- as.data.frame(stats_df_clean)
plotting_df[] <- lapply(plotting_df, function(x) {
  if(is.numeric(x)) return(x)
  factor(x)
})

# some plotting parameters
group_order <- c("Control_1", "Tamoxifen_1", "Control_2", "Tamoxifen_2", "Control_3", "Tamoxifen_3")
plotting_df$GroupByTime <- factor(plotting_df$GroupByTime, levels = group_order)

custom_colors2 <- c("Control_1" = "#CDA0A0", "Tamoxifen_1" = "#4C0000", 
                    "Control_2" = "#eac984", "Tamoxifen_2" = "#7C5403", 
                    "Control_3" = "#97f0da", "Tamoxifen_3" = "#037D5F")

# Morphology score -2.5 to 2.5; branchLengthSkel_4 0 to 75; sphericity 0.05 to 0.3; branchEndNodesPerBranch 0 to 20
custom_xlim <- c(-2.5, 2.5)
custom_ylim <- c(-2.5, 2.5)
VIP <- "Morphology_score"
#"Morphology_score", "sphericity", "branchLengthSkel_4", "branchEndNodesPerBranch", "Timepoint_Group")

plot_violin <- ggbetweenstats(
  data = plotting_df,
  x = GroupByTime,
  y = Morphology_score,
  pairwise.display = "none", 
  title = paste(VIP, " ", "by Timepoint and Group"),
  results.subtitle = FALSE,
  xlab = "Timepoint and Group",
  ylab = VIP,
  show.title = FALSE,  # Hide statistics
  group.order = group_order,  # Set the order of groups
  messages = "hide", # Hide p-values
  ylim = custom_ylim,
  mean.color = "red",
  remove.outliers = TRUE,
  outlier.shape = NA,
) + ggplot2::scale_color_manual(values = custom_colors2) + coord_cartesian(ylim = custom_ylim)

plot_density <- ggplot(plotting_df, aes(x = Morphology_score, fill = GroupByTime)) +
  geom_density(alpha = 0.5, position = "identity") +
  scale_fill_manual(values = custom_colors2) +
  labs(
    title = paste("Density Plot of ", VIP, " by Timepoint and Group"),
    x = VIP,
    y = "Density"
  ) + xlim(custom_xlim)  # Set custom X-axis limits

# Save the plot as a PDF
width_density <- 360*2 / 72
height_density <- 360 / 72
width_violin <- 360*2 / 72
height_violin <- 360 / 72

ggsave(paste0("Density_", VIP, ".pdf"), plot = plot_density, width = width_density, height = height_density)
ggsave(paste0("Violin_", VIP, ".pdf"), plot = plot_violin, width = width_violin, height = height_violin)

# Statistics Tamoxifen vs Control -------------------------------------------------------------
setwd("longDTA_MMQT_Analysis")
stats_df_clean <- read_excel("stats_df_clean.xlsx")
stats_df[] <- lapply(stats_df, function(x) {
  if(is.numeric(x)) return(x)
  factor(x)
})

# Define the comparisons
comparisons <- list(
  c("Control_1", "Tamoxifen_1"),
  c("Control_2", "Tamoxifen_2"),
  c("Control_3", "Tamoxifen_3"),
  c("Control_1", "Control_2"),
  c("Control_2", "Control_3"),
  c("Tamoxifen_1", "Tamoxifen_2"),
  c("Tamoxifen_2", "Tamoxifen_3")
)
# Get the numeric column names
col_names <- names(stats_df)[sapply(stats_df_clean, is.numeric)]
# Initialize an empty dataframe to store results
results_df <- data.frame(matrix(NA, nrow = length(comparisons), ncol = length(col_names)))
rownames(results_df) <- sapply(comparisons, paste, collapse = " vs ")
colnames(results_df) <- col_names
# Function to perform the test
perform_test <- function(levels, col_name) {
  group1 <- stats_df_clean[[col_name]][stats_df_clean$GroupByTime == levels[1]]
  group2 <- stats_df_clean[[col_name]][stats_df_clean$GroupByTime == levels[2]]
  test <- t.test(group1, group2)
  return(test$p.value)
}
# Apply the function to each comparison and each numeric column
for (col_name in col_names) {
  for (i in 1:length(comparisons)) {
    results_df[i, col_name] <- perform_test(comparisons[[i]], col_name)
  }
}
# Save the dataframe as an Excel file
write.xlsx(results_df, "longDTA_Welch.xlsx", rowNames = TRUE)

# Summary Statistics ------------------------------------------------------
setwd("longDTA_MMQT_Analysis")

project <- "longDTA"
stats_df_clean <- read_excel("stats_df_clean.xlsx")

stats_df_2<- as.data.frame(stats_df_clean)
stats_df_2[] <- lapply(stats_df_2, function(x) {
  if(is.numeric(x)) return(x)
  factor(x)
})

# Function to calculate summary statistics and Shapiro-Wilk test
summarize_data <- function(x) {
  data.frame(
    Mean = mean(x, na.rm = TRUE),
    Median = median(x, na.rm = TRUE),
    SD = sd(x, na.rm = TRUE),
    N = length(x),
    Shapiro_Wilk = shapiro.test(x)$p.value
  )
}

# Initialize an empty list to store the results
results <- list()
# Get the group names
groups <- unique(stats_df_2$GroupByTime)

# Get the variable names
variables <-  names(stats_df_2)[sapply(stats_df_2, is.numeric)]

# Loop over the groups
for (g in groups) {
  # Subset the dataframe for the current group
  df_group <- stats_df_2[stats_df_2$GroupByTime == g, ]
  
  # Loop over the variables
  for (v in variables) {
    # Calculate the summary statistics for the current variable
    stats <- summarize_data(df_group[[v]])
    
    # Add the group and variable names to the results
    stats$Group <- g
    stats$Variable <- v
    
    # Add the results to the list
    results <- append(results, list(stats))
  }
}

# Combine the results into a single dataframe
result_flat <- do.call(rbind, results)

# Reshape the dataframe
result_wide <- result_flat %>%
  pivot_wider(
    id_cols = "Variable",
    names_from = "Group",
    values_from = c("Mean", "Median", "SD", "N", "Shapiro_Wilk")
  )

# Write the reshaped dataframe to an Excel file
write.xlsx(result_wide, paste("Summary Statistics_", project, ".xlsx",sep=""),colNames = T,rowNames = F,keepNA = F)

# GLM Statistics - ----------------------------------------------------------
setwd("longDTA_MMQT_Analysis")
stats_df_clean <- read_excel("stats_df_clean.xlsx")
stats_df_clean$Group <- as.factor(stats_df_clean$Group)

# Fit models 
model_All_Morphology_score  <- lmer(Morphology_score ~ GroupByTime + (1|AnimalID), data = stats_df_clean)
model_All_sphericity  <- lmer(sphericity ~ GroupByTime + (1|AnimalID), data = stats_df_clean)
model_All_branchLengthSkel_4  <- lmer(branchLengthSkel_4 ~ GroupByTime + (1|AnimalID), data = stats_df_clean)
model_All_branchEndNodesPerBranch  <- lmer(branchEndNodesPerBranch ~ GroupByTime + (1|AnimalID), data = stats_df_clean)

# Export Pairwise comparisons from models
models_list <- list(model_All_Morphology_score,model_All_sphericity, model_All_branchLengthSkel_4, model_All_branchEndNodesPerBranch)
model_names <- c("model_All_Morphology_score", "model_All_sphericity", "model_All_branchLengthSkel_4", "model_All_branchEndNodesPerBranch")

# loop over the models
for(i in seq_along(models_list)) {
  # calculate the estimated marginal means
  emm <- emmeans(models_list[[i]], ~ GroupByTime)
  # perform pairwise comparisons without adjustment
  pairwise_no_adjust <- pairs(emm, adjust = "none")
  # convert to data frame
  pairwise_df <- as.data.frame(pairwise_no_adjust)
  # write to Excel file
  write_xlsx(pairwise_df, paste0(model_names[i], "_pairwise_comparisons.xlsx"))
}

# Checking Normality ------------------------------------------------------

# Q-Q plot of residuals
qqnorm(residuals(model_All_Morphology_score))
qqline(residuals(model_All_Morphology_score))