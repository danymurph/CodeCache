# Load necessary libraries
library(GEOquery)
library(ggplot2)
library(dplyr)
library(tidyr)
library(pheatmap)
library(limma)
library(caret)
library(e1071)      # For SVM
library(factoextra) # For PCA visualization

# --- Data Loading and Preprocessing ---

# Load the dataset
gse <- getGEO("GSE52553", GSEMatrix = TRUE)
expr_data <- exprs(gse[[1]])  # Assuming only one ExpressionSet in gse

# Get sample metadata
sample_metadata_raw <- pData(gse[[1]])

# Extract Patient_ID from 'title' using string manipulation
extract_patient_id <- function(title) {
  parts <- strsplit(as.character(title), "_")[[1]]
  return(parts[3])
}
patient_ids <- sapply(sample_metadata_raw$title, function(x) extract_patient_id(x))

# Create sample_metadata data frame
sample_metadata <- data.frame(
  Sample_ID = sample_metadata_raw$geo_accession,
  Patient_ID = patient_ids,
  Disease_State = tolower(sample_metadata_raw$`subject status:ch1`),
  Agent = tolower(sample_metadata_raw$`treated with:ch1`)
)

# Ensure the sample IDs match between expression data and metadata
expr_data <- expr_data[, sample_metadata$Sample_ID]

# Log2 transformation
expression_matrix_log <- log2(expr_data + 1)  # Add 1 to avoid log(0)

# Check for missing values
if (sum(is.na(expression_matrix_log)) > 0) {
  expression_matrix_log <- na.omit(expression_matrix_log)
}

# --- Outlier Detection and Removal ---

# Perform PCA
pca_res <- prcomp(t(expression_matrix_log), scale. = TRUE)

# Assign Sample_ID as rownames to pca_res$x
rownames(pca_res$x) <- sample_metadata$Sample_ID

# Create PCA DataFrame
pca_df <- data.frame(
  PC1 = pca_res$x[, 1],
  PC2 = pca_res$x[, 2],
  Sample_ID = sample_metadata$Sample_ID,
  Patient_ID = sample_metadata$Patient_ID,
  Disease_State = sample_metadata$Disease_State,
  Agent = sample_metadata$Agent
)

# Plot PCA to visually identify outliers
ggplot(pca_df, aes(x = PC1, y = PC2, color = Disease_State, shape = Agent)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "PCA of Samples Before Outlier Removal", x = "PC1", y = "PC2")

# Identify outliers based on PCA distances with a stricter threshold
pca_distances <- sqrt(rowSums(pca_res$x[, 1:2]^2))
outlier_threshold <- mean(pca_distances) + 1.5 * sd(pca_distances)  # Stricter threshold for outliers
outliers <- rownames(pca_res$x)[pca_distances > outlier_threshold]
cat("Outliers identified:", outliers, "\n")

# Hierarchical Clustering before removing the outlier
distance_matrix <- dist(t(expression_matrix_log))
hc <- hclust(distance_matrix, method = "average")
plot(hc, labels = paste(sample_metadata$Sample_ID, sample_metadata$Gender, sep = "_"), main = "Hierarchical Clustering Dendrogram Before Outlier Removal")

# Get Patient IDs of outliers
outlier_patient_ids <- unique(sample_metadata$Patient_ID[sample_metadata$Sample_ID %in% outliers])
cat("Outlier Patient IDs:", outlier_patient_ids, "\n")

# Identify all samples corresponding to outlier Patient IDs (outliers and their matching samples)
samples_to_remove <- sample_metadata$Sample_ID[sample_metadata$Patient_ID %in% outlier_patient_ids]
cat("Samples to remove (outliers and their matches):", samples_to_remove, "\n")

# Visual proof of outliers and their matches
pca_df$Outlier <- ifelse(pca_df$Sample_ID %in% samples_to_remove, "Outlier/Match", "Non-Outlier")
ggplot(pca_df, aes(x = PC1, y = PC2, color = Outlier)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "Outlier Detection via PCA", x = "PC1", y = "PC2")

# Remove the identified samples from the expression matrix and metadata
expression_matrix_log <- expression_matrix_log[, !colnames(expression_matrix_log) %in% samples_to_remove]
sample_metadata <- sample_metadata[!sample_metadata$Sample_ID %in% samples_to_remove, ]

# Perform PCA again after removing the outliers
pca_res <- prcomp(t(expression_matrix_log), scale. = TRUE)

# Plot PCA after removing the outliers
pca_df <- data.frame(
  PC1 = pca_res$x[, 1],
  PC2 = pca_res$x[, 2],
  sample_metadata
)

ggplot(pca_df, aes(x = PC1, y = PC2, color = Disease_State, shape = Agent)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "PCA After Removing Outliers", x = "PC1", y = "PC2")

# Hierarchical Clustering after removing the outliers
distance_matrix <- dist(t(expression_matrix_log))
hc <- hclust(distance_matrix, method = "average")
plot(hc, labels = paste(sample_metadata$Sample_ID, sample_metadata$Gender, sep = "_"), main = "Hierarchical Clustering After Removing Outliers")


# --- Gene Filtering ---

# Remove genes with low expression across samples
threshold <- 5  # Expression threshold
keep_genes <- rowSums(expression_matrix_log > threshold) >= (ncol(expression_matrix_log) * 0.5)
expression_matrix_filtered <- expression_matrix_log[keep_genes, ]
cat("Number of genes after filtering:", nrow(expression_matrix_filtered), "\n")

# --- Feature Selection Using One-Way ANOVA with Three Levels ---

# Create a combined factor with four levels (two control types)
sample_metadata$Group <- with(sample_metadata, ifelse(
  Disease_State == "alcoholic" & Agent == "75mm ethanol for 24h", "Alcoholic_Ethanol",
  ifelse(Disease_State == "alcoholic" & Agent == "none (untreated control)", "Alcoholic_Untreated",
         ifelse(Disease_State == "non-alcoholic (control)" & Agent == "75mm ethanol for 24h", "Control_Ethanol",
                "Control_Untreated"
         ))))

# Check group counts
group_counts <- table(sample_metadata$Group)
print(group_counts)

# Split samples by group
group_indices <- split(1:ncol(expression_matrix_filtered), sample_metadata$Group)

# Update the one-way ANOVA function to account for four levels
aov.all.genes <- function(x, s1, s2, s3, s4) {
  x1 <- as.numeric(x[s1])
  x2 <- as.numeric(x[s2])
  x3 <- as.numeric(x[s3])
  x4 <- as.numeric(x[s4])
  fac <- c(rep("Group1", length(x1)), rep("Group2", length(x2)),
           rep("Group3", length(x3)), rep("Group4", length(x4)))
  a.dat <- data.frame(factor = as.factor(fac), express = c(x1, x2, x3, x4))
  p.out <- summary(aov(express ~ factor, data = a.dat))[[1]][1, 5]  # p-value
  return(p.out)
}

# Run the ANOVA with the updated groups
aov.run <- apply(expression_matrix_filtered, 1, aov.all.genes,
                 s1 = group_indices$Alcoholic_Ethanol,
                 s2 = group_indices$Alcoholic_Untreated,
                 s3 = group_indices$Control_Ethanol,
                 s4 = group_indices$Control_Untreated)

# Adjust p-values for multiple testing
adjusted_pvals <- p.adjust(aov.run, method = "fdr")

# Threshold for significance
p_value_threshold <- 0.05

# Get significant genes
significant_genes <- names(adjusted_pvals)[adjusted_pvals < p_value_threshold]
cat("Number of significant genes after ANOVA:", length(significant_genes), "\n")

# --- Plotting Histogram of Adjusted P-values ---

ggplot(data.frame(Adjusted_PValue = adjusted_pvals), aes(x = Adjusted_PValue)) +
  geom_histogram(bins = 50, fill = "blue", color = "black") +
  theme_minimal() +
  labs(title = "Adjusted P-values from One-Way ANOVA", x = "Adjusted P-value", y = "Count")

# --- Dimensionality Reduction and Visualization ---

# Subset the data by significant genes
expression_matrix_sig <- expression_matrix_filtered[significant_genes, ]

# Perform PCA with significant genes
pca_sig <- prcomp(t(expression_matrix_sig), scale. = TRUE)

# Assign Sample_ID as rownames to pca_sig$x
rownames(pca_sig$x) <- sample_metadata$Sample_ID

# Create PCA data frame
pca_sig_df <- data.frame(
  PC1 = pca_sig$x[, 1],
  PC2 = pca_sig$x[, 2],
  sample_metadata
)

# Plot PCA using significant genes
ggplot(pca_sig_df, aes(x = PC1, y = PC2, color = Group, shape = Agent)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "PCA Using Significant Genes", x = "PC1", y = "PC2") +
  guides(shape = guide_legend(title = "Agent"))

# Hierarchical Clustering with significant genes
distance_matrix_sig <- dist(t(expression_matrix_sig))
hc_sig <- hclust(distance_matrix_sig, method = "average")
plot(hc_sig, labels = paste(sample_metadata$Sample_ID, sample_metadata$Group, sep = "_"), 
     main = "Hierarchical Clustering with Significant Genes")

# --- Classification using SVM ---

set.seed(123)  # For reproducibility

# Split data into training and testing sets
train_indices <- createDataPartition(sample_metadata$Group, p = 0.7, list = FALSE)
train_samples <- sample_metadata$Sample_ID[train_indices]
test_samples <- sample_metadata$Sample_ID[-train_indices]

train_data <- t(expression_matrix_sig[, train_samples])
test_data <- t(expression_matrix_sig[, test_samples])
train_labels <- sample_metadata$Group[train_indices]
test_labels <- sample_metadata$Group[-train_indices]

# Train SVM classifier
svm_model <- svm(train_data, as.factor(train_labels), kernel = "linear", probability = TRUE)

# Predict on test data
svm_predictions <- predict(svm_model, test_data)

# Create a confusion matrix for SVM results
confusion <- confusionMatrix(svm_predictions, as.factor(test_labels))

# Print confusion matrix
print(confusion)

# Extract PCA scores for test samples
# Ensure that 'pca_sig$x' has rownames as 'Sample_ID's
# and that 'test_samples' are present in 'pca_sig$x'
if(!all(test_samples %in% rownames(pca_sig$x))){
  stop("Some test samples are not present in the PCA results.")
}

pca_scores_test <- pca_sig$x[test_samples, c("PC1", "PC2")]

# Create a data frame for plotting classification results
classification_results <- data.frame(
  PC1 = pca_scores_test[, "PC1"],
  PC2 = pca_scores_test[, "PC2"],
  Actual = test_labels,
  Predicted = svm_predictions,
  Agent = sample_metadata$Agent[-train_indices]
)

# Verify the classification_results data frame
head(classification_results)

# Ensure 'Actual' and 'Predicted' are factors
classification_results$Actual <- factor(classification_results$Actual)
classification_results$Predicted <- factor(classification_results$Predicted)

# Assign unique shapes to 'Actual' classes
# Define a vector of shapes; adjust as needed
actual_shapes <- c("Alcoholic_Ethanol" = 16, 
                   "Alcoholic_Untreated" = 17, 
                   "Control_Untreated" = 15,
                   "Control_Ethanol" = 18)

# Assign unique colors to 'Predicted' classes
predicted_colors <- c("Alcoholic_Ethanol" = "red", 
                      "Alcoholic_Untreated" = "blue", 
                      "Control_Untreated" = "green",
                      "Control_Ethanol" = "purple")

# Plot SVM classification results
ggplot(classification_results, aes(x = PC1, y = PC2, color = Predicted, shape = Actual)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal() +
  labs(title = "SVM Classification Results", x = "PC1", y = "PC2") +
  guides(color = guide_legend(title = "Predicted Class"),
         shape = guide_legend(title = "Actual Class")) +
  scale_color_manual(values = predicted_colors) +
  scale_shape_manual(values = actual_shapes) +
  theme(legend.position = "right")

# --- Identify Top Discriminant Genes ---

# Calculate F-statistics for each gene
f_stats <- apply(expression_matrix_sig, 1, function(x) {
  model <- aov(x ~ sample_metadata$Group)
  summary(model)[[1]][1, "F value"]
})

# Rank genes based on F-statistics
gene_ranking <- sort(f_stats, decreasing = TRUE)

# Identify top 5 genes in positive direction (high F-statistics)
top_genes_positive <- names(gene_ranking)[1:5]

# Identify top 5 genes in negative direction (low F-statistics)
top_genes_negative <- names(gene_ranking)[(length(gene_ranking) - 4):length(gene_ranking)]

# Combine the top genes from positive and negative directions
top_genes <- c(top_genes_positive, top_genes_negative)

# Print top genes
cat("Top 5 positive and negative discriminant genes based on F-statistics:\n")
print(top_genes)

# Write the top 10 genes to a text file in the specified directory
write.table(top_genes, file = "~/Desktop/DAV_proj/top_10_genes.txt", 
            quote = FALSE, row.names = FALSE, col.names = FALSE)

cat("Top 10 genes (5 positive, 5 negative) have been written to '~/Desktop/DAV_proj/top_10_genes.txt'\n")

# --- End of Script ---
