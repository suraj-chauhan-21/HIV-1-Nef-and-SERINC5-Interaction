#Merging File 
# -------------------- Step 0: Load Required Libraries --------------------
library(dplyr)    # For data manipulation
library(tibble)   # For working with tidy data frames

# -------------------- Step 1: Set the Path to HTSeq-Count Files --------------------
# Update the path to your HTSeq files
htseq_path <- "C:/Users/suraj/Documents/cell_lines/cell lines data"

# -------------------- Step 2: Read Metadata --------------------
# 'infect.txt' contains columns: SRR_ID, Infectivity, Cell_Line
meta <- read.table(
  file = file.path(htseq_path, "infect.txt"),
  header = FALSE, sep = "\t", stringsAsFactors = FALSE,
  col.names = c("SRR_ID", "Infectivity", "Cell_Line")
)

# -------------------- Step 3: Get List of HTSeq Files --------------------
# List all files starting with SRR and ending with .htseq
htseq_files <- list.files(
  path = htseq_path,
  pattern = "^SRR.*\\.htseq$",
  full.names = TRUE
)

# -------------------- Step 4: Initialize List to Store Each SRR File --------------------
htseq_list <- list()

# -------------------- Step 5: Read HTSeq Files and Add to List --------------------
for (file in htseq_files) {
  srr_id <- sub(".htseq$", "", basename(file))  # Extract SRR ID
  data <- read.table(
    file, header = FALSE, sep = "\t", stringsAsFactors = FALSE,
    col.names = c("GeneID", "Count")
  )
  colnames(data)[2] <- srr_id                   # Rename column to SRR ID
  htseq_list[[srr_id]] <- data                  # Add to list
}

# -------------------- Step 6: Merge All Data by GeneID --------------------
merged_data <- Reduce(function(x, y) merge(x, y, by = "GeneID"), htseq_list)

# -------------------- Step 7: Keep Only Relevant SRR Columns Present in Metadata --------------------
available_ids <- intersect(meta$SRR_ID, colnames(merged_data))  # Common SRR IDs
meta <- meta[meta$SRR_ID %in% available_ids, ]                   # Filter metadata
merged_data <- merged_data[, c("GeneID", available_ids)]         # Reorder columns in merged data

# -------------------- Step 8: Prepare Annotation Rows (Metadata Rows) --------------------
srr_id_row      <- c("SRR ID",      meta$SRR_ID)
cell_line_row   <- c("Cell Line",   meta$Cell_Line)
infectivity_row <- c("Infectivity", meta$Infectivity)

# -------------------- Step 9: Combine Annotation and Expression Data --------------------
# Convert gene expression data into character matrix for compatibility
data_matrix <- apply(merged_data, 1, as.character)

# Transpose matrix so that each gene becomes a row and SRR IDs become columns
final_output <- rbind(
  srr_id_row,
  cell_line_row,
  infectivity_row,
  t(data_matrix)   # t() transposes rows ↔ columns
)

# -------------------- Step 10: Write the Final Output to a CSV File --------------------
output_path <- file.path(htseq_path, "merged_htseq_data.csv")
write.table(
  final_output, file = output_path, sep = ",",
  row.names = FALSE, col.names = FALSE, quote = FALSE
)

cat("✅ Merged file successfully written to:\n", output_path, "\n")

# 📦 Load required library
library(dplyr)

# 📁 Define file paths
htseq_dir   <- "C:/Users/suraj/Documents/cell_lines/cell lines data"     # Directory containing SRR htseq-count files
infect_file <- file.path(htseq_dir, "infect.txt")                        # Metadata file
output_file <- file.path(htseq_dir, "merged_data.csv")                   # Output file path

# 📄 Step 1: Read metadata (infect.txt)
infect_data <- read.table(infect_file, header = FALSE, sep = "", stringsAsFactors = FALSE)
colnames(infect_data) <- c("SRR_ID", "Infectivity_Ratio", "Cell_Line_Name")

# 🧾 Step 2: List all htseq-count files starting with "SRR"
srr_files <- list.files(htseq_dir, pattern = "^SRR.*", full.names = TRUE)

# 🧬 Step 3: Initialize container to hold merged data
merged_data <- NULL

# 🔁 Step 4: Read each file and merge counts by Gene_ID
for (file in srr_files) {
  
  # 🔍 Extract SRR ID (used as column name)
  srr_id <- basename(file)
  
  # 🧪 Read gene count file safely
  data <- tryCatch({
    read.table(file, header = FALSE, sep = "\t", stringsAsFactors = FALSE, fill = TRUE, comment.char = "#")
  }, error = function(e) {
    cat("⚠️ Error reading:", file, "\n", e$message, "\n")
    return(NULL)
  })
  
  # 🚫 Skip if reading failed
  if (is.null(data)) next
  
  # ✅ Ensure file has exactly 2 columns: Gene_ID and Count
  if (ncol(data) != 2) {
    cat("⚠️ Skipping", file, "- Incorrect column count\n")
    next
  }
  
  # 🏷 Assign column names
  colnames(data) <- c("Gene_ID", srr_id)
  
  # 🔗 Merge by Gene_ID
  if (is.null(merged_data)) {
    merged_data <- data
  } else {
    merged_data <- merge(merged_data, data, by = "Gene_ID", all = TRUE)
  }
}

# 🧹 Step 5: Replace NAs with 0 (assuming missing genes have zero counts)
merged_data[is.na(merged_data)] <- 0

# 🧠 Step 6: Extract SRR IDs (column names except first one)
srr_ids <- colnames(merged_data)[-1]

# 🧬 Step 7: Add annotation rows (cell line & infectivity)
cell_line_names <- sapply(srr_ids, function(id) {
  idx <- match(id, infect_data$SRR_ID)
  if (!is.na(idx)) return(infect_data$Cell_Line_Name[idx]) else return(NA)
})

infectivity_ratios <- sapply(srr_ids, function(id) {
  idx <- match(id, infect_data$SRR_ID)
  if (!is.na(idx)) return(infect_data$Infectivity_Ratio[idx]) else return(NA)
})

# 📊 Step 8: Combine metadata rows and expression matrix
final_data <- rbind(
  c("SRR_ID", srr_ids),
  c("Cell_Line_Name", cell_line_names),
  c("Infectivity_Ratio", infectivity_ratios),
  merged_data
)

# 💾 Step 9: Export final merged file to CSV
write.csv(final_data, output_file, row.names = FALSE, quote = FALSE)
cat("✅ Merged data saved to:", output_file, "\n")

# 📦 Load required library
library(ggplot2)

# 📄 Step 1: Read infectivity data
infect_data <- read.table("infect.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(infect_data) <- c("Sample", "Infectivity", "CellLine")

# 🔽 Step 2: Sort data by Infectivity (highest to lowest)
infect_data <- infect_data[order(infect_data$Infectivity, decreasing = TRUE), ]

# 🎨 Step 3: Assign custom colors based on Cell Line categories
infect_data$Color <- "gray"  # Default color

# Group 1: T-cell lineage
infect_data$Color[infect_data$CellLine %in% c("JURKAT E6.1", "Jurkat tag")] <- "navy"

# Group 2: B-cell lineage
infect_data$Color[infect_data$CellLine %in% c("bl41", "RAMOS", "CEM A 301", "CEM SS", "HSB2")] <- "darkgreen"

# Group 3: Others / Fibroblast / Non-hematopoietic
infect_data$Color[infect_data$CellLine %in% c("C8166", "WI38", "DAUDI", "MT4", "IMR90", "CEM X 174", "HT1080")] <- "firebrick"

# 📊 Step 4: Plot horizontal bar graph using ggplot2
ggplot(infect_data, aes(x = reorder(CellLine, Infectivity), y = Infectivity, fill = Color)) +
  
  geom_bar(stat = "identity", width = 0.7, color = "black", show.legend = FALSE) +  # Border for aesthetics
  
  scale_fill_identity() +  # Use manually assigned colors directly
  
  coord_flip() +  # Flip coordinates to make horizontal bars
  
  theme_minimal(base_size = 14) +  # Clean theme with larger base font
  
  labs(
    title = "Infectivity Ratio (Nef+ / Nef-) by Cell Line",
    x = "",
    y = "Infectivity Ratio"
  ) +
  
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16, color = "steelblue"),
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(size = 12),
    axis.title.y = element_text(size = 14, face = "bold"),
    panel.grid.major.y = element_blank(),  # Cleaner y-axis
    panel.grid.minor = element_blank()
  )















# 📦 Load required libraries
library(ggplot2)
library(gridExtra)
library(tidyr)
library(dplyr)

# 📄 Step 1: Read the merged expression data from the new directory
data <- read.csv("C:/Users/suraj/Documents/cell_lines/cell lines data/merged_data.csv", stringsAsFactors = FALSE)

# 🛠 Step 2: Rename the first column to "Gene" if it's not already
colnames(data)[1] <- "Gene"

# 🔄 Step 3: Reshape to long format — one row per (Gene × CellLine)
long_data <- pivot_longer(
  data,
  cols = -Gene,
  names_to = "CellLine",
  values_to = "Expression"
)

# 🧹 Step 4: Clean data and apply log transformation
long_data <- long_data %>%
  mutate(Expression = as.numeric(Expression)) %>%
  filter(!is.na(Expression) & is.finite(Expression)) %>%
  mutate(LogExpression = log10(Expression + 1))

# 🎨 Step 5: Generate histograms per cell line
plot_list <- long_data %>%
  split(.$CellLine) %>%
  lapply(function(df) {
    ggplot(df, aes(x = LogExpression)) +
      geom_histogram(binwidth = 0.25, fill = "#69b3a2", color = "black", size = 0.2) +
      labs(
        title = df$CellLine[1],
        x = "Log₁₀(Expression + 1)",
        y = "Gene Count"
      ) +
      theme_minimal(base_size = 12) +
      theme(
        plot.title = element_text(size = 13, face = "bold", hjust = 0.5, color = "steelblue"),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 11, face = "bold")
      ) +
      scale_x_continuous(breaks = seq(0, 6, by = 1), limits = c(0, 6)) +
      coord_cartesian(ylim = c(0, 50000))
  })

# 🧱 Step 6: Display all histograms in a grid layout (3 per row)
do.call(grid.arrange, c(plot_list, ncol = 3))





# 📦 Load required libraries
library(ggplot2)
library(ggpubr)

# 📄 Step 1: Read the merged expression data
df <- read.csv("C:/Users/suraj/Documents/cell_lines/cell lines data/merged_data.csv",
               header = TRUE, stringsAsFactors = FALSE)

# 🧬 Step 2: Extract gene expression for SERINC5 (ENSG00000164300)
# Row matching the SERINC5 gene ID, excluding first column
gene_row <- df[df$Gene == "ENSG00000164300", -1]
gene_counts <- suppressWarnings(as.numeric(gene_row)) # Convert to numeric safely

# 📊 Step 3: Get total reads from row 4 for RPM normalization
total_reads <- suppressWarnings(as.numeric(df[4, -1]))

# 🔬 Step 4: Calculate RPM (Reads Per Million)
rpm <- (gene_counts / total_reads) * 1e6

# 🧮 Step 5: Normalize RPM to a 0–35 scale for plotting
normalized_rpm <- (rpm / max(rpm, na.rm = TRUE)) * 35

# ⚙️ Step 6: Get Infectivity ratio from row 3
infectivity_ratio <- suppressWarnings(as.numeric(df[3, -1]))

# 📋 Step 7: Combine into one dataframe
plot_data <- data.frame(
  Infectivity_Ratio = infectivity_ratio,
  Normalized_RPM = normalized_rpm
)

# 🧼 Step 8: Clean up NA or infinite values
plot_data <- na.omit(plot_data)
plot_data <- plot_data[is.finite(plot_data$Infectivity_Ratio) & is.finite(plot_data$Normalized_RPM), ]

# 🎨 Step 9: Add expression-based color categories
plot_data$Color <- ifelse(
  plot_data$Normalized_RPM > 25, "high",
  ifelse(plot_data$Normalized_RPM > 7, "medium", "low")
)

# 🌟 Step 10: Final Scatter Plot with correlation line
ggplot(plot_data, aes(x = Infectivity_Ratio, y = Normalized_RPM, color = Color)) +
  geom_point(size = 4) +  # Larger points for better visibility
  scale_color_manual(values = c("low" = "red", "medium" = "orange", "high" = "blue")) +
  geom_smooth(method = "lm", se = FALSE, color = "black", size = 0.5) +  # Trend line
  stat_cor(
    method = "pearson",
    aes(label = paste0("italic(r)==", ..r.., "*','~~italic(R)^2==", ..rr.., "*','~~italic(P)==", ..p..)),
    parse = TRUE,
    label.x = 2, label.y = 33,  # Adjust based on max RPM scale
    size = 4.5, color = "black"
  ) +
  labs(
    x = expression("Nef"^"+"*"/Nef"^"-"* " Infectivity Ratio"),
    y = "SERINC5 Expression (Normalized RPM)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    legend.position = "none",
    axis.title = element_text(size = 13, face = "bold"),
    axis.text = element_text(size = 11),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black")
  )







# 🧠 WHY: Shapiro-Wilk test tells us if gene expression data in each cell line follows a normal distribution (bell-shaped curve).
# This helps decide whether to use parametric or non-parametric tests later.

# Load necessary package
library(dplyr)

# Step 1: Load the file with gene expression data
csv_file <- "C:/Users/suraj/Documents/cell_lines/cell lines data/gene id and cell line.csv"
data_raw <- read.csv(csv_file, header = FALSE, stringsAsFactors = FALSE)
# ---

# Step 2: Extract column names (cell line names) from row 1, except the first column (GeneID)
sample_ids <- as.character(data_raw[1, -1])
# ---

# Step 3: Remove metadata rows and rename columns
gene_expr <- data_raw[-(1:3), ]
colnames(gene_expr) <- c("GeneID", sample_ids)
# ---

# Step 4: Make sure all gene expression numbers are actually numeric
gene_expr[,-1] <- suppressWarnings(sapply(gene_expr[,-1], as.numeric))
gene_expr <- as.data.frame(gene_expr, stringsAsFactors = FALSE)
# ---

# Step 5: Define function to do Shapiro-Wilk safely
shapiro_test_safe <- function(values) {
  values <- na.omit(values) 
  if (length(values) > 5000) {
    values <- sample(values, 5000)
  }
  if (length(values) >= 3 && length(unique(values)) > 1) {
    return(shapiro.test(values))
  } else {
    return(list(statistic = NA, p.value = NA)) 
  }
}
# ---

# Step 6: Run the test on all cell lines
normality_results <- lapply(sample_ids, function(samp) shapiro_test_safe(gene_expr[[samp]]))
names(normality_results) <- sample_ids
# ---

# Step 7: Organize results nicely
final_results <- data.frame(
  W_statistic = sapply(normality_results, function(res) round(res$statistic, 4)),
  p_value = sapply(normality_results, function(res) signif(res$p.value, 3)),
  stringsAsFactors = FALSE
)
print(final_results)
# -----------------------------------------------------------------------------------------------------------------------


# 🧠 WHY: We use the Wilcoxon test when the data isn't normally distributed. It compares the values against their own median.

# Reusing previous data loading steps...

# Step 1: Define function to do Wilcoxon safely
wilcoxon_test_safe <- function(values) {
  values <- na.omit(values)
  if (length(values) > 5000) {
    values <- sample(values, 5000)
  }
  if (length(values) >= 3 && length(unique(values)) > 1) {
    return(wilcox.test(values, mu = median(values), exact = FALSE, correct = TRUE))
  } else {
    return(list(statistic = NA, p.value = NA))
  }
}
# ---

# Step 2: Apply it on all samples
wilcoxon_results <- lapply(sample_ids, function(samp) wilcoxon_test_safe(gene_expr[[samp]]))
names(wilcoxon_results) <- sample_ids
# ---

# Step 3: Make a table
final_results <- data.frame(
  W_statistic = sapply(wilcoxon_results, function(res) round(res$statistic, 4)),
  p_value = sapply(wilcoxon_results, function(res) signif(res$p.value, 5)),
  stringsAsFactors = FALSE
)
print(final_results)
# ----------------------------------------------------------------------------------------------------------------------


# 🧠 WHY: Partial correlation removes the effect of other variables so we see how strongly two cell lines are connected without noise.

library(ppcor)
library(ggcorrplot)

# Step 1: Load the cleaned data
data <- read.csv("C:/Users/suraj/Documents/cell_lines/cell lines data/gene id and cell line.csv", stringsAsFactors = FALSE)
data <- data[-c(1:3), ]
colnames(data)[1] <- "Gene_ID"
data[, -1] <- apply(data[, -1], 2, function(x) as.numeric(as.character(x)))
# ---

# Step 2: Compute partial correlation
partial_corr_matrix <- pcor(data[, -1])$estimate
# ---

# Step 3: Plot heatmap
corr_plot <- ggcorrplot(partial_corr_matrix, hc.order = TRUE, type = "lower",
                        lab = TRUE, lab_size = 3,
                        colors = c("blue", "white", "red"),
                        title = "Partial Correlation Matrix of Cell Lines")
png("partial_correlation_plot.png", width = 1000, height = 800)
print(corr_plot)dev.off() 
# ----------------------------------------------------------------------------------------------------------------------


# 🧠 WHY: This is simple Pearson correlation, useful to check how closely two cell lines behave based on expression.

library(ggcorrplot)

data <- read.csv("C:/Users/suraj/Documents/cell_lines/cell lines data/gene id and cell line.csv", stringsAsFactors = FALSE)
data <- data[-c(1:3), ]
colnames(data)[1] <- "Gene_ID"
data[, -1] <- apply(data[, -1], 2, function(x) as.numeric(as.character(x)))

corr_matrix <- cor(data[, -1], use = "complete.obs", method = "pearson")

corr_plot <- ggcorrplot(corr_matrix, hc.order = TRUE, type = "lower",
                        lab = TRUE, lab_size = 3,
                        colors = c("blue", "white", "red"),
                        title = "Complete Correlation Matrix of Cell Lines")
png("complete_correlation_plot.png", width = 1000, height = 800)
print(corr_plot)
dev.off()




# ----------------------------------------------------------------------------------------------------------------------


# 🧠 WHY: Skewness tells us if the data is lopsided (left or right).
# Kurtosis tells us if the data has heavy tails (extreme values).

library(e1071)
file_path <- "C:/Users/suraj/Documents/cell_lines/cell lines data/gene id and cell line.csv"
df <- read.csv(file_path, stringsAsFactors = FALSE)
df_clean <- df[-c(1:3), ]

for (col in names(df_clean)[-1]) {
  df_clean[[col]] <- as.numeric(df_clean[[col]])
}
skewness_vals <- sapply(df_clean[-1], skewness, na.rm = TRUE)
kurtosis_vals <- sapply(df_clean[-1], kurtosis, na.rm = TRUE)

summary_stats <- cbind(Skewness = skewness_vals, Kurtosis = kurtosis_vals)
summary_stats <- as.data.frame(summary_stats)
summary_stats$Cell_Line <- rownames(summary_stats)

# Plotting
ggplot(summary_stats, aes(x = Cell_Line, y = Skewness)) +
  geom_bar(stat = "identity", fill = "tomato") +
  labs(title = "Skewness of Gene Expression", x = "Cell Line", y = "Skewness") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(summary_stats, aes(x = Cell_Line, y = Kurtosis)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  labs(title = "Kurtosis of Gene Expression", x = "Cell Line", y = "Kurtosis") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
# -----------------------------------------------------------------------------------------------------------------
# Set working directory
setwd("C:/Users/suraj/Documents/cell_lines/cell lines data")

# Load required packages
library(org.Hs.eg.db)
library(dplyr)

# Read in data from multiple files and merge into a single dataframe
final_df <- read.delim("SRR2166626.htseq", header = FALSE)
infect_df <- read.delim("infect.txt", header = FALSE)

# Create filenames for other data files
file_names <- paste(infect_df$V1, ".htseq", sep = "")

# Loop through the files, reading each one into a new column of the final dataframe
for (file_name in file_names) {
  if (file_name != "SRR2166626.htseq") {
    temp_df <- read.delim(file_name, header = FALSE)
    final_df <- cbind(final_df, temp_df$V2)
  }
}

# Map ENSEMBL IDs to gene symbols
final_df$GeneID <- mapIds(org.Hs.eg.db, keys = final_df$V1, keytype = "ENSEMBL", column = "SYMBOL")
final_df$V1 <- final_df$GeneID
final_df$GeneID <- NULL

# Rename columns using sample names from infect.txt
colnames(final_df) <- c("GeneID", infect_df$V3)

# Keep only rows with gene data (removing summary rows)
final_df <- final_df[1:58302,]

# Define a normalization function (TPM-style normalization)
normalize <- function(x) {
  return(x / sum(x) * 1000000)
}

# Apply normalization to all samples (excluding GeneID)
normalized_df <- as.data.frame(lapply(final_df[, 2:16], normalize))

# Add GeneID column back
normalized_df <- cbind(normalized_df, GeneID = final_df$GeneID)

# Reorder columns so GeneID is first
normalized_df <- normalized_df[, c(16, 1:15)]

# View final dataframes (optional, can be commented out)
View(final_df)
View(normalized_df)

# --- Plot 1: Scatter plot of SERINC5 expression vs infectivity ratio ---
serinc5_expression <- as.vector(na.omit(as.numeric(unlist(normalized_df[normalized_df$GeneID == "SERINC5", ]))))
inf_ratio <- infect_df$V2

# Color coding based on infectivity levels
colors <- ifelse(inf_ratio < 9, "red", ifelse(inf_ratio >= 5 & inf_ratio < 30, "green", "blue"))

# Plot with regression line
plot(inf_ratio, serinc5_expression, pch = 16, col = colors, 
     main = "SERINC5 Expression vs Infectivity Ratio", 
     xlab = "Infectivity Ratio (Nef+/Nef-)", ylab = "SERINC5 Expression (normalized)")
reg <- lm(serinc5_expression ~ inf_ratio)
abline(reg)

# Correlation (Kendall)
cor(inf_ratio, serinc5_expression, method = "kendall")

# --- Plot 2: Horizontal barplot of infectivity ratios ---
x <- unlist(infect_df$V2)  # Infectivity ratios
y <- unlist(infect_df$V3)  # Sample names (Nef status)

# Set plot margins
par(mar = c(5, 8, 4, 2) + 0.1)

# Assign color by category
color <- rep(c("red", "green", "blue"), c(8, 5, 2))

# Draw horizontal barplot
barplot(x, names.arg = y, horiz = TRUE, col = color, las = 1,
        xlab = "Nef+/Nef- Infectivity Ratio", main = "Infectivity by Nef Status")



# Set working directory
setwd("C:/Users/suraj/Documents/cell_lines/cell lines data")

# Load required packages
library(org.Hs.eg.db)
library(dplyr)

# Read infect.txt properly (assumes tab or space separated with 3 columns)
infect_df <- read.delim("infect.txt", header = FALSE, sep = "", stringsAsFactors = FALSE)

# Read the first file to initialize the dataframe
final_df <- read.delim(paste0(trimws(infect_df$V1[1]), ".htseq"), header = FALSE)

# Loop to read the rest of the files and add columns
for (i in 2:nrow(infect_df)) {
  file_name <- paste0(trimws(infect_df$V1[i]), ".htseq")
  temp_df <- read.delim(file_name, header = FALSE)
  final_df <- cbind(final_df, temp_df$V2)
}

# Map ENSEMBL IDs to gene names
final_df$GeneID <- mapIds(org.Hs.eg.db, keys = final_df$V1, keytype = "ENSEMBL", column = "SYMBOL")
final_df$V1 <- final_df$GeneID
final_df$GeneID <- NULL

# Assign column names: first column = GeneID, others = sample names
colnames(final_df) <- c("GeneID", trimws(infect_df$V3))

# Remove extra rows (metadata), keep only genes
final_df <- final_df[1:58302, ]

# Define normalization function
normalize <- function(x) {
  return(x / sum(x) * 1000000)
}

# Normalize data columns
normalized_df <- as.data.frame(lapply(final_df[, 2:ncol(final_df)], normalize))

# Add back GeneID and reorder columns
normalized_df <- cbind(GeneID = final_df$GeneID, normalized_df)

# View normalized data
View(normalized_df)

# Plot SERINC5 expression vs infectivity
serinc5_expression <- as.vector(na.omit(as.numeric(unlist(normalized_df[normalized_df$GeneID == "SERINC5", 2:ncol(normalized_df)]))))
inf_ratio <- as.numeric(infect_df$V2)
colors <- ifelse(inf_ratio < 9, "red", ifelse(inf_ratio >= 5 & inf_ratio < 30, "green", "blue"))

plot(inf_ratio, serinc5_expression, pch = 16, col = colors,
     xlab = "Infectivity Ratio (Nef+/Nef-)",
     ylab = "SERINC5 Normalized Expression",
     main = "SERINC5 vs. Infectivity")

# Regression line and correlation
reg <- lm(serinc5_expression ~ inf_ratio)
abline(reg, col = "black")
print(cor(inf_ratio, serinc5_expression, method = "kendall"))

# Horizontal barplot for infectivity ratios
x <- unlist(infect_df$V2)
y <- unlist(infect_df$V3)
par(mar = c(5, 8, 4, 2) + 0.1)
color <- rep(c("red", "green", "blue"), c(8, 5, 2))  # adjust based on your group sizes
barplot(x, names.arg = y, horiz = TRUE, col = color, las = 1,
        xlab = "Nef+/Nef- Infectivity Ratio",
        main = "Infectivity Across Cell Lines")



# Spearman correlation test

# Set the gene of interest
gene_of_interest <- "SERINC5"

# Remove rows with missing or empty gene IDs
normalized_df <- subset(normalized_df, !is.na(GeneID) & GeneID != "")

# Create an empty data frame to store correlations
cor.df <- data.frame(name1=NA, name2=NA, correl=NA)

# Loop over all genes except the gene of interest and calculate Spearman correlation coefficient
for(j in 1:nrow(normalized_df)) {
  if (normalized_df[j, "GeneID"] != gene_of_interest) {
    
    # Expression values of the gene of interest
    v1 <- as.numeric(normalized_df[normalized_df$GeneID == gene_of_interest, 2:ncol(normalized_df)])
    
    # Expression values of the current gene
    v2 <- as.numeric(normalized_df[j, 2:ncol(normalized_df)])
    
    # Calculate Spearman correlation
    correl <- cor(v1, v2, method = "spearman")
    
    name1 <- gene_of_interest
    name2 <- normalized_df[j, "GeneID"]
    
    # Add the correlation to the data frame
    dftemp <- data.frame(name1, name2, correl)
    cor.df <- rbind(cor.df, dftemp)
  }
}

# Remove rows with missing values
cor.df <- na.omit(cor.df)

# Sort the data frame by descending correlation coefficient
library(dplyr)
sorted_cordf <- cor.df %>% arrange(desc(correl))

# View the sorted data frame
View(sorted_cordf)




# Get the gene expression data for CD3G and SERINC5
cd3g_expression <- as.vector(na.omit(as.numeric(unlist(normalized_df[normalized_df$GeneID == "CD3G", 2:ncol(normalized_df)]))))
serinc5_expression <- as.vector(na.omit(as.numeric(unlist(normalized_df[normalized_df$GeneID == "SERINC5", 2:ncol(normalized_df)]))))

# Run the correlation test between gene expression of CD3G and SERINC5
cor_test <- cor.test(serinc5_expression, cd3g_expression, method = "spearman")
print(cor_test)

# Correlation coefficient and R-squared
correlation <- cor(cd3g_expression, serinc5_expression)
r_squared <- correlation^2

# Get the infection ratio data
inf_ratio <- infect_df$V2

# Correlation between CD3G expression and Infectivity ratio
cor_test1 <- cor.test(inf_ratio, cd3g_expression, method = "spearman")
p_val <- cor_test1$p.value
print(paste("P-value for CD3G:", p_val))

# Normality tests
shapiro.test(cd3g_expression)
shapiro.test(serinc5_expression)

# QQ plot for CD3G
qqnorm(cd3g_expression)
qqline(cd3g_expression)

# Scatter plot: CD3G (y) vs SERINC5 (x)
colors <- ifelse(serinc5_expression < 100, "red",
                 ifelse(serinc5_expression >= 300 & serinc5_expression < 400, "blue", "green"))

plot(serinc5_expression, cd3g_expression, pch = 16, col = colors,
     xlab = "SERINC5 expression", ylab = "CD3G expression", main = "CD3G vs SERINC5")

# Fit regression line
reg <- lm(cd3g_expression ~ serinc5_expression)
abline(reg, col = "blue")
print(summary(reg))

# Scatter plot: Infectivity ratio (x) vs CD3G (y)
colors <- ifelse(inf_ratio < 9, "red",
                 ifelse(inf_ratio >= 5 & inf_ratio < 30, "green", "blue"))

plot(inf_ratio, cd3g_expression, pch = 16, col = colors,
     xlab = "Infectivity Ratio", ylab = "CD3G Expression", main = "CD3G vs Infectivity Ratio")

# Fit regression line
reg <- lm(cd3g_expression ~ inf_ratio)
abline(reg, col = "blue")
print(summary(reg))

# Horizontal bar plot of CD3G expression
x <- as.numeric(na.omit(normalized_df[normalized_df$GeneID == "CD3G", 2:ncol(normalized_df)]))
y <- as.character(infect_df$V3)

# Set margins to make y-axis names visible
par(mar = c(5, 8, 4, 2) + 0.1)

# Assign colors based on infectivity ratio
colors <- ifelse(inf_ratio >= 50, "blue", "red")

# Create barplot
barplot(x, names.arg = y, horiz = TRUE, col = colors, las = 1,
        main = "Differential expression of CD3G",
        xlab = "CD3G gene expression")

# Set working directory
setwd("C:/Users/suraj/Documents/cell_lines/cell lines data")

# Load required packages
library(org.Hs.eg.db)
library(dplyr)

# Read infect.txt properly (assumes tab or space separated with 3 columns)
infect_df <- read.delim("infect.txt", header = FALSE, sep = "", stringsAsFactors = FALSE)

# Read the first file to initialize the dataframe
final_df <- read.delim(paste0(trimws(infect_df$V1[1]), ".htseq"), header = FALSE)

# Loop to read the rest of the files and add columns
for (i in 2:nrow(infect_df)) {
  file_name <- paste0(trimws(infect_df$V1[i]), ".htseq")
  temp_df <- read.delim(file_name, header = FALSE)
  final_df <- cbind(final_df, temp_df$V2)
}

# Map ENSEMBL IDs to gene names
final_df$GeneID <- mapIds(org.Hs.eg.db, keys = final_df$V1, keytype = "ENSEMBL", column = "SYMBOL")
final_df$V1 <- final_df$GeneID
final_df$GeneID <- NULL

# Assign column names: first column = GeneID, others = sample names
colnames(final_df) <- c("GeneID", trimws(infect_df$V3))

# Remove extra rows (metadata), keep only genes
final_df <- final_df[1:58302, ]

# Define normalization function
normalize <- function(x) {
  return(x / sum(x) * 1000000)
}

# Normalize data columns
normalized_df <- as.data.frame(lapply(final_df[, 2:ncol(final_df)], normalize))

# Add back GeneID and reorder columns
normalized_df <- cbind(GeneID = final_df$GeneID, normalized_df)

# View normalized data
View(normalized_df)

# Plot SERINC5 expression vs infectivity
serinc5_expression <- as.vector(na.omit(as.numeric(unlist(normalized_df[normalized_df$GeneID == "SERINC5", 2:ncol(normalized_df)]))))
inf_ratio <- as.numeric(infect_df$V2)
colors <- ifelse(inf_ratio < 9, "red", ifelse(inf_ratio >= 5 & inf_ratio < 30, "green", "blue"))

plot(inf_ratio, serinc5_expression, pch = 16, col = colors,
     xlab = "Infectivity Ratio (Nef+/Nef-)",
     ylab = "SERINC5 Normalized Expression",
     main = "SERINC5 vs. Infectivity")

# Regression line and correlation
reg <- lm(serinc5_expression ~ inf_ratio)
abline(reg, col = "black")
print(cor(inf_ratio, serinc5_expression, method = "kendall"))

# Horizontal barplot for infectivity ratios
x <- unlist(infect_df$V2)
y <- unlist(infect_df$V3)
par(mar = c(5, 8, 4, 2) + 0.1)
color <- rep(c("red", "green", "blue"), c(8, 5, 2))  # adjust based on your group sizes
barplot(x, names.arg = y, horiz = TRUE, col = color, las = 1,
        xlab = "Nef+/Nef- Infectivity Ratio",
        main = "Infectivity Across Cell Lines")



# Spearman correlation test

# Set the gene of interest
gene_of_interest <- "SERINC5"

# Remove rows with missing or empty gene IDs
normalized_df <- subset(normalized_df, !is.na(GeneID) & GeneID != "")

# Create an empty data frame to store correlations
cor.df <- data.frame(name1=NA, name2=NA, correl=NA)

# Loop over all genes except the gene of interest and calculate Spearman correlation coefficient
for(j in 1:nrow(normalized_df)) {
  if (normalized_df[j, "GeneID"] != gene_of_interest) {
    
    # Expression values of the gene of interest
    v1 <- as.numeric(normalized_df[normalized_df$GeneID == gene_of_interest, 2:ncol(normalized_df)])
    
    # Expression values of the current gene
    v2 <- as.numeric(normalized_df[j, 2:ncol(normalized_df)])
    
    # Calculate Spearman correlation
    correl <- cor(v1, v2, method = "spearman")
    
    name1 <- gene_of_interest
    name2 <- normalized_df[j, "GeneID"]
    
    # Add the correlation to the data frame
    dftemp <- data.frame(name1, name2, correl)
    cor.df <- rbind(cor.df, dftemp)
  }
}

# Remove rows with missing values
cor.df <- na.omit(cor.df)

# Sort the data frame by descending correlation coefficient
library(dplyr)
sorted_cordf <- cor.df %>% arrange(desc(correl))

# View the sorted data frame
View(sorted_cordf)




# Get the gene expression data for CD3G and SERINC5
cd3g_expression <- as.vector(na.omit(as.numeric(unlist(normalized_df[normalized_df$GeneID == "CD3G", 2:ncol(normalized_df)]))))
serinc5_expression <- as.vector(na.omit(as.numeric(unlist(normalized_df[normalized_df$GeneID == "SERINC5", 2:ncol(normalized_df)]))))

# Run the correlation test between gene expression of CD3G and SERINC5
cor_test <- cor.test(serinc5_expression, cd3g_expression, method = "spearman")
print(cor_test)

# Correlation coefficient and R-squared
correlation <- cor(cd3g_expression, serinc5_expression)
r_squared <- correlation^2

# Get the infection ratio data
inf_ratio <- infect_df$V2

# Correlation between CD3G expression and Infectivity ratio
cor_test1 <- cor.test(inf_ratio, cd3g_expression, method = "spearman")
p_val <- cor_test1$p.value
print(paste("P-value for CD3G:", p_val))

# Normality tests
shapiro.test(cd3g_expression)
shapiro.test(serinc5_expression)

# QQ plot for CD3G
qqnorm(cd3g_expression)
qqline(cd3g_expression)

# Scatter plot: CD3G (y) vs SERINC5 (x)
colors <- ifelse(serinc5_expression < 100, "red",
                 ifelse(serinc5_expression >= 300 & serinc5_expression < 400, "blue", "green"))

plot(serinc5_expression, cd3g_expression, pch = 16, col = colors,
     xlab = "SERINC5 expression", ylab = "CD3G expression", main = "CD3G vs SERINC5")

# Fit regression line
reg <- lm(cd3g_expression ~ serinc5_expression)
abline(reg, col = "blue")
print(summary(reg))

# Scatter plot: Infectivity ratio (x) vs CD3G (y)
colors <- ifelse(inf_ratio < 9, "red",
                 ifelse(inf_ratio >= 5 & inf_ratio < 30, "green", "blue"))

plot(inf_ratio, cd3g_expression, pch = 16, col = colors,
     xlab = "Infectivity Ratio", ylab = "CD3G Expression", main = "CD3G vs Infectivity Ratio")

# Fit regression line
reg <- lm(cd3g_expression ~ inf_ratio)
abline(reg, col = "blue")
print(summary(reg))

# Horizontal bar plot of CD3G expression
x <- as.numeric(na.omit(normalized_df[normalized_df$GeneID == "CD3G", 2:ncol(normalized_df)]))
y <- as.character(infect_df$V3)

# Set margins to make y-axis names visible
par(mar = c(5, 8, 4, 2) + 0.1)

# Assign colors based on infectivity ratio
colors <- ifelse(inf_ratio >= 50, "blue", "red")

# Create barplot
barplot(x, names.arg = y, horiz = TRUE, col = colors, las = 1,
        main = "Differential expression of CD3G",
        xlab = "CD3G gene expression")



# Load library
library(pheatmap)

# Subset for a few genes of interest
genes_of_interest <- c("CD3G", "SERINC5", "CD4", "CD8A", "IFNG")  # add/remove as needed
heatmap_data <- normalized_df[normalized_df$GeneID %in% genes_of_interest, ]

# Remove GeneID column, convert to matrix, and scale rows
rownames(heatmap_data) <- heatmap_data$GeneID
heatmap_matrix <- as.matrix(heatmap_data[, -1])
heatmap_matrix <- t(scale(t(heatmap_matrix)))  # row-wise Z-score normalization

# Plot heatmap
pheatmap(heatmap_matrix,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         fontsize_row = 10,
         fontsize_col = 8,
         main = "Heatmap of Selected Gene Expressions")





# Subset for selected genes
genes_of_interest <- c("CD3G", "SERINC5", "CD4", "CD8A", "IFNG")
corr_data <- normalized_df[normalized_df$GeneID %in% genes_of_interest, ]
rownames(corr_data) <- corr_data$GeneID
corr_matrix <- cor(t(corr_data[, -1]), method = "spearman")  # rows = genes

# Correlation heatmap
library(corrplot)
corrplot(corr_matrix, method = "color", type = "upper",
         addCoef.col = "black", tl.col = "black",
         title = "Spearman Correlation Between Genes",
         mar = c(0,0,2,0))

# Load library
library(ppcor)

# Build dataframe with all variables
pcor_data <- data.frame(
  CD3G = cd3g_expression,
  SERINC5 = serinc5_expression,
  INFECTIVITY = inf_ratio
)

# Perform partial correlation
pcor_result <- pcor(pcor_data, method = "spearman")
print(pcor_result)

# Extract CD3G vs SERINC5 controlling for INFECTIVITY
cat("Partial Correlation (CD3G vs SERINC5 | INFECTIVITY):\n")
print(pcor_result$estimate["CD3G", "SERINC5"])
print(paste("p-value:", pcor_result$p.value["CD3G", "SERINC5"]))



# Load required packages
library(ggplot2)
library(ggrepel)
library(tidyverse)

# Transpose expression data for PCA: genes as columns, samples as rows
expr_matrix <- t(as.matrix(normalized_df[, -1]))
colnames(expr_matrix) <- normalized_df$GeneID
rownames(expr_matrix) <- infect_df$V3  # Sample names

# Select top variable genes
var_genes <- apply(expr_matrix, 2, var)
top_genes <- names(sort(var_genes, decreasing = TRUE))[1:10]  # adjust this if needed
pca_input <- expr_matrix[, top_genes]

# Perform PCA
pca_result <- prcomp(pca_input, center = TRUE, scale. = TRUE)

# Create a data frame for PCA scores (samples)
scores_df <- as.data.frame(pca_result$x)
scores_df$Sample <- rownames(scores_df)
scores_df$Infectivity <- infect_df$V2  # attach infectivity ratio if needed

# Create a data frame for PCA loadings (genes)
loadings_df <- as.data.frame(pca_result$rotation)
loadings_df$Gene <- rownames(loadings_df)

# Biplot with ggplot2
ggplot(scores_df, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = Infectivity), size = 3) +
  geom_text_repel(aes(label = Sample), size = 3, max.overlaps = 10) +
  geom_segment(data = loadings_df,
               aes(x = 0, y = 0, xend = PC1 * 5, yend = PC2 * 5),
               arrow = arrow(length = unit(0.2, "cm")), color = "gray30") +
  geom_text_repel(data = loadings_df,
                  aes(x = PC1 * 5, y = PC2 * 5, label = Gene),
                  size = 3, color = "black") +
  theme_minimal() +
  labs(title = "ggplot2 PCA Biplot of Top Genes",
       x = paste0("PC1 (", round(summary(pca_result)$importance[2, 1] * 100, 1), "%)"),
       y = paste0("PC2 (", round(summary(pca_result)$importance[2, 2] * 100, 1), "%)")) +
  scale_color_gradient(low = "blue", high = "red") +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5))
library(ggplot2)



# Add cell type annotations if you have them
infect_df$CellType <- c("T", "T", "T", "T", "T", "T", "T", "T", "B", "B", "Mono", "Mono", "NK")  # Example

# Get CD3G expression vector
cd3g_expr <- as.numeric(normalized_df[normalized_df$GeneID == "CD3G", -1])

# Create data frame
cd3g_df <- data.frame(Expression = cd3g_expr,
                      CellLine = infect_df$V3,
                      CellType = infect_df$CellType)

ggplot(cd3g_df, aes(x = CellType, y = Expression, fill = CellType)) +
  geom_violin(trim = FALSE, alpha = 0.4) +
  geom_jitter(width = 0.2, size = 2) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "white") +
  labs(title = "CD3G Expression Across Cell Types",
       y = "Normalized Expression", x = "Cell Type") +
  theme_minimal()
library(ggplot2)

# Add cell type annotations if you have them
infect_df$CellType <- c("T", "T", "T", "T", "T", "T", "T", "T", "B", "B", "Mono", "Mono", "NK")  # Example

# Get CD3G expression vector
cd3g_expr <- as.numeric(normalized_df[normalized_df$GeneID == "CD3G", -1])

# Create data frame
cd3g_df <- data.frame(Expression = cd3g_expr,
                      CellLine = infect_df$V3,
                      CellType = infect_df$CellType)

ggplot(cd3g_df, aes(x = CellType, y = Expression, fill = CellType)) +
  geom_violin(trim = FALSE, alpha = 0.4) +
  geom_jitter(width = 0.2, size = 2) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "white") +
  labs(title = "CD3G Expression Across Cell Types",
       y = "Normalized Expression", x = "Cell Type") +
  theme_minimal()




