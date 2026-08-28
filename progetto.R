## ----message=FALSE, warning=FALSE---------------------------------------------------------------------------------------------
# setwd("path/to/your/project/folder") # Uncomment if needed

# Install necessary libraries if they are missing (conditional to prevent knitting errors)
if (!requireNamespace("minfi", quietly = TRUE)) BiocManager::install("minfi", ask = FALSE)
if (!requireNamespace("IlluminaHumanMethylation450kmanifest", quietly = TRUE)) BiocManager::install("IlluminaHumanMethylation450kmanifest", ask = FALSE)
if (!requireNamespace("IlluminaHumanMethylation450kanno.ilmn12.hg19", quietly = TRUE)) BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19", ask = FALSE)
if (!requireNamespace("factoextra", quietly = TRUE)) install.packages("factoextra", repos = "https://cloud.r-project.org")
if (!requireNamespace("qqman", quietly = TRUE)) install.packages("qqman", repos = "https://cloud.r-project.org")
if (!requireNamespace("gplots", quietly = TRUE)) install.packages("gplots", repos = "https://cloud.r-project.org")

# Load the required libraries
library(minfi)
library(IlluminaHumanMethylation450kmanifest)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(factoextra)
library(qqman)
library(gplots)

# Load the cleaned Illumina 450K manifest file containing probe annotations
baseDir <- "Input_Data"

load("Illumina450Manifest_clean.RData")


# Read the sample sheet containing metadata and sample information
targets <- read.metharray.sheet(baseDir)

# Creates a RGChannelSet which stores the raw Red and Green light intensities for each physical Address on the microarray
RGset <- read.metharray.exp(targets = targets)

# T explore the RGset object:
# RGset


## -----------------------------------------------------------------------------------------------------------------------------
# Create a dataframe containing the Red fluorescence intensities extracted from RGset
Red <- data.frame(getRed(RGset))
# Check the dimension and the first rows of the object
dim(Red)
head(Red)

# Create a dataframe containing the Green fluorescence intensities extracted from RGset
Green <- data.frame(getGreen(RGset))
# Check the dimension and the first rows of the object
dim(Green)
head(Green)


## -----------------------------------------------------------------------------------------------------------------------------
# Extract the Red fluorescence values for the address assigned to Group 3
Red['67790505', ]

# Extract the Green fluorescence values for the address assigned to Group 3
Green['67790505', ]

# Check in the manifest whether the selected address corresponds to a Type I or Type II probe
Type <- Illumina450Manifest_clean[
  Illumina450Manifest_clean$AddressA_ID == 67790505 | Illumina450Manifest_clean$AddressB_ID == 67790505, 
  'Infinium_Design_Type'
]
Type

# Check in the manifest the color channel of the selected probe
Color <- Illumina450Manifest_clean[
  Illumina450Manifest_clean$AddressA_ID == 67790505 | Illumina450Manifest_clean$AddressB_ID == 67790505,
  'Color_Channel'
]
Color

# Create the final table with Red fluorescence, Green fluorescence, probe type and color channel
data <- data.frame(
  Red = as.numeric(Red['67790505', ]),
  Green = as.numeric(Green['67790505', ]),
  Type = Type,
  Color = Color
)

# Assign sample names as row names of the table
row.names(data) <- colnames(Red)

# Print the final table for Step 3
data


## -----------------------------------------------------------------------------------------------------------------------------
MSet.raw <- preprocessRaw(RGset)


## ----warning=FALSE, message=FALSE---------------------------------------------------------------------------------------------
# Estimate sample-specific QC medians (M and U) and plot them to identify outlier samples
qc <- getQC(MSet.raw)
plotQC(qc)

# Plot the intensity of negative control probes to evaluate the background noise level across samples
controlStripPlot(RGset, control="NEGATIVE")

# Calculate detection p-values for each probe across all samples
detP <- detectionP(RGset)

# Flag probes as 'TRUE' if their p-value is above 0.05, meaning they are NOT statistically significant
failed <- detP>0.05

# Examine the distribution of p-values and look at the logical matrix preview
summary(detP)
head(failed)

# Count the number of passed (FALSE) and failed (TRUE) probes for each sample
output <- apply(failed, 2, table)
output


## -----------------------------------------------------------------------------------------------------------------------------
# extracting Beta and M values
beta <- getBeta(MSet.raw)
M <- getM(MSet.raw)

# Create a data.frame
beta_df <- data.frame(beta)

#filtering the columns of the data frame based on the information contained in the dataframe 'targets'.
#control group (CTRL)
beta_df_control <- beta_df[, targets$Group=='CTRL']

# checking how many samples and probes we have
dim(beta_df_control)

#disease group (DIS)
beta_df_disease <- beta_df[, targets$Group=='DIS']

# checking how many samples and probes we have
dim(beta_df_disease)

# Calculate the mean Beta-value for each probe across samples within each group and plotting
mean_of_beta_c <- apply(beta_df_control, 1, mean, na.rm=TRUE)
mean_of_beta_d <- apply(beta_df_disease, 1, mean, na.rm=TRUE)
mean_of_beta_c_clean <- na.omit(mean_of_beta_c)
mean_of_beta_d_clean <- na.omit(mean_of_beta_d)

#density check
density_control <- density(mean_of_beta_c_clean, na.rm = TRUE)
density_disease <- density(mean_of_beta_d_clean, na.rm = TRUE)

#plot of the density of mean methylation values, dividing the samples in CTRL and DIS
plot(density_control, col='darkgreen')
lines(density_disease, col='purple')
legend("topright", legend = c("CTRL", "DIS"), col = c("darkgreen","purple"), lty = 1)


## ----fig.width=12, fig.height=8, warning=FALSE, message=FALSE-----------------------------------------------------------------
# Subset the manifest into Type I (dfI) and Type II (dfII) probes based on their chemistry
dfI <- Illumina450Manifest_clean[Illumina450Manifest_clean$Infinium_Design_Type == "I", ]
dfI <- droplevels(dfI)
dfII <- Illumina450Manifest_clean[Illumina450Manifest_clean$Infinium_Design_Type == "II", ]
dfII <- droplevels(dfII)

# Separation raw beta values according to probe chemistry (Type I vs Type II)
beta_I <- beta[rownames(beta) %in% dfI$IlmnID, ]
beta_II <- beta[rownames(beta) %in% dfII$IlmnID, ]

# Calculate mean, standard deviation and their densities for raw data
mean_of_beta_I <- apply(beta_I, 1, mean, na.rm = TRUE)
mean_of_beta_II <- apply(beta_II, 1, mean, na.rm = TRUE)
d_mean_of_beta_I <- density(mean_of_beta_I, na.rm = TRUE)
d_mean_of_beta_II <- density(mean_of_beta_II, na.rm = TRUE)

sd_of_beta_I <- apply(beta_I, 1, sd, na.rm = TRUE)
sd_of_beta_II <- apply(beta_II, 1, sd, na.rm = TRUE)
d_sd_of_beta_I <- density(sd_of_beta_I, na.rm = TRUE) 
d_sd_of_beta_II <- density(sd_of_beta_II, na.rm = TRUE)
 
# Apply Quantile Normalization to correct technical biases and standardize distributions
preprocessQuantile_results <- preprocessQuantile(RGset)

# Extract normalized beta values from their GenomicRatioSet
beta_preprocessQuantile <- getBeta(preprocessQuantile_results)

# Separate normalized beta values according to probe chemistry
beta_preprocessQuantile_I <- beta_preprocessQuantile[rownames(beta_preprocessQuantile) %in% dfI$IlmnID, ]
beta_preprocessQuantile_II <- beta_preprocessQuantile[rownames(beta_preprocessQuantile) %in% dfII$IlmnID, ]

# Calculate mean, standard deviation and their densities for normalized data
mean_of_beta_preprocessQuantile_I <- apply(beta_preprocessQuantile_I, 1, mean, na.rm = TRUE)
mean_of_beta_preprocessQuantile_II <- apply(beta_preprocessQuantile_II, 1, mean, na.rm = TRUE)
d_mean_of_beta_preprocessQuantile_I <- density(mean_of_beta_preprocessQuantile_I, na.rm = TRUE)
d_mean_of_beta_preprocessQuantile_II <- density(mean_of_beta_preprocessQuantile_II, na.rm = TRUE)

sd_of_beta_preprocessQuantile_I <- apply(beta_preprocessQuantile_I, 1, sd, na.rm = TRUE)
sd_of_beta_preprocessQuantile_II <- apply(beta_preprocessQuantile_II, 1, sd, na.rm = TRUE)
d_sd_of_beta_preprocessQuantile_I <- density(sd_of_beta_preprocessQuantile_I, na.rm = TRUE)
d_sd_of_beta_preprocessQuantile_II <- density(sd_of_beta_preprocessQuantile_II, na.rm = TRUE)

# Convert the 'Group' column to a factor for grouping and coloring purposes
targets$Group <- as.factor(targets$Group)

# Define a custom color palette for the groups (CTRL and DIS)
palette(c("lightpink", "lightgreen"))

# Set up the layout: 2 rows (Raw vs Normalized) and 3 columns (Mean, SD, Boxplot)
par(mfrow = c(2, 3))

# Panel 1: Density of raw mean beta values
plot(d_mean_of_beta_I, col = "darkgreen", main = "Raw Beta Mean", xlim = c(0, 1), ylim = c(0, 5))
lines(d_mean_of_beta_II, col = "purple")
legend("top", legend = c("Type I", "Type II"), col = c("darkgreen", "purple"), lty = 1, cex = 0.8)


# Panel 2: Density of raw standard deviation beta values
plot(d_sd_of_beta_I, col = "darkgreen", main = "Raw Beta SD", xlim = c(0, 0.6), ylim = c(0, 60))
lines(d_sd_of_beta_II, col = "purple")
legend("topright", legend = c("Type I", "Type II"), col = c("darkgreen", "purple"), lty = 1, cex = 0.8)

# Panel 3: Boxplot of raw beta values (colored by CTRL and DIS groups)
boxplot(beta, col = targets$Group, main = "Boxplot of raw Beta values", las = 2, cex.axis = 0.4)
legend("bottomleft", legend = levels(targets$Group), fill = palette(), cex = 0.8)

# Panel 4: Density of normalized mean beta values
plot(d_mean_of_beta_preprocessQuantile_I, col = "darkgreen", main = "Quantile Beta Mean", xlim = c(0, 1), ylim = c(0, 5))
lines(d_mean_of_beta_preprocessQuantile_II, col = "purple")
legend("top", legend = c("Type I", "Type II"), col = c("darkgreen", "purple"), lty = 1, cex = 0.8)


# Panel 5: Density of normalized standard deviation beta values
plot(d_sd_of_beta_preprocessQuantile_I, col = "darkgreen", main = "Quantile Beta SD", xlim = c(0, 0.6), ylim = c(0, 60))
lines(d_sd_of_beta_preprocessQuantile_II, col = "purple")
legend("topright", legend = c("Type I", "Type II"), col = c("darkgreen", "purple"), lty = 1, cex = 0.8)

# Panel 6: Boxplot of normalized beta values (colored by CTLR and DIS groups)
boxplot(beta_preprocessQuantile, col = targets$Group, main = "Boxplot of normalized Beta values", 
las = 2, cex.axis = 0.8)
legend("bottomleft", legend = levels(targets$Group), fill = palette(), cex = 0.8)

# Reset par setting default (1 row, 1 column) to avoid affecting subsequent plots in the markdown 
par(mfrow = c(1, 1))


## ----message=FALSE, warning=FALSE---------------------------------------------------------------------------------------------
# Perform a PCA on the matrix of normalized beta values generated in step 7 using the function prcomp().
pca_results <- prcomp(t(beta_preprocessQuantile), scale = TRUE)

# We use the fviz_eig function to create the scree plot
fviz_eig(pca_results, addlabels = TRUE, xlab = 'Principal Component', ylab = '% of variance', barfill = "lightblue", barcolor = "lightblue")

# PCA plot per Group
targets$Group <- as.factor(targets$Group)
group_colors <- c("green", "red")
plot(pca_results$x[, 1], pca_results$x[, 2], cex = 1, pch = 19, col = group_colors[targets$Group], xlim = c(-700, 700), ylim = c(-700, 700), xlab = "PC1", ylab = "PC2", main = 'PCA (Groups)')
text(pca_results$x[, 1], pca_results$x[, 2], labels = rownames(pca_results$x), cex = 0.4, pos = 3)
legend("bottomright", legend = levels(targets$Group), col = group_colors, pch = 19, cex = 1.0)

# PCA plot per Sex
targets$Sex <- as.factor(targets$Sex)
sex_colors <- c("purple", "orange")
plot(pca_results$x[, 1], pca_results$x[, 2], cex = 1, pch = 19, col = sex_colors[targets$Sex], xlim = c(-700, 700), ylim = c(-700, 700), xlab = "PC1", ylab = "PC2", main = 'PCA (Sex)')
text(pca_results$x[, 1], pca_results$x[, 2], labels = rownames(pca_results$x), cex = 0.4, pos = 3)
legend("bottomright", legend = levels(targets$Sex), col = sex_colors, pch = 19, cex = 1.0)

# PCA plot colored by batch
targets$Slide <- as.factor(targets$Slide)
batch_colors <- rainbow(nlevels(targets$Slide))
plot(pca_results$x[, 1], pca_results$x[, 2], cex = 1, pch = 19, col = batch_colors[targets$Slide], xlim = c(-700, 700), ylim = c(-700, 700), xlab = "PC1", ylab = "PC2", main = "PCA (Batch)")
text(pca_results$x[, 1], pca_results$x[, 2], labels = rownames(pca_results$x), cex = 0.4, pos = 3)
legend("bottomright", legend = levels(targets$Slide), col = batch_colors, pch = 19, cex = 1.0)


## -----------------------------------------------------------------------------------------------------------------------------
# Prepare the normalized matrix and check dimensions
beta_preprocessQuantile_step9 <- beta_preprocessQuantile
dim(beta_preprocessQuantile_step9)

# Define the categorical phenotype groups
targets$Group <- as.factor(targets$Group)
group_step9 <- targets$Group

# Safety checks before testing
stopifnot(ncol(beta_preprocessQuantile_step9) == nrow(targets))
stopifnot(length(group_step9) == ncol(beta_preprocessQuantile_step9))

# Remove probes containing missing values
beta_preprocessQuantile_step9_clean <- beta_preprocessQuantile_step9[complete.cases(beta_preprocessQuantile_step9), ]
dim(beta_preprocessQuantile_step9_clean)

# Identify differentially methylated probes
dmp_preprocessQuantile <- dmpFinder(beta_preprocessQuantile_step9_clean, pheno = group_step9, type = "categorical")

# Order the table by raw p-value and extract identifiers
dmp_preprocessQuantile <- dmp_preprocessQuantile[order(dmp_preprocessQuantile$pval), ]
dmp_preprocessQuantile$CpG <- rownames(dmp_preprocessQuantile)

# Calculate group specific means
mean_beta_CTRL_step9 <- rowMeans(beta_preprocessQuantile_step9_clean[, group_step9 == "CTRL", drop = FALSE], na.rm = TRUE)
mean_beta_DIS_step9  <- rowMeans(beta_preprocessQuantile_step9_clean[, group_step9 == "DIS", drop = FALSE], na.rm = TRUE)

# Append mean profiles to the output dataset
dmp_preprocessQuantile$mean_beta_CTRL <- mean_beta_CTRL_step9[dmp_preprocessQuantile$CpG]
dmp_preprocessQuantile$mean_beta_DIS  <- mean_beta_DIS_step9[dmp_preprocessQuantile$CpG]

# Calculate biological effect direction (Delta Beta)
dmp_preprocessQuantile$deltaBeta <- dmp_preprocessQuantile$mean_beta_DIS - dmp_preprocessQuantile$mean_beta_CTRL

# Optimize column sequence layout
reordered_cols <- c("CpG", "mean_beta_CTRL", "mean_beta_DIS", "deltaBeta")
dmp_preprocessQuantile <- dmp_preprocessQuantile[, c(reordered_cols, setdiff(colnames(dmp_preprocessQuantile), reordered_cols))]

# Isolate and print the top 20 significant features
top20_dmp_preprocessQuantile <- head(dmp_preprocessQuantile, 20)
top20_dmp_preprocessQuantile


## -----------------------------------------------------------------------------------------------------------------------------
# Calculate Bonferroni and Benjamini-Hochberg (BH) adjusted p-values
dmp_preprocessQuantile$pval_Bonferroni <- p.adjust(dmp_preprocessQuantile$pval, method = "bonferroni")
dmp_preprocessQuantile$pval_BH <- p.adjust(dmp_preprocessQuantile$pval, method = "BH")

# Set the significance threshold
threshold <- 0.05

# Count the number of significant probes for each method
sig_nominal <- sum(dmp_preprocessQuantile$pval < threshold)
sig_bonferroni <- sum(dmp_preprocessQuantile$pval_Bonferroni < threshold)
sig_BH <- sum(dmp_preprocessQuantile$pval_BH < threshold)

# Construct a summary table comparing the stringency of different multiple testing corrections
multiple_testing_summary <- data.frame(
  `Statistical Approach` = c("Nominal p-value (Unadjusted)", "Bonferroni Correction (FWER)", "Benjamini-Hochberg (BH/FDR)"),
  `Significant CpG Probes` = c(sig_nominal, sig_bonferroni, sig_BH),
  check.names = FALSE
)

# Display the final summary of differentially methylated loci across significance thresholds
multiple_testing_summary


## ----volcano-and-manhattan, fig.width=10, fig.height=5------------------------------------------------------------------------
# PRODUCE VOLCANO PLOT

# Isolate the normalized Beta values for each group and calculate their average profiles per probe
mean_ctrl <- apply(beta_preprocessQuantile_step9_clean[, targets$Group == "CTRL"], 1, mean, na.rm=TRUE)
mean_dis  <- apply(beta_preprocessQuantile_step9_clean[, targets$Group == "DIS"], 1, mean, na.rm=TRUE)

# Calculate the difference between the average methylation values (biological effect size)
delta <- mean_dis - mean_ctrl

# Match the probe order with the dmpFinder statistical output to ensure rows align correctly
delta_ordered <- delta[rownames(dmp_preprocessQuantile)]

# Create a dataframe with two columns: one for the delta values and the other for the -log10 of p-values
Volc_df <- data.frame(Delta = delta_ordered, NegLogP = -log10(dmp_preprocessQuantile$pval))

# Set layout to 1 row and 2 columns so the plots are side by side
par(mfrow = c(1, 2))

# Generate the base scatter plot for the Volcano Plot (cloud of grey dots)
plot(Volc_df$Delta, Volc_df$NegLogP, pch = 16, cex = 0.5, col = "grey70",
     xlab = "Delta Beta (DIS - CTRL)", ylab = "-log10(p-value)", main = "Volcano Plot")

# Add a horizontal significance threshold line at p-value = 0.05 and vertical threshold lines indicating a 10% difference in methylation (biological effect)
abline(h = -log10(0.05), col = "red", lty = 2)     # Significance threshold (p-value = 0.05)
abline(v = c(-0.1, 0.1), col = "blue", lty = 2)    # Biological threshold (Delta Beta = ±10%)

# Extract and isolate the statistically and biologically significant results and Highlight the significant probes on the plot by overlaying them in red
significant_features <- dmp_preprocessQuantile$pval < 0.05 & abs(delta_ordered) > 0.1
points(Volc_df$Delta[significant_features], Volc_df$NegLogP[significant_features], pch = 16, cex = 0.6, col = "red")


# PRODUCE MANHATTAN PLOT

# Merge statistical p-values with genomic coordinates from the manifest file
dmp_preprocessQuantile$IlmnID <- rownames(dmp_preprocessQuantile)
input_Manhattan <- merge(dmp_preprocessQuantile[, c("IlmnID", "pval")], Illumina450Manifest_clean[, c("IlmnID", "CHR", "MAPINFO")], by="IlmnID")

# Format and convert chromosomes (1-22, X, Y) into numbers for plot compatibility
order_chr <- c(1:22, "X", "Y")
input_Manhattan$CHR <- factor(input_Manhattan$CHR, levels=order_chr)
input_Manhattan$CHR <- as.numeric(input_Manhattan$CHR)

# Generate the Manhattan Plot and add a red significance threshold line
manhattan(input_Manhattan, snp="IlmnID", chr="CHR", bp="MAPINFO", p="pval",
          col=c("lightblue", "blue", "darkred", "green", "orange", "magenta", "darkgreen", "yellow"), 
          main="Manhattan Plot", suggestiveline = FALSE, genomewideline = FALSE)
abline(h = -log10(0.05), col = "red", lty = 2)

# Reset layout to default
par(mfrow = c (1,1))


## ----top100-heatmaps, fig.width=10, fig.height=7------------------------------------------------------------------------------
# 1. We extract the IDs of the top 100 most significant probes
top100_probes <- dmp_preprocessQuantile$CpG[1:100]

# 2. We subset the normalized beta matrix to keep ONLY these 100 probes
input_heatmap <- as.matrix(beta_preprocessQuantile_step9_clean[top100_probes, ])

# 3. We create dynamic color bars based on the targets dataframe (Group: CTRL vs DIS)
colorbar_group <- ifelse(targets$Group == "CTRL", "green", "orange")

# Distance between clusters is found with average linkage mapping (GROUPS)
heatmap.2(input_heatmap,
          col = colorRampPalette(c("navy", "white", "firebrick3"))(100),
          Rowv = TRUE,
          Colv = TRUE,
          hclustfun = function(x) hclust(x, method = 'average'),
          dendrogram = "both",
          key = TRUE,
          ColSideColors = colorbar_group,
          density.info = "none",
          trace = "none",
          scale = "none",
          symm = FALSE,
          main = "Average Linkage - Groups",
          cexCol = 0.8,
          labRow = FALSE) # Hides probe IDs to make the heatmap cleaner

# 4. We perform the analysis comparing sample using sex groups
# Robust check to handle both "Male"/"Female" and "M"/"F" formats in targets$Sex
colorbar_sex <- ifelse(targets$Sex %in% c("Male", "M"), "cyan",
                       ifelse(targets$Sex %in% c("Female", "F"), "pink", "grey"))

# Distance between clusters is found with average linkage mapping (SEX)
# Using the exact same color ramp palette to ensure scientific consistency
heatmap.2(input_heatmap,
          col = colorRampPalette(c("navy", "white", "firebrick3"))(100),
          Rowv = TRUE,
          Colv = TRUE,
          hclustfun = function(x) hclust(x, method = 'average'),
          dendrogram = "both",
          key = TRUE,
          ColSideColors = colorbar_sex,
          density.info = "none",
          trace = "none",
          scale = "none",
          symm = FALSE,
          main = "Average Linkage - Sex",
          cexCol = 0.8,
          labRow = FALSE)

