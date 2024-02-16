# TMAP


# Data Processing Pipeline for Gene Expression Analysis

This README outlines the steps required for processing gene expression data, from downloading the data to applying batch normalization. The pipeline is designed for use with various platforms and includes specific instructions for handling different types of annotation files.

## Step 1: Downloading Data

1. **Download the gene expression data using the `GEOparser` function.**
   - Ensure you input the correct series name. Some series IDs are part of super series; in such cases, open the super series page to locate the correct Series and Platform IDs.
   - Example: `GEOparser("GSE46903", "GPL6947")`

## Step 2: Converting Platform IDs to Gene Names in Expression Matrix

### 2.1 For Affymetrix Annotation File

- **Input:** Affymetrix platform ID annotation file.
- **Function:** `rearrange_Affy("GPL6244.csv")`
- **Note:** This process will return an error if the GEO soft file lacks an `ID_REF` instance. This issue will be addressed in future updates.

The function creates a new file named `GPL6244_new.csv`, which contains the gene name as the first column and the gene platform ID as the second column. This file should be used as input for converting platform ID to gene ID in the downloaded expression matrix file.

- **Next Step:** Use `ANNOT("GSE42495-GPL6244.csv", "GPL6244_new.csv").add_geneName()` to add gene names to your expression matrix.

### 2.2 For Platforms Other Than Affymetrix

1. Open the downloaded GPL platform annotation file.
2. Keep only the ID and gene name columns, removing all others.
3. Run the following steps to add gene names: `ANNOT("GSE46903-GPL6947.csv", "GPL6947.csv").add_geneName()`

This creates a new file named `GSE46903-GPL6947withGeneName.csv`, with gene names as row names.

## Step 3: Concatenating Gene Expression Matrix Files

- Add a folder name containing all the gene expression matrix CSV files with gene name as the first column.
  - Example: `concat("Macrophages")`

## Step 4: Removing Background Effect

- Use `BACKEFF("Combined.csv", output="Background_effect_removed.csv")` to remove the background effect from your data.

## Step 5: Applying Batch Normalization

Batch normalization can be applied using either quantile or lowess normalization, or both.

- **For Quantile Normalization:** `BATCH_NORM("Background_effect_removed.csv", output="Normalized_quantile.csv", method="quantile")`
- **For Lowess Normalization:** `BATCH_NORM("Normalized_quantile.csv", output="Normalized_lowess.csv", method="lowess")`

**Note:** Lowess normalization may not be suitable for data containing treatment and control groups together in one file, as it may force the expression values to align to the regression slope constructed by both datasets.

## Additional Steps

Following the main processing pipeline, additional steps may include background correction, normalization, feature selection, dimensionality reduction, clustering, marker identification, annotation validation, and differential expression analysis.

## Usage Notes

- Replace example file names and IDs with those relevant to your data.
- Ensure correct handling of specific platform annotation files for accurate gene name mapping.
