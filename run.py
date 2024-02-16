# function cal
from TMAP import *

#source ~/ml_plat/bin/activate


"""

#Step 1. Downloadingg data

# Make sure you add correct series name. Some series IDs are super series. Open the super seires page and look for correct Series and platform IDs
#GEOparser("GSE46903", "GPL6947")


#Step 2. Converting platform IDs to gene names in experession matrix


#2.1. For affymetrics annotation file

#Input Affymetrics platform id annotation file as input

rearrange_Affy("GPL6244.csv")

Note: It will return Error if the the GEO soft file lacks ID_REF instance. This will be fixed later.

#This will create a new file names GPL6244_new.csv. Above function pulls platform gene ID column and gene name colum, splits and removes any other information and creates a new file GPL6244_new.csv which has gene name as first column and gene platform id as second column. Use this as input when converting paltform id to gene ID in the downloaded expression matrix file using steps below.

ANNOT("GSE42495-GPL6244.csv", "GPL6244_new.csv").add_geneName()


#2.2. for platforms other than Affy,
1. open the downloaded GPL platform annotation file
2. keep ID and gene name columns and remove all other columns.
3. Run steps shown below

ANNOT("GSE46903-GPL6947.csv", "GPL6947.csv").add_geneName()

This will create a new file named GSE46903-GPL6947withGeneName.csv. This has gene names as row names.


#Step 3. Concatenating gene expression matrix files that has gene name as first column.

Add a folder name that contains all the gene expression matrix csv files that has gene name as first column, as shown below

concat("Macrophages")


#Step 4: Remove background effect

BACKEFF("Combined.csv", output="Background_effect_removed.csv")


#Step 5: Apply batch normalization

Apply Batch normalization using quantile or lowess normalization or both.



To apply quantile normalization use the code below
BATCH_NORM("Background_effect_removed.csv", output="Normalized_quantile.csv", method="quantile")

Then to apply lowess normalization use the code below
BATCH_NORM("Normalized_quantile.csv", output="Normalized_lowess.csv", method="lowess")



NOTE: VVI
Lowes normalization results may not be useful when running it on data that contains treatment and control grouped together in one file. Lowess normalized may forces the expression values to align to the regression slope constructed byboth data.



"""

#BACKEFF("Final_Merged_with_CSC_Endothelial.csv", output="Background_effect_removed.csv")
#BATCH_NORM("Background_effect_removed.csv", output="Normalized_quantile.csv", method="quantile")
#BATCH_NORM("Normalized_quantile.csv", output="Normalized_lowess.csv", method="lowess")


#concat("merge_input")



#BACKEFF("Final_Merged_with_CSC_Endothelial_1.csv", output="Background_effect_removed.csv")

#BATCH_NORM("Background_effect_removed.csv", output="Normalized_quantile.csv", method="quantile")
#BATCH_NORM("Normalized_quantile.csv", output="Normalized_lowess.csv", method="lowess")


#BATCH_NORM("Final_Merged_with_CSC_Endothelial_1.csv", output="Normalized_quantile.csv", method="quantile")

#feature_importance_measure(gene_expression_file="Normalized_quantile.csv", target_file="annotation_final_CSC_Endothelial.csv", top_diff_genes_file="features.txt", n_features_to_select=100, output="selected_genes.csv")




#stuff to do
# 1. Apply background corrections - Done
# 2. Normalization - Done
# 2. Apply feature selection to select n number of genes
# 3. perform dimensionality reduction
# 4. cluster the samples using reduced dimensions
# 5. report important markers for clusters
# 6. check how well cluster satisfies the assigned annotations
# 7. Run differential expression between the clusters
# 8.


#GEOparser("GSE183590", "GPL18573")
#ANNOT("GSE183590-GPL18573.csv", "GPL18573.csv").add_geneName()


#rearrange_ILMN_special("GPL6884.csv")





#rearrange_Rosetta("GPL2029.csv")
#ANNOT("GSE3299-GPL2029.csv", "GPL2029_new.csv").add_geneName()

#print("GSE20364")
#GEOparser("GSE20364", "GPL10054")
#rearrange_Affy("GPL96.csv")
#ANNOT("GSE68465-GPL96.csv", "GPL96_new.csv").add_geneName()

#concat("merge")
