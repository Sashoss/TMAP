# First, install the necessary packages if they are not already installed

# Load the packages
import scipy.stats as stats
import statsmodels.stats.multitest as multitest
import pandas as pd
import statsmodels.formula.api as smf
import statsmodels.stats.multitest as smm
import numpy as np
import statsmodels.api as sm

class DIFF_EXP:
    def __init__(self, filename, output, method="linear_model"):
    
        self.exprs_df = pd.read_csv(filename, index_col=0, header=0)
        
        if method=="linear_model":
            expression_result = self.linear_model_diff_exp()
            expression_result_bayes = self.apply_empirical_bayes(expression_result, "padj")
            expression_result.to_csv(output)
            

    def apply_contrast_fit(self, exprs_df, contrast_vectors):
        """Applies contrast fit to the gene expression data using the given contrast vectors.

        Args:
            gene_expression_data: A 2D numpy array with rows representing genes and columns representing samples.
            contrast_vectors: A list of 1D numpy arrays representing the contrast vectors.

        Returns:
            A 2D numpy array with the same shape as the input gene expression data, but with the    contrast fit applied.
            
        # Example usage
        df = pd.DataFrame({'gene': ['A', 'B'], 'expression': [[1, 2, 3], [4, 5, 6]]})
        contrast_vectors = [np.array([1, -1, 0]), np.array([0, 1, -1])]
        df_contrast_fitted = apply_contrast_fit(df, 'expression', contrast_vectors)
        print(df_contrast_fitted)
        """
        # Compute the contrasts for each gene
        
        gene_expression_data = exprs_df.values

        contrasts = []
        for gene_expression in gene_expression_data:
            gene_contrasts = []
            for contrast_vector in contrast_vectors:
                gene_contrast = np.dot(gene_expression, contrast_vector)
                gene_contrasts.append(gene_contrast)
            contrasts.append(gene_contrasts)

        # Apply the contrast fit to the gene expression data
        contrast_fitted_data = gene_expression_data - np.array(contrasts).T
        df[gene_expression_column] = contrast_fitted_data
        return df
    
    def linear_model_diff_exp(self):
        # add a column group that contains group category of samples to compare
        
        self.exprs_df["group"] = pd.to_numeric(self.exprs_df["group"])
        
        groups = self.exprs_df["group"]
        exprs = self.exprs_df.drop("group", axis=1)

        
        #mean_expr_by_group = self.exprs_df.groupby("group", axis=1).mean()
        
        #mean_expr_by_group = exprs_df.groupby('group', as_index=False).mean()
        
        #mean_expr_by_group = exprs_df.groupby("group", axis=1).mean()

        # Initialize an empty list to store the results
        results = []

        # Iterate over the rows of the dataframe (i.e., the genes)
        for gene in exprs.index:
            # Fit a linear model to the data for the current gene with group as a predictor variable and mean expression level as the response variable
            """
            # Make sure to set reference group as 0 in the dataframe for this to use 0 group as reference
            """
            lm = smf.gls(f"exprs['{gene}'] ~ C(groups, Treatment(reference=0))", data=exprs).fit()
            # Extract the log fold change, p value, and B value from the model summary
            logfc = lm.params[1]
            pval = lm.pvalues[1]
            bval = lm.params[0]
            # Calculate the adjusted p value using the Benjamini-Hochberg method
            padj = smm.multipletests(pval, method="fdr_bh")[1][0]
            # Append the results to the list
            results.append([gene, logfc, pval, padj, bval])

        # Convert the list of results to a pandas dataframe
        df_results = pd.DataFrame(results, columns=["gene", "logfc", "pval", "padj", "bval", "eb"])
        
        differentially_expressed_genes = df_results.index[padj < 0.05]
        return df_results, differentially_expressed_genes
    
    

        


        
    
    
    
