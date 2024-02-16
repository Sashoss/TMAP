
import pandas as pd
from sklearn.decomposition import FactorAnalysis, PCA

from sklearn_extensions.fastica_ import XCA


class REDUCE:
    def __init__(self, df)
    
        print("stuff")
    

    def apply_factor_analysis(self, df, gene_expression_column, n_factors):
        """Applies factor analysis to the gene expression data in the given DataFrame.

        Args:
            df: A Pandas DataFrame with a column containing the gene expression data.
            gene_expression_column: The name of the column containing the gene expression data.
            n_factors: The number of factors to extract.

        Returns:
            A tuple containing the factors and the BIC score for the fitted model.
            
        # Example usage
        df = pd.DataFrame({'gene': ['A', 'B', 'C', 'D'], 'expression': [[1, 2, 3], [4, 5, 6], [7, 8, 9], [10, 11, 12]]})
        factors, bic = apply_factor_analysis(df, 'expression', 2)
        print(factors)
        """
        # Extract the gene expression data as a 2D numpy array
        gene_expression_data = df[gene_expression_column].values

        # Fit a factor analysis model to the data
        model = FactorAnalysis(n_components=n_factors)
        model.fit(gene_expression_data)

        # Extract the factors and BIC score from the model
        factors = model.components_
        bic = model.bic(gene_expression_data)

        return factors, bic
        
        
    def apply_xca(df, gene_expression_column, n_components):
        """Applies extreme component analysis to the gene expression data in the given DataFrame.

        Args:
            df: A Pandas DataFrame with a column containing the gene expression data.
            gene_expression_column: The name of the column containing the gene expression data.
            n_components: The number of components to extract.

        Returns:
            A tuple containing the components and the BIC score for the fitted model.
            
        
        # Example usage
        df = pd.DataFrame({'gene': ['A', 'B', 'C', 'D'], 'expression': [[1, 2, 3], [4, 5, 6], [7, 8, 9], [10, 11, 12]]})
        components, bic = apply_xca(df, 'expression', 2)
        print(components)
        """
        # Extract the gene expression data as a 2D numpy array
        gene_expression_data = df[gene_expression_column].values

        # Fit an XCA model to the data
        model = XCA(n_components=n_components)
        model.fit(gene_expression_data)

        # Extract the components and BIC score from the model
        components = model.components_
        bic = model.bic(gene_expression_data)

        return components, bic
        
    def apply_pca(df, gene_expression_column, n_components):
        """Applies principal component analysis to the gene expression data in the given DataFrame.

        Args:
            df: A Pandas DataFrame with a column containing the gene expression data.
            gene_expression_column: The name of the column containing the gene expression data.
            n_components: The number of components to extract.

        Returns:
            A tuple containing the components and the explained variance ratios for the fitted model.
            
        df = pd.DataFrame({'gene': ['A', 'B', 'C', 'D'], 'expression': [[1, 2, 3], [4, 5, 6], [7, 8, 9], [10, 11, 12]]})
        components, explained_variance_ratios = apply_pca(df, 'expression', 2)
        print(components)
        """
        
        # Extract the gene expression data as a 2D numpy array
        gene_expression_data = df[gene_expression_column].values

        # Fit a PCA model to the data
        model = PCA(n_components=n_components)
        model.fit(gene_expression_data)

        # Extract the components and explained variance ratios from the model
        components = model.components_
        explained_variance_ratios = model.explained_variance_ratio_

        return components, explained_variance_ratios






