import pandas as pd
from sklearn.feature_selection import RFE
from sklearn.linear_model import LinearRegression


def select_important_genes(gene_expression_file, target_file, top_diff_genes_file, n_features_to_select, output):
    """Identifies the most important genes from the given gene expression data.

    Args:
        df: A Pandas DataFrame with columns containing the gene expression data and the target column.
        gene_expression_columns: The names of the columns containing the gene expression data.
        target_column: The name of the column containing the target variable.
        n_features_to_select: The number of features (genes) to select.

    Returns:
        A list of the most important genes.
    
    RFE citation: Guyon, I., Weston, J., Barnhill, S., & Vapnik, V., “Gene selection for cancer classification using support vector machines”, Mach. Learn., 46(1-3), 389–422, 2002.
    """
    
    # Extract the gene expression data and target variable as numpy arrays
    gene_expression_df = pd.read_csv(gene_expression_file, index_col=0, header=0).T
    fdata = open(top_diff_genes_file, "r")
    store_diffGene = list()
    for i in fdata:
        gname = i.strip()
        store_diffGene.append(gname)
        
    gene_expression_df = gene_expression_df[store_diffGene]
    
    #print(gene_expression_df.shape, "shape")
    
    #print(gene_expression_df[gene_expression_df.isnull().any(axis=1)], "NA row")
    target_df = pd.read_csv(target_file, index_col=0, header=0)
    X = gene_expression_df.values
    y = target_df.values
    
    #print(X.shape, y.shape)
    # Use RFE to identify the most important genes
    model = LinearRegression()
    rfe = RFE(model, n_features_to_select=n_features_to_select)
    rfe.fit(X, y)
    
    print(rfe.ranking_)

    # Extract the indices of the selected genes
    gene_indices = rfe.get_support(indices=True)

    # Extract the names of the selected genes from the DataFrame
    selected_genes = gene_expression_df.T.iloc[gene_indices].index
    
    print(selected_genes)
    selected_genes.to_csv(output)

    return True
    
def feature_importance_measure(gene_expression_file, target_file, top_diff_genes_file, n_features_to_select, output):
    import eli5
    from eli5.sklearn import PermutationImportance
    
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.ensemble import GradientBoostingClassifier
    
    
    gene_expression_df = pd.read_csv(gene_expression_file, index_col=0, header=0).T
    
    
    fdata = open(top_diff_genes_file, "r")
    store_diffGene = list()
    for i in fdata:
        gname = i.strip()
        store_diffGene.append(gname)
        
    gene_expression_df = gene_expression_df[store_diffGene]
    print(gene_expression_df)
    target_df = pd.read_csv(target_file, index_col=0, header=0)
    X = gene_expression_df.values
    y = target_df.values
    y = y.reshape(-1)
    
    model = GradientBoostingClassifier(n_estimators=1000, learning_rate=0.01, random_state=0)
    
    #model = RandomForestClassifier()
    model.fit(X, y)

    # Compute permutation importance
    permuter = PermutationImportance(model, random_state=42)
    permuter.fit(X, y)
    importances = permuter.feature_importances_
    print(importances)
    return True

