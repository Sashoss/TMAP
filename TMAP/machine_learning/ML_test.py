




import pandas as pd
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split

def train_model_with_pca(df, gene_expression_columns, target_column, n_components):
    """Trains a logistic regression model to classify cells as either endothelial or non-endothelial, using PCA to reduce the dimensionality of the data.

    Args:
        df: A Pandas DataFrame with columns containing the gene expression data and the target column.
        gene_expression_columns: The names of the columns containing the gene expression data.
        target_column: The name of the column containing the target variable (either 'endothelial' or 'non-endothelial').
        n_components: The number of components to keep after performing PCA.

    Returns:
        The trained logistic regression model.
        
    # Example usage
    df = pd.read_csv('microarray_samples.csv')
    model = train_model_with_pca(df, ['gene1', 'gene2', 'gene3'], 'cell_type', n_components=2)
    """
    
    # Extract the gene expression data and target variable as numpy arrays
    X = df[gene_expression_columns].values
    y = df[target_column].values

    # Split the data into training and testing sets
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

    # Perform PCA on the training data
    pca = PCA(n_components=n_components)
    X_train_pca = pca.fit_transform(X_train)

    # Transform the testing data using the same PCA object
    X_test_pca = pca.transform(X_test)

    # Train the logistic regression model
    model = LogisticRegression(solver='lbfgs')
    model.fit(X_train_pca, y_train)

    return model


