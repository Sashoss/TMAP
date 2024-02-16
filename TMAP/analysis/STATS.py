


import numpy as np
from scipy.linalg import eigvalsh

        
from sklearn.utils import resample
from sklearn.linear_model import lasso_path
from sklearn.metrics import mutual_info_score

from sklearn.covariance import GraphicalLasso

"""
Algorithm: Estimating a sparse gene network from gene expression data using factor analysis and graphical Lasso with bootstrap-based threshold selection

Input:
- Gene expression matrix X of size (n_samples, n_genes)
- Number of factors k to select
- Regularization parameter alpha for Graphical Lasso
- Number of bootstrap samples n_bootstraps for threshold selection

Output:
- Important nodes and edges in the sparse gene network

1. Preprocessing: Center and scale the gene expression matrix X to have zero mean and unit variance.

2. Compute the covariance matrix: Compute the sample covariance matrix of X using the formula Cov(X) = (X^T X)/(n_samples - 1).

3. Factor analysis: Apply factor analysis to select the optimal number of factors using BIC. The steps are as follows:
   a. Fit a factor analysis model to the covariance matrix using the specified number of factors k.
   b. Compute the log-likelihood of the model on the data.
   c. Compute the BIC for the model using the formula BIC = -2*log_likelihood + k*(k+1)*log(n_genes)/2.
   d. Repeat steps a-c for different values of k, and choose the value of k that minimizes the BIC.

4. Create new covariance matrix: Using the optimal number of factors selected in step 3, create a new covariance matrix using the formula Cov_new = L L^T + D, where L is the factor loadings matrix and D is a diagonal matrix with the unique variances of the factors.

5. Apply graphical Lasso with bootstrap-based threshold selection: Apply graphical Lasso to the new covariance matrix using the specified regularization parameter alpha. Then, estimate the optimal threshold for the sparse precision matrix using bootstrap-based stability selection. The steps are as follows:
   a. Fit the Graphical Lasso estimator to the new covariance matrix using the specified value of alpha.
   b. Compute the absolute values of the precision matrix.
   c. Generate n_bootstraps bootstrap samples from the data, and for each sample, fit the Graphical Lasso estimator to the sample and compute the stability scores of the edges as the sum of the mutual information values across the bootstrap samples.
   d. Compute the optimal threshold using stability selection, and obtain the important nodes and edges in the sparse gene network.
   e. pull the important nodes and edges using the optimal threshold.

6. Output the important nodes and edges in the sparse gene network.

Note: This algorithm is just an example, and the specific details and implementation may depend on the data and the specific methods used. It is important to carefully validate and evaluate the results of any analysis before drawing conclusions or making scientific claims.

"""


class STAT_MODELS:

    def __init__(self, filename):
        self.filenam= filename

    def factor_analysis(self, cov_mat):
        n, _ = cov_mat.shape

        # Calculate the eigenvalues and eigenvectors of the covariance matrix
        eigvals, eigvecs = eigvalsh(cov_mat)

        # Sort the eigenvalues and eigenvectors in descending order
        idx = np.argsort(eigvals)[::-1]
        eigvals = eigvals[idx]
        eigvecs = eigvecs[:, idx]

        # Calculate the BIC for different number of factors
        bic = []
        for k in range(1, n+1):
            # Calculate the factor loadings for k factors
            loadings = eigvecs[:, :k] * np.sqrt(eigvals[:k])

            # Calculate the residual covariance matrix
            residual = cov_mat - loadings @ loadings.T

            # Calculate the log-likelihood of the data
            log_likelihood = -n/2 * (np.sum(np.log(eigvals[:k])) + np.trace(np.linalg.inv(residual) @ (n-k) * np.eye(n)))

            # Calculate the number of parameters
            num_params = k*n - k*(k+1)/2

            # Calculate the BIC for k factors
            bic.append(-2*log_likelihood + num_params*np.log(n))

        # Find the optimal number of factors based on the minimum BIC
        opt_k = np.argmin(bic) + 1

        # Calculate the factor loadings and the residual covariance matrix for the optimal number of factors
        opt_loadings = eigvecs[:, :opt_k] * np.sqrt(eigvals[:opt_k])
        opt_residual = cov_mat - opt_loadings @ opt_loadings.T

        # Return the factor loadings and the BIC for different number of factors
        
        return opt_loadings, opt_residual
    
    def convert_to_covariance_matrix(self, loadings, residual):
        n, k = loadings.shape

        # Calculate the diagonal matrix of residual variances
        residual_variances = np.diag(residual)

        # Calculate the matrix of specific factors
        specific_factors = np.diag(residual_variances)

        # Calculate the new covariance matrix
        cov_new = loadings @ loadings.T + specific_factors

        return cov_new


    def stability_selection(self, cov_mat, n_samples, n_bootstraps, alpha):
        # Compute the mutual information between pairs of variables
        
        """
        Here, we first compute the mutual information matrix between pairs of variables, which is a measure of the statistical dependence between variables. We then use bootstrapping to generate n_bootstraps samples of size n_samples from the data and estimate the precision matrix for each sample using Lasso with the regularization parameter alpha. We then compute the stability scores of the edges as the sum of the mutual information values for the edges across the samples.

        Finally, we use the stability selection method to estimate the best threshold for the sparse precision matrix based on the stability scores. The threshold is chosen as the minimum value that selects at least one edge and that is above a decreasing sequence of percentiles of the stability scores. The best threshold is returned as the threshold with the highest number of selected edges.
        """
        
        mi_mat = np.zeros_like(cov_mat)
        for i in range(cov_mat.shape[0]):
            for j in range(i, cov_mat.shape[0]):
                mi = mutual_info_score(cov_mat[:, i], cov_mat[:, j])
                mi_mat[i, j] = mi
                mi_mat[j, i] = mi

        # Initialize the stability scores
        stability_scores = np.zeros_like(cov_mat)

        # Generate bootstrap samples and estimate the precision matrix
        for i in range(n_bootstraps):
            sample_indices = resample(range(cov_mat.shape[0]), n_samples)
            sample_cov_mat = cov_mat[sample_indices][:, sample_indices]
            sample_mi_mat = mi_mat[sample_indices][:, sample_indices]
            alphas, coefs, _ = lasso_path(sample_cov_mat, alpha=alpha, fit_intercept=False, copy_X=True)
            for j, alpha in enumerate(alphas):
                nonzero_indices = np.where(coefs[:, j] != 0)[0]
                for k in range(len(nonzero_indices)):
                    for l in range(k, len(nonzero_indices)):
                        stability_scores[nonzero_indices[k], nonzero_indices[l]] += sample_mi_mat[nonzero_indices[k], nonzero_indices[l]]

        # Compute the final stability scores
        stability_scores /= n_bootstraps

        # Compute the optimal threshold using the stability selection method
        threshold_scores = []
        for i in range(1, cov_mat.shape[0]):
            threshold = np.percentile(stability_scores, 100*(1-i/cov_mat.shape[0]))
            selected_indices = np.where(stability_scores > threshold)
            n_selected = len(selected_indices[0])
            if n_selected > 0:
                threshold_scores.append((threshold, n_selected))

        if len(threshold_scores) == 0:
            return None

        best_threshold, _ = max(threshold_scores, key=lambda x: x[1])
        return best_threshold


    def graphical_lasso(self, cov_new, threshold):
        
        # Convert the covariance matrix to a precision matrix using Graphical Lasso
        """
        Here, we first fit the Graphical Lasso estimator to the new covariance matrix cov_new, which converts it to a sparse precision matrix prec_mat. We then take the absolute values of the precision matrix and threshold it to obtain a sparse matrix sparse_prec_mat, where only the elements with absolute values greater than the threshold are set to 1.

        We extract the indices of the non-zero elements of the sparse precision matrix using the np.where() function, and use them to obtain the indices of the important nodes and edges in the network. We combine the two sets of indices using the union operator | to obtain the set of important nodes, and zip the two sets of indices to obtain the set of important edges.

        Note that the threshold value is chosen based on your requirements and the characteristics of your data. A higher threshold will result in a sparser network with fewer edges, while a lower threshold will result in a denser network with more edges. You may need to experiment with different threshold values to find a good balance between sparsity and connectivity in your network.

        """
        
        # Estimate the best threshold using stability selection
        n_samples = 100 # set the number of samples to use in bootstrapping
        n_bootstraps = 1000 # set the number of bootstraps to use in stability selection
        alpha = 0.01 # set the regularization parameter for Lasso
        best_threshold = stability_selection(cov_mat, n_samples, n_bootstraps, alpha)

        # Apply Graphical Lasso with the estimated threshold
        graphical_lasso = GraphicalLasso(alpha=best_threshold)
        graphical_lasso.fit(cov_mat)
        prec_mat = graphical_lasso.precision_

        # Get the absolute values of the precision matrix
        abs_prec_mat = np.abs(prec_mat)

        # Extract the indices of the non-zero elements of the sparse precision matrix
        nonzero_indices = np.where(abs_prec_mat > best_threshold)

        # Get the important nodes and edges from the sparse precision matrix
        important_nodes = list(set(nonzero_indices[0]) | set(nonzero_indices[1]))
        important_edges = [(i, j) for i, j in zip(nonzero_indices[0], nonzero_indices[1])]

        print("Important nodes:", important_nodes)
        print("Important edges:", important_edges)
        
        return important_edges, important_nodes
    



    

