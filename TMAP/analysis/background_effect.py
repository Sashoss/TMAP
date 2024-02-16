# First, install the necessary packages if they are not already installed

# Load the packages

import pandas as pd
import numpy as np
from sklearn.model_selection import GridSearchCV
from tqdm import tqdm


class BACKEFF:
    def __init__(self, filename, method="normexp", output=None):
        self.exprs_df = pd.read_csv(filename, index_col=0, header=0)
        exprs_df_normalized = self.normexp_convolution(self.exprs_df)
        if output != None:
            exprs_df_normalized.to_csv(output)
        
    
    def normexp_convolution(self, data):
        """
        Removes the background effect in microarray gene expression data samples using normal-exponential (normexp) convolution.
    
        Parameters:
        - data (pandas.DataFrame): The microarray gene expression data.
        - sigma (float): The standard deviation parameter of the normal distribution.
        - tau (float): The decay constant of the exponential distribution.
    
        Returns:
        - pandas.DataFrame: The microarray gene expression data with the background effect removed.
        
        # Load the microarray gene expression data into a pandas.DataFrame
        data = pd.read_csv("gene_expression_data.csv")

        # Remove the background effect using normexp convolution with sigma = 5 and tau = 0.1
        convolved_data = normexp_convolution(data, sigma=5, tau=0.1)

        """
        
        print("\nRunning background correction of %s samples\n" %(len(data.columns)))
        # Create a pandas.DataFrame to store the convolved data
        
    
        # Iterate over the columns (samples) in the data
        
        alpha_vals = [1, 2, 3, 4, 5]
        beta_vals = [0.1, 0.2, 0.3, 0.4, 0.5]
        least_mse = 10000
        cutoff_crossed = None
        
        pbar = tqdm(total=25*len(data.columns))
        
        for alpha in alpha_vals:
            for beta in beta_vals:
                convolved_data = data.copy()
                
                for col in data.columns:
                    # Extract the expression values for the current sample
                    mean = data[col].mean()
                    std = data[col].std()
                    convolved_data[col] =  alpha * (data[col] - mean) + (1 - alpha) * (beta * np.exp(-(data[col] - mean)**2 / (2 * std**2)))
                    pbar.update(1)
                
                #print(convolved_data)
                #check the mean squared difference between original and convolved data
                mse = ((convolved_data - data) ** 2).mean().mean()
                #print(sigma, tau, mse)
                
                if mse < least_mse:
                    least_mse_conv_data = convolved_data
                    least_mse = mse
                    
        pbar.close()
        cutoff_crossed = least_mse_conv_data
    
        return cutoff_crossed
        

    






        
    
    
    
