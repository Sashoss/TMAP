# First, install the necessary packages if they are not already installed

# Load the packages

import pandas as pd
import statsmodels.formula.api as smf

import numpy as np


class DIFF_EXP:
    def __init__(self, filename, output=None, method="linear_model"):
        self.df = pd.read_csv(filename, index_col=0, header=0)
        df_corrected = self.batch_effect()
        if output != None:
            df_corrected.to_csv(output)
            
        
    def batch_effect_lm(self):
        # Assume that the dataframe is called 'df' and has multiple columns with gene expression values,
        # and that the column names contain information about the batch that each sample belongs to

        # Create a new column in the dataframe with the batch information
        self.df["batch"] = self.df.columns.str.extract(r'(batch\d+)', expand=False)
        # Reshape the dataframe from wide to long format
        df_long = pd.melt(self.df, id_vars="batch", var_name="gene", value_name="expression")
        # Fit a linear model to the data with gene expression as the response variable and batch as a predictor
        lm = smf.ols("expression ~ batch", data=df_long).fit()
        
        df_long["batch_effect"] = df_long["expression"] - lm.predict(df_long[["batch"]])
        # Reshape the dataframe back to wide format
        df_corrected = df_long.pivot(index="gene", columns="batch", values="batch_effect")
        return df_corrected

        
    
    
    
