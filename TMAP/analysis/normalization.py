import pandas as pd
import sys
from statsmodels.nonparametric.smoothers_lowess import lowess
from tqdm import tqdm


class BATCH_NORM:
    def __init__(self, filename, method="lowess", output=None):
        self.exprs_df = pd.read_csv(filename, index_col=0, header=0)
        if method == "lowess":
            exprs_df_normalized = self.lowess_normalize(self.exprs_df)
        elif method=="quantile":
            exprs_df_normalized = self.quantile_normalize(self.exprs_df)
        else:
            sys.exit("Incorrect normalization method called. Please use method == 'lowess' or method=='quantile'\n")
            
        if output != None:
            exprs_df_normalized.to_csv(output)
                        
    def progress_bar_wrapper(self, pbar):
        def wrapper(x):
            # Your code here
        
            # Update the progress bar
            pbar.update(1)
        
            return x
        return wrapper
        
    def lowess_normalize(self, df):
        print("Running Lowess normalization for %s genes from %s samples" %(len(df.index), len(df.columns)))
        df = df.T
        df1 = df.reset_index(drop=True)
        df_normalized = df1.copy()
        pbar = tqdm(total=len(df.columns))
        for column in df1.columns:
            print(df_normalized[column].shape, df1[column].shape, df1.index.shape)
            df_normalized[column] = lowess(df1[column], df1.index, frac=1.0/20, it=3)[:,1]
            pbar.update(1)
            
        pbar.close()
        df_normalized.index = df.index
        df_normalized.columns = df.columns
        return df_normalized.T

    def quantile_normalize(self, df):
        rank_mean = df.stack().groupby(df.rank(method='first').stack().astype(int)).mean()
        return df.rank(method='min').stack().astype(int).map(rank_mean).unstack()


