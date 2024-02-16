import GEOparse
import pandas as pd
import logging

class GEOparser:
    def __init__(self, series_id, platform):
        self.id = series_id
        self.platform=platform
        self.pull_data()
        
    def pull_data(self):

        #logging.getLogger("GEOparse").setLevel(logging.ERROR)
        gse = GEOparse.get_GEO(self.id, annotate_gpl=True, destdir="./")
        gpl = GEOparse.get_GEO(self.platform, destdir="./")
        
        ma_data = gse.pivot_samples("VALUE")
        ma_data.to_csv(self.id+"-"+self.platform+".csv")
        
        ma_data_gpl = gse.gpls[self.platform].table

        ma_data_gpl.to_csv(self.platform+".csv")
        
    
    
        

