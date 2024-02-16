
import pandas as pd
import numpy as np

class ANNOT:
    def __init__(self, expr_filename=None, platform_filename=None, remove_redundant_genename=True):
        self.expr_filename = expr_filename
        self.platform_filename = platform_filename
        self.remove_redundant_genename = remove_redundant_genename
        
        
    def read_exprs(self):
        fdata = open(self.expr_filename).readlines()
        store_exprs = dict()
        redundant_genename_test = list()
        for i, lines in enumerate(fdata):
            lobj = [x.strip() for x in lines.split(",")]
            if i == 0:
                headerObj = lobj[1:]
                continue
                
            #geneName = self.annot[lobj[0]]
            #print(list(self.annot.keys()))
            
            try:
                geneName = self.annot[lobj[0]].upper()
            except:
                #print(lobj[0])
                continue
                
            if self.remove_redundant_genename:
                try:
                    redundant_genename_test.index(geneName)
                except:
                    store_exprs[geneName] = {sampleName:exp for sampleName, exp in zip(headerObj, lobj[1:])}
                    redundant_genename_test.append(geneName)

        
        retStore = pd.DataFrame.from_dict(store_exprs, orient='index')
        #print(store_exprs)
        retStore.to_csv(self.expr_filename.split(".")[0] + "_withGeneName.csv")
        
        return True
        
    def read_platform(self):
        fdata = open(self.platform_filename).readlines()
        self.annot = dict()
        for i, lines in enumerate(fdata):
            lobj = lines.split(",")
            if lobj[0].strip() == "":
                continue
                
            geneName = lobj[0].strip()
            if "/" in geneName:
                genes = [x.strip() for x in geneName.split("/") if x.strip() != ""]
            elif "," in geneName:
                genes = [x.strip() for x in geneName.split(",") if x.strip() != ""]
            else:
                genes = [geneName]
                
            for gene in genes:
                self.annot[lobj[1].strip()] = gene
            
        return True
        
        
    def add_geneName(self):
        self.exprs_df = pd.read_csv(self.expr_filename, index_col=0, header=0)
        self.plat = pd.read_csv(self.platform_filename, index_col=0, header=0)
        self.read_platform()
        self.read_exprs()



def rearrange_ILMN_special(filename):
    df = pd.read_csv(filename, index_col=0)

    df = df[["ID", "Symbol"]]
    gcol = "Symbol"

    rows = df.index
    with open(filename.split(".")[0] + "_new.csv", "w") as fp:
        fp.write("ID, gene\n")
        for head in rows:

            rowObj = df.loc[head]
            id = rowObj["ID"]
            geneName = rowObj[[gcol]].values[0]
            if pd.isna(geneName):
                continue

            if "///" in geneName:
                geneName = geneName.split("///")[0]

            elif "//" in geneName:
                geneName = geneName.split("//")[0]

            fp.write("%s, %s\n" % (geneName, id))

    return True
        
        
def rearrange_Affy(filename):
    df = pd.read_csv(filename, index_col=0)
    if "Gene Symbol" in df.columns.to_list():
        df = df[["ID", "Gene Symbol"]]
        gcol = "Gene Symbol"
    elif "gene_assignment" in df.columns.to_list():
        df = df[["ID", "gene_assignment"]]
        gcol = "gene_assignment"
        

    rows = df.index
    with open(filename.split(".")[0] + "_new.csv", "w") as fp:
        fp.write("ID, gene\n")
        for head in rows:

            rowObj = df.loc[head]
            id = rowObj["ID"]
            geneName = rowObj[[gcol]].values[0]
            if pd.isna(geneName):
                continue
            
            if "///" in geneName:
                geneName = geneName.split("///")[0]
            
            elif "//" in geneName:
                geneName = geneName.split("//")[0]
                
            fp.write("%s, %s\n" %(geneName, id))
            
   
    return True
    
def rearrange_Rosetta(filename):
    fdata = open(filename)
    with open(filename.split(".")[0] + "_new.csv", "w") as fp:
        fp.write("ID,gene\n")
        for i, lines in enumerate(fdata):
            lobj = lines.split(",")
            fp.write("%s,%s\n" %(lobj[3], lobj[1]))
            
    return True


def concat(folderName):
    import glob
    
    files = glob.glob(folderName + "/*.csv")
    store_pd = dict()
    with open("Sample_id_annotations.csv", "w") as fp:
        fp.write("Sample ID,Series Filename\n")
        for filename in files:
            store_pd[filename] = pd.read_csv(filename, index_col=0, header=0)
            for sname in store_pd[filename].columns:
                fp.write("%s,%s\n" %(sname, filename.split("/")[-1].rstrip(".csv")))
        
    result = pd.concat(list(store_pd.values()), axis=1, join="inner")
    result.to_csv("Combined.csv")
    print(result.shape)
    return True
    
        
