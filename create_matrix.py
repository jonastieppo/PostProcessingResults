#Given a nxn array, write the bmatrix in Tex format
# %%
import numpy as np
import pandas as pd
class NumptArrrayToBmatrix():

    def read_out_arquive(self, file):
        self.information=pd.read_csv(file,sep='\t',header=None)
    
    def print_outarquive(self):
        print(self.information)

    def WriteTexBmatrixFromNpArray(self,array):

        f = open("tex_arquive.tex", "w")
        f.write("\\begin{bmatrix}\n")
        def TransformArray():
            for each_row in range(6):
                for each_column in range(6):
                    if each_column<5:
                        f.write("{:8.2E} & ".format(array[each_row,each_column]))
                    else: 
                        f.write("{:8.2E} \\\\ \n".format(array[each_row,each_column]))

        TransformArray()
        f.write("\\end{bmatrix}\n")
        f.close()

    def WriteTexBmatrix(self):
        
        f = open("tex_arquive.tex", "w")
        f.write("\\begin{bmatrix}\n")
        def TransformArray():
            for each_row in self.information.index:
                for each_column in self.information.columns:
                    if each_column<max(self.information.columns):
                        f.write("{:8.2E} & ".format(self.information.iloc[each_row,each_column]))
                    else: 
                        f.write("{:8.2E} \\\\ \n".format(self.information.iloc[each_row,each_column]))

        TransformArray()
        f.write("\\end{bmatrix}\n")
        f.close()

# %%
# Class  = NumptArrrayToBmatrix()
# # %%
# import os
# curren_folder = os.getcwd()
# folder = curren_folder+r"\Plain Weave matrices\MBC"
# stiffness_arquive = folder+r"\stiffness.out"
# Class.read_out_arquive(stiffness_arquive)
# # %%
# Class.print_outarquive()
# Class.WriteTexBmatrix()
# # %%
# stiffness_arquive = folder+r"\compliance.out"
# Class.read_out_arquive(stiffness_arquive)
# Class.print_outarquive()
# Class.WriteTexBmatrix()
# %%
