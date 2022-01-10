
#%%
from postprocessing_files import *
import os
import unittest
# %%
curren_folder = os.getcwd()

# DO NOT FORGET TO PUT IN VOIGT NOTATION

results = os.getcwd()+r"\MBC"
# stiffness_arquive = fr"D:\Masters\Masters\Thesis\Text\Rato\Capitulo de Livro - Jonas dissertação\Visual Studio Code\PostProcessingResults\teste\stiffness.out"
ClassInitiation = AnisoIndex()
ClassInitiation.HomogenizeResults(results_folder=results)
ClassInitiation.ExportStiffness(folder="teste")
ClassInitiation.CalculateIndicators("stiffness.out")


Au = ClassInitiation.aniso_index_user_material_simmetric
Az = ClassInitiation.ZenerIndex
Fr = ClassInitiation.Frobenius

# print("Norms: \n\n ++++++++++++")
# ClassInitiation.PrintPretty(Au)
# ClassInitiation.PrintPretty(Az)
# print(Fr)

# %%

matrix = ClassInitiation.stiffness_matrix_simmetric
K=ClassInitiation.BulkModulusTensor(matrix)
print("K invariants:")
PrintPretty(K.I1)
PrintPretty(K.I2)
PrintPretty(K.I3)
L=ClassInitiation.DeviatoryModulusTensor(matrix)
print("L invariants:")
PrintPretty(L.I1)
PrintPretty(L.I2)
PrintPretty(L.I3)


def WriteTexBmatrix(matrix):
    
    f = open("tex_arquive.tex", "w")
    f.write("\\begin{bmatrix}\n")
    def TransformArray():
        for each_row in range(matrix.shape[0]):
            for each_column in range(matrix.shape[1]):
                if each_column<matrix.shape[1]-1:
                    # f.write("{:8.2E} & ".format(matrix[each_row,each_column]))
                    f.write("{:8.2f} & ".format(matrix[each_row,each_column]))
                else: 
                    # f.write("{:8.2E} \\\\ \n".format(matrix[each_row,each_column]))
                    f.write("{:8.2f} \\\\ \n".format(matrix[each_row,each_column]))

    TransformArray()
    f.write("\\end{bmatrix}\n")
    f.close()

C_new = TensorContracion(matrix)
C_new = C_new.C_contract
WriteTexBmatrix(C_new)


# WriteTexBmatrix(ClassInitiation.L)
# s_new = np.linalg.inv(c_new)

# print(s_new)
# WriteTexBmatrix(s_new)
# %%
from CLT import*

C_new=GeneralOrtho()
S_new = np.linalg.inv(C_new)
print(S_new)
WriteTexBmatrix(S_new)
# %%

S_CLT=S_anisotropic(Eng)
C_CLT = np.linalg.inv(S_CLT)
WriteTexBmatrix(S_CLT) 

# ABD=ABD_from_homogenization(Eng)

# WriteTexBmatrix(ABD)
# %%
'''
Calculating Invariants
'''


Invar2 = Invariants()
ClassInitiation.PrincipalDirectionsOfAnisotropy(ClassInitiation.stiffness_matrix_simmetric)
Invar2.matrix = ClassInitiation.K.matrix
Invar2.CalculateInvariants()
Invar2.PrinpipalAxys()
# print("Eig vec")
# print(Invar2.eig_val)
PrincipalTrans=Invar2.Principal

WriteTexBmatrix(PrincipalTrans)

# %%
'''
Comparing Tranformation Matrix
'''
ClassInitiation.PrincipalDirectionsOfAnisotropy(ClassInitiation.experimental_matrix_stiffness)
Invar = Invariants()
Invar.matrix = ClassInitiation.K.matrix
Invar.CalculateInvariants()
# Invar.PrinpipalAxys()

T1 = Invar.eig_vec

Invar2 = Invariants()
ClassInitiation.PrincipalDirectionsOfAnisotropy(ClassInitiation.stiffness_matrix_simmetric)
Invar2.matrix = ClassInitiation.K.matrix
Invar2.CalculateInvariants()
Invar2.PrinpipalAxys()
# Invar2.TransformationMatrix()
T2 = Invar2.T

# %%
Compare = CompareTransformationMatrix(T1,T2)

Compare.compare
# %%
Contraction = TensorContracion(tensor=ClassInitiation.stiffness_matrix_simmetric)
ContractedMatrix=Contraction.C_contract
WriteTexBmatrix(ContractedMatrix)