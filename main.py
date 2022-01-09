
#%%
from postprocessing_files import *
import os
# %%
curren_folder = os.getcwd()

# DO NOT FORGET TO PUT IN VOIGT NOTATION

stiffness_arquive = "stiffess_PBC.out"
ClassInitiation = AnisoIndex(stiffness_arquive)
# A=ClassInitiation.CalculateCalculateEngineeringConstants()

# Au = ClassInitiation.aniso_index_user_material_simmetric
# Az = ClassInitiation.ZenerIndex
# Fr = ClassInitiation.Frobenius

# print("Norms: \n\n ++++++++++++")
# ClassInitiation.PrintPretty(Au)
# ClassInitiation.PrintPretty(Az)
# print(Fr)

# %%

# 
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

# ClassInitiation.PrincipalDirectionsOfAnisotropy(matrix)
# ClassInitiation.K

EngDictPlain = {

    'E1':24.8E3,
    'E2':24.8E3,
    'E3':8.5E3,
    'G12':4.2E3,
    'G23':4.2E3,
    'G13':4.2E3,
    'v12':0.1,
    'v13':0.28,
    'v23':0.28
}



# EngDictTwill = {

#     'E1':19.2e3,
#     'E2':19.2e3,
#     'E3':10.92e3,
#     'G12':3.92e3,
#     'G23':3.78e3,
#     'G13':3.78e3,
#     'v12':0.13,
#     'v13':0.305,
#     'v23':0.305
# }

# EngDictSatin = {

#     'E1':25.6e3,
#     'E2':25.6e3,
#     'E3':15.65e3,
#     'G12':5.67e3,
#     'G23':5.42e3,
#     'G13':5.42e3,
#     'v12':0.13,
#     'v13':0.283,
#     'v23':0.283
# }
# ClassInitiation.SetExperimentalProperties(EngDictPlain)
# # ClassInitiation.ZenerIndexExperimental
# matrix = ClassInitiation.experimental_matrix_stiffness
# K=ClassInitiation.BulkModulusTensor(matrix)
# print("K invariants:")
# PrintPretty(K.I1)
# PrintPretty(K.I2)
# PrintPretty(K.I3)
# L=ClassInitiation.DeviatoryModulusTensor(matrix)
# print("L invariants:")
# PrintPretty(L.I1)
# PrintPretty(L.I2)
# PrintPretty(L.I3)

'''
CLT Comparison
'''
# from CLT import*

# Eight_harness_experimental_CLT()

matrix = ClassInitiation.stiffness_matrix_simmetric

ClassContraction = TensorContracion(matrix)
c_new=ClassContraction.C_contract

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