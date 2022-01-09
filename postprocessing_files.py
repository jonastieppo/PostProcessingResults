# Module to calculte the Kelving-Voigt bulk and shear modulus, by an input matrix,  then 
# calculate the general anysotropy index
# %%
from create_matrix import NumptArrrayToBmatrix
import numpy as np
import math as m
# from PostProcessing import CalculateCalculateEngineeringConstants
# from principal_directions_anisotropy import PrincipalAnisotropyDiretions
class Invariants():

    def __init__(self):

        self.matrix = np.zeros((3,3))
        self.eig_vec = None
        self.eig_val = None

    def CalculateI1(self):
        Matrix = self.matrix

        e,v = np.linalg.eig(Matrix)
        self.eig_vec = v
        self.eig_val = e

        print(v)
        self.I1=sum(e)

    def CalculateI2(self):
        
        Matrix = self.matrix
        e,v = np.linalg.eig(Matrix)

        self.I2=e[0]*e[1]+e[0]*e[2]+e[1]*e[2]

    def CalculateI3(self):

        Matrix = self.matrix
        e,v = np.linalg.eig(Matrix)

        self.I3=e[0]*e[1]*e[2]

    def CalculateInvariants(self):

        self.CalculateI1()
        self.CalculateI2()
        self.CalculateI3()

    def OrderEigenVectors(self):
        '''
        Method that will ordenated the eigen values
        '''
        #biding eigvectors and eigen values:

        zipped_values = zip(self.eig_val,self.eig_vec)

        #Sorting the eigvalues and eigvectors
        print("+++++++++++++++++\nEigen values and Eigenvection before ordening:")
        print("Eigen Values:")
        print(self.eig_val)
        print("Eigen Vectors:")
        print(self.eig_vec)
        sorted_e_val_e_vec = sorted(zipped_values,reverse=True)
        
        self.eig_val,self.eig_vec = zip(*sorted_e_val_e_vec)
        self.eig_val = list(self.eig_val)
        self.eig_vec = list(self.eig_vec)


    def TransformationMatrix(self):
        '''
        Defines the transformation matrix
        '''

        self.T = np.array([self.eig_vec[0],
                           self.eig_vec[1],
                           self.eig_vec[2]])

        # self.T=self.T.transpose()

    def PrinpipalAxys(self):
        '''
        Defines the principal axys of a matrix
        '''
        
        self.OrderEigenVectors()
        self.TransformationMatrix()
        self.Principal = np.matmul(np.matmul(self.T,self.matrix),self.T.transpose())

        return self.Principal


class CompareTransformationMatrix():

    def __init__(self,T1,T2):

        self.compare = np.multiply(T1,T2)

        for i in range(3):
            for j in range(3):
                self.compare[i,j] = np.rad2deg(np.arccos(self.compare[i,j]))


class TensorContracion():

    def __init__(self,tensor):
        
        self.C_contract = np.zeros((3,3))

        I = [0,1,5]
        a=0
        for i in  I:
            b=0
            for j in I:
                self.C_contract[a,b] = tensor[i,j]-tensor[i,2]*tensor[2,j]/tensor[2,2]
                b=b+1
            a=a+1
        
class AnisoIndex(NumptArrrayToBmatrix):


    def __init__(self,file):

        self.read_matrix(file)
        self.WriteVoigtNotationOfTheMatrix()
        self.OrthoMaterial()
        self.IsoTransMaterial()
        self.IsoMaterial()
        
        # Ortho
        self.Kv_ortho,self.Gv_ortho=self.Voigt(self.ortho_matrix_stiffness)
        self.Kr_ortho,self.Gr_ortho=self.Reuss(self.ortho_matrix_compliance)
        self.aniso_ortho=self.AnisoIndexCalculation(self.Gv_ortho,self.Gr_ortho,self.Kv_ortho,self.Kr_ortho)

        #Transversaly Isotropic

        self.Kv_trans_iso,self.Gv_trans_iso=self.Voigt(self.iso_trans_matrix_stiffness)
        self.Kr_trans_iso,self.Gr_trans_iso=self.Reuss(self.iso_trans_matrix_compliance)
        self.aniso_trans_iso=self.AnisoIndexCalculation(self.Gv_trans_iso,self.Gr_trans_iso,self.Kv_trans_iso,self.Kr_trans_iso)

        #Isotropic

        self.Kv_iso,self.Gv_iso=self.Voigt(self.iso_matrix_stiffness)
        self.Kr_iso,self.Gr_iso=self.Reuss(self.iso_matrix_compliance)
        self.aniso_iso=self.AnisoIndexCalculation(self.Gv_iso,self.Gr_iso,self.Kv_iso,self.Kr_iso)

        Kv,Gv = self.Voigt(self.stiffness_matrix)
        Kr,Gr = self.Reuss(self.compliance_matrix)

        self.aniso_index_user_material = self.AnisoIndexCalculation(Gv,Gr,Kv,Kr)

        # Calculating simmetric matrix

        Kv,Gv = self.Voigt(self.stiffness_matrix_simmetric)
        Kr,Gr = self.Reuss(self.compliance_matrix_simmetric)
        self.Kv = Kv
        self.Gv = Gv
        self.Kr = Kr
        self.Gr = Gr
        self.aniso_index_user_material_simmetric = self.AnisoIndexCalculation(Gv,Gr,Kv,Kr)
        self.ZenerIndexCalculation()
        self.FrobenniusNorm()
        # 

    def WriteTexBmatrixFromNpArray(self,array,name="tex_arquive"):
        return super().WriteTexBmatrixFromNpArray(array,name=name)

    def PrintPretty(string, format="{:8.3}"):
        print(f"{format}".format(string))

    def CalculateCalculateEngineeringConstants(self):
        '''
        Method to calculate the Engi;neering Constants
        '''
        # df = pd.DataFrame(ComplianceMatrix)

        ComplianceMatrix = self.compliance_matrix_original

        E1 = 1/ComplianceMatrix[0,0]
        E2 = 1/ComplianceMatrix[1,1]
        E3 = 1/ComplianceMatrix[2,2]
        G12 = 1/ComplianceMatrix[3,3]
        G23 = 1/ComplianceMatrix[4,4]
        G31 = 1/ComplianceMatrix[5,5]

        v21 = -E2*ComplianceMatrix[0,1]
        v31 = -E3*ComplianceMatrix[0,2]
        eta1_12 = G12*ComplianceMatrix[0,3]
        eta_23 = G23*ComplianceMatrix[0,4]
        eta1_31 = G31*ComplianceMatrix[0,5]

        v12 = -E1*ComplianceMatrix[1,0]
        v32 = -E3*ComplianceMatrix[1,2]
        eta2_12 = G12*ComplianceMatrix[1,3]
        eta2_23 = G23*ComplianceMatrix[1,4]
        eta2_31 = G31*ComplianceMatrix[1,5]

        v13 = -E1*ComplianceMatrix[2,0]
        v23 = -E2*ComplianceMatrix[2,1]
        eta3_12 = G12*ComplianceMatrix[2,3]
        eta3_23 = G23*ComplianceMatrix[2,4]
        eta3_31 = G31*ComplianceMatrix[2,5]

        eta12_1 = E1*ComplianceMatrix[3,0]
        eta12_2 = E2*ComplianceMatrix[3,1]
        eta12_3 = E3*ComplianceMatrix[3,2]
        mi12_23 = G23*ComplianceMatrix[3,4]
        mi12_31 = G31*ComplianceMatrix[3,5]

        eta23_1 = E1*ComplianceMatrix[4,0]
        eta23_2 = E2*ComplianceMatrix[4,1]
        eta23_3 = E3*ComplianceMatrix[4,2]
        mi23_12 = G12*ComplianceMatrix[4,3]
        mi23_31 = G23*ComplianceMatrix[4,5]

        eta31_1 = E1*ComplianceMatrix[5,0]
        eta31_2 = E2*ComplianceMatrix[5,1]
        eta31_3 = E3*ComplianceMatrix[5,2]
        mu31_12 = G12*ComplianceMatrix[5,3]
        mi31_23 = G23*ComplianceMatrix[5,4]

        EngineeringConstants = {
            "E1": E1,
            "E2": E2,
            "E3" : E3,
            "G12" : G12,
            "G23" : G23,
            "G31" : G31,
            "v21" : v21,
            "v31" : v31,
            "eta1_12" : eta1_12,
            "eta_23" : eta_23,
            "eta1_31" : eta1_31,
            "v12" : v12,
            'v32' : v32,
            'eta2_12' : eta2_12,
            'eta2_23' : eta2_23,
            'eta2_31' : eta2_31,
            'v13' : v13,
            'v23' : v23,
            'eta3_12' : eta3_12,
            'eta3_23' : eta3_23,
            'eta3_31' : eta3_31,
            'eta12_1' : eta12_1,
            'eta12_2' : eta12_2,
            'eta12_3' : eta12_3,
            'mi12_23' : mi12_23,
            'mi12_31' : mi12_31,
            'eta23_1' : eta23_1,
            'eta23_2' : eta23_2,
            'eta23_3' : eta23_3,
            'mi23_12' : mi23_12,
            'mi23_31' : mi23_31,
            'eta31_1' : eta31_1,
            'eta31_2' : eta31_2,
            'eta31_3' : eta31_3,
            'mu31_12' : mu31_12,
            'mi31_23' : mi31_23
        }
        
        self.EngineerinConstants = EngineeringConstants
        return EngineeringConstants


    def TransformToNumpy(self):
        self.stiffness_matrix = np.zeros((6,6))
        for each_row in self.information.index:
            for each_column in self.information.columns:
                self.stiffness_matrix[each_row,each_column] = self.information.iloc[each_row,each_column]


        # self.stiffness_matrix_simmetric = self.TurnSimm(self.stiffness_matrix)

        self.compliance_matrix = np.linalg.inv(self.stiffness_matrix)

        A = np.zeros((6,6))

        for i in range(6):
            for j in range(6):
                A[i,j]=self.compliance_matrix[i,j]

        self.compliance_matrix_original = A
        self.stiffness_matrix_original =np.linalg.inv(A)

        # self.WriteVoigtNotationOfTheMatrix()
        self.stiffness_matrix = np.linalg.inv(self.compliance_matrix)
        self.stiffness_matrix_simmetric = self.TurnSimm(self.stiffness_matrix)
        self.compliance_matrix_simmetric = np.linalg.inv(self.stiffness_matrix_simmetric)
        
    def read_matrix(self,file):
        # Reads a matrix and stores as an 6x6 numpy array
        ClassInit = NumptArrrayToBmatrix()
        ClassInit.read_out_arquive(file)
        self.information=ClassInit.information

        self.TransformToNumpy()

    def WriteVoigtNotationOfTheMatrix(self):
        '''
        Write the .out results in the voigt notation
        '''

        C = self.stiffness_matrix_simmetric 
        C_4_16 = [i for i in C[3,0:6]]
        C_6_16 = [i for i in C[5,0:6]]

        C[3,0:6]=C_6_16
        C[5,0:6]=C_4_16

        # Swapping fourth column with sixth column
        temp = [i for i in C[0:6,3]]
        C[0:6,3]=[i for i in C[0:6,5]]
        C[0:6,5]=temp

        # Swapping fourth line with fifth line
        temp = [i for i in C[3,0:6]]
        C[3,0:6]=[i for i in C[4,0:6]]
        C[4,0:6]=temp

        # Swapping fourth column with fifth column
        temp = [i for i in C[0:6,3]]
        C[0:6,3]=[i for i in C[0:6,4]]
        C[0:6,4]=temp

        self.stiffness_matrix_simmetric = C
        self.compliance_matrix_simmetric = np.linalg.inv(self.stiffness_matrix_simmetric)

        return None


    def Voigt(self,matrix):
        '''
        Depicts the linear relation between isotropic shear Modulus and Bulk Modulus
        '''

        C11 = matrix[0,0]
        C22 = matrix[1,1]
        C33 = matrix[1,1]
        C12 = matrix[0,1]
        C23 = matrix[1,2]
        C13 = matrix[0,2]
        C44 = matrix[3,3]
        C55 = matrix[4,4]
        C66 = matrix[5,5]

        Kv = ((C11+C22+C33)+2*(C12+C23+C13))/9
        Gv = ((C11+C22+C33)-(C12+C23+C13)+3*(C44+C55+C66))/15

        return Kv,Gv


    def IsoMaterial(self):

        E1 = 100000
        E2 = 100000
        E3 = 100000
        v13 = 0.35
        v12 = 0.35
        v23 = 0.35
        G12 = E1/(2*(1+v12))
        G23 = G12
        G13 = G12


        S11 = 1/E1
        S22 = 1/E2
        S33 = 1/E3
        S44 = 1/G12
        S55 = 1/G23
        S66 =1/G13
        S12 = -v12/E1
        S13 = -v13/E1
        S23 = -v23/E2

        self.iso_matrix_compliance = np.zeros((6,6))

        self.iso_matrix_compliance[0,0] = S11
        self.iso_matrix_compliance[1,1] = S22
        self.iso_matrix_compliance[2,2] = S33
        self.iso_matrix_compliance[3,3] = S44
        self.iso_matrix_compliance[4,4] = S55
        self.iso_matrix_compliance[5,5] = S66
        self.iso_matrix_compliance[0,1] = S12
        self.iso_matrix_compliance[0,2] = S13
        self.iso_matrix_compliance[1,2] = S23

        self.iso_matrix_compliance[1,0] = S12
        self.iso_matrix_compliance[2,0] = S13
        self.iso_matrix_compliance[2,1] = S23

        self.iso_matrix_stiffness = np.linalg.inv(self.iso_matrix_compliance)


    def IsoTransMaterial(self):

        E1 = 24.8E3
        E3 = 8.5E3
        v13 = 0.28
        v12 = 0.1
        G13 = 4.2E3

        S11 = 1/E1
        S22 = 1/E1
        S33 = 1/E3
        S44 = 2*(1+v12)/E1
        S55 = 1/G13
        S66 =1/G13
        S12 = -v12/E1
        S13 = -v13/E3
        S23 = -v13/E3

        self.iso_trans_matrix_compliance = np.zeros((6,6))

        self.iso_trans_matrix_compliance[0,0] = S11
        self.iso_trans_matrix_compliance[1,1] = S22
        self.iso_trans_matrix_compliance[2,2] = S33
        self.iso_trans_matrix_compliance[3,3] = S44
        self.iso_trans_matrix_compliance[4,4] = S55
        self.iso_trans_matrix_compliance[5,5] = S66
        self.iso_trans_matrix_compliance[0,1] = S12
        self.iso_trans_matrix_compliance[0,2] = S13
        self.iso_trans_matrix_compliance[1,2] = S23

        self.iso_trans_matrix_compliance[1,0] = S12
        self.iso_trans_matrix_compliance[2,0] = S13
        self.iso_trans_matrix_compliance[2,1] = S23

        self.iso_trans_matrix_stiffness = np.linalg.inv(self.iso_trans_matrix_compliance)

    def PlainWeaveSida(self):

        E1 = 24.8E3
        E2 = 24.8E3
        E3 = 8.5E3
        v13 = 0.28
        v12 = 0.1
        v23 = 0.28
        G12 = 6.5E3
        G23 = 4.2E3
        G13 = 4.2E3


        S11 = 1/E1
        S22 = 1/E2
        S33 = 1/E3
        S44 = 1/G12
        S55 = 1/G23
        S66 =1/G13
        S12 = -v12/E1
        S13 = -v13/E1
        S23 = -v23/E2


        self.plainweave_matrix_compliance = np.zeros((6,6))

        self.plainweave_matrix_compliance[0,0] = S11
        self.plainweave_matrix_compliance[1,1] = S22
        self.plainweave_matrix_compliance[2,2] = S33
        self.plainweave_matrix_compliance[3,3] = S44
        self.plainweave_matrix_compliance[4,4] = S55
        self.plainweave_matrix_compliance[5,5] = S66
        self.plainweave_matrix_compliance[0,1] = S12
        self.plainweave_matrix_compliance[0,2] = S13
        self.plainweave_matrix_compliance[1,2] = S23

        self.plainweave_matrix_compliance[1,0] = S12
        self.plainweave_matrix_compliance[2,0] = S13
        self.plainweave_matrix_compliance[2,1] = S23

        self.plainweave_matrix_stiffness = np.linalg.inv(self.plainweave_matrix_compliance)

        self.stiffness_matrix_simmetric = self.plainweave_matrix_stiffness

        Kv,Gv = self.Voigt(self.plainweave_matrix_stiffness)
        Kr,Gr = self.Reuss(self.plainweave_matrix_compliance)

        self.scida_plainweave_aniso_index = self.AnisoIndexCalculation(Gv,Gr,Kv,Kr)

    def OrthoMaterial(self):

        E1 = 10000
        E2 = 5000
        E3 = 1000
        v13 = 0.3
        v12 = 0.35
        v23 = 0.32
        G12 = 900
        G23 = 800
        G13 = 700


        S11 = 1/E1
        S22 = 1/E2
        S33 = 1/E3
        S44 = 1/G12
        S55 = 1/G23
        S66 =1/G13
        S12 = -v12/E1
        S13 = -v13/E1
        S23 = -v23/E2

        self.ortho_matrix_compliance = np.zeros((6,6))

        self.ortho_matrix_compliance[0,0] = S11
        self.ortho_matrix_compliance[1,1] = S22
        self.ortho_matrix_compliance[2,2] = S33
        self.ortho_matrix_compliance[3,3] = S44
        self.ortho_matrix_compliance[4,4] = S55
        self.ortho_matrix_compliance[5,5] = S66
        self.ortho_matrix_compliance[0,1] = S12
        self.ortho_matrix_compliance[0,2] = S13
        self.ortho_matrix_compliance[1,2] = S23

        self.ortho_matrix_compliance[1,0] = S12
        self.ortho_matrix_compliance[2,0] = S13
        self.ortho_matrix_compliance[2,1] = S23

        self.ortho_matrix_stiffness = np.linalg.inv(self.ortho_matrix_compliance)

    def Reuss(self,matrix):
        '''
        Reuss relations for K and G
        '''

        S11 = matrix[0,0]
        S22 = matrix[1,1]
        S33 = matrix[2,2]
        S44 = matrix[3,3]
        S55 = matrix[4,4]
        S66 = matrix[5,5]
        S12 = matrix[0,1]
        S13 = matrix[0,2]
        S23 = matrix[1,2]

        Kr = 1/((S11+S22+S33)+2*(S12+S23+S13))
        Gr = 4*(S11+S22+S33)-4*(S12+S23+S13)+3*(S44+S55+S66)
        Gr = 15/Gr

        return Kr,Gr

    def TurnSimm(self,matrix):
        '''
        Turns a matrix simmetric
        '''
        return (matrix+np.transpose(matrix))/2

    def AnisoIndexCalculation(self,Gv,Gr,Kv,Kr):

        return 5*Gv/Gr+Kv/Kr-6

    def PrincipalDirectionsOfAnisotropy(self,matrix):
        '''
        Returns the eign values in a descending order
        '''

        # For volumetric part

        self.K = self.BulkModulusTensor(matrix)
        self.K.eig_val,self.K.eig_vec = np.linalg.eig(self.K.matrix)

        ev_list = list(zip(self.K.eig_val, self.K.eig_vec))

        self.ev_list = ev_list
        try:
            ev_list.sort(reverse=True)
        except:
            pass

        e,v=zip(*ev_list)
        v = list(v)
        e = list (e)

        v = np.array([v[0]])
        self.K.eig_vec = v # sorting the eign values in descending order
        self.K.eig_val = e

        # For the deviatoric Part

        self.L = self.DeviatoryModulusTensor(matrix)
        self.L.eig_val,self.L.eig_vec = np.linalg.eig(self.L.matrix)

        ev_list = list(zip(self.L.eig_val, self.L.eig_vec))

        self.ev_list = ev_list
        try:
            ev_list.sort(reverse=True)
        except:
            pass

        e,v=zip(*ev_list)
        v = list(v)
        e = list (e)
        
        self.L.eig_vec = v # sorting the eign values in descending order
        self.L.eig_val = e

        return None

    def SetTransformationMatrix(self):
        '''
        Builds the transformation matrix, ordering the eingvector
        '''
        v1 = self.K.eig_vec[0].transpose()
        v2 = self.K.eig_vec[1].transpose()
        v3 = self.K.eig_vec[2].transpose()

        self.TransformationMatrix = np.array([v1,v2,v3])

    # def WriteVoigtNotationOfTheMatrix(self,alternative_matrix='None'):
    #     '''
    #     Writes the voigt notation of the read matrix from Finite Element Analysis
    #     '''
    #     if alternative_matrix=='None':
    #         S = self.compliance_matrix
    #     else:
    #         S = alternative_matrix
    #     # S = np.random.randint(100, size=(6, 6))

    #     # print("Old C")
    #     # print(S)
        
    #     self.old_compliance = np.zeros((6,6))

    #     for i in range(6):
    #         for j in range(6):
    #             self.old_compliance[i,j] = S[i,j]

    #     trans1 = np.zeros((6,6))
    #     trans2 = np.zeros((6,6))

    #     def ChangePositions(a,b,c,d):
    #         '''
    #         Change the position of the element C[a,b] with the element
    #         C[c,d]
    #         '''
    #         aux = S[c,d]
    #         # print(aux)
    #         # print(C[a,b])
    #         S[c,d]=S[a,b]
    #         S[a,b]=aux
        
    #     # Changing xy->xz

    #     lines_to_change = [0,1,2]

    #     for line in lines_to_change:
    #         ChangePositions(line,5,line,3)

    #     columns_to_change = [0,1,2,3,4,5]

    #     for column in columns_to_change:
    #         ChangePositions(3,column,5,column)

    #     ChangePositions(3,5,3,3)
    #     ChangePositions(5,5,5,3)

    #     for i in range(6):
    #         for j in range(6):
    #             trans1[i,j]=S[i,j]

    #     # Changing xz->yz

    #     lines_to_change = [0,1,2]

    #     for line in lines_to_change:
    #         ChangePositions(line,4,line,3)
        
    #     columns_to_change = [0,1,2,3,4,5]

    #     for column in columns_to_change:
    #         ChangePositions(3,column,4,column)
        
    #     ChangePositions(3,4,3,3)
    #     ChangePositions(4,4,4,3)

    #     for i in range(6):
    #         for j in range(6):
    #             trans2[i,j]=S[i,j]

    #     self.compliance_matrix = S
    #     # print("xy->xz")
    #     # print(trans1)
    #     # print("xz->yz")
    #     # print(trans2)

    def IndexTranformation(self,i,j):

        # Swapping index if necessary: 2,1 ->1,2

        if i>j:
            aux=j
            j=i
            i=aux

        Struct = {
                '00':0,
                '11':1,
                '22':2,
                '01':3,
                '12':4,
                '02':5,
        }

        return Struct[f"{i}{j}"]
    
    def BulkModulusTensor(self,matrix):

        K = Invariants()
        K_matrix = K.matrix

        # M = self.ortho_matrix_stiffness
        M = matrix

        for i in range(3):
            for j in range(i,3):
                first_index = self.IndexTranformation(i,j)
                # print(f"(i,j)={i},{j}->{first_index}")
                K_matrix[i,j] = M[first_index,0]+M[first_index,1]+M[first_index,2]
                K_matrix[j,i]=K_matrix[i,j]
                   
        K.matrix = K_matrix
        K.CalculateInvariants()

        return K

    def DeviatoryModulusTensor(self,matrix):

        '''
        Matrix is a 6x6 matrix 
        '''

        L = Invariants()
        L_matrix = L.matrix
        M = matrix

        for i in range(3):
            for j in range(i,3):
                A = []
                for k in range(3):
                    m = self.IndexTranformation(i,k) #transforming using voigt
                    n =  self.IndexTranformation(j,k) #transforming using voigt
                    A.append(M[m,n])
                # print(f"(i,j)={i},{j}->{first_index}")
                L_matrix[i,j] = sum(A)
                L_matrix[j,i]=L_matrix[i,j]

        L.matrix = L_matrix
        L.CalculateInvariants()

        return L

    def ZenerIndexCalculation(self,**kwargs):
        '''
        Calculates the zener index of a matrix
        '''
        C44 = self.stiffness_matrix_simmetric[3,3]
        C11 = self.stiffness_matrix_simmetric[0,0]
        C12 = self.stiffness_matrix_simmetric[0,1]

        if kwargs.get('experimental_matrix'):

            C44 = self.experimental_matrix_stiffness[3,3]
            C11 = self.experimental_matrix_stiffness[0,0]
            C12 = self.experimental_matrix_stiffness[0,1]
            self.ZenerIndexExperimental= 2*C44/(C11-C12)

        self.ZenerIndex = 2*C44/(C11-C12)

    def ClosestIsotropicDistance(self):

        eig,vec = np.linalg.eig(self.stiffness_matrix_simmetric)
        # print(eig)
        # i = np.array(eig)

        i = np.array([1,2,3,4])

        # print(type(i))
        n = i.shape[0]

        def Product(members):
            prod=1

            for element in members:
                prod = prod*element
            return prod

        def Summatory(members):
            sum=0
            for element in members:
                sum = sum+element

            return sum
        
        def NormalizedEignVector(eig_i):
            print(eig_i)
            return 2*self.Kr*eig_i

        def ProductMember(i):
            Product_member = []
            for element in range(i.shape[0]):
                Product_member.append(NormalizedEignVector(i[element]))

            return Product_member

        def Summatory_expression(eig_i,i):
            n = i.shape[0]

            LambdaNormal = (NormalizedEignVector(eig_i))**(-n)
            print(LambdaNormal)
            Produtory = Product(ProductMember(i))

            return (np.log(LambdaNormal*Produtory))**2

        def SummatoryMembers(i):
            member_i=[]
            n = i.shape[0]
            # Summatory_expression = np.log()
            for each_member in range(n):
                eig_i = each_member
                member_i.append(Summatory_expression(eig_i,i))

            return member_i

        expression = 1/n*(Summatory(SummatoryMembers(i)))**(0.5)

        return (Summatory(SummatoryMembers(i)))

    def FrobenniusNorm(self,**kwargs):

        x = self.stiffness_matrix_simmetric
        Fr = np.linalg.norm(x,ord='fro')

        if kwargs.get('experimental_matrix'):
            x = self.experimental_matrix_stiffness
            Fr = np.linalg.norm(x,ord='fro')
            self.Fr_experimental =  Fr

        self.Frobenius = Fr
        return Fr

    def CompositeStiffnessMatrix(self,StrandProps,ResinProps):
        '''
        Method to create the stiffness matrix for the strand and for the matrix
        of a composite. StrandProps and ResinProps will be both dictionaries

        The StransProps and ResinProps will have followin kwords:

        E1 v12 G12
        E2 v23 G23
        E3 v13 G13

        '''

        E1 = StrandProps['E1']
        E2 = StrandProps['E2']
        E3 = StrandProps['E3']
        G12 = StrandProps['G12']
        G23 = StrandProps['G23']
        G13 = StrandProps['G13']
        v12 = StrandProps['v12']
        v13 = StrandProps['v13']
        v23 = StrandProps['v23']

        Em = ResinProps['E1']
        Vm = ResinProps['v12']

        S11 = 1/E1
        S22 = 1/E2
        S33 = 1/E3
        S44 = 1/G12
        S55 = 1/G23
        S66 =1/G13
        S12 = -v12/E1
        S13 = -v13/E1
        S23 = -v23/E2

        # Strand Compliance Matrix

        self.strand_matrix_compliance = np.zeros((6,6))

        self.strand_matrix_compliance[0,0] = S11
        self.strand_matrix_compliance[1,1] = S22
        self.strand_matrix_compliance[2,2] = S33
        self.strand_matrix_compliance[3,3] = S44
        self.strand_matrix_compliance[4,4] = S55
        self.strand_matrix_compliance[5,5] = S66
        self.strand_matrix_compliance[0,1] = S12
        self.strand_matrix_compliance[0,2] = S13
        self.strand_matrix_compliance[1,2] = S23

        self.strand_matrix_compliance[1,0] = S12
        self.strand_matrix_compliance[2,0] = S13
        self.strand_matrix_compliance[2,1] = S23

        self.strand_matrix_stiffness = np.linalg.inv(self.strand_matrix_compliance)
        
        # Resin Matrix Stiffness

        self.resin_matrix_compliance = np.zeros((6,6))
        self.resin_matrix_compliance[0,0]=1
        self.resin_matrix_compliance[0,1]=-Vm
        self.resin_matrix_compliance[0,2]=-Vm
        self.resin_matrix_compliance[1,1]=1
        self.resin_matrix_compliance[2,2]=1
        self.resin_matrix_compliance[3,3]=1+Vm
        self.resin_matrix_compliance[4,4]=self.resin_matrix_compliance[3,3]
        self.resin_matrix_compliance[5,5]=self.resin_matrix_compliance[3,3]
        self.resin_matrix_compliance[1,2]=Vm

        for i in range(6):
            for j in range(i):
                self.resin_matrix_compliance[i,j]=self.resin_matrix_compliance[j,i]


        self.resin_matrix_compliance = Em*self.resin_matrix_compliance

        self.resin_matrix_stiffness = np.linalg.inv(self.resin_matrix_compliance)

    def AverageStiffness(self):
        '''
        This method will calculate the average stiffness of a strand
        '''
        C_s_0 = self.strand_matrix_stiffness
        # print(C_s_0)
        self.TrasnformationTMaterial(90)
        C_s_90 = self.stiffness_matrix_rotate
        # print(C_s_90)
        self.Avg_Strand_Stiffness = (C_s_0+C_s_90)/2

    def TrasnformationTMaterial(self,theta):
        '''
        Matrix to tranform a material matrix 
        '''
        theta = m.radians(theta)
        T11 = (np.cos(theta))**2
        T12 = np.sin(theta)**2
        T13 = 0 
        T14 = 0
        T15 = 0
        T16 = -np.sin(2*theta)
        T21 = (np.sin(theta))**2
        T22 = (np.cos(theta))**2
        T23 = 0
        T24 = 0
        T25 = 0
        T26 = np.sin(2*theta)
        T31 = 0
        T32 = 0
        T33 = 1
        T34 = 0
        T35 = 0
        T36 = 0
        T41 = 0
        T42 = 0
        T43 = 0
        T44 = np.cos(theta)
        T45 = np.sin(theta)
        T46 = 0
        T51 = 0
        T52 = 0
        T53 = 0
        T54 = -np.sin(theta)
        T55 = np.cos(theta)
        T56 = 0
        T61 = np.sin(theta)*np.cos(theta)
        T62 = -np.cos(theta)*np.cos(theta)
        T63 = 0
        T64 = 0
        T65 = 0
        T66 = (np.cos(theta))**2-(np.sin(theta))**2 
        
        T = np.array([[T11,T12,T13,T14,T15,T16],
                     [T21,T22,T23,T24,T25,T26],
                     [T31,T32,T33,T34,T35,T36],
                     [T41,T42,T43,T44,T45,T46],
                     [T51,T52,T53,T54,T55,T56],
                     [T61,T62,T63,T64,T65,T66]])

        self.T = T
        # teste = np.random.randint(100, size=(6, 6))
        # trans1 = np.zeros((6,6))
        # for i in range(6):
        #     for j in range(6):
        #         trans1[i,j]=teste[i,j]
        # print(trans1)
        self.stiffness_matrix_rotate = np.matmul((np.matmul(T,self.strand_matrix_stiffness)),T.transpose())
        # print(self.stiffness_matrix)

    def CorrectStiffness(self,factor):
        '''
        Method to correct the stiffness matrix by a factor

        Plain weave ---- 1.173
        Eight Harness Satin ---- 1.001
        Twill-weave ---------  1.017

        '''

        #old procedure

        self.stiffness_matrix_simmetric = self.stiffness_matrix_simmetric*factor
        self.stiffness_matrix =  self.stiffness_matrix*factor
        self.ortho_matrix_compliance = np.linalg.inv(self.stiffness_matrix)
        self.compliance_matrix_simmetric = np.linalg.inv(self.stiffness_matrix_simmetric)
        self.stiffness_matrix_original = self.stiffness_matrix_original*factor
        self.compliance_matrix_original = np.linalg.inv(self.stiffness_matrix_original)
            
    def CorrectStiffnessV2(self,Vg_1,Vg_2):
        '''
        Method to correct the stiffness matrix
        '''
        # self.PlainWeaveSida()
        self.AverageStiffness()

        fact = Vg_1-Vg_2
        self.stiffness_matrix_simmetric = self.stiffness_matrix_simmetric+(self.Avg_Strand_Stiffness-self.resin_matrix_stiffness)*fact
        self.compliance_matrix_simmetric = np.linalg.inv(self.stiffness_matrix_simmetric)
        #recalculating compliance

        self.compliance_matrix_simmetric

    def SetExperimentalProperties(self,EngDict):
        '''
        Method to set the current Stiffness Matrix considering experimental Results.

        The EngDict will have followin kwords:

        E1 v12 G12
        E2 v23 G23
        E3 v13 G13

        '''

        E1 = EngDict['E1']
        E2 = EngDict['E2']
        E3 = EngDict['E3']
        G12 = EngDict['G12']
        G23 = EngDict['G23']
        G13 = EngDict['G13']
        v12 = EngDict['v12']
        v13 = EngDict['v13']
        v23 = EngDict['v23']

        S11 = 1/E1
        S22 = 1/E2
        S33 = 1/E3
        S44 = 1/G12
        S55 = 1/G23
        S66 =1/G13
        S12 = -v12/E1
        S13 = -v13/E1
        S23 = -v23/E2

        self.experimental_matrix_compliance = np.zeros((6,6))

        self.experimental_matrix_compliance[0,0] = S11
        self.experimental_matrix_compliance[1,1] = S22
        self.experimental_matrix_compliance[2,2] = S33
        self.experimental_matrix_compliance[3,3] = S44
        self.experimental_matrix_compliance[4,4] = S55
        self.experimental_matrix_compliance[5,5] = S66
        self.experimental_matrix_compliance[0,1] = S12
        self.experimental_matrix_compliance[0,2] = S13
        self.experimental_matrix_compliance[1,2] = S23

        self.experimental_matrix_compliance[1,0] = S12
        self.experimental_matrix_compliance[2,0] = S13
        self.experimental_matrix_compliance[2,1] = S23

        self.experimental_matrix_stiffness = np.linalg.inv(self.experimental_matrix_compliance)

        self.Kv_experimental,self.Gv_experimental=self.Voigt(self.experimental_matrix_stiffness)
        self.Kr_experimental,self.Gr_experimental=self.Reuss(self.experimental_matrix_compliance)
        #performing parameters calculation for the experimental results:
        self.aniso_experimental=self.AnisoIndexCalculation(self.Gv_experimental,self.Gr_experimental,self.Kv_experimental,self.Kr_experimental)
        self.ZenerIndexCalculation(experimental_matrix=True)
        self.FrobenniusNorm(experimental_matrix=True)

# %%
if __name__=="__main__":

    import os

    curren_folder = os.getcwd()

    def PrintPretty(string):
        print("{:8.3}".format(string))

    # weaves = ["Plain Weave matrices",]
    # folder = curren_folder+r"\Plain Weave matrices\UDBC"
    folder = curren_folder
    stiffness_arquive = folder+r"\stiffness_unidirectional.out"
    ClassInitiation = AnisoIndex(stiffness_arquive)

    # EglassVynil = {

    #     'E1':57.5E3,
    #     'E2':18.8E3,
    #     'E3':18.8E3,
    #     'G12':7.44E3,
    #     'G23':7.26E3,
    #     'G13':7.26E3,
    #     'v12':0.25,
    #     'v13':0.29,
    #     'v23':0.29
    # }

    # VynilEster = {

    #     'E1':3.4E3,
    #     'v12':0.35
    # }

    # ClassInitiation.CompositeStiffnessMatrix(EglassVynil,VynilEster)
    ClassInitiation.CorrectStiffness(factor=1)
    ClassInitiation.CalculateCalculateEngineeringConstants()
    # Matrix = ClassInitiation.stiffness_matrix_simmetric
    # ClassInitiation.CorrectStiffness(factor=1.1)
    A=ClassInitiation.CalculateCalculateEngineeringConstants()
    list_to_show = zip(['E1','E2','G12','v12'],[A['E1'],A['E2'],A['G12'],A['v12'],A['v12']])
    for label,item in list_to_show:
        if label in ['E1','E2','G12']: 
            item = item*0.001 #transforming to GPa
        print("{},{:1.4}".format(label,item))



    Au = ClassInitiation.aniso_index_user_material_simmetric
    Az = ClassInitiation.ZenerIndex
    Fr = ClassInitiation.Frobenius

    # print("Norms: \n\n ++++++++++++")
    # PrintPretty(Au)
    # PrintPretty(Az)
    # print(Fr)

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
# %%
