# Module to calculte the Kelving-Voigt bulk and shear modulus, by an input matrix,  then 
# calculate the general anysotropy index
# %%
from create_matrix import NumptArrrayToBmatrix
import numpy as np
import math as m
from PostProcessing_homot_tool import PostProcessing
import os
from CLT import*
# from kelvin_reuss import WriteTexBmatrix
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
        
class AnisoIndex(NumptArrrayToBmatrix,PostProcessing):

    
    def __init__(self,file=False):
        self.name = "SUBC"
        pass

    def CalculateIndicators(self,file):
        '''
        Read Stiffness.out file, and calculate indicators
        '''


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

    def CheckIfFolderExist(self,folder):
        '''
        Check if the folder exists. If don't, create it
        '''
        CHECK_FOLDER = os.path.isdir(folder)
        
        if not CHECK_FOLDER:
            os.makedirs(folder)


    def HomogenizeResults(self,results_folder = 'Default',BoundaryConditions='ALL'):
        self.PostProcessFiniteElement = PostProcessing(results_folder,BoundaryConditions)
        return None

    def CalculateEngConstants(self):
        self.PostProcessFiniteElement.CalculateEngConstants()
        self.EngineeringConstants = self.PostProcessFiniteElement.EngineeringConstants
        return None


    def ShowEngineeringConstants(self, Property='ALL'):
        return self.PostProcessFiniteElement.ShowEngineeringConstants(Property=Property)

    def ExportStiffness(self, folder, name="stiffness"):
        self.CheckIfFolderExist(folder)
        return self.PostProcessFiniteElement.ExportStiffness(folder,name)

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

        C_original = self.stiffness_matrix_original 

        C_4_16 = [i for i in C_original[3,0:6]]
        C_6_16 = [i for i in C_original[5,0:6]]

        C_original[3,0:6]=C_6_16
        C_original[5,0:6]=C_4_16

        # Swapping fourth column with sixth column
        temp = [i for i in C_original[0:6,3]]
        C_original[0:6,3]=[i for i in C_original[0:6,5]]
        C_original[0:6,5]=temp

        # Swapping fourth line with fifth line
        temp = [i for i in C_original[3,0:6]]
        C_original[3,0:6]=[i for i in C_original[4,0:6]]
        C_original[4,0:6]=temp

        # Swapping fourth column with fifth column
        temp = [i for i in C_original[0:6,3]]
        C_original[0:6,3]=[i for i in C_original[0:6,4]]
        C_original[0:6,4]=temp

        self.stiffness_matrix_original_voigt = C_original

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
    
    def BulkModulusTensor(self):

        K = Invariants()
        K_matrix = K.matrix

        # M = self.ortho_matrix_stiffness
        M = self.stiffness_matrix_simmetric

        for i in range(3):
            for j in range(i,3):
                first_index = self.IndexTranformation(i,j)
                # print(f"(i,j)={i},{j}->{first_index}")
                K_matrix[i,j] = M[first_index,0]+M[first_index,1]+M[first_index,2]
                K_matrix[j,i]=K_matrix[i,j]
                   
        K.matrix = K_matrix
        K.CalculateInvariants()

        return K

    def DeviatoryModulusTensor(self):

        '''
        Matrix is a 6x6 matrix 
        '''

        L = Invariants()
        L_matrix = L.matrix
        M = self.stiffness_matrix_simmetric

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

class ComparetoCLT():
    '''
    Class to execute the comparison within the CLT
    '''
    def __init__(self,matrix):

        self.Contraction = TensorContracion(matrix)
        self.C_contracted = self.Contraction.C_contract

        self.S_contracted = np.linalg.inv(self.C_contracted)
        self.Calculate_eng()
        self.name = "SUBC"


    def Calculate_eng(self):
        '''
        Method to calculate de engineering constants from the contrated matrix
        '''

        S = self.S_contracted

        E11 = 1/S[0,0]
        E22 = 1/S[1,1]
        G12 = 1/S[2,2]
        v12 = -S[0,1]/S[1,1]

        self.EngineeringConst = {
            "E1":E11,
            "E2":E22,
            "G12":G12,
            "v12":v12
        }


    def Calculate_CLT_eng(self,fiber,matrix,h):
        '''
        Calculates the engineering constants based on Barberos approach
        '''
        return compare_to_experiments(fiber,matrix,h)

class GenerateTablesAndMatrix3D(NumptArrrayToBmatrix):

    def __init__(self,SUBC,PBC,MBC,textile):
        

        self.SUBC = SUBC
        self.PBC = PBC
        self.MBC = MBC
        self.textile = textile
        self.WriteCaptionNames()





    def DeleteTexFiles(self):
        list_files = os.listdir(os.getcwd())
        tex_files = [a for a in list_files if a.endswith(".tex")]
        [os.remove(file) for file in tex_files]
        return "Deleted All *.tex files"
        
    def WriteCaptionNames(self):

        if self.textile == "Plain-weave-Test":
            self.caption_name = r"{\rightTwillOneTwoMiddle}"
        
        if self.textile == "plainWeave":
            self.caption_name = r"{\plainWeaveMiddle}"
        if self.textile == "twoTwoTwillWeave":
            self.caption_name = r"{\twoTwoTwillWeaveMiddle}"
        if self.textile == "fiveHarnessSatinWeave":
            self.caption_name = r"{\fiveHarnessSatinWeave}"
        if self.textile == "eightHarnessSatinWeave":
            self.caption_name = r"{\eightHarnessSatinWeave}"
        if self.textile == "basketWeave":
            self.caption_name = r"{\basketWeave}"
        if self.textile == "leftTwillOneTwo":
            self.caption_name = r"{\leftTwillOneTwoMiddle}"
        if self.textile == "rightTwillOneTwo":
            self.caption_name = r"{\rightTwillOneTwoMiddle}"


    def Write3dMatrix(self,tridimensionalObj):
        self.bmatrix_name = self.textile+"_"+tridimensionalObj.name
        super().WriteTexBmatrixFromNpArray(tridimensionalObj.stiffness_matrix_original_voigt,self.bmatrix_name)
        self.GenerateEquationforMatrix(tridimensionalObj)

    def OpenTextFile(self,file):
        with open(f"{file}.tex", encoding='utf-8') as texfile:
            self.texcontent = texfile.read()
            

    def BindAllFiles(self,separation="\n\n"):
        '''
        Method to bind all files generated. The easiest way is to generate a file, and then, automatically add to a variable, and after generate the file
        '''
        tex_all_files_content = fr"\section{{Resultados da an??lise 3D para o {self.caption_name}}}"
        self.Write3dMatrix(self.SUBC)
        self.OpenTextFile(f"b_matrix_{self.bmatrix_name}")
        tex_all_files_content += separation+self.texcontent
        self.Write3dMatrix(self.PBC)
        self.OpenTextFile(f"b_matrix_{self.bmatrix_name}")
        tex_all_files_content += separation+self.texcontent
        self.Write3dMatrix(self.MBC)
        self.OpenTextFile(f"b_matrix_{self.bmatrix_name}")
        tex_all_files_content += separation+self.texcontent

        self.GenerateTableEngConst()
        self.OpenTextFile("3dTable")
        tex_all_files_content += separation+self.texcontent


        self.GenerateTableInvar()
        self.OpenTextFile("3dTableInvariants")
        tex_all_files_content+=separation+self.texcontent

        self.GenerateTableIndicators()
        self.OpenTextFile("3dTableIndicators")
        tex_all_files_content+=separation+self.texcontent

        with open(f"general_tex_{self.textile}.text", 'w', encoding='utf-8') as newfile:
            newfile.write(tex_all_files_content)

        return 


    def GenerateEquationforMatrix(self,tridimensionalObj):
        '''
        \begin{equation}
            {\stiffnessPlane}_{\text{\twoTwoTwillWeaveShort}}^{\proposedMBC} =
            \begin{bmatrix} \label{eq:2_d_twill_weave_MBC}
        21915.77 & 5121.02 & 5027.82 \\
        5121.02 & 22076.02 & 5033.24 \\
        5027.82 & 5033.24 & 10078.52 \\
        \end{bmatrix} 
        \end{equation}

        '''
        name=self.bmatrix_name + ".tex"
        with open(name) as bmatrix:
            b_matrix_text=bmatrix.read()

        if tridimensionalObj.name=="SUBC":
            BC_name = "SUBC"
        if tridimensionalObj.name=="PBC":
            BC_name = "PBC"
        if tridimensionalObj.name=="MBC":
            BC_name = "proposedMBC"


        l1 = r"\begin{equation}"+"\n" 
        # COM NOME DO TEXTIL
        # l2 = r"{\stiffnessHomogenizedVoigt}_{\text{"+"\\"+f"{self.textile}"+r"}}^{"+"\\"+f"{BC_name}"+r"} ="+"\n"
        # SEM NOME DO TEXTIL
        l2 = r"{\stiffnessHomogenizedVoigt}^{"+"\\"+f"{BC_name}"+r"} ="+"\n"
        l3 = r"\end{equation}"+"\n"
        f = open(f"b_matrix_{self.bmatrix_name}.tex",'w')
        f.write("".join([l1,l2,b_matrix_text,l3]))
        f.close()

    # def LatexTableKeyWords(self,**kwargs):
    #     '''
    #     Create code snnipets
    #     '''
    #     if kwargs.get("begin"):
    #         return ""

    def FormatTwoDecimals(self,number):
        return "{:.2f}".format(number)

    def GenerateTableEngConst(self):
        '''
        \begin{table}[h]
        \centering
        \caption{Effective material properties comparison for E-glass/vinylester plain-weave composite.}
        \begin{tabular}{ccccccccccc}
        \hline
                                        b.c & $E_1$ & $E_2$ & $E_3$  & $G_{12}$ &  $G_{23}$ & $G_{31}$ & $\nu_{13}$ & $\nu_{23}$  & $\nu_{12}$\\
        \hline
        %INSERT DATA HERE
        SUBC                            & 25.18 & 25.81 & 11.06 & 4.74 & 2.70 & 2.97 &  0.42 & 0.35 & 0.13\\
        PBC                             & 24.66 & 25.12 & 10.87 & 4.84 & 2.76 & 2.76 &  0.41 & 0.41 & 0.13\\
        MBC$^*$                             & 24.78 & 25.52 & 10.92 & 4.84 & 2.44 & 2.88 &  0.41 & 0.41 & 0.13\\
        \hline
        \end{tabular} \label{tab:plain_weave_results_compare}
        \end{table}
        '''
        
        Eng = self.ReturnEngineeringConstant3D
        l1 = "\\begin{table}[H]\n"
        l2 = "\centering\n"
        l3 = f"\\caption{{Compara????o de propriedadades efetivas para a trama {self.caption_name} }}\n"
        l4 = "\\begin{tabular}{ccccccccccc}\n"
        l5 = "\hline\n"
        l6 = "                                c.c & $E_1$ & $E_2$ & $E_3$  & $G_{12}$ &  $G_{23}$ & $G_{31}$ & $\\nu_{13}$ & $\\nu_{23}$  & $\\nu_{12}$"+r"\\"+ "\n"
        l7 = "\hline\n"
        l8 = "%\n"
        # Just loop trought the boundary conditions
        # BC
        l9 = []
        BC = [self.SUBC,self.PBC,self.MBC]
        print(self.SUBC)
        for each_bc in BC:
            if each_bc.name=="SUBC":
                BC_name = "SUBC"
            if each_bc.name=="PBC":
                BC_name = "PBC"
            if each_bc.name=="MBC":
                BC_name = "proposedMBCText"

            self.tridimensionalObj = each_bc
            l9.append(fr"\{BC_name}                            & {Eng('E1')} & {Eng('E2')} & {Eng('E3')} & {Eng('G12')} & {Eng('G23')} & {Eng('G31')} &  {Eng('v31',factor=1)} & {Eng('v23',factor=1)} & {Eng('v12',factor=1)}\\"+"\n")

        l9 = "".join(l9)
        l10 = "\hline\n"
        l11 = f"\end{{tabular}} \label{{tab:engconst_3d{self.textile}_{BC_name}}}\n"
        l12 = "\end{table}"

        f = open("3dTable.tex","+w",encoding='utf-8')
        f.write("".join([l1,l2,l3,l4,l5,l6,l7,l8,l9,l10,l11,l12]))
        f.close()
        return None


    def GenerateTableInvar(self):
        '''
        \begin{table}[H]
        \centering
        \caption{Invariantes de $\bulkModulus$ e $\deviatoricModulus$ para o \plainWeaveMiddle}
        \begin{tabular}{@{}lllllllll@{}}
        \toprule
        c.c.      & $I^K_1$ & $I^K_2$ & $I^K_3$ & $I^L_1$ & $I^L_2$ & $I^L_3$  \\ \midrule
        \SUBC       & $1.08 \times 10^5$    & $3.82 \times 10^9$    & $4.36 \times 10^{13}$   & $9.22 \times 10^4$  & $2.76 \times 10^9$      & $2.65 \times 10^{13}$      \\
        \PBC       & $1.08 \times 10^5$    & $3.79 \times 10^9$    & $4.31 \times 10^{13}$   & $9.09 \times 10^4$  & $2.68 \times 10^9$      & $2.55 \times 10^{13}$   \\
        \proposedMBCText       & $1.08 \times 10^5$    & $3.80 \times 10^9$    & $4.34 \times 10^{13}$   & $9.08 \times 10^4$  & $2.68 \times 10^9$      & $2.55 \times 10^{13}$ \\ \bottomrule
        \end{tabular} \label{tab:ivariants_plain_weave}
        \end{table}
        '''
        l1 = "\\begin{table}[H]\n"
        l2 = "\centering\n"
        l3 = r"\caption{Invariantes de $\bulkModulus$ e $\deviatoricModulus$ para o " + f"{self.caption_name}" + "}\n"
        l4 = "\\begin{tabular}{@{}lllllllll@{}}\n"
        l5 = r"\toprule"+"\n"
        l6 = r" c.c.      & $I^K_1$ & $I^K_2$ & $I^K_3$ & $I^L_1$ & $I^L_2$ & $I^L_3$  \\ \midrule"+"\n"

        l7 = []
        BC = [self.SUBC,self.PBC,self.MBC]
        print(self.SUBC)
        for each_bc in BC:
            if each_bc.name=="SUBC":
                BC_name = "SUBC"
            if each_bc.name=="PBC":
                BC_name = "PBC"
            if each_bc.name=="MBC":
                BC_name = "proposedMBCText"

            self.tridimensionalObj = each_bc
            I = self.ReturnInvariant()
            l7.append(fr"\{BC_name}                            & ${I[0]}$    & ${I[1]}$    & ${I[2]}$   & ${I[3]}$  & ${I[4]}$      & ${I[5]}$      \\"+"\n")
        
        l7 = "".join(l7)
        l8 = r"\bottomrule"+"\n"
        # l9 = r"\end{tabular} \label{tab:ivariants_plain_weave}"+"\n"
        l9 = f"\end{{tabular}} \label{{tab:invariants3d_{self.textile}_{BC_name}}}\n"
        l10 = r"\end{table}"+"\n"

        f = open('3dTableInvariants.tex','w',encoding='utf-8')
        f.write("".join([l1,l2,l3,l4,l5,l6,l7,l8,l9,l10]))
        f.close()
        return None


    def GenerateTableIndicators(self):
        '''

        '''

        l1 = "\\begin{table}[H]\n"
        l2 = "\centering\n"
        l3 = r"\caption{Indicadores anisotr??picos para o "+f"{self.caption_name}"+r" analisado}"+"\n"
        l4 = r"\begin{tabular}{cccc}"+"\n"
        l5 = r"\toprule"+"\n"
        l6 = r"c.c. &  $\universalAnisotropyIndex$ & $\zenerAnisotropyIndex$ & $\frobenius$ \\ \midrule"+"\n"
        
        l7 = []
        BC = [self.SUBC,self.PBC,self.MBC]
        print(self.SUBC)
        for each_bc in BC:
            if each_bc.name=="SUBC":
                BC_name = "SUBC"
            if each_bc.name=="PBC":
                BC_name = "PBC"
            if each_bc.name=="MBC":
                BC_name = "proposedMBCText"

            self.tridimensionalObj = each_bc

            Au,Az,Fr = self.ReturnIndicator()
        
            l7.append(fr"\{BC_name}    &  {Au} & {Az} & {Fr}" +r"\\ "+"\n")

        l7 = "".join(l7)
        l8 = r"\bottomrule"+"\n"
        l9 = r"\end{tabular}"+r"\label{"+fr"tab:indicator_{self.textile}_{BC_name}"+"}"+"\n"
        l10 = r"\end{table}"+"\n"

        f = open('3dTableIndicators.tex','w',encoding='utf-8')
        f.write("".join([l1,l2,l3,l4,l5,l6,l7,l8,l9,l10]))
        f.close()

        return None

    def ReturnIndicator(self):

        Au = self.tridimensionalObj.aniso_index_user_material_simmetric
        Az = self.tridimensionalObj.ZenerIndex
        Fr = self.tridimensionalObj.Frobenius

        return "{:1.2f}".format(Au),"{:1.2f}".format(Az),int(Fr)

        
    def ReturnPrettyFormatExponencial(self,number):
        A="{:1.2E}".format(number).split("E")[0]
        B="{:1.2E}".format(number).split("E")[1]
        B = int(B)

        return A+r" \times"+" 10^{"+f"{B}"+"}"
    

    def ReturnInvariant(self):
        K = self.tridimensionalObj.BulkModulusTensor()
        L = self.tridimensionalObj.DeviatoryModulusTensor()
        R = self.ReturnPrettyFormatExponencial
        return [R(K.I1),R(K.I2),R(K.I3),R(L.I1),R(L.I2),R(L.I3)]

    def ReturnEngineeringConstant3D(self,Prop,factor=0.001):
        '''
        factor = correction of GPa
        '''
        # must defint 3d object
        return self.FormatTwoDecimals(self.tridimensionalObj.PostProcessFiniteElement.EngineeringConstants[f'{Prop}']*factor)

class GenerateTablesAndMatrix2D(GenerateTablesAndMatrix3D, NumptArrrayToBmatrix):

    def __init__(self, SUBC, PBC, MBC, textile):
        # self.WriteCaptionNames()
        super().__init__(SUBC, PBC, MBC, textile)

    def WriteCaptionNames(self):
        return super().WriteCaptionNames()
        
    def BindAllFiles(self, separation="\n\n"):
        '''
        Method to bind all files generated. The easiest way is to generate a file, and then, automatically add to a variable, and after generate the file
        '''
        tex_all_files_content = fr"\section{{Resultados da an??lise 2D para o {self.caption_name}}}"
        self.WriteBmatrix2D(self.SUBC)
        self.OpenTextFile(f"b_matrix_{self.bmatrix_name}")
        tex_all_files_content += separation+self.texcontent
        self.WriteBmatrix2D(self.PBC)
        self.OpenTextFile(f"b_matrix_{self.bmatrix_name}")
        tex_all_files_content += separation+self.texcontent
        self.WriteBmatrix2D(self.MBC)
        self.OpenTextFile(f"b_matrix_{self.bmatrix_name}")
        tex_all_files_content += separation+self.texcontent

        self.GenerateTableEngConst()
        self.OpenTextFile("2dEngineeringConst")
        tex_all_files_content += separation+self.texcontent

        with open(f"general_tex_{self.textile}_2d.text", 'w', encoding='utf-8') as newfile:
            newfile.write(tex_all_files_content)

    def WriteBmatrix2D(self,bidimensionalOBJ):
        self.bmatrix_name=self.textile+"_"+bidimensionalOBJ.name+"_2d"
        super().WriteTexBmatrixFromNpArray(bidimensionalOBJ.C_contracted,self.bmatrix_name)
        self.GenerateEquationforMatrix(bidimensionalOBJ)

    def GenerateEquationforMatrix(self,bidimensionalOBJ):
        '''
        \begin{equation}
            {\stiffnessPlane}_{\text{\twoTwoTwillWeaveShort}}^{\proposedMBC} =
            \begin{bmatrix} \label{eq:2_d_twill_weave_MBC}
        21915.77 & 5121.02 & 5027.82 \\
        5121.02 & 22076.02 & 5033.24 \\
        5027.82 & 5033.24 & 10078.52 \\
        \end{bmatrix} 
        \end{equation}

        '''
        name=self.bmatrix_name + ".tex"
        with open(name) as bmatrix:
            b_matrix_text=bmatrix.read()

        if bidimensionalOBJ.name=="SUBC":
            BC_name = "SUBC"
        if bidimensionalOBJ.name=="PBC":
            BC_name = "PBC"
        if bidimensionalOBJ.name=="MBC":
            BC_name = "proposedMBC"

        # b_matrix_text="\n".join(b_matrix_text)

        l1 = r"\begin{equation}"+"\n" 
        # COM NOME DO TEXTIL
        # l2 = r"{\stiffnessPlane}_{\text{"+"\\"+f"{self.textile}"+r"}}^{"+"\\"+f"{BC_name}"+r"} ="+"\n"
        # SEM NOME DO TEXTIL
        l2 = r"{\stiffnessPlane}_{\text{"+"\\"+f"{self.textile}"+r"}}^{"+"\\"+f"{BC_name}"+r"} ="+"\n"
        l3 = r"\end{equation}"+"\n"
        f = open(f"b_matrix_{self.bmatrix_name}.tex",'w')
        f.write("".join([l1,l2,b_matrix_text,l3]))
        f.close()

    def GenerateOrganizedResults(self):
        '''
        Method to generate a giant text file with orna
        '''
        pass
    
    def GenerateTableEngConst(self):
        '''
        \begin{table}[H]
        \centering
        \caption{Constantes de engenharia da rigidez reduzida para o \twoTwoTwillWeaveMiddle}
        \begin{tabular}{@{}lllllllll@{}}
        % 
        \toprule
        c.c.      & $E_{11}$ & $E_{22}$ & $G_{12}$ & $\nu_{12}$ \\ \midrule
        \SUBC       & 20.96    & 20.33    & 3.46     & 0.22      \\
        \PBC        & 20.56    & 20.76    & 3.49     & 0.23      \\ 
        \proposedMBCText        & 20.66    & 20.87    & 3.49     & 0.23  \\ \midrule
        \ROM        & 18.76    & 18.97    & 2.48    &  0.12   \\
        \bottomrule
        \end{tabular}
        \label{ec2}
        \end{table}
        '''
        l1 = "\\begin{table}[H]\n"
        l2 = r"\centering"+"\n"
        l3 = r"\caption{Constantes de engenharia da rigidez reduzida para o padr??o {"+"\\"+f"{self.textile}"+r"}}"+"\n"
        l4 = r"\begin{tabular}{@{}lllllllll@{}}"+"\n"
        l5 = r"%"+"\n"
        l6 = r"\toprule"+"\n"
        l7 = r"c.c.      & $E_{11}$ & $E_{22}$ & $G_{12}$ & $\nu_{12}$ \\ \midrule"+"\n"

        l8 = []
        BC = [self.SUBC,self.PBC,self.MBC]

        for each_bc in BC:
            if each_bc.name=="SUBC":
                BC_name = "SUBC"
            if each_bc.name=="PBC":
                BC_name = "PBC"
            if each_bc.name=="MBC":
                BC_name = "proposedMBCText"

            self.bidimensionalOBJ = each_bc

            E = self.ReturnEngConst2D

            l8.append(fr"\{BC_name}       & {E('E1')}     & {E('E2')}     & {E('G12')}      & {E('v12',factor=1)}      \\"+"\n")

        l8 = "".join(l8)
        l9 = r"\bottomrule"+"\n"
        l10 = f"\end{{tabular}} \label{{tab:eng_const_2d_{self.caption_name}_{BC_name}}}\n"
        l11 = r"\end{table}"+"\n"

        f = open('2dEngineeringConst.tex','w', encoding='utf-8')
        f.write("".join([l1,l2,l3,l4,l5,l6,l7,l8,l9,l10,l11]))
        f.close()

        return None
        
    def FormatTwoDecimals(self, number):
        return super().FormatTwoDecimals(number)

    def ReturnEngConst2D(self,EngConst,factor=0.001):
        return self.FormatTwoDecimals(self.bidimensionalOBJ.EngineeringConst[EngConst]*factor)

    
# %%
