
#%%
from postprocessing_files import *
import os
import unittest
from graphcreator import CreateGraphs
# %%


# Falta criar os arquivos pro caso bidimensional
# Materials

E_glass_Vinylester = {'E1':73E3,
                  'E2':73E3,
                  'G12':30.4E3,
                  'G23':30.4E3,
                  'v12':0.2,
                  'v23':0.2,
                  'vf':0.5535,
                  'Name':'E-glass'
                    }


VinilEsterDerakane = {'E1':3.4E3,
                      'v12':0.34,
                      'Name':r'Vinylester'
                        }

EnginerringPropertiesResults = {}
# DONT CHANGE THE ORDER
textiles = ["plainWeave","twoTwoTwillWeave","fiveHarnessSatinWeave","eightHarnessSatinWeave","basketWeave","leftTwillOneTwo","rightTwillOneTwo"]
# textiles = ["plainWeave"]
# stiffness_arquive = fr"D:\Masters\Masters\Thesis\Text\Rato\Capitulo de Livro - Jonas dissertação\Visual Studio Code\PostProcessingResults\teste\stiffness.out"
def SetResultsFolder(textile):
  return os.getcwd()+fr"\{textile}"

def SetStiffnessOutPut(results):
  return results+r"\stiffness.out"

# def ChooseTextilesName

for each_textile in textiles:
  results = SetResultsFolder(each_textile)+r"\UDBC"
  SUBC = AnisoIndex()
  SUBC.name ="SUBC"
  SUBC.HomogenizeResults(results_folder=results)
  SUBC.ExportStiffness(folder=results)
  SUBC.CalculateIndicators(SetStiffnessOutPut(results))
  SUBC.CalculateEngConstants()
  EnginerringPropertiesResults[f'{each_textile}_{SUBC.name}']=SUBC.EngineeringConstants

  matrix = SUBC.stiffness_matrix_simmetric
  SUBC_2D = ComparetoCLT(matrix)
  SUBC_2D.name = "SUBC"
  SUBC_2D.Calculate_CLT_eng(E_glass_Vinylester,VinilEsterDerakane,10)
  EnginerringPropertiesResults[f'{each_textile}_{SUBC.name}_2D']=SUBC_2D.EngineeringConst

  results = SetResultsFolder(each_textile)+r"\PBC"
  PBC = AnisoIndex()
  PBC.name ="PBC"
  PBC.HomogenizeResults(results_folder=results)
  PBC.ExportStiffness(folder=results)
  PBC.CalculateIndicators(SetStiffnessOutPut(results))
  PBC.CalculateEngConstants()
  EnginerringPropertiesResults[f'{each_textile}_{PBC.name}']=PBC.EngineeringConstants


  matrix = PBC.stiffness_matrix_simmetric
  PBC_2D = ComparetoCLT(matrix)
  PBC_2D.name = "PBC"
  CLT_Eng = PBC_2D.Calculate_CLT_eng(E_glass_Vinylester,VinilEsterDerakane,10)
  EnginerringPropertiesResults[f'{each_textile}_{PBC_2D.name}_2D']=PBC_2D.EngineeringConst

  results = SetResultsFolder(each_textile)+r"\MBC"
  MBC = AnisoIndex()
  MBC.name ="MBC"
  MBC.HomogenizeResults(results_folder=results)
  MBC.ExportStiffness(folder=results)
  MBC.CalculateIndicators(SetStiffnessOutPut(results))
  MBC.CalculateEngConstants()
  EnginerringPropertiesResults[f'{each_textile}_{MBC.name}']=MBC.EngineeringConstants

  matrix = MBC.stiffness_matrix_simmetric
  MBC_2D = ComparetoCLT(matrix)
  MBC_2D.name = "MBC"
  MBC_2D.Calculate_CLT_eng(E_glass_Vinylester,VinilEsterDerakane,10)
  EnginerringPropertiesResults[f'{each_textile}_{MBC_2D.name}_2D']=MBC_2D.EngineeringConst

  Results = GenerateTablesAndMatrix3D(SUBC,PBC,MBC,each_textile)

  Results.BindAllFiles()
  Results.DeleteTexFiles()

  Results2D = GenerateTablesAndMatrix2D(SUBC_2D,PBC_2D,MBC_2D,each_textile)
  # Results2D.WriteBmatrix2D(SUBC_2D)
  Results2D.BindAllFiles()

  Results2D.GenerateTableEngConst()

# %%
Plotter = CreateGraphs()
Plotter.Eng3D(EnginerringPropertiesResults,"E1")
Plotter.CreateLatexFigure()
Plotter.Eng3D(EnginerringPropertiesResults,"E2")
Plotter.CreateLatexFigure()
Plotter.Eng3D(EnginerringPropertiesResults,"G12")
Plotter.CreateLatexFigure()
Plotter.Eng3D(EnginerringPropertiesResults,"E3")
Plotter.CreateLatexFigure()
Plotter.Eng3D(EnginerringPropertiesResults,"G31")
Plotter.CreateLatexFigure()

Plotter.Eng2D(EnginerringPropertiesResults,"E1")
Plotter.CreateLatexFigure()
Plotter.Eng2D(EnginerringPropertiesResults,"E2")
Plotter.CreateLatexFigure()
Plotter.Eng2D(EnginerringPropertiesResults,"G12")
Plotter.CreateLatexFigure()
# %%
def OpenTexAndTakeContect(file):
  with open(file,encoding='utf-8') as tex_ready:
    tex_content = tex_ready.read()
  return tex_content

def CreateGiantTex():

    list_files = os.listdir(os.getcwd())
    tex_files = [a for a in list_files if a.endswith(".text")]
    string = [OpenTexAndTakeContect(file) for file in tex_files]
    
    with open("GiantTex.tex",'w',encoding='utf-8') as newfile:
      for each_file in string:
        newfile.write(each_file)
    return None

CreateGiantTex()
# %%
# def DeleteTexFiles():
#   list_files = os.listdir(os.getcwd())
#   tex_files = [a for a in list_files if a.endswith(".tex")]
#   removal=[os.remove(file) for file in tex_files]
#   return "Deleted All .tex files"
# %%
