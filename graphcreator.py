
# %%

import numpy as np
import matplotlib.pyplot as plt
import re
EnginerringPropertiesResultsS={"dictkeySUBC":{"E1":10},
                              "dictkeyMBC":{"E1":20},
                              "dictkeyPBC":{"E1":30},
                              "dictkeySUBC2D":{"E1":10},
                              "dictkeyMBC2D":{"E1":20},
                              "dictkeyPC2D":{"E1":30},
                                }

class CreateGraphsTest:

    def returnGroupsPosition(self,dictkey):
        '''
        Returns the user-defined array position for the group
        '''
        if re.search("dictkey",dictkey):
            return 0

    def Eng3D(self,EnginerringPropertiesResults,property,testedCases=1):
 
        # set width of bars
        barWidth = 0.25
        
        propSUBC = [None]*testedCases
        propPBC = [None]*testedCases
        propMBC = [None]*testedCases



        for each_key in EnginerringPropertiesResults.keys():
            if not re.search("2D",each_key):
                if re.search("SUBC",each_key):
                    propSUBC[self.returnGroupsPosition(each_key)]=EnginerringPropertiesResults[each_key][property]
                if re.search("PBC",each_key):
                    propPBC[self.returnGroupsPosition(each_key)]=EnginerringPropertiesResults[each_key][property]
                if re.search("MBC",each_key):
                    propMBC[self.returnGroupsPosition(each_key)]=EnginerringPropertiesResults[each_key][property]

        
        # Set position of bar on X axis
        r1 = np.arange(len(propSUBC))
        r2 = [x + barWidth for x in r1]
        r3 = [x + barWidth for x in r2]
        
        # Make the plot
        fig, ax = plt.subplots()  # a figure with a single Axeses
        

        ax.bar(r1, propSUBC, color='#526760', width=barWidth, edgecolor='white', label='SUBC')
        ax.bar(r2, propPBC, color='#9ED8DB', width=barWidth, edgecolor='white', label='PBC')
        ax.bar(r3, propMBC, color='#1D3354', width=barWidth, edgecolor='white', label='MBC$^{\mho}$')
        
        # Add xticks on the middle of the group bars
        ax.set_xlabel('Grupo', fontweight='bold')
        ax.set_ylabel(f'{property}', fontweight='bold')
        # plt.xticks([r + barWidth for r in range(len(propSUBC))], ['A', 'B', 'C', 'D', 'E','F','G'])
        ax.set_xticks([r + barWidth for r in range(len(propSUBC))])
        ax.set_xticklabels(['A'])
        
        # Create legend & Show graphic
        ax.legend(bbox_to_anchor=(1.05, 1.05))
        # plt.show()
        self.figure_name = f"ENG_{property}.png"
        fig.set_size_inches(6.2, 4)
        fig.savefig(self.figure_name)
        print(fig.get_figwidth())
        print(fig.get_figheight())
        self.property = property


    def CreateLatexFigure(self):
        '''
            \begin{figure}[H]
                \centering
                    \centering
                    \includegraphics[width=0.8\textwidth]{fig/fig11.jpg}
                    \caption{Geometrias do RVE (sem resina) gerada para cada padrão da \Fig{fig:weave_patterns}. }
                    \label{fig:texgen_model}
            \end{figure}
        '''
        l1 = r"\begin{figure}[H]"+"\n"
        l2 = r"    \centering"+"\n"
        l3 = r"\includegraphics[width=0.8\textwidth]{fig/"+f"{self.figure_name}"+r"}"+"\n"
        l4 = fr"\caption{{Comparação de propriedades efetivas para {self.property}}}"+"\n"
        l5 = fr"        \label{{fig:comparacao_{self.property}}}"+"\n"
        l6 = r"\end{figure}"+"\n"

        with open(f"comparacao_{self.property}.tex",'w',encoding='utf-8') as nf:
            nf.write("".join([l1,l2,l3,l4,l5,l6]))

        return "Done image tex creating"

# Plotter = CreateGraphsTest()
# Plotter.Eng3D(EnginerringPropertiesResultsS,"E1")
# Plotter.CreateLatexFigure()

# %%
class CreateGraphs():

    def __init__(self):

        self.CLT_calculatedProperties = {'E1': 25910.600153707783,
                                        'E2': 25910.600153707783,
                                        'G12': 4018.437296404772,
                                        'v12': 0.09759649285067286}

    def returnGroupsPosition(self,dictkey):
        '''
        Returns the user-defined array position for the group
        '''
        if re.search("plainWeave",dictkey):
            return 0
        if re.search("twoTwoTwillWeave",dictkey):
            return 1
        if re.search("fiveHarnessSatinWeave",dictkey):
            return 2
        if re.search("eightHarnessSatinWeave",dictkey):
            return 3
        if re.search("basketWeave",dictkey):
            return 4
        if re.search("leftTwillOneTwo",dictkey):
            return 5
        if re.search("rightTwillOneTwo",dictkey):
            return 6

    def ReturnNiceLatex(self,Property):
        '''
        Just to return a nice latex
        '''
        if Property=="E1":
            return r"$E_{11}$"
        if Property=="G12":
            return r"$G_{12}$"
        if Property=="E3":
            return r"$E_{33}$"
        if Property=="E2":
            return r"$E_{22}$"
        if Property=="G23":
            return r"$G_{23}$"
        if Property=="v13":
            return r"$nu_{13}$"
        if Property=="v23":
            return r"$nu_{23}$"
        if Property=="v12":
            return r"$nu_{12}$"
        if Property=="G31":
            return r"$G_{31}$"
        
        return Property


    def Eng2D(self,EnginerringPropertiesResults,property,testedCases=7):
 
        # set width of bars
        barWidth = 0.25
        
        propSUBC = [None]*testedCases
        propPBC = [None]*testedCases
        propMBC = [None]*testedCases

        for each_key in EnginerringPropertiesResults.keys():
            if not re.search("2D",each_key):
                if re.search("SUBC",each_key):
                    propSUBC[self.returnGroupsPosition(each_key)]=EnginerringPropertiesResults[each_key][property]
                if re.search("PBC",each_key):
                    propPBC[self.returnGroupsPosition(each_key)]=EnginerringPropertiesResults[each_key][property]
                if re.search("MBC",each_key):
                    propMBC[self.returnGroupsPosition(each_key)]=EnginerringPropertiesResults[each_key][property]

        
        # Set position of bar on X axis
        r1 = np.arange(len(propSUBC))
        r2 = [x + barWidth for x in r1]
        r3 = [x + barWidth for x in r2]
        
        # Make the plot
        fig, ax = plt.subplots()  # a figure with a single Axeses

        ax.bar(r1, propSUBC, color='#526760', width=barWidth, edgecolor='white', label='SUBC')
        ax.bar(r2, propPBC, color='#9ED8DB', width=barWidth, edgecolor='white', label='PBC')
        ax.bar(r3, propMBC, color='#1D3354', width=barWidth, edgecolor='white', label='MBC$^{\mho}$')

        if property=="E1":
            point1 = [0,self.CLT_calculatedProperties[property]]
            point2 =[barWidth*len(propSUBC)*3,self.CLT_calculatedProperties[property]]
            
            ax.plot([point1[0],point2[0]],[point1[1],point2[1]],label="CLT")

        if property=="E2":
            point1 = [0,self.CLT_calculatedProperties[property]]
            point2 =[barWidth*len(propSUBC)*3,self.CLT_calculatedProperties[property]]
            
            ax.plot([point1[0],point2[0]],[point1[1],point2[1]],label="CLT")

        if property=="G12":
            point1 = [0,self.CLT_calculatedProperties[property]]
            point2 =[barWidth*len(propSUBC)*3,self.CLT_calculatedProperties[property]]
            
            ax.plot([point1[0],point2[0]],[point1[1],point2[1]],label="CLT")
        
        # Add xticks on the middle of the group bars
        ax.set_xlabel('Grupo', fontweight='bold')
        ax.set_ylabel(f'{self.ReturnNiceLatex(property)}')
        # plt.xticks([r + barWidth for r in range(len(propSUBC))], ['A', 'B', 'C', 'D', 'E','F','G'])
        ax.set_xticks([r + barWidth for r in range(len(propSUBC))])
        ax.set_xticklabels(['A', 'B', 'C', 'D', 'E','F','G'])
        
        # Create legend & Show graphic
        ax.legend()
        # plt.show()
        self.figure_name = f"ENG2D_{property}.png"
        fig.set_size_inches(6.3, 4)
        fig.savefig(self.figure_name)
        self.property = property



    def Eng3D(self,EnginerringPropertiesResults,property,testedCases=7):
 
        # set width of bars
        barWidth = 0.25
        
        propSUBC = [None]*testedCases
        propPBC = [None]*testedCases
        propMBC = [None]*testedCases

        for each_key in EnginerringPropertiesResults.keys():
            if not re.search("2D",each_key):
                if re.search("SUBC",each_key):
                    propSUBC[self.returnGroupsPosition(each_key)]=EnginerringPropertiesResults[each_key][property]
                if re.search("PBC",each_key):
                    propPBC[self.returnGroupsPosition(each_key)]=EnginerringPropertiesResults[each_key][property]
                if re.search("MBC",each_key):
                    propMBC[self.returnGroupsPosition(each_key)]=EnginerringPropertiesResults[each_key][property]

        
        # Set position of bar on X axis
        r1 = np.arange(len(propSUBC))
        r2 = [x + barWidth for x in r1]
        r3 = [x + barWidth for x in r2]
        
        # Make the plot
        fig, ax = plt.subplots()  # a figure with a single Axeses

        ax.bar(r1, propSUBC, color='#526760', width=barWidth, edgecolor='white', label='SUBC')
        ax.bar(r2, propPBC, color='#9ED8DB', width=barWidth, edgecolor='white', label='PBC')
        ax.bar(r3, propMBC, color='#1D3354', width=barWidth, edgecolor='white', label='MBC$^{\mho}$')
        
        # Add xticks on the middle of the group bars
        ax.set_xlabel('Grupo', fontweight='bold')
        ax.set_ylabel(f'{self.ReturnNiceLatex(property)}')
        # plt.xticks([r + barWidth for r in range(len(propSUBC))], ['A', 'B', 'C', 'D', 'E','F','G'])
        ax.set_xticks([r + barWidth for r in range(len(propSUBC))])
        ax.set_xticklabels(['A', 'B', 'C', 'D', 'E','F','G'])

        if property=="E1":
            point1 = [0,self.CLT_calculatedProperties[property]]
            point2 =[barWidth*len(propSUBC)*3,self.CLT_calculatedProperties[property]]
            
            ax.plot([point1[0],point2[0]],[point1[1],point2[1]],label="CLT")

        if property=="E2":
            point1 = [0,self.CLT_calculatedProperties[property]]
            point2 =[barWidth*len(propSUBC)*3,self.CLT_calculatedProperties[property]]
            
            ax.plot([point1[0],point2[0]],[point1[1],point2[1]],label="CLT")

        if property=="G12":
            point1 = [0,self.CLT_calculatedProperties[property]]
            point2 =[barWidth*len(propSUBC)*3,self.CLT_calculatedProperties[property]]
            
            ax.plot([point1[0],point2[0]],[point1[1],point2[1]],label="CLT")


        # Create legend & Show graphic
        ax.legend()
        # plt.show()
        self.figure_name = f"ENG3D_{property}.png"
        # fig.set_size_inches(7, 4)
        fig.set_size_inches(6.3, 4)
        fig.savefig(self.figure_name)
        self.property = property




    def CreateLatexFigure(self):
        '''
            \begin{figure}[H]
                \centering
                    \centering
                    \includegraphics[width=0.8\textwidth]{fig/fig11.jpg}
                    \caption{Geometrias do RVE (sem resina) gerada para cada padrão da \Fig{fig:weave_patterns}. }
                    \label{fig:texgen_model}
            \end{figure}
        '''
        l1 = r"\begin{figure}[H]"+"\n"
        l2 = r"    \centering"+"\n"
        l3 = r"\includegraphics[width=0.8\textwidth]{fig/"+f"{self.figure_name}"+r"}"+"\n"
        l4 = fr"\caption{{Comparação de propriedades efetivas para {self.property}}}"+"\n"
        l5 = fr"        \label{{fig:comparacao_{self.property}}}"+"\n"
        l6 = r"\end{figure}"+"\n"

        with open(f"comparacao_{self.figure_name}.tex",'w',encoding='utf-8') as nf:
            nf.write("".join([l1,l2,l3,l4,l5,l6]))

        return "Done image tex creating"
