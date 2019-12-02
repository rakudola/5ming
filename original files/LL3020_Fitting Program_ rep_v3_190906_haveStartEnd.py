import numpy as np
import pylab
import csv
from scipy.optimize import curve_fit
from scipy.stats import t
import matplotlib.pyplot as plt

#basic information --set up variable 
#method infodisgunsh
denatureMethod="CD" # choose 'CD' for chemical denature, 'HD' for heat denature
CD_xaxis_FinalGdmClConcentration_value=[0,0.43,0.87,1.3,1.74,2.17,2.61,3.04,3.48]
CD_conc_upbound=3.48
CD_conc_lowbound=0
CD_conc_range=CD_conc_upbound-CD_conc_lowbound
HD_xasis_normalizedTemp_value=[0,0.129277567,0.224334601,0.338403042,0.429657795,0.539923954,0.646387833,0.775665399,0.882129278,1]
HD_Temp_upbound=63
HD_Temp_lowbound=36.7
HD_temp_range=HD_Temp_upbound-HD_Temp_lowbound
CI_cutoff=0.3

#run info
descriptionOfTheRun="TTR-control"# what ever desciption can distinguish 
numOfReplicate="1" # number of replicate used to generate 1 fit 
plotColor="2"#1=blue,2=organe,3=purple,4=brown

#file selection
infile_locationAndName="G:\BYU OneDrive_G\OneDrive - BYU Office 365\BYU\===Price lab===\Data\Folding Project\LL3020_TTR_full Range_\\LL3020_pythonIn.csv"
outfile_locationAndName="G:\BYU OneDrive_G\OneDrive - BYU Office 365\BYU\===Price lab===\Data\Folding Project\LL3020_TTR_full Range_\\LL3020_fitoutput_modOnly.csv"
outImage_location="G:\Plot for fitting\LL3020\\"

#info for each peptide from infile
rowLoc=0
IDLoc=1
NameLoc=2
PeptideLoc=3
ReporterLoc=4
StartLoc=5
EndLoc=6

#rep1 info
rep1_notebookCode="LL3020" #the notebook code of the run
rep1_numPtsLocation=7# type in the number get from excel column to python raw translator
rep1_Statrt=8 # type in the number get from excel column to python raw translator
rep1_End=16# type in the number get from excel column to python raw translator

"""#rep2 info
rep2_notebookCode="LL1191" #the notebook code of the run
rep2_numPtsLocation=16# type in the number get from excel column to python raw translator
rep2_Statrt=17 # type in the number get from excel column to python raw translator
rep2_End=26# type in the number get from excel column to python raw translator"""

"""#rep3_info
rep3_notebookCode="LL1101" #the notebook code of the run
rep3_numPtsLocation=1# type in the number get from excel column to python raw translator
rep3_Statrt=5 # type in the number get from excel column to python raw translator
rep3_End=6# type in the number get from excel column to python raw translator"""

"""#rep4_info
rep4_notebookCode="LL1101" #the notebook code of the run
rep4_numPtsLocation=1# type in the number get from excel column to python raw translator
rep4_Statrt=5 # type in the number get from excel column to python raw translator
rep4_End=6# type in the number get from excel column to python raw translator"""


#plot color Code, usually don't need to change
color1_dataAndFit="#1f77b4"#blue
color1_Chalf="#aec7e8"
color1_CI="#c6dbef"


color2_dataAndFit="#ff7f0e"#orgnae
color2_Chalf="#ffbb78"
color2_CI="#fdd0a2"

color3_dataAndFit="#5254a3"#purple
color3_Chalf="#9e9ac8"
color3_CI="#dadaeb"

color4_dataAndFit="#8c6d31"#brown
color4_Chalf="#e7ba52"
color4_CI="#e7cb94"




def sigmoid(x,B,A,Chalf, b):#the fitting equation 
     y =B+((A-B)/(1+np.exp((-1/b)*(Chalf-x))))
     return y

infile=open(infile_locationAndName,"r")#"r"=reading mode
reader=csv.reader(infile)#create reader object
header=next(reader)





#base on replication numbers
if numOfReplicate=="1":
    xdata1=[]
    finalData=[] 
    for i in range(rep1_Statrt,rep1_End+1):#5 is the data start, 15 is the end of data+1
        flow=float(header[i])#making the input value to number
        xdata1.append(flow)
    for row in reader:
        # input y value
        ydata1=[]
        for i in range(rep1_Statrt,rep1_End+1):# it means row[5] to row[15]
            flow=float(row[i])#making the input value to number
            ydata1.append(flow)
        
        #baisc information input 
        title=row[rowLoc]
        ID=row[IDLoc]
        proteinName=row[NameLoc]
        Peptide=row[PeptideLoc]
        Start=row[StartLoc]
        End=row[EndLoc]
        numOfReporter=row[ReporterLoc]
        numOf_rep1_pts=row[rep1_numPtsLocation]
        
        #baisc information ouput 
        tempRow=[title]#[0]
        tempRow.append(ID)#[1]
        tempRow.append(proteinName)#[2]
        tempRow.append(Peptide)#[3]
        tempRow.append(Start)#[4]
        tempRow.append(End)#[5]
        tempRow.append(numOfReporter)#[6]
        tempRow.append(numOf_rep1_pts)#[7]
    
        try:#if statment, when it is not error   
            #normalized the data
            normalized_ydata1=[]
            normalized_xdata1=[]
                
            #normalizing rep1 dataset
            if denatureMethod=="CD":#first point is 0, last point is 1                
                for element in ydata1: # read from the begining, find first data !=0, make it A
                    if element !=0:
                        A=element# set 0
                        break
                for element in reversed(ydata1): #read from the end, find first data !=0, make it B
                    if element !=0:
                        B=element#set 1
                        #print ("A",A)
                        break
                for i in range(len(ydata1)):
                    if ydata1[i] !=0:
                        element=(ydata1[i]-A)/(B-A) #3 for each data !=0, normalized_ydata=(data-A)/(B-A)
                        normalized_ydata1.append(element)
                        normalized_xdata1.append(xdata1[i]) #4 for each normalized_ydata, find the xdata, make it normalized_xdata
                    else:        
                        normalized_ydata1.append("error")
                        normalized_xdata1.append(xdata1[i])  
                        
            if denatureMethod=="HD":#first point is 1, last point is 0 
                for element in ydata1: #1 read from the begining, find first data !=0, make it B
                    if element !=0:
                        B=element# set 1
                        break
                for element in reversed(ydata1): #2 read from the end, find first data !=0, make it A
                    if element !=0:
                        A=element#set 0
                        break
                for i in range(len(ydata1)):
                    if ydata1[i] !=0:
                        element=(ydata1[i]-A)/(B-A) #3 for each data !=0, normalized_ydata=(data-A)/(B-A)
                        normalized_ydata1.append(element)
                        normalized_xdata1.append(xdata1[i]) #4 for each normalized_ydata, find the xdata, make it normalized_xdata
                    else:        
                        normalized_ydata1.append("error")
                        normalized_xdata1.append(xdata1[i])        
                        
            #Set data set for fit and plot
            xplot1=[]
            yplot1=[]
            
            for i in range (len(normalized_ydata1)):   
                if normalized_ydata1[i]!="error":
                    yplot1.append(normalized_ydata1[i])
                    xplot1.append(normalized_xdata1[i])                    
            
            #fit Rep1 to s curve
            rep1_fit_result_output=[]
            popt1, pcov1 = curve_fit(sigmoid, xplot1, yplot1)
            rep1_fit_result_output.extend(popt1)#extend merge two list to one big list  
                                                #popt has all the parameter (variable) in regression equation (x,B, A, Chalf, b)
                                                                                                        #row [-, 0, 1, 2,    3] in  rep1_fit_result_output   
            if denatureMethod=="CD":
                x1 = np.linspace(0, 4, 50)
            if denatureMethod=="HD": 
                x1 = np.linspace(0, 1.2, 50)
            y1= sigmoid(x1, *popt1)# *popt split the two variable 
            rep1_Chalf=rep1_fit_result_output[2]
            rep1_b=rep1_fit_result_output[3]
            
            #calculate AL confidence interval
            numOfPoint_rep1_repA=float(numOf_rep1_pts)
            sumPoints_AL=numOfPoint_rep1_repA
            tempRow.append(sumPoints_AL)#[6]
            rep1_Standard_error=np.sqrt(np.diag(pcov1))
            rep1_fit_result_output.extend(rep1_Standard_error)#[4,5,6,7]=[Berror, Aerror, Chalferror, berror]
            rep1_ConfidenceIntervrep1_chalf=t.ppf(.975, (sumPoints_AL-1))* rep1_fit_result_output[6]/np.sqrt(sumPoints_AL)
            rep1_fit_result_output.append(rep1_ConfidenceIntervrep1_chalf)#[8]
            if denatureMethod=="CD":
                CI_ratioTOrange=rep1_ConfidenceIntervrep1_chalf/CD_conc_range
                rep1_fit_result_output.append(CI_ratioTOrange)#[9]
            if denatureMethod=="HD":
                CI_ratioTOrange=rep1_ConfidenceIntervrep1_chalf/HD_temp_range
                rep1_fit_result_output.append(CI_ratioTOrange)#[9]
            rep1_CI_Chalf_lowbound=rep1_Chalf-rep1_fit_result_output[8]
            rep1_CI_Chalf_upbound=rep1_Chalf+rep1_fit_result_output[8]
            rep1_fit_result_output.append(rep1_CI_Chalf_lowbound)#[10]
            rep1_fit_result_output.append(rep1_CI_Chalf_upbound)#[11]
            rep1_ConfidenceIntervrep1_b=t.ppf(.975, (sumPoints_AL-1))* rep1_fit_result_output[7]/np.sqrt(sumPoints_AL)
            rep1_fit_result_output.append(rep1_ConfidenceIntervrep1_b)#[12]
            rep1_CI_b_lowbound=rep1_b-rep1_fit_result_output[12]
            rep1_CI_b_upbound=rep1_b+rep1_fit_result_output[12]
            rep1_fit_result_output.append(rep1_CI_b_lowbound)#[13]
            rep1_fit_result_output.append(rep1_CI_b_upbound)#[14]
            
            #compute AL r squared
            residuals=(yplot1)-sigmoid((xplot1),*popt1)
            ss_res=np.sum(residuals**2)
            ss_tot=np.sum(((yplot1)-np.mean(yplot1))**2)
            r_squared1=1-(ss_res/ss_tot)
            rep1_fit_result_output.append(r_squared1)#[15]; append add something to the list   
            if denatureMethod=="HD":
                Chalf1_temp=(rep1_Chalf*HD_temp_range)+HD_Temp_lowbound
                rep1_fit_result_output.append(Chalf1_temp)#[16]
            if denatureMethod=="CD":
                Chalf1_normalized=rep1_Chalf/CD_conc_range
                rep1_fit_result_output.append(Chalf1_normalized)#[16]
                

            
            tempRow.extend(rep1_fit_result_output)
            tempRow.extend(normalized_ydata1)

        
            #generate the figure
            if CI_ratioTOrange<CI_cutoff:
                graphName=descriptionOfTheRun+"_"+"row#"+title+" (has CI)"+"\n"+ID +"||"+ proteinName+"\n"+Peptide
            else:
                graphName=descriptionOfTheRun+"_"+"row#"+title+"\n"+ID +"||"+ proteinName+"\n"+Peptide
                
            #plot color
            if plotColor=="1":#choose from blue, orange, purple, brown
                plt.plot(xplot1, yplot1,color=color1_dataAndFit,ls=':',marker='o',label='rep1_'+rep1_notebookCode)
                #plt.plot(xplot2, yplot2,color=color1_dataAndFit,ls=':',marker='^',label='rep1_1057')
                pylab.plot(x1,y1,color1_dataAndFit, label='fit')
                pylab.axvline(rep1_Chalf,color=color1_Chalf,ls='-',label='C1/2')#chalf
                if CI_ratioTOrange<CI_cutoff:
                    plt.axvspan((rep1_Chalf-rep1_fit_result_output[8]), (rep1_Chalf+rep1_fit_result_output[8]),color=color1_CI)#CI
            if plotColor=="2":#choose from blue, orange, purple, brown
                plt.plot(xplot1, yplot1,color=color2_dataAndFit,ls=':',marker='o',label='rep1_'+rep1_notebookCode)
                #plt.plot(xplot2, yplot2,color=color2_dataAndFit,ls=':',marker='^',label='rep1_1057')
                pylab.plot(x1,y1,color2_dataAndFit, label='fit')
                pylab.axvline(rep1_Chalf,color=color2_Chalf,ls='-',label='C1/2')#chalf
                if CI_ratioTOrange<CI_cutoff:
                    plt.axvspan((rep1_Chalf-rep1_fit_result_output[8]), (rep1_Chalf+rep1_fit_result_output[8]),color=color2_CI)#CI
            if plotColor=="3":#choose from blue, orange, purple, brown
                plt.plot(xplot1, yplot1,color=color3_dataAndFit,ls=':',marker='o',label='rep1_'+rep1_notebookCode)
                #plt.plot(xplot2, yplot2,color=color3_dataAndFit,ls=':',marker='^',label='rep1_1057')
                pylab.plot(x1,y1,color3_dataAndFit, label='fit')
                pylab.axvline(rep1_Chalf,color=color3_Chalf,ls='-',label='C1/2')#chalf
                if CI_ratioTOrange<CI_cutoff:
                    plt.axvspan((rep1_Chalf-rep1_fit_result_output[8]), (rep1_Chalf+rep1_fit_result_output[8]),color=color3_CI)#CI
            if plotColor=="4":#choose from blue, orange, purple, brown
                plt.plot(xplot1, yplot1,color=color4_dataAndFit,ls=':',marker='o',label='rep1_'+rep1_notebookCode)
                #plt.plot(xplot2, yplot2,color=color4_dataAndFit,ls=':',marker='^',label='rep1_1057')
                pylab.plot(x1,y1,color4_dataAndFit, label='fit')
                pylab.axvline(rep1_Chalf,color=color4_Chalf,ls='-',label='C1/2')#chalf
                if CI_ratioTOrange<CI_cutoff:
                    plt.axvspan((rep1_Chalf-rep1_fit_result_output[8]), (rep1_Chalf+rep1_fit_result_output[8]),color=color4_CI)#CI
            
            #plot label
            if denatureMethod=="CD":
                pylab.xlabel('denaturant concentration')
                pylab.ylabel('% labeled')
            if denatureMethod=="HD": 
                pylab.xlabel('normalized temperature')
                pylab.ylabel('% soluble')
            #pylab.ylim(0.3, 2.2)
            pylab.legend(loc='best')
            pylab.title(graphName,fontsize=10)
            pylab.savefig(outImage_location+title+"_"+proteinName+".jpg")
            pylab.clf()#clearfig, start a new piece , make sure do it before or after a set of instruction
    
        except :#define the error
            print (title)         
    
        finalData.append(tempRow)  
        
        


    #output the data file        
    outfile=open(outfile_locationAndName,"w",newline="")        
    writer=csv.writer(outfile)#create write object
    newheader=["row","protein ID","ProteinNAme","Peptide","Start","End","numOfReporter","numOf_rep1_pts","sumPoints"]
    newheader_rep1_fitOutPut_CD=["B","A","Chalf","b","B_err","A_err","Chalf_err","b_err","CI_Chalf","CIratioTOrange","CI_Chalf_low","CI_Chalf_up","CI_b","CI_b_low","CI_b_up","r_square","Chalf_normalized"]
    newheader_rep1_fitOutPut_HD=["B","A","Chalf","b","B_err","A_err","Chalf_err","b_err","CI_Chalf","CIratioTOrange","CI_Chalf_low","CI_Chalf_up","CI_b","CI_b_low","CI_b_up","r_square","Chalf_temp"]
    if denatureMethod=="CD":
        newheader.extend(newheader_rep1_fitOutPut_CD)
        normalized_header1=CD_xaxis_FinalGdmClConcentration_value
        newheader.extend(normalized_header1)
    if denatureMethod=="HD": 
        newheader.extend(newheader_rep1_fitOutPut_HD)
        normalized_header1=HD_xasis_normalizedTemp_value    
        newheader.extend(normalized_header1)
    
    writer.writerow(newheader)
    for addvariable in finalData:
        writer.writerow(addvariable)
    outfile.close()
    infile.close()

    
