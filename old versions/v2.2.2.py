import numpy as np
import pylab
import csv
from scipy.optimize import curve_fit
from scipy.stats import t
import matplotlib.pyplot as plt

# basic information --set up variable
# method infodisgunsh
denatureMethod = "HD" # choose 'CD' for chemical denature, 'HD' for heat denature
CD_xaxis_FinalGdmClConcentration_value = [0, 0.43, 0.87, 1.3, 1.74, 2.17, 2.61, 3.04, 3.48]
CD_conc_upbound = 3.48
CD_conc_lowbound = 0
CD_conc_range = CD_conc_upbound - CD_conc_lowbound
HD_xasis_normalizedTemp_value = [0, 0.129277567, 0.224334601, 0.338403042, 0.429657795, 0.539923954, 0.646387833, 0.775665399, 0.882129278, 1]
HD_Temp_upbound = 63
HD_Temp_lowbound = 36.7
HD_temp_range = HD_Temp_upbound - HD_Temp_lowbound
CI_cutoff = 0.3

# run info
descriptionOfTheRun = "Heat.Denature High Protein Diet" # what ever desciption can distinguish
numOfCondition = "2" # number of replicate used to generate 1 fit
condition1_plotColor = "1" # 1=blue,2=orange,3=purple,4=brown
condition2_plotColor = "2" # 1=blue,2=orange,3=purple,4=brown

# hard-coded dictionaries; change to dynamic input later
condition = {
    1: {'description' : 'Ad Libitum', 'plot color' : "1",       # condition 1
        1 : {'notebook code' : '1085',                          # condition 1, replicate 1
             'numPtsLocation' : 4,                              # input csv column for numPtsLocation
             'Start' : 5,                                       # input csv column for Start
             'End' : 14,                                        # input csv column for End
             }},                     
    2: {'description' : 'Dietary restriction', 'plot color' : 2}}

# file selection
infile_locationAndName = "../test_multiInput.csv"
outfile_locationAndName = "../test_multiOutput_2.csv"
outImage_location = "../testoutput2_2//"

# info for each peptide from infile
rowLoc = 0
IDLoc = 1
NameLoc = 2
PeptideLoc = 3
# ReporterLoc=4

# condition1 info
condition1_description = "Ad Libitum" # the notebook code of the run
# condition1 rep1
condi1_rep1_notebookcode = "1085"
condi1_rep1_numPtsLocation = 4 # type in the number get from excel column to python raw translator
condi1_rep1_Start = 5 # type in the number get from excel column to python raw translator
condi1_rep1_End = 14 # type in the number get from excel column to python raw translator
# condition1 rep2
condi1_rep2_notebookcode = "1091"
condi1_rep2_numPtsLocation = 15 # type in the number get from excel column to python raw translator
condi1_rep2_Start = 16 # type in the number get from excel column to python raw translator
condi1_rep2_End = 25 # type in the number get from excel column to python raw translator
# condition1 rep3
condi1_rep3_notebookcode = "1105"
condi1_rep3_numPtsLocation = 26 # type in the number get from excel column to python raw translator
condi1_rep3_Start = 27 # type in the number get from excel column to python raw translator
condi1_rep3_End = 36 # type in the number get from excel column to python raw translator
# condition1 rep4
condi1_rep4_notebookcode = "1107"
condi1_rep4_numPtsLocation = 37 # type in the number get from excel column to python raw translator
condi1_rep4_Start = 38 # type in the number get from excel column to python raw translator
condi1_rep4_End = 47 # type in the number get from excel column to python raw translator

#condition2 info
condition2_description="Dietary restriction" #the notebook code of the run
#condition2 rep1
condi2_rep1_notebookcode="1085"
condi2_rep1_numPtsLocation=48# type in the number get from excel column to python raw translator
condi2_rep1_Start=49 # type in the number get from excel column to python raw translator
condi2_rep1_End=58# type in the number get from excel column to python raw translator
#condition2 rep2
condi2_rep2_notebookcode="1091"
condi2_rep2_numPtsLocation=59# type in the number get from excel column to python raw translator
condi2_rep2_Start=60 # type in the number get from excel column to python raw translator
condi2_rep2_End=69# type in the number get from excel column to python raw translator
#condition2 rep3
condi2_rep3_notebookcode="1105"
condi2_rep3_numPtsLocation=70# type in the number get from excel column to python raw translator
condi2_rep3_Start=71 # type in the number get from excel column to python raw translator
condi2_rep3_End=80# type in the number get from excel column to python raw translator
#condition2 rep4
condi2_rep4_notebookcode="1107"
condi2_rep4_numPtsLocation=81# type in the number get from excel column to python raw translator
condi2_rep4_Start=82 # type in the number get from excel column to python raw translator
condi2_rep4_End=91# type in the number get from excel column to python raw translator



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

def sigmoid(x, B, A, Chalf, b): # the fitting equation
    """Fitting equation. Fits data to sigmoid curve.
    
    Returns y
    """
    y = B + ((A - B) / (1 + np.exp((-1 / b) * (Chalf - x))))
    return y

def read_xdata(xdata, start, end):    # line 195
    """Reads x data from csv input file"""
    #condition 1 info
    for i in range(start, end + 1): # 5 is the data start, 15 is the end of data+1
        flow = float(header[i]) # making the input value to number
        xdata.append(flow)
        
def read_ydata(ydata, start, end):
    """Reads y data from csv input file"""
    for i in range(start, end + 1): # it means row[5] to row[15]
        flow = float(row[i]) # making the input value to number
        ydata.append(flow)
        
def normalize_data(xdata, ydata, normalized_xdata, normalized_ydata, tempRow):
    """Normalizes data"""
    # normalizing condition1 dataset
    if denatureMethod == "CD": # first point is 0, last point is 1      
        # condition1 rep1          
        for element in ydata: # read from the begining, find first data !=0, make it A
            if element != 0:
                A = element # set 0
                break
        for element in reversed(ydata): # read from the end, find first data !=0, make it B
            if element != 0:
                B = element # set 1
                break
        for i in range(len(ydata)):
            if ydata[i] != 0:
                element = (ydata[i] - A) / (B - A) # 3 for each data !=0, normalized_ydata=(data-A)/(B-A)
                normalized_ydata.append(element)
                normalized_xdata.append(xdata[i]) # 4 for each normalized_ydata, find the xdata, make it normalized_xdata
            else:        
                normalized_ydata.append("error")
                normalized_xdata.append(xdata[i])  
                
        tempRow.extend(normalized_ydata)
        
    elif denatureMethod == "HD": # first point is 1, last point is 0 
        # condition1 rep1       
        for element in ydata: # 1 read from the begining, find first data !=0, make it B
            if element != 0:
                B = element # set 1
                break
        for element in reversed(ydata): # 2 read from the end, find first data !=0, make it A
            if element != 0:
                A = element # set 0
                break
        for i in range(len(ydata)):
            if ydata[i] != 0:
                element = (ydata[i] - A) / (B - A) #3 for each data !=0, normalized_ydata=(data-A)/(B-A)
                normalized_ydata.append(element)
                normalized_xdata.append(xdata[i]) #4 for each normalized_ydata, find the xdata, make it normalized_xdata
            else:        
                normalized_ydata.append("error")
                normalized_xdata.append(xdata[i])
                
        tempRow.extend(normalized_ydata)
        
def plot_data(xplot, yplot, normalized_xdata, normalized_ydata):
    """Begins compiling data for plotting"""
    for i in range (len(normalized_ydata)):   
        if normalized_ydata[i] != "error":
            yplot.append(normalized_ydata[i])
            xplot.append(normalized_xdata[i])

infile = open(infile_locationAndName, "r") # "r"=reading mode
reader = csv.reader(infile) # create reader object
header = next(reader)


print(infile_locationAndName)


#base on replication numbers

if numOfCondition=="2":
    condition[1][1]['xdata'] = []
    xdata1_2=[]
    xdata1_3=[]
    xdata1_4=[]

    xdata2_1=[]
    xdata2_2=[]
    xdata2_3=[]
    xdata2_4=[]

    finalData=[]
    
    #condition 1 info
    read_xdata(condition[1][1]['xdata'], condition[1][1]['Start'], condition[1][1]['End'])
    read_xdata(xdata1_2, condi1_rep2_Start, condi1_rep2_End)
    read_xdata(xdata1_3, condi1_rep3_Start, condi1_rep3_End)
    read_xdata(xdata1_4, condi1_rep4_Start, condi1_rep4_End)

    #condition 2 info
    read_xdata(xdata2_1, condi2_rep1_Start, condi2_rep1_End)
    read_xdata(xdata2_2, condi2_rep2_Start, condi2_rep2_End)
    read_xdata(xdata2_3, condi2_rep3_Start, condi2_rep3_End)
    read_xdata(xdata2_4, condi2_rep4_Start, condi2_rep4_End)

    for row in reader:
        try:      
            #condition 1 input y value
            condition[1][1]['ydata'] = []
            ydata1_2=[]
            ydata1_3=[]
            ydata1_4=[]
            
            read_ydata(condition[1][1]['ydata'], condition[1][1]['Start'], condition[1][1]['End'])
            read_ydata(ydata1_2, condi1_rep2_Start, condi1_rep2_End)
            read_ydata(ydata1_3, condi1_rep3_Start, condi1_rep3_End)
            read_ydata(ydata1_4, condi1_rep4_Start, condi1_rep4_End)
    
            #condition 2 input y value
            ydata2_1=[]
            ydata2_2=[]
            ydata2_3=[]
            ydata2_4=[]
            
            read_ydata(ydata2_1, condi2_rep1_Start, condi2_rep1_End)
            read_ydata(ydata2_2, condi2_rep2_Start, condi2_rep2_End)
            read_ydata(ydata2_3, condi2_rep3_Start, condi2_rep3_End)
            read_ydata(ydata2_4, condi2_rep4_Start, condi2_rep4_End)
        except:
            pass
                  
        #basic information input 
        title = row[rowLoc]
        ID = row[IDLoc]
        proteinName = row[NameLoc]
        Peptide = row[PeptideLoc]
        
        # numOfReporter = row[ReporterLoc]
        numOf_condition1_1_pts = row[condition[1][1]['numPtsLocation']]
        numOf_condition1_2_pts = row[condi1_rep2_numPtsLocation]
        numOf_condition1_3_pts = row[condi1_rep3_numPtsLocation]
        numOf_condition1_4_pts = row[condi1_rep4_numPtsLocation]
        numOf_condition2_1_pts = row[condi2_rep1_numPtsLocation]
        numOf_condition2_2_pts = row[condi2_rep2_numPtsLocation]
        numOf_condition2_3_pts = row[condi2_rep3_numPtsLocation]
        numOf_condition2_4_pts = row[condi2_rep4_numPtsLocation]
        
        # baisc information ouput 
        tempRow = [title] # [0]
        tempRow.append(ID) # [1]
        tempRow.append(proteinName) # [2]
        tempRow.append(Peptide) # [3]
        # tempRow.append(numOfReporter)#[4]
        
    
        try: # if statment, when it is not error   
            # normalized the data
            condition[1][1]['normalized_ydata'] = []
            condition[1][1]['normalized_xdata'] = []
            normalized_ydata1_2=[]
            normalized_xdata1_2=[]
            normalized_ydata1_3=[]
            normalized_xdata1_3=[]
            normalized_ydata1_4=[]
            normalized_xdata1_4=[]
            
            normalized_ydata2_1=[]
            normalized_xdata2_1=[]
            normalized_ydata2_2=[]
            normalized_xdata2_2=[]
            normalized_ydata2_3=[]
            normalized_xdata2_3=[]
            normalized_ydata2_4=[]
            normalized_xdata2_4=[]
                
            # normalizing condition1 dataset
            normalize_data(condition[1][1]['xdata'], condition[1][1]['ydata'], condition[1][1]['normalized_xdata'], condition[1][1]['normalized_ydata'], tempRow)
            normalize_data(xdata1_2, ydata1_2, normalized_xdata1_2, normalized_ydata1_2, tempRow)
            normalize_data(xdata1_3, ydata1_3, normalized_xdata1_3, normalized_ydata1_3, tempRow)
            normalize_data(xdata1_4, ydata1_4, normalized_xdata1_4, normalized_ydata1_4, tempRow)
            
            #normalizing condition2 dataset
            normalize_data(xdata2_1, ydata2_1, normalized_xdata2_1, normalized_ydata2_1, tempRow)
            normalize_data(xdata2_2, ydata2_2, normalized_xdata2_2, normalized_ydata2_2, tempRow)
            normalize_data(xdata2_3, ydata2_3, normalized_xdata2_3, normalized_ydata2_3, tempRow)
            normalize_data(xdata2_4, ydata2_4, normalized_xdata2_4, normalized_ydata2_4, tempRow)
            
            #Set data set for fit and plot
            condition[1][1]['xplot'] = []
            condition[1][1]['yplot'] = []
            xplot1_2=[]
            yplot1_2=[]
            xplot1_3=[]
            yplot1_3=[]
            xplot1_4=[]
            yplot1_4=[]
            
            xplot2_1=[]
            yplot2_1=[]
            xplot2_2=[]
            yplot2_2=[]
            xplot2_3=[]
            yplot2_3=[]
            xplot2_4=[]
            yplot2_4=[]
            
            #condition1
            plot_data(condition[1][1]['xplot'], condition[1][1]['yplot'], condition[1][1]['normalized_xdata'], condition[1][1]['normalized_ydata'])
            """for i in range (len(normalized_ydata1_1)):   
                if normalized_ydata1_1[i]!="error":
                    condition[1][1]['yplot'].append(normalized_ydata1_1[i])
                    condition[1][1]['xplot'].append(normalized_xdata1_1[i])"""
            
            for i in range (len(normalized_ydata1_2)):   
                if normalized_ydata1_2[i]!="error":
                    yplot1_2.append(normalized_ydata1_2[i])
                    xplot1_2.append(normalized_xdata1_2[i])  

            for i in range (len(normalized_ydata1_3)):   
                if normalized_ydata1_3[i]!="error":
                    yplot1_3.append(normalized_ydata1_3[i])
                    xplot1_3.append(normalized_xdata1_3[i]) 

            for i in range (len(normalized_ydata1_4)):   
                if normalized_ydata1_4[i]!="error":
                    yplot1_4.append(normalized_ydata1_4[i])
                    xplot1_4.append(normalized_xdata1_4[i])  
                    
            #condition2
            for i in range (len(normalized_ydata2_1)):   
                if normalized_ydata2_1[i]!="error":
                    yplot2_1.append(normalized_ydata2_1[i])
                    xplot2_1.append(normalized_xdata2_1[i])                    
            
            for i in range (len(normalized_ydata2_2)):   
                if normalized_ydata2_2[i]!="error":
                    yplot2_2.append(normalized_ydata2_2[i])
                    xplot2_2.append(normalized_xdata2_2[i]) 

            for i in range (len(normalized_ydata2_3)):   
                if normalized_ydata2_3[i]!="error":
                    yplot2_3.append(normalized_ydata2_3[i])
                    xplot2_3.append(normalized_xdata2_3[i]) 
            
            for i in range (len(normalized_ydata2_4)):   
                if normalized_ydata2_4[i]!="error":
                    yplot2_4.append(normalized_ydata2_4[i])
                    xplot2_4.append(normalized_xdata2_4[i]) 
 
                   
            # fit condition1 to s curve
            condition1_fit_result_output = []
            popt1, pcov1 = curve_fit(sigmoid, condition[1][1]['xplot'] + xplot1_2 + xplot1_3 + xplot1_4, condition[1][1]['yplot'] + yplot1_2 + yplot1_3 + yplot1_4)
            condition1_fit_result_output.extend(popt1)#extend merge two list to one big list  
                                                #popt has all the parameter (variable) in regression equation (x,B, A, Chalf, b)                                                                                            #row [-, 0, 1, 2,    3] in  condition1_fit_result_output   
            if denatureMethod == "CD":
                x1 = np.linspace(0, 4, 50)
            if denatureMethod == "HD": 
                x1 = np.linspace(0, 1.2, 50)
            y1 = sigmoid(x1, *popt1)# *popt split the two variable 
            condition1_Chalf = condition1_fit_result_output[2]
            condition1_b = condition1_fit_result_output[3]
            
         
            # calculate condition1 confidence interval
            numOfPoint_condition1_repA = float(numOf_condition1_1_pts)
            numOfPoint_condition1_repB = float(numOf_condition1_2_pts)
            numOfPoint_condition1_repC = float(numOf_condition1_3_pts)
            numOfPoint_condition1_repD = float(numOf_condition1_4_pts)
            sumPoints_condition1 = numOfPoint_condition1_repA + numOfPoint_condition1_repB + numOfPoint_condition1_repC + numOfPoint_condition1_repD
            tempRow.append(numOf_condition1_1_pts) # [5]
            tempRow.append(numOf_condition1_2_pts) # [5]
            tempRow.append(numOf_condition1_3_pts) # [5]
            tempRow.append(numOf_condition1_4_pts) # [5]
            tempRow.append(sumPoints_condition1)   # [6]

            condition1_Standard_error = np.sqrt(np.diag(pcov1))
            condition1_fit_result_output.extend(condition1_Standard_error)#[4,5,6,7]=[Berror, Aerror, Chalferror, berror]
            condition1_ConfidenceIntervcondition1_chalf = t.ppf(.975, (sumPoints_condition1-1)) * condition1_fit_result_output[6] / np.sqrt(sumPoints_condition1)
            condition1_fit_result_output.append(condition1_ConfidenceIntervcondition1_chalf)#[8]
            if denatureMethod=="CD":
                CI_ratioTOrange1=condition1_ConfidenceIntervcondition1_chalf/CD_conc_range
                condition1_fit_result_output.append(CI_ratioTOrange1)#[9]
            if denatureMethod=="HD":
                CI_ratioTOrange1 = condition1_ConfidenceIntervcondition1_chalf/1
                condition1_fit_result_output.append(CI_ratioTOrange1)#[9]
            condition1_CI_Chalf_lowbound = condition1_Chalf - condition1_fit_result_output[8]
            condition1_CI_Chalf_upbound = condition1_Chalf + condition1_fit_result_output[8]
            condition1_fit_result_output.append(condition1_CI_Chalf_lowbound)#[10]
            condition1_fit_result_output.append(condition1_CI_Chalf_upbound)#[11]
            condition1_ConfidenceIntervcondition1_b=t.ppf(.975, (sumPoints_condition1-1))* condition1_fit_result_output[7]/np.sqrt(sumPoints_condition1)
            condition1_fit_result_output.append(condition1_ConfidenceIntervcondition1_b)#[12]
            condition1_CI_b_lowbound=condition1_b-condition1_fit_result_output[12]
            condition1_CI_b_upbound=condition1_b+condition1_fit_result_output[12]
            condition1_fit_result_output.append(condition1_CI_b_lowbound)#[13]
            condition1_fit_result_output.append(condition1_CI_b_upbound)#[14]
            
            #compute condition1 r squared
            residuals=(condition[1][1]['yplot']+yplot1_2+yplot1_3+yplot1_4)-sigmoid((condition[1][1]['xplot']+xplot1_2+xplot1_3+xplot1_4),*popt1)
            ss_res=np.sum(residuals**2)
            ss_tot=np.sum(((condition[1][1]['yplot']+yplot1_2+yplot1_3+yplot1_4)-np.mean(condition[1][1]['yplot']+yplot1_2+yplot1_3+yplot1_4))**2)
            r_squared1=1-(ss_res/ss_tot)
            condition1_fit_result_output.append(r_squared1)#[15]; append add something to the list   
            if denatureMethod=="HD":
                Chalf1_temp=(condition1_Chalf*HD_temp_range)+HD_Temp_lowbound
                condition1_fit_result_output.append(Chalf1_temp)#[16]
            if denatureMethod=="CD":
                Chalf1_normalized=condition1_Chalf/CD_conc_range
                condition1_fit_result_output.append(Chalf1_normalized)#[16]
                
            
            
            tempRow.extend(condition1_fit_result_output)

            
             #fit condition2 to s curve
            condition2_fit_result_output=[]
            popt2, pcov2 = curve_fit(sigmoid, xplot2_1+xplot2_2+xplot2_3+xplot2_4, yplot2_1+yplot2_2+yplot2_3+yplot2_4)
            condition2_fit_result_output.extend(popt2)#extend merge two list to one big list  
                                                #popt has all the parameter (variable) in regression equation (x,B, A, Chalf, b)
                                                                                           #row [-, 0, 1, 2,    3] in  condition1_fit_result_output   
           
            if denatureMethod=="CD":
                x2 = np.linspace(0, 4, 50)
            if denatureMethod=="HD": 
                x2 = np.linspace(0, 1.2, 50)
            y2= sigmoid(x2, *popt2)# *popt split the two variable 
            condition2_Chalf=condition2_fit_result_output[2]
            condition2_b=condition2_fit_result_output[3]
            
            
            #calculate condition2 confidence interval
            numOfPoint_condition2_repA=float(numOf_condition2_1_pts)
            numOfPoint_condition2_repB=float(numOf_condition2_2_pts)
            numOfPoint_condition2_repC=float(numOf_condition2_3_pts)
            numOfPoint_condition2_repD=float(numOf_condition2_4_pts)
            sumPoints_condition2=numOfPoint_condition2_repA+numOfPoint_condition2_repB+numOfPoint_condition2_repC+numOfPoint_condition2_repD
            tempRow.append(numOf_condition2_1_pts)#[6]
            tempRow.append(numOf_condition2_2_pts)#[6]
            tempRow.append(numOf_condition2_3_pts)#[6]
            tempRow.append(numOf_condition2_4_pts)#[6]
            tempRow.append(sumPoints_condition2)#[6]

            condition2_Standard_error=np.sqrt(np.diag(pcov2))
            condition2_fit_result_output.extend(condition2_Standard_error)#[4,5,6,7]=[Berror, Aerror, Chalferror, berror]
            condition2_ConfidenceIntervcondition2_chalf=t.ppf(.975, (sumPoints_condition2-1))* condition2_fit_result_output[6]/np.sqrt(sumPoints_condition2)
            condition2_fit_result_output.append(condition2_ConfidenceIntervcondition2_chalf)#[8]
            if denatureMethod=="CD":
                CI_ratioTOrange2=condition2_ConfidenceIntervcondition2_chalf/CD_conc_range
                condition2_fit_result_output.append(CI_ratioTOrange2)#[9]
            if denatureMethod=="HD":
                CI_ratioTOrange2=condition2_ConfidenceIntervcondition2_chalf/1
                condition2_fit_result_output.append(CI_ratioTOrange2)#[9]
            condition2_CI_Chalf_lowbound=condition2_Chalf-condition2_fit_result_output[8]
            condition2_CI_Chalf_upbound=condition2_Chalf+condition2_fit_result_output[8]
            condition2_fit_result_output.append(condition2_CI_Chalf_lowbound)#[10]
            condition2_fit_result_output.append(condition2_CI_Chalf_upbound)#[11]
            condition2_ConfidenceIntervcondition2_b=t.ppf(.975, (sumPoints_condition2-1))* condition2_fit_result_output[7]/np.sqrt(sumPoints_condition2)
            condition2_fit_result_output.append(condition2_ConfidenceIntervcondition2_b)#[12]
            condition2_CI_b_lowbound=condition2_b-condition2_fit_result_output[12]
            condition2_CI_b_upbound=condition2_b+condition2_fit_result_output[12]
            condition2_fit_result_output.append(condition2_CI_b_lowbound)#[13]
            condition2_fit_result_output.append(condition2_CI_b_upbound)#[14]
            
            #compute condition2 r squared
            residuals2=(yplot2_1+yplot2_2+yplot2_3+yplot2_4)-sigmoid((xplot2_1+xplot2_2+xplot2_3+xplot2_4),*popt2)
            ss_res2=np.sum(residuals2**2)
            ss_tot2=np.sum(((yplot2_1+yplot2_2+yplot2_3+yplot2_4)-np.mean(yplot2_1+yplot2_2+yplot2_3+yplot2_4))**2)
            r_squared2=1-(ss_res2/ss_tot2)
            condition2_fit_result_output.append(r_squared2)#[15]; append add something to the list   
            if denatureMethod=="HD":
                Chalf2_temp=(condition2_Chalf*HD_temp_range)+HD_Temp_lowbound
                condition2_fit_result_output.append(Chalf2_temp)#[16]
            if denatureMethod=="CD":
                Chalf2_normalized=condition2_Chalf/CD_conc_range
                condition2_fit_result_output.append(Chalf2_normalized)#[16]
                

            tempRow.extend(condition2_fit_result_output)



            #delta calculation 
            if denatureMethod=="CD":
                delta_Chalf=condition2_Chalf-condition1_Chalf
                delta_Chalf_percent=delta_Chalf/condition1_Chalf
                delta_b=condition2_b-condition1_b
                delta_b_percent=delta_b/condition1_b
            if denatureMethod=="HD":
                delta_Chalf=Chalf2_temp-Chalf1_temp
                delta_Chalf_percent=delta_Chalf/Chalf1_temp
                delta_b=condition2_b-condition1_b
                delta_b_percent=delta_b/condition1_b

            tempRow.append(delta_Chalf)
            tempRow.append(delta_Chalf_percent)
            tempRow.append(delta_b)
            tempRow.append(delta_b_percent)

            """#CI_chalf_significant=IF(AL>DR,IF(AND(AL>Dru,DR<Ald),"significant",""),if(AND(AL<DRd,DR>Alu),"significant","")
            if(condition1_Chalf>condition2_Chalf):
                if(condition1_Chalf>condition2_CI_Chalf_upbound and condition2_Chalf<condition1_CI_Chalf_lowbound):
                    CI_chalf_significant="significant"
            elif(condition1_Chalf<condition2_CI_Chalf_lowbound and condition2_Chalf>condition1_CI_Chalf_upbound):
                CI_chalf_significant = "significant"
            else:
                CI_chalf_significant = ""
            tempRow.append(CI_chalf_significant)

            #CI_chalf_overlap=IF(ALC>DRC,IF(ald<DRu,"overlap",""),IF(DRd<ALu,"overlap","")))
            if(condition1_Chalf>condition2_Chalf):
                if(condition1_CI_Chalf_lowbound<condition2_CI_Chalf_upbound):
                    CI_chlaf_overlap="overlap"
            elif(condition2_CI_Chalf_lowbound<condition1_CI_Chalf_upbound):
                CI_chlaf_overlap = "overlap"
            else:
                CI_chlaf_overlap = ""
            tempRow.append(CI_chlaf_overlap)"""

            
            
            
            #generate the figure
            if CI_ratioTOrange1<CI_cutoff and CI_ratioTOrange2<CI_cutoff:
                graphName=descriptionOfTheRun+"_"+"row#"+title+" (has CI)"+"\n"+ID +"||"+ proteinName+"\n"+Peptide
            else:
                graphName=descriptionOfTheRun+"_"+"row#"+title+"\n"+ID +"||"+ proteinName+"\n"+Peptide
                
            #condition1_plot color
            if condition1_plotColor=="1":#choose from blue, orange, purple, brown
                plt.plot(condition[1][1]['xplot'], condition[1][1]['yplot'],color=color1_dataAndFit,ls=':',marker='o',label=condition[1][1]['notebook code'])
                plt.plot(xplot1_2, yplot1_2,color=color1_dataAndFit,ls=':',marker='^',label=condi1_rep2_notebookcode)
                plt.plot(xplot1_3, yplot1_3,color=color1_dataAndFit,ls=':',marker='s',label=condi1_rep3_notebookcode)
                plt.plot(xplot1_4, yplot1_4,color=color1_dataAndFit,ls=':',marker='*',label=condi1_rep4_notebookcode)
                pylab.plot(x1,y1,color1_dataAndFit, label=condition[1]['description']+' fit')
                pylab.axvline(condition1_Chalf,color=color1_Chalf,ls='-',label='C1/2')#chalf
                if CI_ratioTOrange1<CI_cutoff:
                    plt.axvspan((condition1_Chalf-condition1_fit_result_output[8]), (condition1_Chalf+condition1_fit_result_output[8]),color=color1_CI)#CI
            if condition1_plotColor=="2":#choose from blue, orange, purple, brown
                plt.plot(condition[1][1]['xplot'], condition[1][1]['yplot'],color=color2_dataAndFit,ls=':',marker='o',label=condition[1][1]['notebook code'])
                plt.plot(xplot1_2, yplot1_2,color=color2_dataAndFit,ls=':',marker='^',label=condi1_rep2_notebookcode)
                plt.plot(xplot1_3, yplot1_3,color=color2_dataAndFit,ls=':',marker='s',label=condi1_rep3_notebookcode)
                plt.plot(xplot1_4, yplot1_4,color=color2_dataAndFit,ls=':',marker='*',label=condi1_rep4_notebookcode)
                pylab.plot(x1,y1,color2_dataAndFit, label=condition[1]['description']+' fit')
                pylab.axvline(condition1_Chalf,color=color2_Chalf,ls='-',label='C1/2')#chalf
                if CI_ratioTOrange1<CI_cutoff:
                    plt.axvspan((condition1_Chalf-condition1_fit_result_output[8]), (condition1_Chalf+condition1_fit_result_output[8]),color=color2_CI)#CI
            if condition1_plotColor=="3":#choose from blue, orange, purple, brown
                plt.plot(condition[1][1]['xplot'], condition[1][1]['yplot'],color=color3_dataAndFit,ls=':',marker='o',label=condition[1][1]['notebook code'])
                plt.plot(xplot1_2, yplot1_2,color=color3_dataAndFit,ls=':',marker='^',label=condi1_rep2_notebookcode)
                plt.plot(xplot1_3, yplot1_3,color=color3_dataAndFit,ls=':',marker='s',label=condi1_rep3_notebookcode)
                plt.plot(xplot1_4, yplot1_4,color=color3_dataAndFit,ls=':',marker='*',label=condi1_rep4_notebookcode)
                pylab.plot(x1,y1,color3_dataAndFit, label=condition[1]['description']+' fit')
                pylab.axvline(condition1_Chalf,color=color3_Chalf,ls='-',label='C1/2')#chalf
                if CI_ratioTOrange1<CI_cutoff:
                    plt.axvspan((condition1_Chalf-condition1_fit_result_output[8]), (condition1_Chalf+condition1_fit_result_output[8]),color=color3_CI)#CI
            if condition1_plotColor=="4":#choose from blue, orange, purple, brown
                plt.plot(condition[1][1]['xplot'], condition[1][1]['yplot'],color=color4_dataAndFit,ls=':',marker='o',label=condition[1][1]['notebook code'])
                plt.plot(xplot1_2, yplot1_2,color=color4_dataAndFit,ls=':',marker='^',label=condi1_rep2_notebookcode)
                plt.plot(xplot1_3, yplot1_3,color=color4_dataAndFit,ls=':',marker='s',label=condi1_rep3_notebookcode)
                plt.plot(xplot1_4, yplot1_4,color=color4_dataAndFit,ls=':',marker='*',label=condi1_rep4_notebookcode)
                pylab.plot(x1,y1,color4_dataAndFit, label=condition[1]['description']+' fit')
                pylab.axvline(condition1_Chalf,color=color4_Chalf,ls='-',label='C1/2')#chalf
                if CI_ratioTOrange1<CI_cutoff:
                    plt.axvspan((condition1_Chalf-condition1_fit_result_output[8]), (condition1_Chalf+condition1_fit_result_output[8]),color=color4_CI)#CI

            #condition2_plot color
            if condition2_plotColor=="1":#choose from blue, orange, purple, brown
                plt.plot(xplot2_1, yplot2_1,color=color1_dataAndFit,ls=':',marker='o',label=condi2_rep1_notebookcode)
                plt.plot(xplot2_2, yplot2_2,color=color1_dataAndFit,ls=':',marker='^',label=condi2_rep2_notebookcode)
                plt.plot(xplot2_3, yplot2_3,color=color1_dataAndFit,ls=':',marker='s',label=condi2_rep3_notebookcode)
                plt.plot(xplot2_4, yplot2_4,color=color1_dataAndFit,ls=':',marker='*',label=condi2_rep4_notebookcode)
                pylab.plot(x2,y2,color1_dataAndFit, label=condition2_description+' fit')
                pylab.axvline(condition2_Chalf,color=color1_Chalf,ls='-',label='C1/2')#chalf
                if CI_ratioTOrange2<CI_cutoff:
                    plt.axvspan((condition2_Chalf-condition2_fit_result_output[8]), (condition2_Chalf+condition2_fit_result_output[8]),color=color1_CI)#CI
            if condition2_plotColor=="2":#choose from blue, orange, purple, brown
                plt.plot(xplot2_1, yplot2_1,color=color2_dataAndFit,ls=':',marker='o',label=condi2_rep1_notebookcode)
                plt.plot(xplot2_2, yplot2_2,color=color2_dataAndFit,ls=':',marker='^',label=condi2_rep2_notebookcode)
                plt.plot(xplot2_3, yplot2_3,color=color2_dataAndFit,ls=':',marker='s',label=condi2_rep3_notebookcode)
                plt.plot(xplot2_4, yplot2_4,color=color2_dataAndFit,ls=':',marker='*',label=condi2_rep4_notebookcode)
                pylab.plot(x2,y2,color2_dataAndFit, label=condition2_description+' fit')
                pylab.axvline(condition2_Chalf,color=color2_Chalf,ls='-',label='C1/2')#chalf
                if CI_ratioTOrange2<CI_cutoff:
                    plt.axvspan((condition2_Chalf-condition2_fit_result_output[8]), (condition2_Chalf+condition2_fit_result_output[8]),color=color2_CI)#CI
            if condition2_plotColor=="3":#choose from blue, orange, purple, brown
                plt.plot(xplot2_1, yplot2_1,color=color3_dataAndFit,ls=':',marker='o',label=condi2_rep1_notebookcode)
                plt.plot(xplot2_2, yplot2_2,color=color3_dataAndFit,ls=':',marker='^',label=condi2_rep2_notebookcode)
                plt.plot(xplot2_3, yplot2_3,color=color3_dataAndFit,ls=':',marker='s',label=condi2_rep3_notebookcode)
                plt.plot(xplot2_4, yplot2_4,color=color3_dataAndFit,ls=':',marker='*',label=condi2_rep4_notebookcode)
                pylab.plot(x2,y2,color3_dataAndFit, label=condition2_description+' fit')
                pylab.axvline(condition2_Chalf,color=color3_Chalf,ls='-',label='C1/2')#chalf
                if CI_ratioTOrange2<CI_cutoff:
                    plt.axvspan((condition2_Chalf-condition2_fit_result_output[8]), (condition2_Chalf+condition2_fit_result_output[8]),color=color3_CI)#CI
            if condition2_plotColor=="4":#choose from blue, orange, purple, brown
                plt.plot(xplot2_1, yplot2_1,color=color4_dataAndFit,ls=':',marker='o',label=condi2_rep1_notebookcode)
                plt.plot(xplot2_2, yplot2_2,color=color4_dataAndFit,ls=':',marker='^',label=condi2_rep2_notebookcode)
                plt.plot(xplot2_3, yplot2_3,color=color4_dataAndFit,ls=':',marker='s',label=condi2_rep3_notebookcode)
                plt.plot(xplot2_4, yplot2_4,color=color4_dataAndFit,ls=':',marker='*',label=condi2_rep4_notebookcode)
                pylab.plot(x2,y2,color4_dataAndFit, label=condition2_description+' fit')
                pylab.axvline(condition2_Chalf,color=color4_Chalf,ls='-',label='C1/2')#chalf
                if CI_ratioTOrange2<CI_cutoff:
                    plt.axvspan((condition2_Chalf-condition2_fit_result_output[8]), (condition2_Chalf+condition2_fit_result_output[8]),color=color4_CI)#CI
            
            #plot label
            if denatureMethod=="CD":
                pylab.xlabel('denaturant concentration')
                pylab.ylabel('% labeled')
            if denatureMethod=="HD": 
                pylab.xlabel('normalized temperature')
                pylab.ylabel('% soluble')
            #pylab.ylim(0.3, 2.2)
            pylab.legend(loc='best',fontsize = 'x-small')
            pylab.title(graphName,fontsize=10)
            pylab.savefig(outImage_location+title+"_"+proteinName+".jpg")
            pylab.clf()#clearfig, start a new piece , make sure do it before or after a set of instruction
    
        except :#define the error
            print (title)         
    
        finalData.append(tempRow)      
    #output the data file        
    outfile=open(outfile_locationAndName,"w",newline="")        
    writer=csv.writer(outfile)#create write object
    newheader=["row","protein ID","ProteinNAme","Peptide"]
    newheader_condition1_fitOutPut_CD=["numOf_condition1_1_pts","numOf_condition1_2_pts","numOf_condition1_3_pts","numOf_condition1_4_pts","sumPoints_condition1","B","A","condition1_Chalf","b","B_err","A_err","Chalf_err","b_err","CI_Chalf","CIratioTOrange","CI_Chalf_low","CI_Chalf_up","CI_b","CI_b_low","CI_b_up","r_square","condition1_Chalf_normalized"]
    newheader_condition1_fitOutPut_HD=["numOf_condition1_1_pts","numOf_condition1_2_pts","numOf_condition1_3_pts","numOf_condition1_4_pts","sumPoints_condition1","B","A","condition1_Chalf","b","B_err","A_err","Chalf_err","b_err","CI_Chalf","CIratioTOrange","CI_Chalf_low","CI_Chalf_up","CI_b","CI_b_low","CI_b_up","r_square","condition1_Chalf_temp"]
    
    newheader_condition2_fitOutPut_CD=["numOf_condition2_1_pts","numOf_condition2_2_pts","numOf_condition2_3_pts","numOf_condition2_4_pts","sumPoints_condition2","B","A","condition2_Chalf","b","B_err","A_err","Chalf_err","b_err","CI_Chalf","CIratioTOrange","CI_Chalf_low","CI_Chalf_up","CI_b","CI_b_low","CI_b_up","r_square","condition2_Chalf_normalized","dChalf","dChalf_%","db","db%","CI_Chalf_significant","CI_chalf_overlap"]
    newheader_condition2_fitOutPut_HD=["numOf_condition2_1_pts","numOf_condition2_2_pts","numOf_condition2_3_pts","numOf_condition2_4_pts","sumPoints_condition2","B","A","condition2_Chalf","b","B_err","A_err","Chalf_err","b_err","CI_Chalf","CIratioTOrange","CI_Chalf_low","CI_Chalf_up","CI_b","CI_b_low","CI_b_up","r_square","condition2_Chalf_temp","dChalf","dChalf_%","db","db%","CI_Chalf_significant","CI_chalf_overlap"]
    
    if denatureMethod=="CD":
 
        #condition 1 normalized header
        normalized_header1_1=CD_xaxis_FinalGdmClConcentration_value
        newheader.extend(normalized_header1_1)
        normalized_header1_2=CD_xaxis_FinalGdmClConcentration_value
        newheader.extend(normalized_header1_2)
        normalized_header1_3=CD_xaxis_FinalGdmClConcentration_value
        newheader.extend(normalized_header1_3)
        normalized_header1_4=CD_xaxis_FinalGdmClConcentration_value
        newheader.extend(normalized_header1_4)
        
        #condition 2 normalized header
        normalized_header2_1=CD_xaxis_FinalGdmClConcentration_value
        newheader.extend(normalized_header2_1)
        normalized_header2_2=CD_xaxis_FinalGdmClConcentration_value
        newheader.extend(normalized_header2_2)
        normalized_header2_3=CD_xaxis_FinalGdmClConcentration_value
        newheader.extend(normalized_header2_3)
        normalized_header2_4=CD_xaxis_FinalGdmClConcentration_value
        newheader.extend(normalized_header2_4)
        #fit header
        newheader.extend(newheader_condition1_fitOutPut_CD)
        newheader.extend(newheader_condition2_fitOutPut_CD)
        
    elif denatureMethod=="HD": 
        #condition 1 normalized header        
        normalized_header1_1=HD_xasis_normalizedTemp_value    
        newheader.extend(normalized_header1_1)
        normalized_header1_2=HD_xasis_normalizedTemp_value    
        newheader.extend(normalized_header1_2)
        normalized_header1_3=HD_xasis_normalizedTemp_value    
        newheader.extend(normalized_header1_3)   
        normalized_header1_4=HD_xasis_normalizedTemp_value    
        newheader.extend(normalized_header1_4)          
        #condition 2 normalized header
        normalized_header2_1=HD_xasis_normalizedTemp_value    
        newheader.extend(normalized_header2_1)
        normalized_header2_2=HD_xasis_normalizedTemp_value    
        newheader.extend(normalized_header2_2)
        normalized_header2_3=HD_xasis_normalizedTemp_value    
        newheader.extend(normalized_header2_3)
        normalized_header2_4=HD_xasis_normalizedTemp_value    
        newheader.extend(normalized_header2_4)
        #fit header
        newheader.extend(newheader_condition1_fitOutPut_HD)
        newheader.extend(newheader_condition2_fitOutPut_HD)    

    writer.writerow(newheader)
    for addvariable in finalData:
        writer.writerow(addvariable)
    outfile.close()
    infile.close()

    
