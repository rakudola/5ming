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

# dictionary for condition/replicate information
condition = {
    # run information
    'conditions' : 2, # number of conditions used to generate 1 fit
    'run desc' : 'Heat.Denature High Protein Diet', # what ever desciption can distinguish
    'denature method' : 'HD', # choose 'CD' for chemical denature, 'HD' for heat denature
    
    # condition 1 information
    1: {'description' : 'Ad Libitum', 'plot color' : 1, 'replicates' : 4,
        # condition 1, replicate 1
        1 : {'notebook code' : '1085',                                      # condition 1, replicate 1
             'numPtsLocation' : 4,                                          # input csv column for numPtsLocation
             'Start' : 5,                                                   # input csv column for Start
             'End' : 14,                                                    # input csv column for End
             'marker' : 'o',                                                # marker for plot
             },
        # condition 1, replicate 2
        2 : {'notebook code' : '1091',                                      # condition 1, replicate 2
             'numPtsLocation' : 15,                                         # input csv column for numPtsLocation
             'Start' : 16,                                                  # input csv column for Start
             'End' : 25,                                                    # input csv column for End
             'marker' : '^',                                                # marker for plot
            },
        # condition 1, replicate 3
        3 : {'notebook code' : '1105',                                      # condition 1, replicate 3
             'numPtsLocation' : 26,                                         # input csv column for numPtsLocation
             'Start' : 27,                                                  # input csv column for Start
             'End' : 36,                                                    # input csv column for End
             'marker' : 's',                                                # marker for plot
            },
        # condition 1, replicate 4
        4 : {'notebook code' : '1107',                                      # condition 1, replicate 4
             'numPtsLocation' : 37,                                         # input csv column for numPtsLocation
             'Start' : 38,                                                  # input csv column for Start
             'End' : 47,                                                    # input csv column for End
             'marker' : '*',                                                # marker for plot
            },
        },     

    # condition 2 information
    2: {'description' : 'Dietary restriction', 'plot color' : 2, 'replicates' : 4,   # condition 1
        # condition 2, replicate 1
        1 : {'notebook code' : '1085',                                      # condition 2, replicate 1
             'numPtsLocation' : 48,                                         # input csv column for numPtsLocation
             'Start' : 49,                                                  # input csv column for Start
             'End' : 58,                                                    # input csv column for End
             'marker' : 'o',                                                # marker for plot
             },
        # condition 2, replicate 2
        2 : {'notebook code' : '1091',                                      # condition 2, replicate 2
             'numPtsLocation' : 59,                                         # input csv column for numPtsLocation
             'Start' : 60,                                                  # input csv column for Start
             'End' : 69,                                                    # input csv column for End
             'marker' : '^',                                                # marker for plot
            },
        #condition 2, replicate 3
        3 : {'notebook code' : '1105',                                      # condition 2, replicate 3
             'numPtsLocation' : 70,                                         # input csv column for numPtsLocation
             'Start' : 71,                                                  # input csv column for Start
             'End' : 80,                                                    # input csv column for End
             'marker' : 's',                                                # marker for plot
            },
        # condition 2, replicate 4
        4 : {'notebook code' : '1107',                                      # condition 2, replicate 4
             'numPtsLocation' : 81,                                         # input csv column for numPtsLocation
             'Start' : 82,                                                  # input csv column for Start
             'End' : 91,                                                    # input csv column for End
             'marker' : '*',                                                # marker for plot
             },
        },
}
             

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

# plot color Code, usually don't need to change
# 'data' is for data and fit
colors = {
    1: {'data' : '#1f77b4', 'Chalf' : '#aec7e8', 'CI' : '#c6dbef'},     # blue
    2: {'data' : '#ff7f0e', 'Chalf' : '#ffbb78', 'CI' : '#fdd0a2'},     # orange
    3: {'data' : '#5254a3', 'Chalf' : '#9e9ac8', 'CI' : '#dadaeb'},     # purple
    4: {'data' : '#8c6d31', 'Chalf' : '#e7ba52', 'CI' : '#e7cb94'}}     # brown

def sigmoid(x, B, A, Chalf, b): # the fitting equation
    """Fitting equation. Fits data to sigmoid curve.
    
    Returns y
    """
    y = B + ((A - B) / (1 + np.exp((-1 / b) * (Chalf - x))))
    return y

def read_xdata(con, rep):
    """Reads x data from csv input file"""
    condition[con][rep]['xdata'] = []
    for i in range(condition[con][rep]['Start'], condition[con][rep]['End'] + 1): # data start to data end + 1
        flow = float(header[i]) # making the input value to number
        condition[con][rep]['xdata'].append(flow)
        
def read_ydata(con, rep):
    """Reads y data from csv input file"""
    condition[con][rep]['ydata'] = []
    for i in range(condition[con][rep]['Start'], condition[con][rep]['End'] + 1): # data start to data end + 1
        flow = float(row[i]) # making the input value to number
        condition[con][rep]['ydata'].append(flow)
        
def normalize_data(con, rep):
    """Normalizes data"""
    condition[con][rep]['normalized_ydata'] = []
    condition[con][rep]['normalized_xdata'] = []
    # normalizing condition1 dataset
    if denatureMethod == "CD": # first point is 0, last point is 1      
        # condition1 rep1          
        for element in condition[con][rep]['ydata']: # read from the begining, find first data !=0, make it A
            if element != 0:
                A = element # set 0
                break
        for element in reversed(condition[con][rep]['ydata']): # read from the end, find first data !=0, make it B
            if element != 0:
                B = element # set 1
                break
        for i in range(len(condition[con][rep]['ydata'])):
            if condition[con][rep]['ydata'][i] != 0:
                element = (condition[con][rep]['ydata'][i] - A) / (B - A) # 3 for each data !=0, normalized_ydata=(data-A)/(B-A)
                condition[con][rep]['normalized_ydata'].append(element)
                condition[con][rep]['normalized_xdata'].append(condition[con][rep]['xdata'][i]) # 4 for each normalized_ydata, find the xdata, make it normalized_xdata
            else:        
                condition[con][rep]['normalized_ydata'].append("error")
                condition[con][rep]['normalized_xdata'].append(condition[con][rep]['xdata'][i])  
                
        tempRow.extend(condition[con][rep]['normalized_ydata'])
        
    elif denatureMethod == "HD": # first point is 1, last point is 0 
        # condition1 rep1       
        for element in condition[con][rep]['ydata']: # 1 read from the begining, find first data !=0, make it B
            if element != 0:
                B = element # set 1
                break
        for element in reversed(condition[con][rep]['ydata']): # 2 read from the end, find first data !=0, make it A
            if element != 0:
                A = element # set 0
                break
        for i in range(len(condition[con][rep]['ydata'])):
            if condition[con][rep]['ydata'][i] != 0:
                element = (condition[con][rep]['ydata'][i] - A) / (B - A) #3 for each data !=0, normalized_ydata=(data-A)/(B-A)
                condition[con][rep]['normalized_ydata'].append(element)
                condition[con][rep]['normalized_xdata'].append(condition[con][rep]['xdata'][i]) #4 for each normalized_ydata, find the xdata, make it normalized_xdata
            else:        
                condition[con][rep]['normalized_ydata'].append("error")
                condition[con][rep]['normalized_xdata'].append(condition[con][rep]['xdata'][i])
                
        tempRow.extend(condition[con][rep]['normalized_ydata'])
        
def plot_data(con, rep):
    """Begins compiling data for plotting"""
    condition[con][rep]['xplot'] = []
    condition[con][rep]['yplot'] = []
    for i in range (len(condition[con][rep]['normalized_ydata'])):
        if condition[con][rep]['normalized_ydata'][i] != "error":
            condition[con][rep]['yplot'].append(condition[con][rep]['normalized_ydata'][i])
            condition[con][rep]['xplot'].append(condition[con][rep]['normalized_xdata'][i])
            
# below function is giant wip and has not been tested!!!
def plot_color(col, con):
    for x in range(1, condition[con]['replicates']):
        plt.plot(condition[con][x]['xplot'], condition[con][x]['yplot'], color = colors[col]['data'],
                 ls = ':', marker = condition[con][x]['marker'],
                 label = condition[con][x]['notebook code'])
        # plt.plot(xplot1_2, yplot1_2, color = colors[col]['data'], ls = ':', marker = '^', label = condi1_rep2_notebookcode)
        # plt.plot(xplot1_3, yplot1_3,color=colors[col]['data'],ls=':',marker='s',label=condi1_rep3_notebookcode)
        # plt.plot(xplot1_4, yplot1_4,color=colors[col]['data'],ls=':',marker='*',label=condi1_rep4_notebookcode)
        pylab.plot(condition[1]['x'], condition[1]['y'], colors[col]['data'], label = condition[con]['description'] + ' fit')
        pylab.axvline(condition[1]['C half'], color = colors[col]['Chalf'], ls = '-', label = 'C1/2')#chalf
        if CI_ratioTOrange1 < CI_cutoff:
            plt.axvspan((condition[1]['C half']-condition[1]['fit output'][8]),
                        (condition[1]['C half']+condition[1]['fit output'][8]),
                        color = colors[1]['CI']) # CI
 


infile = open(infile_locationAndName, "r") # "r"=reading mode
reader = csv.reader(infile) # create reader object
header = next(reader)


print(infile_locationAndName)


#base on replication numbers

if numOfCondition=="2":

    finalData = []
    
    for i in range(1, condition['conditions'] + 1):         # for each condition (i)
        for j in range(1, condition[i]['replicates'] + 1):  # for each replicate (j) of condition (i)
            read_xdata(i, j)                                # read xdata from input file for condition i, replicate j

    for row in reader:
        try:
            for i in range(1, condition['conditions'] + 1):         # for each condition (i)
                for j in range(1, condition[i]['replicates'] + 1):  # for each replicate (j) of condition (i)
                    read_ydata(i, j)                                # read xdata from input file for condition i, replicate j
        except:
            pass
                  
        #basic information input 
        title = row[rowLoc]
        ID = row[IDLoc]
        proteinName = row[NameLoc]
        Peptide = row[PeptideLoc]
        
        # numOfReporter = row[ReporterLoc]
        condition[1][1]['num pts'] = row[condition[1][1]['numPtsLocation']]
        condition[1][2]['num pts'] = row[condition[1][2]['numPtsLocation']]
        condition[1][3]['num pts'] = row[condition[1][3]['numPtsLocation']]
        condition[1][4]['num pts'] = row[condition[1][4]['numPtsLocation']]
        
        condition[2][1]['num pts'] = row[condition[2][1]['numPtsLocation']]
        condition[2][2]['num pts'] = row[condition[2][2]['numPtsLocation']]
        condition[2][3]['num pts'] = row[condition[2][3]['numPtsLocation']]
        condition[2][4]['num pts'] = row[condition[2][4]['numPtsLocation']]
        
        # baisc information ouput 
        tempRow = [title] # [0]
        tempRow.append(ID) # [1]
        tempRow.append(proteinName) # [2]
        tempRow.append(Peptide) # [3]
        # tempRow.append(numOfReporter)#[4]
        
    
        try: # if statment, when it is not error   
            # normalized the data

            for i in range(1, condition['conditions'] + 1):
                for j in range(1, condition[i]['replicates'] + 1):
                    normalize_data(i, j)
            
            #Set data set for fit and plot
            for i in range(1, condition['conditions'] + 1):
                for j in range(1, condition[i]['replicates'] + 1):
                    plot_data(i, j)
                   
            # fit condition1 to s curve
            condition[1]['fit output'] = []
            # giant string of  xplots and whatever below just concatenated; write funtion for?
            condition[1]['popt'], condition[1]['pcov'] = curve_fit(sigmoid, condition[1][1]['xplot'] + condition[1][2]['xplot'] + condition[1][3]['xplot'] + condition[1][4]['xplot'], condition[1][1]['yplot'] + condition[1][2]['yplot'] + condition[1][3]['yplot'] + condition[1][4]['yplot'])
            condition[1]['fit output'].extend(condition[1]['popt'])#extend merge two list to one big list  
                                                #popt has all the parameter (variable) in regression equation (x,B, A, Chalf, b)
                                                                                            #row [-, 0, 1, 2,    3] in  condition[1]['fit output']   
            if denatureMethod == "CD":
                condition[1]['x'] = np.linspace(0, 4, 50)
            elif denatureMethod == "HD": 
                condition[1]['x'] = np.linspace(0, 1.2, 50)
            condition[1]['y'] = sigmoid(condition[1]['x'], *condition[1]['popt'])# *popt split the two variable 
            condition[1]['C half'] = condition[1]['fit output'][2]
            condition[1]['b'] = condition[1]['fit output'][3]
            
         
            # calculate condition1 confidence interval
            condition[1][1]['float num pts'] = float(condition[1][1]['num pts'])
            condition[1][2]['float num pts'] = float(condition[1][2]['num pts'])
            condition[1][3]['float num pts'] = float(condition[1][3]['num pts'])
            condition[1][4]['float num pts'] = float(condition[1][4]['num pts'])
            condition[1]['sum pts'] = condition[1][1]['float num pts'] + condition[1][2]['float num pts'] + condition[1][3]['float num pts'] + condition[1][4]['float num pts']
            tempRow.append(condition[1][1]['num pts'])  # [5]
            tempRow.append(condition[1][2]['num pts'])      # [5]
            tempRow.append(condition[1][3]['num pts'])      # [5]
            tempRow.append(condition[1][4]['num pts'])      # [5]
            tempRow.append(condition[1]['sum pts'])     # [6]

            condition[1]['std error'] = np.sqrt(np.diag(condition[1]['pcov']))
            condition[1]['fit output'].extend(condition[1]['std error'])#[4,5,6,7]=[Berror, Aerror, Chalferror, berror]
            condition1_ConfidenceIntervcondition1_chalf = t.ppf(.975, (condition[1]['sum pts']-1)) * condition[1]['fit output'][6] / np.sqrt(condition[1]['sum pts'])
            condition[1]['fit output'].append(condition1_ConfidenceIntervcondition1_chalf)#[8]
            if denatureMethod=="CD":
                CI_ratioTOrange1=condition1_ConfidenceIntervcondition1_chalf/CD_conc_range
                condition[1]['fit output'].append(CI_ratioTOrange1)#[9]
            elif denatureMethod=="HD":
                CI_ratioTOrange1 = condition1_ConfidenceIntervcondition1_chalf/1
                condition[1]['fit output'].append(CI_ratioTOrange1)#[9]
            condition1_CI_Chalf_lowbound = condition[1]['C half'] - condition[1]['fit output'][8]
            condition1_CI_Chalf_upbound = condition[1]['C half'] + condition[1]['fit output'][8]
            condition[1]['fit output'].append(condition1_CI_Chalf_lowbound)#[10]
            condition[1]['fit output'].append(condition1_CI_Chalf_upbound)#[11]
            condition1_ConfidenceIntervcondition1_b=t.ppf(.975, (condition[1]['sum pts']-1))* condition[1]['fit output'][7]/np.sqrt(condition[1]['sum pts'])
            condition[1]['fit output'].append(condition1_ConfidenceIntervcondition1_b)#[12]
            condition1_CI_b_lowbound=condition[1]['b']-condition[1]['fit output'][12]
            condition1_CI_b_upbound=condition[1]['b']+condition[1]['fit output'][12]
            condition[1]['fit output'].append(condition1_CI_b_lowbound)#[13]
            condition[1]['fit output'].append(condition1_CI_b_upbound)#[14]
            
            #compute condition1 r squared
            residuals=(condition[1][1]['yplot']+condition[1][2]['yplot']+condition[1][3]['yplot']+condition[1][4]['yplot'])-sigmoid((condition[1][1]['xplot']+condition[1][2]['xplot']+condition[1][3]['xplot']+condition[1][4]['xplot']),*condition[1]['popt'])
            ss_res=np.sum(residuals**2)
            ss_tot=np.sum(((condition[1][1]['yplot']+condition[1][2]['yplot']+condition[1][3]['yplot']+condition[1][4]['yplot'])-np.mean(condition[1][1]['yplot']+condition[1][2]['yplot']+condition[1][3]['yplot']+condition[1][4]['yplot']))**2)
            r_squared1=1-(ss_res/ss_tot)
            condition[1]['fit output'].append(r_squared1)#[15]; append add something to the list   
            if denatureMethod=="HD":
                Chalf1_temp=(condition[1]['C half']*HD_temp_range)+HD_Temp_lowbound
                condition[1]['fit output'].append(Chalf1_temp)#[16]
            if denatureMethod=="CD":
                Chalf1_normalized=condition[1]['C half']/CD_conc_range
                condition[1]['fit output'].append(Chalf1_normalized)#[16]
                
            
            
            tempRow.extend(condition[1]['fit output'])

            
             #fit condition2 to s curve
            condition2_fit_result_output=[]
            popt2, pcov2 = curve_fit(sigmoid, condition[2][1]['xplot']+condition[2][2]['xplot']+condition[2][3]['xplot']+condition[2][4]['xplot'], condition[2][1]['yplot']+condition[2][2]['yplot']+condition[2][3]['yplot']+condition[2][4]['yplot'])
            condition2_fit_result_output.extend(popt2)#extend merge two list to one big list  
                                                #popt has all the parameter (variable) in regression equation (x,B, A, Chalf, b)
                                                                                           #row [-, 0, 1, 2,    3] in  condition[1]['fit output']   
           
            if denatureMethod=="CD":
                x2 = np.linspace(0, 4, 50)
            if denatureMethod=="HD": 
                x2 = np.linspace(0, 1.2, 50)
            y2= sigmoid(x2, *popt2)# *popt split the two variable 
            condition2_Chalf=condition2_fit_result_output[2]
            condition2_b=condition2_fit_result_output[3]
            
            
            #calculate condition2 confidence interval
            numOfPoint_condition2_repA=float(condition[2][1]['num pts'])
            numOfPoint_condition2_repB=float(condition[2][2]['num pts'])
            numOfPoint_condition2_repC=float(condition[2][3]['num pts'])
            numOfPoint_condition2_repD=float(condition[2][4]['num pts'])
            sumPoints_condition2=numOfPoint_condition2_repA+numOfPoint_condition2_repB+numOfPoint_condition2_repC+numOfPoint_condition2_repD
            tempRow.append(condition[2][1]['num pts'])#[6]
            tempRow.append(condition[2][2]['num pts'])#[6]
            tempRow.append(condition[2][3]['num pts'])#[6]
            tempRow.append(condition[2][4]['num pts'])#[6]
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
            residuals2=(condition[2][1]['yplot']+condition[2][2]['yplot']+condition[2][3]['yplot']+condition[2][4]['yplot'])-sigmoid((condition[2][1]['xplot']+condition[2][2]['xplot']+condition[2][3]['xplot']+condition[2][4]['xplot']),*popt2)
            ss_res2=np.sum(residuals2**2)
            ss_tot2=np.sum(((condition[2][1]['yplot']+condition[2][2]['yplot']+condition[2][3]['yplot']+condition[2][4]['yplot'])-np.mean(condition[2][1]['yplot']+condition[2][2]['yplot']+condition[2][3]['yplot']+condition[2][4]['yplot']))**2)
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
                delta_Chalf=condition2_Chalf-condition[1]['C half']
                delta_Chalf_percent=delta_Chalf/condition[1]['C half']
                delta_b=condition2_b-condition[1]['b']
                delta_b_percent=delta_b/condition[1]['b']
            if denatureMethod=="HD":
                delta_Chalf=Chalf2_temp-Chalf1_temp
                delta_Chalf_percent=delta_Chalf/Chalf1_temp
                delta_b=condition2_b-condition[1]['b']
                delta_b_percent=delta_b/condition[1]['b']

            tempRow.append(delta_Chalf)
            tempRow.append(delta_Chalf_percent)
            tempRow.append(delta_b)
            tempRow.append(delta_b_percent)

            
            
            
            # generate the figure
            if CI_ratioTOrange1 < CI_cutoff and CI_ratioTOrange2 < CI_cutoff:
                graphName=descriptionOfTheRun + "_" + "row#" + title+ " (has CI)" + "\n" + ID + "||" + proteinName + "\n" + Peptide
            else:
                graphName=descriptionOfTheRun + "_" + "row#" + title + "\n" + ID + "||" + proteinName + "\n" + Peptide
                
            # condition1_plot color
            if condition[1]['plot color'] == 1: # choose from blue, orange, purple, brown
                plt.plot(condition[1][1]['xplot'], condition[1][1]['yplot'], color = colors[1]['data'], ls = ':', marker = 'o', label = condition[1][1]['notebook code'])
                plt.plot(condition[1][2]['xplot'], condition[1][2]['yplot'], color = colors[1]['data'], ls = ':', marker = '^', label = condition[1][2]['notebook code'])
                plt.plot(condition[1][3]['xplot'], condition[1][3]['yplot'],color=colors[1]['data'],ls=':',marker='s',label=condition[1][3]['notebook code'])
                plt.plot(condition[1][4]['xplot'], condition[1][4]['yplot'],color=colors[1]['data'],ls=':',marker='*',label=condition[1][4]['notebook code'])
                pylab.plot(condition[1]['x'],condition[1]['y'],colors[1]['data'], label=condition[1]['description']+' fit')
                pylab.axvline(condition[1]['C half'],color=colors[1]['Chalf'],ls='-',label='C1/2')#chalf
                if CI_ratioTOrange1<CI_cutoff:
                    plt.axvspan((condition[1]['C half']-condition[1]['fit output'][8]), (condition[1]['C half']+condition[1]['fit output'][8]),color=colors[1]['CI'])#CI
                    
            elif condition[1]['plot color'] == 2: # choose from blue, orange, purple, brown
                plt.plot(condition[1][1]['xplot'], condition[1][1]['yplot'],color=colors[2]['data'],ls=':',marker='o',label=condition[1][1]['notebook code'])
                plt.plot(condition[1][2]['xplot'], condition[1][2]['yplot'],color=colors[2]['data'],ls=':',marker='^',label = condition[1][2]['notebook code'])
                plt.plot(condition[1][3]['xplot'], condition[1][3]['yplot'],color=colors[2]['data'],ls=':',marker='s',label=condition[1][3]['notebook code'])
                plt.plot(condition[1][4]['xplot'], condition[1][4]['yplot'],color=colors[2]['data'],ls=':',marker='*',label=condition[1][4]['notebook code'])
                pylab.plot(condition[1]['x'],condition[1]['y'],colors[2]['data'], label=condition[1]['description']+' fit')
                pylab.axvline(condition[1]['C half'],color=colors[2]['Chalf'],ls='-',label='C1/2')#chalf
                if CI_ratioTOrange1<CI_cutoff:
                    plt.axvspan((condition[1]['C half']-condition[1]['fit output'][8]), (condition[1]['C half']+condition[1]['fit output'][8]),color=colors[2]['CI'])#CI
                    
            elif condition[1]['plot color'] == 3: # choose from blue, orange, purple, brown
                plt.plot(condition[1][1]['xplot'], condition[1][1]['yplot'],color=colors[3]['data'],ls=':',marker='o',label=condition[1][1]['notebook code'])
                plt.plot(condition[1][2]['xplot'], condition[1][2]['yplot'],color=colors[3]['data'],ls=':',marker='^',label = condition[1][2]['notebook code'])
                plt.plot(condition[1][3]['xplot'], condition[1][3]['yplot'],color=colors[3]['data'],ls=':',marker='s',label=condition[1][3]['notebook code'])
                plt.plot(condition[1][4]['xplot'], condition[1][4]['yplot'],color=colors[3]['data'],ls=':',marker='*',label=condition[1][4]['notebook code'])
                pylab.plot(condition[1]['x'],condition[1]['y'],colors[3]['data'], label=condition[1]['description']+' fit')
                pylab.axvline(condition[1]['C half'],color=colors[3]['Chalf'],ls='-',label='C1/2')#chalf
                if CI_ratioTOrange1<CI_cutoff:
                    plt.axvspan((condition[1]['C half']-condition[1]['fit output'][8]), (condition[1]['C half']+condition[1]['fit output'][8]),color=colors[3]['CI'])#CI
                    
            elif condition[1]['plot color'] == 4:#choose from blue, orange, purple, brown
                plt.plot(condition[1][1]['xplot'], condition[1][1]['yplot'],color=colors[4]['data'],ls=':',marker='o',label=condition[1][1]['notebook code'])
                plt.plot(condition[1][2]['xplot'], condition[1][2]['yplot'],color=colors[4]['data'],ls=':',marker='^',label = condition[1][2]['notebook code'])
                plt.plot(condition[1][3]['xplot'], condition[1][3]['yplot'],color=colors[4]['data'],ls=':',marker='s',label=condition[1][3]['notebook code'])
                plt.plot(condition[1][4]['xplot'], condition[1][4]['yplot'],color=colors[4]['data'],ls=':',marker='*',label=condition[1][4]['notebook code'])
                pylab.plot(condition[1]['x'],condition[1]['y'],colors[4]['data'], label=condition[1]['description']+' fit')
                pylab.axvline(condition[1]['C half'],color=colors[4]['Chalf'],ls='-',label='C1/2')#chalf
                if CI_ratioTOrange1<CI_cutoff:
                    plt.axvspan((condition[1]['C half']-condition[1]['fit output'][8]), (condition[1]['C half']+condition[1]['fit output'][8]),color=colors[4]['CI'])#CI

            #condition2_plot color
            if condition2_plotColor=="1":#choose from blue, orange, purple, brown
                plt.plot(condition[2][1]['xplot'], condition[2][1]['yplot'],color=colors[1]['data'],ls=':',marker='o',label=condition[2][1]['notebook code'])
                plt.plot(condition[2][2]['xplot'], condition[2][2]['yplot'],color=colors[1]['data'],ls=':',marker='^',label=condition[2][2]['notebook code'])
                plt.plot(condition[2][3]['xplot'], condition[2][3]['yplot'],color=colors[1]['data'],ls=':',marker='s',label=condition[2][3]['notebook code'])
                plt.plot(condition[2][4]['xplot'], condition[2][4]['yplot'],color=colors[1]['data'],ls=':',marker='*',label=condition[2][4]['notebook code'])
                pylab.plot(x2,y2,colors[1]['data'], label=condition[2]['description']+' fit')
                pylab.axvline(condition2_Chalf,color=colors[1]['Chalf'],ls='-',label='C1/2')#chalf
                if CI_ratioTOrange2<CI_cutoff:
                    plt.axvspan((condition2_Chalf-condition2_fit_result_output[8]), (condition2_Chalf+condition2_fit_result_output[8]),color=colors[1]['CI'])#CI
            elif condition2_plotColor=="2":#choose from blue, orange, purple, brown
                plt.plot(condition[2][1]['xplot'], condition[2][1]['yplot'],color=colors[2]['data'],ls=':',marker='o',label=condition[2][1]['notebook code'])
                plt.plot(condition[2][2]['xplot'], condition[2][2]['yplot'],color=colors[2]['data'],ls=':',marker='^',label=condition[2][2]['notebook code'])
                plt.plot(condition[2][3]['xplot'], condition[2][3]['yplot'],color=colors[2]['data'],ls=':',marker='s',label=condition[2][3]['notebook code'])
                plt.plot(condition[2][4]['xplot'], condition[2][4]['yplot'],color=colors[2]['data'],ls=':',marker='*',label=condition[2][4]['notebook code'])
                pylab.plot(x2,y2,colors[2]['data'], label=condition[2]['description']+' fit')
                pylab.axvline(condition2_Chalf,color=colors[2]['Chalf'],ls='-',label='C1/2')#chalf
                if CI_ratioTOrange2<CI_cutoff:
                    plt.axvspan((condition2_Chalf-condition2_fit_result_output[8]), (condition2_Chalf+condition2_fit_result_output[8]),color=colors[2]['CI'])#CI
            elif condition2_plotColor=="3":#choose from blue, orange, purple, brown
                plt.plot(condition[2][1]['xplot'], condition[2][1]['yplot'],color=colors[3]['data'],ls=':',marker='o',label=condition[2][1]['notebook code'])
                plt.plot(condition[2][2]['xplot'], condition[2][2]['yplot'],color=colors[3]['data'],ls=':',marker='^',label=condition[2][2]['notebook code'])
                plt.plot(condition[2][3]['xplot'], condition[2][3]['yplot'],color=colors[3]['data'],ls=':',marker='s',label=condition[2][3]['notebook code'])
                plt.plot(condition[2][4]['xplot'], condition[2][4]['yplot'],color=colors[3]['data'],ls=':',marker='*',label=condition[2][4]['notebook code'])
                pylab.plot(x2,y2,colors[3]['data'], label=condition[2]['description']+' fit')
                pylab.axvline(condition2_Chalf,color=colors[3]['Chalf'],ls='-',label='C1/2')#chalf
                if CI_ratioTOrange2<CI_cutoff:
                    plt.axvspan((condition2_Chalf-condition2_fit_result_output[8]), (condition2_Chalf+condition2_fit_result_output[8]),color=colors[3]['CI'])#CI
            elif condition2_plotColor=="4":#choose from blue, orange, purple, brown
                plt.plot(condition[2][1]['xplot'], condition[2][1]['yplot'],color=colors[4]['data'],ls=':',marker='o',label=condition[2][1]['notebook code'])
                plt.plot(condition[2][2]['xplot'], condition[2][2]['yplot'],color=colors[4]['data'],ls=':',marker='^',label=condition[2][2]['notebook code'])
                plt.plot(condition[2][3]['xplot'], condition[2][3]['yplot'],color=colors[4]['data'],ls=':',marker='s',label=condition[2][3]['notebook code'])
                plt.plot(condition[2][4]['xplot'], condition[2][4]['yplot'],color=colors[4]['data'],ls=':',marker='*',label=condition[2][4]['notebook code'])
                pylab.plot(x2,y2,colors[4]['data'], label=condition[2]['description']+' fit')
                pylab.axvline(condition2_Chalf,color=colors[4]['Chalf'],ls='-',label='C1/2')#chalf
                if CI_ratioTOrange2<CI_cutoff:
                    plt.axvspan((condition2_Chalf-condition2_fit_result_output[8]), (condition2_Chalf+condition2_fit_result_output[8]),color=colors[4]['CI'])#CI
            
            #plot label
            if denatureMethod=="CD":
                pylab.xlabel('denaturant concentration')
                pylab.ylabel('% labeled')
            elif denatureMethod=="HD": 
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