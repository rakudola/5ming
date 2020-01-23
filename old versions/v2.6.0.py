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
descriptionOfTheRun = "Heat.Denature High Protein Diet" # whatever desciption can distinguish

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

# plot color Code, usually don't need to change
# 'data' is for data and fit
colors = {
    1: {'data' : '#1f77b4', 'C half' : '#aec7e8', 'CI' : '#c6dbef'},     # blue
    2: {'data' : '#ff7f0e', 'C half' : '#ffbb78', 'CI' : '#fdd0a2'},     # orange
    3: {'data' : '#5254a3', 'C half' : '#9e9ac8', 'CI' : '#dadaeb'},     # purple
    4: {'data' : '#8c6d31', 'C half' : '#e7ba52', 'CI' : '#e7cb94'}}     # brown


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
    """Normalizes x and y data"""
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


def fit_scurve(con):
    """fit condition to s-curve"""
    # fit condition to s curve
    condition[con]['fit output'] = []
    # giant string of  xplots and whatever below just concatenated; write funtion for?
    condition[con]['popt'], condition[con]['pcov'] = curve_fit(sigmoid, condition[con][1]['xplot'] + condition[con][2]['xplot'] + condition[con][3]['xplot'] + condition[con][4]['xplot'], condition[con][1]['yplot'] + condition[con][2]['yplot'] + condition[con][3]['yplot'] + condition[con][4]['yplot'])
    condition[con]['fit output'].extend(condition[con]['popt'])#extend merge two list to one big list  
                                        #popt has all the parameter (variable) in regression equation (x,B, A, Chalf, b)
                                                                                    #row [-, 0, 1, 2,    3] in  condition[1]['fit output']   
    if denatureMethod == "CD":
        condition[con]['x'] = np.linspace(0, 4, 50)
    elif denatureMethod == "HD": 
        condition[con]['x'] = np.linspace(0, 1.2, 50)
    condition[con]['y'] = sigmoid(condition[con]['x'], *condition[con]['popt']) # *popt split the two variable 
    condition[con]['C half'] = condition[con]['fit output'][2]
    condition[con]['b'] = condition[con]['fit output'][3]


def confidence_interval(con):
    """calculate confidence interval for c 1/2 and b"""
    condition[con]['sum pts'] = []
    for j in range(1, condition[con]['replicates'] + 1):
        condition[con][j]['float num pts'] = float(condition[con][j]['num pts'])

    # BELOW LINE NEEDS TO BE FIXED
    condition[con]['sum pts'] = condition[con][1]['float num pts'] + condition[con][2]['float num pts'] + condition[con][3]['float num pts'] + condition[con][4]['float num pts']
    """for j in range(1, condition[con]['replicates'] + 1):
        condition[con]['sum pts'] = condition[con]['sum pts'] + condition[con][j]['float num pts']"""

    for j in range(1, condition[con]['replicates'] + 1):
        tempRow.append(condition[con][j]['num pts']) # [5]

    tempRow.append(condition[con]['sum pts']) # [6]

    condition[con]['std error'] = np.sqrt(np.diag(condition[con]['pcov']))
    condition[con]['fit output'].extend(condition[con]['std error'])#[4,5,6,7]=[Berror, Aerror, Chalferror, berror]
    condition[con]['CI_Chalf'] = t.ppf(.975, (condition[con]['sum pts'] - 1)) * condition[con]['fit output'][6] / np.sqrt(condition[con]['sum pts'])
    condition[con]['fit output'].append(condition[con]['CI_Chalf'])#[8]

    if denatureMethod == "CD":
        condition[con]['ratioTOrange'] = condition[con]['CI_Chalf'] / CD_conc_range
        condition[con]['fit output'].append(condition[con]['ratioTOrange'])#[9]
    elif denatureMethod == "HD":
        condition[con]['ratioTOrange'] = condition[con]['CI_Chalf']/1
        condition[con]['fit output'].append(condition[con]['ratioTOrange'])#[9]

    condition[con]['CI lowbound'] = condition[con]['C half'] - condition[con]['fit output'][8]
    condition[con]['CI upbound'] = condition[con]['C half'] + condition[con]['fit output'][8]
    condition[con]['fit output'].append(condition[con]['CI lowbound'])#[10]
    condition[con]['fit output'].append(condition[con]['CI upbound'])#[11]
    condition[con]['CI_b'] = t.ppf(.975, (condition[con]['sum pts'] - 1)) * condition[con]['fit output'][7] / np.sqrt(condition[con]['sum pts'])
    condition[con]['fit output'].append(condition[con]['CI_b'])#[12]
    condition[con]['CI_b lowbound'] = condition[con]['b'] - condition[con]['fit output'][12]
    condition[con]['CI_b upbound'] = condition[con]['b'] + condition[con]['fit output'][12]
    condition[con]['fit output'].append(condition[con]['CI_b lowbound'])#[13]
    condition[con]['fit output'].append(condition[con]['CI_b upbound'])#[14]


def r_squared(con):
    """calculate r squared"""
    condition[con]['residuals'] = (condition[con][1]['yplot']+condition[con][2]['yplot']+condition[con][3]['yplot']+condition[con][4]['yplot'])-sigmoid((condition[con][1]['xplot']+condition[con][2]['xplot']+condition[con][3]['xplot']+condition[con][4]['xplot']),*condition[con]['popt'])
    condition[con]['ss_res'] = np.sum(condition[con]['residuals']**2)
    condition[con]['ss_tot'] = np.sum(((condition[con][1]['yplot']+condition[con][2]['yplot']+condition[con][3]['yplot']+condition[con][4]['yplot'])-np.mean(condition[con][1]['yplot']+condition[con][2]['yplot']+condition[con][3]['yplot']+condition[con][4]['yplot']))**2)
    condition[con]['r squared'] = 1 - (condition[con]['ss_res'] / condition[con]['ss_tot'])
    condition[con]['fit output'].append(condition[con]['r squared'])#[15]; append add something to the list   
    if denatureMethod == "HD":
        condition[con]['Chalf temp'] = (condition[con]['C half'] * HD_temp_range) + HD_Temp_lowbound
        condition[con]['fit output'].append(condition[con]['Chalf temp'])#[16]
    elif denatureMethod == "CD":
        condition[con]['Chalf normalized'] = condition[con]['C half'] / CD_conc_range
        condition[con]['fit output'].append(condition[con]['Chalf normalized'])#[16]


def plot_color(col, con):
    """plot color"""
    for x in range(1, condition[con]['replicates'] + 1):
        plt.plot(condition[con][x]['xplot'], condition[con][x]['yplot'], color = colors[col]['data'],
                 ls = ':', marker = condition[con][x]['marker'],
                 label = condition[con][x]['notebook code'])

    pylab.plot(condition[con]['x'], condition[con]['y'], colors[col]['data'], label = condition[con]['description'] + ' fit')
    pylab.axvline(condition[con]['C half'], color = colors[col]['C half'], ls = '-', label = 'C1/2') # chalf
    if condition[con]['ratioTOrange'] < CI_cutoff:
        plt.axvspan((condition[con]['C half'] - condition[con]['fit output'][8]),
                    (condition[con]['C half'] + condition[con]['fit output'][8]),
                    color = colors[col]['CI']) # CI


infile = open(infile_locationAndName, "r") # "r" = reading mode
reader = csv.reader(infile) # create reader object
header = next(reader)

print(infile_locationAndName)


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
    """for i in range(1, condition['conditions'] + 1):
        for j in range(1, condition[i]['replicates'] + 1):
            condition[i][j]['num pts'] = row[condition[1][1]['numPtsLocation']]"""
    condition[1][1]['num pts'] = row[condition[1][1]['numPtsLocation']]
    condition[1][2]['num pts'] = row[condition[1][2]['numPtsLocation']]
    condition[1][3]['num pts'] = row[condition[1][3]['numPtsLocation']]
    condition[1][4]['num pts'] = row[condition[1][4]['numPtsLocation']]
    
    condition[2][1]['num pts'] = row[condition[2][1]['numPtsLocation']]
    condition[2][2]['num pts'] = row[condition[2][2]['numPtsLocation']]
    condition[2][3]['num pts'] = row[condition[2][3]['numPtsLocation']]
    condition[2][4]['num pts'] = row[condition[2][4]['numPtsLocation']]
    
    # basic information ouput 
    tempRow = [title] # [0]
    tempRow.append(ID) # [1]
    tempRow.append(proteinName) # [2]
    tempRow.append(Peptide) # [3]
    

    try: # if statment, when it is not error   
        # normalize the data
        for i in range(1, condition['conditions'] + 1):
            for j in range(1, condition[i]['replicates'] + 1):
                normalize_data(i, j)
        
        #Set data set for fit and plot
        for i in range(1, condition['conditions'] + 1):
            for j in range(1, condition[i]['replicates'] + 1):
                plot_data(i, j)

        
        for i in range(1, condition['conditions'] + 1):
            fit_scurve(i)
            confidence_interval(i)
            r_squared(i)
            tempRow.extend(condition[i]['fit output'])


        # delta calculation; need to fix + adjust for diff number conditions
        # 'D' is shorthand for delta
        if denatureMethod == "CD":
            condition['D_C half'] = condition[2]['C half'] - condition[1]['C half']
            condition['D_C half percent'] = condition['D_C half'] / condition[1]['C half']
            condition['D_b'] = condition[2]['b'] - condition[1]['b']
            condition['D_b percent'] = condition['D_b'] / condition[1]['b']
        elif denatureMethod == "HD":
            condition['D_C half'] = condition[2]['Chalf temp'] - condition[1]['Chalf temp']
            condition['D_C half percent'] = condition['D_C half'] / condition[1]['Chalf temp']
            condition['D_b'] = condition[2]['b'] - condition[1]['b']
            condition['D_b percent'] = condition['D_b'] / condition[1]['b']
    
        tempRow.append(condition['D_C half'])
        tempRow.append(condition['D_C half percent'])
        tempRow.append(condition['D_b'])
        tempRow.append(condition['D_b percent'])


        # generate the figure
        if condition[1]['ratioTOrange'] < CI_cutoff and condition[2]['ratioTOrange'] < CI_cutoff:
            graphName=descriptionOfTheRun + "_" + "row#" + title + " (has CI)" + "\n" + ID + "||" + proteinName + "\n" + Peptide
        else:
            graphName=descriptionOfTheRun + "_" + "row#" + title + "\n" + ID + "||" + proteinName + "\n" + Peptide

        for i in range(1, condition['conditions'] + 1):
            plot_color(condition[i]['plot color'], i)
        
        # plot label
        if denatureMethod == "CD":
            pylab.xlabel('denaturant concentration')
            pylab.ylabel('% labeled')
        elif denatureMethod == "HD": 
            pylab.xlabel('normalized temperature')
            pylab.ylabel('% soluble')
        # pylab.ylim(0.3, 2.2)
        pylab.legend(loc='best',fontsize = 'x-small')
        pylab.title(graphName,fontsize=10)
        pylab.savefig(outImage_location+title+"_"+proteinName+".jpg")
        pylab.clf() # clearfig, start a new piece , make sure do it before or after a set of instruction

    except: # define the error
        print (title)

    finalData.append(tempRow)
      
# output the data file        
outfile = open(outfile_locationAndName, "w", newline = "")        
writer = csv.writer(outfile) # create write object
newheader = ["row", "protein ID", "ProteinNAme", "Peptide"]

for i in range(1, condition['conditions'] + 1):
    condition[i]['newheader fit output'] = []
    for j in range(1, condition[i]['replicates'] + 1):
        temp_string = "numOf_condition" + str(i) + "_" + str(j) + "_pts"
        condition[i]['newheader fit output'] += [temp_string]

    temp_string = "sumPoints_condition" + str(i)
    condition[i]['newheader fit output'] += [temp_string]

    condition[i]['newheader fit output'] += ["B", "A"]

    temp_string = "condition" + str(i) + "_Chalf"
    condition[i]['newheader fit output'] += [temp_string]

    condition[i]['newheader fit output'] += ["b", "B_err", "A_err", "Chalf_err", "b_err", "CI_Chalf", "CIratioTOrange", "CI_Chalf_low", "CI_Chalf_up", "CI_b", "CI_b_low", "CI_b_up", "r_square"]

    if denatureMethod == "CD":
        temp_string = "condition" + str(i) + "_Chalf_normalized"
        condition[i]['newheader fit output'] += [temp_string]
    elif denatureMethod == "HD":
        temp_string = "condition" + str(i) + "_Chalf_temp"
        condition[i]['newheader fit output'] += [temp_string]

    if i == (condition['conditions']):
        condition[i]['newheader fit output'] += ["dChalf","dChalf_%","db","db%","CI_Chalf_significant","CI_chalf_overlap"]

if denatureMethod == "CD":
 
    # normalized header
    for i in range(1, condition['conditions'] + 1):
        for j in range(1, condition[i]['replicates'] + 1):
            newheader.extend(CD_xaxis_FinalGdmClConcentration_value)
    
elif denatureMethod == "HD": 
    
    # normalized header
    for i in range(1, condition['conditions'] + 1):
        for j in range(1, condition[i]['replicates'] + 1):
            newheader.extend(HD_xasis_normalizedTemp_value)  

# fit header
for i in range(1, condition['conditions'] + 1):
    newheader.extend(condition[i]['newheader fit output'])

writer.writerow(newheader)
for addvariable in finalData:
    writer.writerow(addvariable)
outfile.close()
infile.close()
