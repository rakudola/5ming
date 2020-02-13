import numpy as np
import pylab
import csv
from scipy.optimize import curve_fit
from scipy.stats import t
import matplotlib.pyplot as plt

# Setting PRINT_ALL to "True" will print some graphs with bad data
# Setting PRINT_ALL to "False" will not print graphs with bad data
PRINT_ALL = False

# dictionary for condition/replicate information
condition = {
    # run information
    'run desc' : 'Heat.Denature High Protein Diet', # what ever desciption can distinguish
    'denature method' : 'HD', # choose 'CD' for chemical denature, 'HD' for heat denature
}
             
# file selection
infile_locationAndName = "../test_multiInput.csv"
outfile_locationAndName = "../test_multiOutput_2.csv"
outImage_location = "../testoutput2_2//"

"""Code below here does not need to be modified for a standard run"""

# basic information --set up variable
# method infodisgunsh
CD_xaxis_FinalGdmClConcentration_value = [0, 0.43, 0.87, 1.3, 1.74, 2.17, 2.61, 3.04, 3.48]
CD_conc_upbound = 3.48
CD_conc_lowbound = 0
CD_conc_range = CD_conc_upbound - CD_conc_lowbound
HD_xasis_normalizedTemp_value = [0, 0.129277567, 0.224334601, 0.338403042, 0.429657795, 0.539923954, 0.646387833, 0.775665399, 0.882129278, 1]
HD_Temp_upbound = 63
HD_Temp_lowbound = 36.7
HD_temp_range = HD_Temp_upbound - HD_Temp_lowbound
CI_cutoff = 0.3

# info for each peptide from infile
rowLoc = 0
IDLoc = 1
NameLoc = 2
PeptideLoc = 3

# plot color code, usually don't need to change
colors = {
    # 'data' is for data and fit
    # numbers are scrambled due to issue with the way the color variable is calculated. Pain
    2: {'data' : '#1f77b4', 'C half' : '#aec7e8', 'CI' : '#c6dbef'},     # blue
    3: {'data' : '#ff7f0e', 'C half' : '#ffbb78', 'CI' : '#fdd0a2'},     # orange
    4: {'data' : '#5254a3', 'C half' : '#9e9ac8', 'CI' : '#dadaeb'},     # purple
    1: {'data' : '#8c6d31', 'C half' : '#e7ba52', 'CI' : '#e7cb94'}}     # brown

plot_markers = { 1 : 'o', 2 : '^', 3 : 's', 4 : '*'}

# DEBUG VAR, IGNORE IF NOT DEBUGGING
PRINT_EXC = True

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
    if condition['denature method'] == "CD": # first point is 0, last point is 1      
        # condition1 rep1          
        for element in condition[con][rep]['ydata']: # read from the begining, find first data !=0, make it A
            if element != 0:
                A = element # set 0
                break
        for element in reversed(condition[con][rep]['ydata']): # read from the end, find first data !=0, make it B
            if element != 0:
                B = element # set 1
                break
            
        #if B > A CD, A > B HD
        # another column (tell is positive or negative slope)
        for i in range(len(condition[con][rep]['ydata'])):
            if condition[con][rep]['ydata'][i] != 0:
                element = (condition[con][rep]['ydata'][i] - A) / (B - A) # 3 for each data !=0, normalized_ydata=(data-A)/(B-A)
                condition[con][rep]['normalized_ydata'].append(element)
                condition[con][rep]['normalized_xdata'].append(condition[con][rep]['xdata'][i]) # 4 for each normalized_ydata, find the xdata, make it normalized_xdata
            else:        
                condition[con][rep]['normalized_ydata'].append("error")
                condition[con][rep]['normalized_xdata'].append(condition[con][rep]['xdata'][i])  
                
        tempRow.extend(condition[con][rep]['normalized_ydata'])
        
    elif condition['denature method'] == "HD": # first point is 1, last point is 0 
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
    condition[con]['all xplot'] = []
    condition[con]['all yplot'] = []
    for j in range(1, condition[i]['replicates'] + 1):
        condition[con]['all xplot'].extend(condition[con][j]['xplot'])
        condition[con]['all yplot'].extend(condition[con][j]['yplot'])
    condition[con]['popt'], condition[con]['pcov'] = curve_fit(sigmoid, condition[con]['all xplot'], condition[con]['all yplot'])
    condition[con]['fit output'].extend(condition[con]['popt'])#extend merge two list to one big list  
                                        #popt has all the parameter (variable) in regression equation (x,B, A, Chalf, b)
                                                                                    #row [-, 0, 1, 2,    3] in  condition[1]['fit output']   
    if condition['denature method'] == "CD":
        condition[con]['x'] = np.linspace(0, 4, 50)
    elif condition['denature method'] == "HD": 
        condition[con]['x'] = np.linspace(0, 1.2, 50)
    condition[con]['y'] = sigmoid(condition[con]['x'], *condition[con]['popt']) # *popt split the two variable 
    condition[con]['C half'] = condition[con]['fit output'][2]
    condition[con]['b'] = condition[con]['fit output'][3]


def confidence_interval(con):
    """calculate confidence interval for c 1/2 and b"""
    condition[con]['sum pts'] = []
    for j in range(1, condition[con]['replicates'] + 1):
        condition[con][j]['float num pts'] = float(condition[con][j]['num pts'])

    TEMP_VAR = 0
    for j in range(1, condition[con]['replicates'] + 1):
        TEMP_VAR += condition[con][j]['float num pts']
    condition[con]['sum pts'] = TEMP_VAR

    for j in range(1, condition[con]['replicates'] + 1):
        tempRow.append(condition[con][j]['num pts']) # [5]

    tempRow.append(condition[con]['sum pts']) # [6]

    condition[con]['std error'] = np.sqrt(np.diag(condition[con]['pcov']))
    condition[con]['fit output'].extend(condition[con]['std error'])#[4,5,6,7]=[Berror, Aerror, Chalferror, berror]
    condition[con]['CI_Chalf'] = t.ppf(.975, (condition[con]['sum pts'] - 1)) * condition[con]['fit output'][6] / np.sqrt(condition[con]['sum pts'])
    condition[con]['fit output'].append(condition[con]['CI_Chalf'])#[8]

    if condition['denature method'] == "CD":
        condition[con]['ratioTOrange'] = condition[con]['CI_Chalf'] / CD_conc_range
        condition[con]['fit output'].append(condition[con]['ratioTOrange'])#[9]
    elif condition['denature method'] == "HD":
        condition[con]['ratioTOrange'] = condition[con]['CI_Chalf'] / 1
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
    # changed below
    condition[con]['residuals'] = (condition[con]['all yplot']) - sigmoid(condition[con]['all xplot'],*condition[con]['popt'])
    condition[con]['ss_res'] = np.sum(condition[con]['residuals']**2)
    # changed below
    condition[con]['ss_tot'] = np.sum(((condition[con]['all yplot']) - np.mean(condition[con]['all yplot'])) ** 2)
    condition[con]['r squared'] = 1 - (condition[con]['ss_res'] / condition[con]['ss_tot'])
    condition[con]['fit output'].append(condition[con]['r squared'])#[15]; append add something to the list   
    if condition['denature method'] == "HD":
        condition[con]['Chalf temp'] = (condition[con]['C half'] * HD_temp_range) + HD_Temp_lowbound
        condition[con]['fit output'].append(condition[con]['Chalf temp'])#[16]
    elif condition['denature method'] == "CD":
        condition[con]['Chalf normalized'] = condition[con]['C half'] / CD_conc_range
        condition[con]['fit output'].append(condition[con]['Chalf normalized'])#[16]


def plot_color(col, con):
    """plot color"""
    for j in range(1, condition[con]['replicates'] + 1):
        plt.plot(condition[con][j]['xplot'], condition[con][j]['yplot'], color = colors[col]['data'],
                 ls = ':', marker = condition[con][j]['marker'],
                 label = condition[con][j]['notebook code'])

    pylab.plot(condition[con]['x'], condition[con]['y'], colors[col]['data'], label = condition[con]['description'] + ' fit')
    pylab.axvline(condition[con]['C half'], color = colors[col]['C half'], ls = '-', label = 'C1/2') # chalf
    if condition[con]['ratioTOrange'] < CI_cutoff:
        plt.axvspan((condition[con]['C half'] - condition[con]['fit output'][8]),
                    (condition[con]['C half'] + condition[con]['fit output'][8]),
                    color = colors[col]['CI']) # CI

infile = open(infile_locationAndName, "r") # "r" = reading mode
reader = csv.reader(infile) # create reader object
header = next(reader)

if (condition['denature method'] == "HD"):
    total_cell = 0
    
    for cell in header:
        total_cell += 1

    condition['conditions'] = 0
    cell_num = 4
    con_num = 0 # condition of current rep
    rep_num = 1 # current rep num
    con_dict = {}
    rep_dict = {}
    con_dict[1] = ""
    num_pts_str = ""

    while (cell_num < total_cell):
        num_pts_str = header[cell_num]
        underscore_loc = header[cell_num].find("_")
        con_letter = num_pts_str[underscore_loc + 1]
        found = False

        for i in con_dict:
            if (con_letter == con_dict[i]):
                found = True
                rep_dict[i] += 1
                rep_num = rep_dict[i]
                con_num = i
                condition[con_num]['replicates'] += 1
                condition[con_num][rep_num] = {}
                break
    
        if (not found):
            con_num += 1
            rep_dict[con_num] = 1
            condition['conditions'] += 1    # total conditions
            con_dict[con_num] = con_letter
            rep_num = 1
            condition[con_num] = {}
            condition[con_num][rep_num] = {}
            condition[con_num]['replicates'] = 1

        
        condition[con_num][rep_num]['notebook code'] = ""
        for num in range(0, underscore_loc):
            condition[con_num][rep_num]['notebook code'] += num_pts_str[num]

        condition[con_num][rep_num]['numPtsLocation'] = cell_num
        condition[con_num][rep_num]['Start'] = cell_num + 1
        condition[con_num][rep_num]['End'] = cell_num + 10
        cell_num += 11

condition[1]['description'] = 'Ad Libitum'
condition[2]['description'] = 'Dietary restriction'

finalData = []

for i in range(1, condition['conditions'] + 1):         # for each condition (i)
    for j in range(1, condition[i]['replicates'] + 1):  # for each replicate (j) of condition (i)
        read_xdata(i, j)                                # read xdata from input file for condition i, replicate j

for row in reader:
    """if statement should force EOF at first blank row.
    if input file has blank lines, you can disable this if/else statement, but
    the program may try to read blank lines after EOF
    """
    if any(x.strip() for x in row):
        try:
            for i in range(1, condition['conditions'] + 1):         # for each condition (i)
                for j in range(1, condition[i]['replicates'] + 1):  # for each replicate (j) of condition (i)
                    read_ydata(i, j)                                # read ydata from input file for condition i, replicate j
        except:
            pass
                  
        # basic information input 
        title = row[rowLoc]
        ID = row[IDLoc]
        proteinName = row[NameLoc]
        Peptide = row[PeptideLoc]
        
        for i in range(1, condition['conditions'] + 1):
            # assigns condition plot color; cycles between 4 colors
            condition[i]['plot color'] = (i % 4) + 1
            for j in range(1, condition[i]['replicates'] + 1):
                # assigns replicate plot marker; cycles between 4 markers
                condition[i][j]['marker'] = plot_markers[(j % 4) + 1]
    
        for i in range(1, condition['conditions'] + 1):
            for j in range(1, condition[i]['replicates'] + 1):
                condition[i][j]['num pts'] = row[condition[i][j]['numPtsLocation']]
        
        # basic information ouput 
        tempRow = [title] # [0]
        tempRow.append(ID) # [1]
        tempRow.append(proteinName) # [2]
        tempRow.append(Peptide) # [3]
        
        # checks if there was an error generating the graph
        error_graph = False
         
        # normalize the data
        try:
            for i in range(1, condition['conditions'] + 1):
                for j in range(1, condition[i]['replicates'] + 1):
                    normalize_data(i, j)
        except:
            if (PRINT_EXC):
                print("row " + title + ": Data could not be normalized")

            if (PRINT_ALL):
                """If PRINT_ALL is set to True, the error will print and
                the program will continue."""
                error_graph = True
        
        #Set data set for fit and plot
        try:
            for i in range(1, condition['conditions'] + 1):
                for j in range(1, condition[i]['replicates'] + 1):
                    plot_data(i, j)
        except:
            if (PRINT_EXC):
                print("row " + title + ": Data could not be compiled for plotting")

            if (PRINT_ALL):
                error_graph = True

        try:
            for i in range(1, condition['conditions'] + 1):
                fit_scurve(i)
                confidence_interval(i)
                r_squared(i)
                tempRow.extend(condition[i]['fit output'])
        except:
            if (PRINT_EXC):
                print("row " + title + ": Statistical data could not be generated")

            if (PRINT_ALL):
                error_graph = True

        if (condition['conditions'] == 2):
            # delta calculation; 'D' is shorthand for delta
            if condition['denature method'] == "CD":
                condition['D_C half'] = condition[2]['C half'] - condition[1]['C half']
                condition['D_C half percent'] = condition['D_C half'] / condition[1]['C half']
                condition['D_b'] = condition[2]['b'] - condition[1]['b']
                condition['D_b percent'] = condition['D_b'] / condition[1]['b']
            elif condition['denature method'] == "HD":
                condition['D_C half'] = condition[2]['Chalf temp'] - condition[1]['Chalf temp']
                condition['D_C half percent'] = condition['D_C half'] / condition[1]['Chalf temp']
                condition['D_b'] = condition[2]['b'] - condition[1]['b']
                condition['D_b percent'] = condition['D_b'] / condition[1]['b']
        
            tempRow.append(condition['D_C half'])
            tempRow.append(condition['D_C half percent'])
            tempRow.append(condition['D_b'])
            tempRow.append(condition['D_b percent'])
    
            # Check if C half is significant
            try:
                if(condition[1]['C half'] > condition[2]['C half']):
                    if(condition[1]['C half'] > condition[2]['CI upbound'] and condition[2]['C half'] < condition[1]['CI lowbound']):
                        CI_chalf_significant = "significant"
                elif(condition[1]['C half'] < condition[2]['CI lowbound'] and condition[2]['C half'] > condition[1]['CI upbound']):
                    CI_chalf_significant = "significant"
                else:
                    CI_chalf_significant = ""
                tempRow.append(CI_chalf_significant)
                
            except:
                CI_chalf_significant = ""
                tempRow.append(CI_chalf_significant)
                print("row " + title + ": Issue with C half significant check")

            # Check if there is a CI overlap
            try:
                if(condition[1]['C half'] > condition[2]['C half']):
                    if(condition[1]['CI lowbound'] < condition[2]['CI upbound']):
                        CI_chalf_overlap = "overlap"
                elif(condition[2]['CI lowbound'] < condition[1]['CI upbound']):
                    CI_chalf_overlap = "overlap"
                else:
                    CI_chalf_overlap = ""
                tempRow.append(CI_chalf_overlap)
                
            except:
                CI_chalf_overlap = ""
                tempRow.append(CI_chalf_overlap)
                print("row " + title + ": Issue with CI overlap check")

        # generate the figure
        try:
            if condition[1]['ratioTOrange'] < CI_cutoff and condition[2]['ratioTOrange'] < CI_cutoff:
                graphName = condition['run desc'] + "_" + "row#" + title + " (has CI)" + "\n" + ID + "||" + proteinName + "\n" + Peptide
            else:
                graphName = condition['run desc'] + "_" + "row#" + title + "\n" + ID + "||" + proteinName + "\n" + Peptide
    
            for i in range(1, condition['conditions'] + 1):
                plot_color(condition[i]['plot color'], i)
            
            # plot label
            if condition['denature method'] == "CD":
                pylab.xlabel('denaturant concentration')
                pylab.ylabel('% labeled')
            elif condition['denature method'] == "HD": 
                pylab.xlabel('normalized temperature')
                pylab.ylabel('% soluble')
            # pylab.ylim(0.3, 2.2)
            pylab.legend(loc='best', fontsize = 'x-small')
            pylab.title(graphName, fontsize = 10)
            
            # generate figure's filename
            figure_filename = outImage_location + title + "_" + proteinName
            """if there was an error in generating the figure, then
            ' (data issue)' will be added to the end of the figure name.
            Only applies if PRINT_ALL = true."""
            if (error_graph):
                figure_filename += " (data issue)"
            figure_filename += ".jpg"
            
            pylab.savefig(figure_filename)
            pylab.clf() # clearfig, start a new piece, make sure do it before or after a set of instruction

        except:
            if (PRINT_EXC):
                print("row " + title + ": Plot could not be generated")
                
            if (PRINT_ALL):
                error_graph = True

        finalData.append(tempRow)

    else:
        break

infile.close()  

# output the data file        
outfile = open(outfile_locationAndName, "w", newline = "")        
writer = csv.writer(outfile) # create write object
newheader = ["row", "protein ID", "ProteinNAme", "Peptide"]

# create headers
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
    condition[i]['newheader fit output'] += ["b", "B_err", "A_err",
             "Chalf_err", "b_err", "CI_Chalf", "CIratioTOrange",
             "CI_Chalf_low", "CI_Chalf_up", "CI_b", "CI_b_low", "CI_b_up",
             "r_square"]

    if condition['denature method'] == "CD":
        temp_string = "condition" + str(i) + "_Chalf_normalized"
        condition[i]['newheader fit output'] += [temp_string]
    elif condition['denature method'] == "HD":
        temp_string = "condition" + str(i) + "_Chalf_temp"
        condition[i]['newheader fit output'] += [temp_string]

    if (condition['conditions'] == 2 and i == 2):
        condition[i]['newheader fit output'] += ["dChalf", "dChalf_%", "db",
                 "db%", "CI_Chalf_significant", "CI_chalf_overlap"]

if condition['denature method'] == "CD":
    # normalized header
    for i in range(1, condition['conditions'] + 1):
        for j in range(1, condition[i]['replicates'] + 1):
            newheader.extend(CD_xaxis_FinalGdmClConcentration_value)
elif condition['denature method'] == "HD":
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