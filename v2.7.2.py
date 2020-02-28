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
    'run desc' : 'testrun_desc', # description of the run
    'denature method' : 'CD', # choose 'CD' for chemical denature, 'HD' for heat denature
    'slope' : "positive", # ignore this; declaring a variable for later use

    # condition 1 description
    1: {'description' : 'Single-Condition Run'},

    # If you have more conditions, you can uncomment the lines below or
    # copy-paste to make your own. Just keep going up one number per line and
    # make sure they're in the same order as your input file.

    2: {'description' : 'put description here'},
    # 3: {'description' : 'put description here'},
    # 4: {'description' : 'put description here'},
}

# file selection
infile_locationAndName = "../CD input/CD_LP_1025+1023+1019_python.csv"
outfile_locationAndName = "../test_CDoutput_3.csv"
outImage_location = "../CD output_3//"

"""
Code below here does not need to be modified for a standard run
"""

# basic information --set up variable
CD_xaxis_FinalGdmClConcentration_value = [0, 0.43, 0.87, 1.3, 1.74, 2.17, 2.61, 3.04, 3.48]
CD_conc_upbound = 3.48
CD_conc_lowbound = 0
CD_conc_range = CD_conc_upbound - CD_conc_lowbound
CI_cutoff = 0.3

# plot color code, usually don't need to change
colors = {
    # 'data' is for data and fit
    2: {'data' : '#1f77b4', 'C half' : '#aec7e8', 'CI' : '#c6dbef'},     # blue
    3: {'data' : '#ff7f0e', 'C half' : '#ffbb78', 'CI' : '#fdd0a2'},     # orange
    4: {'data' : '#5254a3', 'C half' : '#9e9ac8', 'CI' : '#dadaeb'},     # purple
    1: {'data' : '#8c6d31', 'C half' : '#e7ba52', 'CI' : '#e7cb94'}}     # brown

plot_markers = { 1 : 'o', 2 : '^', 3 : 's', 4 : '*'}

def sigmoid(x, B, A, Chalf, b):
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

    # first point is 0, last point is 1
    for element in condition[con][rep]['ydata']: # read from the begining, find first data != 0, make it A
        if element != 0:
            A = element # set 0
            break
    for element in reversed(condition[con][rep]['ydata']): # read from the end, find first data != 0, make it B
        if element != 0:
            B = element # set 1
            break
        
    if (A > B):
        temp_var = A
        A = B
        B = temp_var
        # condition['denature method'] == "CD"
        condition['slope'] = "positive"
    else:
        # condition['denature method'] == 'HD'
        condition['slope'] = "negative"
        
    # if B > A CD, A > B HD; another column (pos/neg slope)
    for i in range(len(condition[con][rep]['ydata'])):
        if condition[con][rep]['ydata'][i] != 0:
            element = (condition[con][rep]['ydata'][i] - A) / (B - A) # 3 for each data != 0, normalized_ydata=(data-A)/(B-A)
            condition[con][rep]['normalized_ydata'].append(element)
            condition[con][rep]['normalized_xdata'].append(condition[con][rep]['xdata'][i]) # 4 for each normalized_ydata, find the xdata, make it normalized_xdata
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
    condition[con]['x'] = np.linspace(0, 4, 50)
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

    condition[con]['ratioTOrange'] = condition[con]['CI_Chalf'] / CD_conc_range
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
    condition[con]['ss_res'] = np.sum(condition[con]['residuals'] ** 2)
    # changed below
    condition[con]['ss_tot'] = np.sum(((condition[con]['all yplot']) - np.mean(condition[con]['all yplot'])) ** 2)
    condition[con]['r squared'] = 1 - (condition[con]['ss_res'] / condition[con]['ss_tot'])
    condition[con]['fit output'].append(condition[con]['r squared'])#[15]; append add something to the list
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

# DEBUG VAR, IGNORE IF NOT DEBUGGING
PRINT_EXC = True

"""
Main code
"""

infile = open(infile_locationAndName, "r") # "r" = reading mode
reader = csv.reader(infile) # create reader object
header = next(reader)

rowLoc = 0
IDLoc = 0
NameLoc = 0
PeptideLoc = 0

# find which columns have what data
if (condition['denature method'] == "CD"):
    condition['conditions'] = 0
    conditions = []
    rep_num = 1
    con_num = 0
    cell_num = 0

    for cell in header:
        if (cell.find("#pt") != -1):
            underscore_loc = cell.find("_")
            notebook_str = ""

            for i in range(0, underscore_loc):
                notebook_str += cell[i]
            
            if (cell[underscore_loc + 1].lower() in conditions):
                rep_num += 1
                con_num = conditions.index(cell[underscore_loc + 1].lower()) + 1
                condition[con_num]['replicates'] += 1
                condition[con_num][rep_num] = {}

            else:
                conditions += cell[underscore_loc + 1].lower()
                con_num += 1
                condition['conditions'] += 1
                condition[con_num] = {}
                condition[con_num]['description'] = "dummy_desc"

                condition[con_num]['replicates'] = 1
                rep_num = 1
                condition[con_num][rep_num] = {}
            
            condition[con_num][rep_num]['numPtsLocation'] = cell_num
            condition[con_num][rep_num]['notebook code'] = notebook_str

        elif (cell == "0"):
            condition[con_num][rep_num]['Start'] = cell_num
        elif (cell == "3.48"):
            condition[con_num][rep_num]['End'] = cell_num
        
        cell_num += 1

    if ("raw#" in header):
        rowLoc = header.index("raw#")
    elif ("row#" in header):
        rowLoc = header.index("row#")
    else:
        rowLoc = 0
        
    IDLoc = header.index("Protein Accession")
    NameLoc = header.index("Description")

    if ("real peptide" in header):
        PeptideLoc = header.index("real peptide")
    elif ("Peptide" in header):
        PeptideLoc = header.index("Peptide")
    else:
        PeptideLoc = 3

finalData = []

for i in range(1, condition['conditions'] + 1):         # for each condition (i)
    for j in range(1, condition[i]['replicates'] + 1):  # for each replicate (j) of condition (i)
        read_xdata(i, j)                                # read xdata from input file for condition i, replicate j

for row in reader:
    """if statement should force EOF at first blank row.
    if input file has blank lines, you can disable this if/else statement, but
    the program may try to read blank lines after EOF.
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
        
        # Set data set for fit and plot
        try:
            for i in range(1, condition['conditions'] + 1):
                for j in range(1, condition[i]['replicates'] + 1):
                    plot_data(i, j)
        except:
            if (PRINT_EXC):
                print("row " + title + ": Data could not be compiled for plotting")

            if (PRINT_ALL):
                error_graph = True

        # fit data to s curve
        try:
            for i in range(1, condition['conditions'] + 1):
                fit_scurve(i)
        except:
            if (PRINT_EXC):
                print("row " + title + ": Data could not be fit")

            if (PRINT_ALL):
                error_graph = True

        # calculate confidence interval
        try:
            for i in range(1, condition['conditions'] + 1):
                confidence_interval(i)
        except:
            if (PRINT_EXC):
                print("row " + title + ": Confidence interval could not be calculated")

            if (PRINT_ALL):
                error_graph = True

        # calculate r squared
        try:
            for i in range(1, condition['conditions'] + 1):
                r_squared(i)
        except:
            if (PRINT_EXC):
                print("row " + title + ": R squared be calculated")

            if (PRINT_ALL):
                error_graph = True
            
        tempRow.extend(condition[i]['fit output'])

        if (condition['conditions'] == 2):
            # delta calculation; 'D' is shorthand for delta
            condition['D_C half'] = condition[2]['C half'] - condition[1]['C half']
            condition['D_C half percent'] = condition['D_C half'] / condition[1]['C half']
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
            if (condition['conditions'] == 2):
                if condition[1]['ratioTOrange'] < CI_cutoff and condition[2]['ratioTOrange'] < CI_cutoff:
                    graphName = condition['run desc'] + "_" + "row#" + title + " (has CI)" + "\n" + ID + "||" + proteinName + "\n" + Peptide
                else:
                    graphName = condition['run desc'] + "_" + "row#" + title + "\n" + ID + "||" + proteinName + "\n" + Peptide
            else:
                if condition[1]['ratioTOrange'] < CI_cutoff:
                    graphName = condition['run desc'] + "_" + "row#" + title + " (has CI)" + "\n" + ID + "||" + proteinName + "\n" + Peptide
                else:
                    graphName = condition['run desc'] + "_" + "row#" + title + "\n" + ID + "||" + proteinName + "\n" + Peptide
    
            for i in range(1, condition['conditions'] + 1):
                plot_color(condition[i]['plot color'], i)
            
            # plot label
            pylab.xlabel('denaturant concentration')
            pylab.ylabel('% labeled')
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

        tempRow.append(condition['slope'])
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

    temp_string = "condition" + str(i) + "_Chalf_normalized"
    condition[i]['newheader fit output'] += [temp_string]

    if (condition['conditions'] == 2 and i == 2):
        condition[i]['newheader fit output'] += ["dChalf", "dChalf_%", "db",
                 "db%", "CI_Chalf_significant", "CI_chalf_overlap"]

# normalized header
for i in range(1, condition['conditions'] + 1):
    for j in range(1, condition[i]['replicates'] + 1):
        newheader.extend(CD_xaxis_FinalGdmClConcentration_value)

# fit header
for i in range(1, condition['conditions'] + 1):
    newheader.extend(condition[i]['newheader fit output'])

writer.writerow(newheader)

for addvariable in finalData:
    writer.writerow(addvariable)
    
outfile.close()