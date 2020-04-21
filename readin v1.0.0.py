"""
author's note: vil schoenheit please step on me
Loc variables are broken:
    -taken from global var value in this code, not updating w/ functions
    -please. please why cant you just work i am so tired i have so many finals
"""
# remove these once file contents finalized
# import numpy as np
# import pylab
import csv
# from scipy.optimize import curve_fit
# from scipy.stats import t
# import matplotlib.pyplot as plt
import file_combine as fc

condition = {
    # run information
    'run desc' : 'testrun_desc', # description of the run
    'denature method' : 'CD', # choose 'CD' for chemical denature, 'HD' for heat denature
    'slope' : "positive", # ignore this; declaring a variable for later use
    'conditions' : 0,
}

conditions = []
header = []
rowLoc = 0
IDLoc = 1
NameLoc = 2
PeptideLoc = 3

def read_masterfile(master_file):
    """Reads in the user input from the master input file

    Returns infile_locationAndName, outfile_locationAndName, outImage_location
    """
    infile = open(master_file, "r")
    reader = csv.reader(infile)
    header = next(reader)
    
    run_desc_index = 0
    out_name_index = 0
    graph_index = 0
    cond_index = 0
    rep_index = 0
    inf_index = 0
    inf_direct = 0
    
    for cell in header:
        if ("run description" in header):
            run_desc_index = header.index("run description")
        if ("outfile name" in header):
            out_name_index = header.index("outfile name")
        if ("graph output folder" in header):
            graph_index = header.index("graph output folder")
        if ("conditions" in header):
            cond_index = header.index("conditions")
        if ("replicates" in header):
            rep_index = header.index("replicates")
        if ("infile directory" in header):
            inf_direct = header.index("infile directory")
        if ("infile name" in header):
            inf_index = header.index("infile name")
    
    conditions = []
    first_row = False
    for row in reader:
        # if statement below should force EOF at first blank row
        if any(x.strip() for x in row):
            if (not first_row):
                condition['run desc'] = row[run_desc_index]
        
                if (row[out_name_index].find(".csv") == -1):
                    outfile_locationAndName = row[out_name_index] + ".csv"
                else:
                    outfile_locationAndName = row[out_name_index]
                print(outfile_locationAndName)
                        
                outImage_location = row[graph_index] + "//"
                print(outImage_location)
            
            # if the replicate cell is blank, the program assumes user input is
            # finished
            if (row[rep_index] == ""):
                break
    
            if (row[cond_index] in conditions):
                pass
            elif (row[cond_index] != ""):
                conditions += row[cond_index]
                
            print(row[rep_index])
    
            # only runs if infile cell isn't blank
            if(row[inf_index] != ""):
                if (row[inf_index].find(".csv") == -1):
                    infile_locationAndName = row[inf_index] + ".csv"
                else:
                    infile_locationAndName = row[inf_index]
                print(infile_locationAndName)
            
            if(row[inf_direct] != ""):
                if (row[inf_direct].find("//") == -1):
                    dir_path = row[inf_direct] + "//"
                else:
                    dir_path = row[inf_direct]
                print(dir_path)
            first_row = True
    
    
    # file combine (file_combine.py)
    first_file = dir_path + "protein-peptides.csv"
    second_file = dir_path + "proteins.csv"
    infile_locationAndName = dir_path + "test_infile.csv"
    fc.file_combine(first_file, second_file, infile_locationAndName)
    
    conditions = []
    
    return(infile_locationAndName, outfile_locationAndName, outImage_location)

def define_headers(infile_locationAndName):
    infile = open(infile_locationAndName, "r") # "r" = reading mode
    reader = csv.reader(infile) # create reader object
    header = next(reader)
    
    # find which columns have what data
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

    elif (condition['denature method'] == "CD"):
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
            
        if ("Protein Accession" in header):
            IDLoc = header.index("Protein Accession")
        elif ("Protein" in header):
            IDLoc = header.index("Protein")
        else:
            IDLoc = 1

        if ("Description" in header):
            NameLoc = header.index("Description")
        else:
            NameLoc = 2
    
        if ("real peptide" in header):
            PeptideLoc = header.index("real peptide")
        elif ("Peptide" in header):
            PeptideLoc = header.index("Peptide")
        else:
            PeptideLoc = 3