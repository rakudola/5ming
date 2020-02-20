import csv

first_file = "../CD input/protein-peptides.csv" # should be protein-peptides file
first_headers = { # holds the names of the headers to look for in first_file
        # 'variable name' (do not change) : "header name"
        'protein_acc' : "Protein Accession", # header for protein accession
        'full_pep' : "Peptide", # heaeder for the full peptide
        'start_index' : "A000", # header for first point
        'end_index' : "A6", # header for last point
        'start_loc' : "Start", # header for peptide start location
        'end_loc' : "End", # header for peptide end location
}

second_file = "../CD input/proteins.csv" # should be proteins file
second_headers = { # holds the names of the headers to look for in second_file
        # 'variable name' (do not change) : "header name"
        'accession' : "Accession", # header for protein accession
        'description' : "Description" # header for protein description
}
output_file = "../CD input/test file combine.csv"

"""
Changing code below here will affect functionality
"""

# INPUT FIRST FILE
dataA = []
FIRST_ID = 1 # first file's protein accession index in headerA

with open(first_file, "r") as first_infile:
    Areader = csv.reader(first_infile) # create reader object
    full_header = next(Areader)
    headerA = ["row#", "Protein", "Description", "real peptide", "#reporter",
               "Start", "End", "#pt", "0", "0.43", "0.87", "1.3", "1.74",
               "2.17", "2.61", "3.04", "3.48"]

    protein_acc = full_header.index(first_headers['protein_acc'])
    full_pep = full_header.index(first_headers['full_pep'])
    start_index = full_header.index(first_headers['start_index'])
    end_index = full_header.index(first_headers['end_index'])
    start_loc = full_header.index(first_headers['start_loc'])
    end_loc = full_header.index(first_headers['end_loc'])

    for row in Areader:
        row_info = []
        row_info += [row[protein_acc]]
        
        PEP_START = 2 # first index in peptide that is actual peptide name
        PEP_END = len(row[full_pep]) - 2 # last index in peptide that is actual peptide name
        full_peptide = ""
        for i in range(PEP_START, PEP_END):
            full_peptide += row[full_pep][i]
        
        row_info += [full_peptide]
        row_info += [""] # num reporters
        row_info += [row[start_loc]]
        row_info += [row[end_loc]]

        points_list = []
        num_pts = 0
        for i in range(start_index, end_index + 1):
            points_list += [row[i]]
            if (row[i] != "0"):
                num_pts += 1

        if (num_pts < 4):
            continue # if there are less than 4 points, skip adding this row

        row_info += [num_pts]
        row_info += points_list

        dataA += [row_info]
    
dataA.sort() # sorts info alphabetically by Protein column

# INPUT SECOND FILE
dataB = []
SECOND_ID = 0 # second file's protein accession index in headerB

with open(second_file, "r") as second_infile:
    Breader=csv.reader(second_infile) # create reader object
    full_header = next(Breader)
    headerB = ["Accession", "Description"]

    protein_acc = full_header.index(second_headers['accession'])
    desc_index = full_header.index(second_headers['description'])

    for row in Breader:
        row_info = []
        row_info += [row[protein_acc]]

        cut_from = row[desc_index].find(" OS=")
        full_desc = ""
        for i in range(0, cut_from):
            full_desc += row[desc_index][i]

        row_info += [full_desc]
        dataB += [row_info]

# USE THE KEY AND COMBINE TWO LIST-WAY2
combineData = []

for i in range(len(dataA)):
    dataA[i].insert(0, str(i + 1))
    for j in range(len(dataB)):
        if dataB[j][SECOND_ID] == dataA[i][FIRST_ID]:
             dataA[i].insert(2, dataB[j][1])



# OUTPUT DATA FILE   
with open(output_file,"w",newline="") as outfile:
    writer = csv.writer(outfile) # create write object
    newheader = headerA
    writer.writerow(newheader)
    for addvariable in dataA:
        writer.writerow(addvariable)