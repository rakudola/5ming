# 5ming
Python code to visually organize PEAKS peptide analysis for protein unfolding

## Current Version
### Warning: currently uploaded version is wrong version haha oops
v2.6.2 (2020/1/24)
- If PRINT_ALL is set to True, graphs that wouldn't have been generated otherwise will now include " (data issue)" in filename
- C half significance and confidence interval overlap code is working and printing to csv output
- One 2-condition-4-replicate-only bug squished!

## Goals
- [x] Allow for multiple condition/replicate analysis without copy-pasting code (c. 2020/1/21)
- [ ] Better documentation (2020/1/21)
- [ ] Move functions not returning anything back into main code for legibility? (2020/1/21)
- [ ] Name this software
- [ ] Modification frequency + cleaning up and consolidating files (2019/9/24)
- [ ] Analyze PEAKS output directly (2019/10/24)
- [ ] Allow for .csv input (in case not using PEAKS) if correct headers (2020/1/17)
- [ ] Option to combine replicates to single curve (2019/10/24)
- [ ] r^2 value (2019/11/14)*
- [ ] C 1/2 value in range - 0 <  C 1/2 < 3.48 (2019/11/14)*
- [ ] C 1/2 vs. location of peptide graphs (2019/11/14)*
- [ ] Match non-unique peptides to their respective proteins (2019/11/26)**
- [ ] Allow run analysis even if missing replicates/conditions (2019/11/26)**

### Footnotes
- *see page 18 of notebook
- **see page 22

## Bugs
- [x] csv.reader reading blank lines after data (s. 2020/1/23)
- [ ] Certain functions only work with 2-condition, 4-replicate runs
