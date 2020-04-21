# C-Half
Software to visually organize PEAKS peptide analysis for protein unfolding (under development)

## Files
- v2.8.1.py: main file
- file combine v1.1.0.py: combines input files
- readin v1.0.0.py: reads in combined input files

## Current Versions
### v2.8.1 (2020/4/21)
- Beginning reimplementation of heat denaturation
- Minor bug fixes

### file combine 1.1.0.py (2020/4/10)
- Supported in main file v2.8.0-current version
- Takes three inputs (first input file, second input file, output file)
- Accessed in main code as fc

### readin 1.0.0.py (2020/4/21)
- Supported in main file v2.8.1-current version
- Takes in master input file and combined files from file combine
- Accessed in main code as ri

## Goals
- [x] Allow for multiple condition/replicate analysis without copy-pasting code (c. 2020/1/21)
- [x] Name this software (c. 2020/2/21)
- [x] Allow for .csv input (in case not using PEAKS) if correct headers (c. 2020/2/28)
- [ ] Modification frequency + cleaning up and consolidating files (2019/9/24)
- [ ] Analyze PEAKS output directly (2019/10/24)
- [ ] Option to combine replicates to single curve (2019/10/24)
- [ ] r^2 value (2019/11/14)*
- [ ] C 1/2 value in range - 0 <  C 1/2 < 3.48 (2019/11/14)*
- [ ] C 1/2 vs. location of peptide graphs (2019/11/14)*
- [ ] Match non-unique peptides to their respective proteins (2019/11/26)**
- [ ] Allow run analysis even if missing replicates/conditions (2019/11/26)**
- [ ] Better documentation (2020/1/21)
- [ ] Use con/rep description from folder name (2020/3/20)
- [ ] Master input file (2020/3/20)
- [ ] Limit graphing to 3 conditions (2020/3/20)
- [ ] Multiple input files - per condition/replicate (2020/3/26)
- [ ] Add heat denature back in + fix temperature normalization -> reconversion (2020/3/27)
- [ ] Combine based on modifications -> 1 curve (2020/3/27)

### Footnotes
- *see page 18 of notebook
- **see page 22
