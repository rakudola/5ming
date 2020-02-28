# Changelog
NOTE: minor updates (change in last number in software number) were largely not recorded before 2019/12/2

## v2.7.0
(2020/2/11)
- Beginning read-in input
- Currently does not scale outside of 2-condition, 4-replicate run from modified input, but development is underway!

### v2.7.1
(2020/2/13)
- Support ending for HD input. Next version will only accept CD input

## v2.6.0
(2019/1/21)
- Program now can scale to different numbers of conditions and replicates save for a handful of bugs
- Prior raw data manipulation and some in-code variable changes still required, but no more copy-pasting functions!

### v2.6.1
(2019/1/21)
- Errors (from the main-body try statement) are now printed to console
- Colors and markers for plot are now chosen automatically; no need to change variables for a normal run
- Added CI_chalf_significant and CI_chlaf_overlap calculations back in as a comment; will need to be debugged
- Code should actually stop at EOF instead of running for a bunch of empty lines
- Plots for some data that previously wouldn't run can now be generated. This can be toggled with the PRINT_ALL variable
- Required user input is now at the top of the code and divided from the rest of the code

### v2.6.2
(2020/1/24)
- If PRINT_ALL is set to True, graphs that wouldn't have been generated otherwise will now include " (data issue)" in filename
- C half significance and confidence interval overlap code is working and printing to csv output
- One 2-condition-4-replicate-only bug squished!

### v2.6.3
(2020/1/30)
-Squished a few more bugs that weren't allowing for dynamic number of conditions and replicates

### v2.6.4
(2020/1/30)
- Finished debugging most of the 2-condition-4-replicate-only bugs. Still need to adjust C 1/2 significance + CI overlap, but those should (finally) be the last of these specific bugs

### v2.6.5
(2020/1/30)
- Delta calculations, C 1/2 significance, and CI overlap are reserved for 2-condition runs
- Program is finally fully functional with different number of conditions and replicates!!

### v2.6.6
(2020/1/30)
- Changed how plot generation handles error
- Program now closes infile after finished reading from it (rather than at the end of the program)

## v2.5.0
2019/12/5
- Rolled back to v2.3.3 due to issue with plot_math function
- Breaking up plot_math into smaller functions during implementation
- New fit_scurve function

### v2.5.1
2019/12/6
- New confidence_interval function

### v2.5.2
2019/12/10
- New r_squared function
- Removed some extraneous code

### v2.5.3
(2019/12/10)
- plot_color function now operational
- Added whitespace for legibility, changed 'Chalf' vars to 'C half'

### v2.5.4
(2019/1/17)
- Added delta_calculation function; only works for 2 conditions at the moment

### v2.5.5
(2019/1/21)
- Moved delta_calculation function back into body; was not functioning properly as a separate function (still only works with 2 conditions for the moment)
- Condensed normalized header printing (now allows for different numbers of conditions and replcates)

## v2.4.0 (BUG)
2019/12/2
- plot_math functional, but only works for 2 conditions w/ 4 replicates each right now
- I spent a week debugging plot_math. It gets its own update.

## v2.3.0
2019/11/18
- Iterative functions now working
- Condition 2 for current input file data now in dictionary implementation (program still only compatible with input file Lavender provided)

## v2.2.0
2019/11/11
- Fixed list scope bug

## v.2.1.0
2019/11/5
- Started framework for dynamic number of conditions and replicates

## v2.0.1
2019/10/31
- Changed base program to Lavender's 2-condition, 4-replicate file
- Began implementing functions to condense code
- Still only compatible with input file Lavender provided

## v1.0.1
2019/9/18
- Original code from Lavender; 1 condition only
### v1.0.2
2019/10/10
- Added exception to y-value input (for row in reader) to get around a ValueError preventing graphs from printing with deleted columns in input file
- Minor convention updates
