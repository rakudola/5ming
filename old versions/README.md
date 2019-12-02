# Changelog
NOTE: minor update (change in last number in software number) changes were largely not recorded before 2019/12/2

## v2.4.0
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
