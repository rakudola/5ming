# Changelog
NOTE: minor updates (change in last number in software number) were largely not recorded before 2019/12/2

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

### v2.5.
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
