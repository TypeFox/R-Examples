# An example rule file

# Testing the example database db.csv

# Rules regarding ELEVATION
is.integer(ELEVATION)  # is right datatype
is_within_range(ELEVATION, -1, 8)  # is between min max
-2 < ELEVATION  #
ELEVATION < 9  #
-2 < ELEVATION & ELEVATION < 9  # alternative form of range testing

is.character(LOCATION)  #
LOCATION %in% LETTERS[1:10]  #

is.character(GENUS)  #
is_proper_name(GENUS) 
