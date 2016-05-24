# An example rule file

# Testing the example database db.csv

# Rules regarding ELEVATION
sapply(ELEVATION, is.integer)  # is right datatype
is_within_range(ELEVATION, -1, 8)  # is between min max
sapply(LOCATION, is.character)  #
LOCATION %in% LETTERS[1:10]  #


sapply(GENUS, is.character)  #
is_proper_name(GENUS) 
