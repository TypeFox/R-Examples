######################################################################
# mylevels() returns levels if given a factor, otherwise 0 
######################################################################
mylevels <- function(x) if (is.factor(x)) levels(x) else 0

