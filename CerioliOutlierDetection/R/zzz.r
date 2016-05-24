# Miscellaneous utility functions 
# for the 'CerioliOutlierDetection' package
#
# Christopher G. Green
# 2014

# Fraction of the sample used to compute the MCD in
# the maximum breakdown point case
max.bdp.mcd.alpha <- function(n,p) floor((n+p+1)/2)/n
