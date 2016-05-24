#########################################################################################
# A replication of MatLab repmat function!
#
#
# Source:
# R FOR OCTAVE USERS
#    version 0.4
#    Copyright (C) 2001 Robin Hankin
# 	 http://cran.r-project.org/doc/contrib/R-and-octave.txt
#########################################################################################
repmat <- function(a,n,m) {kronecker(matrix(1,n,m),a)}