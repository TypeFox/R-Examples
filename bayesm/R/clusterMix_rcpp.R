clusterMix=function(zdraw,cutoff=.9,SILENT=FALSE,nprint=BayesmConstant.nprint){
#
#
# revision history:
#   written by p. rossi 9/05
#
# purpose: cluster observations based on draws of indicators of 
#   normal mixture components
#
# arguments:
#   zdraw is a R x nobs matrix of draws of indicators (typically output from rnmixGibbs)
#   the rth row of zdraw contains rth draw of indicators for each observations
#   each element of zdraw takes on up to p values for up to p groups. The maximum
#   number of groups is nobs.  Typically, however, the number of groups will be small
#   and equal to the number of components used in the normal mixture fit.
#
#   cutoff is a cutoff used in determining one clustering scheme it must be 
#   a number between .5 and 1.
# 
#   nprint - print every nprint'th draw
#
# output:
#   two clustering schemes each with a vector of length nobs which gives the assignment
#   of each observation to a cluster
#
#   clustera (finds zdraw with similarity matrix closest to posterior mean of similarity)
#   clusterb (finds clustering scheme by assigning ones if posterior mean of similarity matrix
#             > cutoff and computing associated z )
#
# define needed functions
#
# ------------------------------------------------------------------------------------------   

#
# check arguments
#
if(missing(zdraw)) {pandterm("Requires zdraw argument -- R x n matrix of indicator draws")}
if(nprint<0) {pandterm('nprint must be an integer greater than or equal to 0')}
#
# check validity of zdraw rows -- must be integers in the range 1:nobs
#
nobs=ncol(zdraw)
R=nrow(zdraw)
if(sum(zdraw %in% (1:nobs)) < ncol(zdraw)*nrow(zdraw))
   {pandterm("Bad zdraw argument -- all elements must be integers in 1:nobs")}
cat("Table of zdraw values pooled over all rows",fill=TRUE)
print(table(zdraw))
#
# check validity of cuttoff
if(cutoff > 1 || cutoff < .5) {pandterm(paste("cutoff invalid, = ",cutoff))}

###################################################################
# Keunwoo Kim
# 10/06/2014
###################################################################
out=clusterMix_rcpp_loop(zdraw, cutoff, SILENT, nprint)
###################################################################

return(list(clustera=as.vector(out$clustera),clusterb=as.vector(out$clusterb)))
}