near3d<-function(x,pt,fr=.8,m){
# determine which values in x are near pt
# based on fr * cov.mve
#
# x is assumed to be an n by p matrix
# pt is a vector of length p (a point in p-space).
# m is cov.mve(x) computed by runm3d
#

if(!is.matrix(x))stop("Data are not stored in a matrix.")
dis<-sqrt(mahalanobis(x,pt,m$cov))
dflag<-dis < fr
dflag
}
