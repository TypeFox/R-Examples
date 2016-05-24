ogk<-function(x,sigmamu=taulc,v=gkcov,n.iter=1,beta=.9,...){
#
# Compute robust (weighted) covariance matrix in Maronna and Zamar
# (2002, Technometrics, eq. 7).
#
# x is an n by p matrix
# n.iter number of iterations. 1 seems to be best
# sigmamu is any user supplied function having the form
#   sigmamu(x,mu.too=F) and which computes a robust measure of
#   of dispersion if mu.too=F. If mu.too=T, it returns
#   a robust measure of location as well.
# v is any robust covariance
#
if(!is.matrix(x))stop("x should be a matrix")
x<-elimna(x)  # remove any rows with missing data
temp<-ogk.pairwise(x,sigmamu=sigmamu,v=v,n.iter=n.iter,beta=beta,...)
list(center=temp$wcenter,cov=temp$wcovmat)
}

