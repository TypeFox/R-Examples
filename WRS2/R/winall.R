winall<-function(m,tr=.2){
#
#    Compute the Winsorized correlation and covariance matrix for the
#    data in the n by p matrix m.
#
#    This function also returns the two-sided significance level
#
if(is.data.frame(m))m=as.matrix(m)
if(!is.matrix(m))stop("The data must be stored in a n by p matrix")
wcor<-matrix(1,ncol(m),ncol(m))
wcov<-matrix(0,ncol(m),ncol(m))
siglevel<-matrix(NA,ncol(m),ncol(m))
for (i in 1:ncol(m)){
ip<-i
for (j in ip:ncol(m)){
val<-wincor(m[,i],m[,j],tr)
wcor[i,j]<-val$cor
wcor[j,i]<-wcor[i,j]
if(i==j)wcor[i,j]<-1
wcov[i,j]<-val$cov
wcov[j,i]<-wcov[i,j]
if(i!=j){
siglevel[i,j]<-val$siglevel
siglevel[j,i]<-siglevel[i,j]
}
}
}
list(cor=wcor,cov=wcov,p.values=siglevel)
}

