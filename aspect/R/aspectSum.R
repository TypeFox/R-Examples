#sum of the r^pow elements of the correlation matrix 

aspectSum<-function(r,pow = 1) {
m<-dim(r)[1]
list(f=sum(r^pow),g=pow*r^(pow-1))
}

