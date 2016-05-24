basis<-function(x,nb){
  B<-outer(x,2*pi*(1:nb-1),FUN="*")
  B<-cos(B)
  return(B)
}

vbetafun<-function(smthpar,nb,varmult){
  vbeta<-smthpar^abs(1:nb-1)
  vbeta<-varmult^2*vbeta/sum(vbeta)
  vbeta<-diag(vbeta)
  return(vbeta)
}

# vbetafun<-function(lambda=2,m=5,nb,varmult){
#   vbeta<-(m+lambda*(1:nb-1))^-m
#   vbeta<-varmult^2*vbeta/sum(vbeta)
#   vbeta<-diag(vbeta)
#   return(vbeta)
# }
