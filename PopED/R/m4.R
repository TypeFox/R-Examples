## Function translated automatically using 'matlab.to.r()'
## Author: Andrew Hooker

m4 <- function(mv,n){
#
# size: (samps per individ^2 x samps per individ^2)
  ns=n^2
  mm4=zeros(ns,ns)
  for(ct1 in 1:n){
    for(ct2 in 1:n){
      mm4[((ct1-1)*n+1):((ct1-1)*n+n),((ct2-1)*n+1):((ct2-1)*n+n)]= mv[ct1,ct2]*mv*2
   }
  }
  ret=mm4
return( ret ) 
}

