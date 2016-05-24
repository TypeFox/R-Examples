summary.monte1<-function(object, digits=3, ...){
z<-object
nvar<-z$nvar


obs.skew.matrix<-obs.kurt.matrix<-matrix(99,1,nvar)                                  
                

  obs.cormat  <- cor(z$data)
  obs.skewvec <- apply(z$data, 2, skew)
  obs.kurtvec <- apply(z$data, 2, kurt)


ans <- z["call"]


cat("\nCall monte1:","\n")
print(z$call)

cat("\nSeed = ",z$seed,"\n")



cat("\n\n\nExpected correlation matrix:","\n")
print(round(z$cormat,digits))

cat("\n\nObserved correlation matrix:","\n")
print(round(obs.cormat,digits))


cat("\nExpected indicator skewness:")
    print(round(z$skewvec, digits))
   
cat("\nObserved indicator skewness:","\n")
print(round(obs.skewvec,digits))

cat("\n\nExpected indicator kurtosis:")
    print(round(z$kurtvec,digits))

cat("\nObserved indicator kurtosis:","\n")
 print(round(obs.kurtvec,digits))


 
 ans$seed<-z$seed
 ans$nvar <- nvar
 ans$cormat  <- z$cormat
 ans$obs.cormat<-obs.cormat
 ans$skewvec <- z$skewvec
 ans$obs.skewvec <-obs.skewvec
 ans$kurtvec <- z$kurtvec
 ans$obs.kurtvec <- obs.kurtvec
 
 
 class(ans)<-"summary.monte1"
invisible(ans)
}
