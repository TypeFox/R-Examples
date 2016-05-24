print.summary.monte1<-function(x,digits=3, ...){
z<-x
nvar<-z$nvar



cat("\nCall monte1:","\n")
print(z$call)

cat("\nSeed = ",z$seed,"\n")



cat("\n\n\nExpected correlation matrix:","\n")
print(round(z$cormat,digits))

cat("\n\nObserved correlation matrix:","\n")
print(round(z$obs.cormat,digits))


cat("\nExpected indicator skewness:")
    print(round(z$skewvec, digits))
   
cat("\nObserved indicator skewness:","\n")
print(round(z$obs.skewvec,digits))

cat("\n\nExpected indicator kurtosis:")
    print(round(z$kurtvec,digits))

cat("\nObserved indicator kurtosis:","\n")
 print(round(z$obs.kurtvec,digits))


}
