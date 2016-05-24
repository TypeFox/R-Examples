`print.poly3test` <-
function(x, digits=4, ...)
{


if(x$method=="BP")
 {methvar<-"Bailer-Portier"}
if(x$method=="BW")
 {methvar<-"Bieler-Williams"}
if(x$method=="ADD1")
 {methvar<-"Add-1"}
if(x$method=="ADD2")
 {methvar<-"Add-2"}


#  Table of sample estimates

sample.estimate<-rbind(x$sample.estimate$Y, x$sample.estimate$n, x$sample.estimate$nadj, x$sample.estimate$estimate)
rownames(sample.estimate)<-c("x","n", "adjusted n", "adjusted estimate")

cat("Sample estimates using poly-",x$k,"-adjustment","\n")
print(round(sample.estimate, digits=digits))
cat(" ","\n")

# The contrast matrix:

cat("Contrast matrix:","\n")
print(x$cmat)
cat(" ","\n")

if(x$dist=="MVN")
{
# The p-value of the maximum test
cat("Union-Intersection test using", methvar," variance estimator:","\n")
cat("P-value of the maximum test:","\n")
print(round(x$pval, digits=digits))
cat(" ","\n")
}
else{
cat("Local Wald test using",methvar," variance estimator:","\n")
}


# A table of testresults

testresult <- cbind(x$estimate, x$teststat, x$p.val.adj)

colnames(testresult) <- c("estimate","testat","p.val.adj") 

print(round(testresult, digits=digits))

cat(" ","\n")

invisible(x)

}

