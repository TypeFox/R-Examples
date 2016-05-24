`print.binomRDtest` <-
function(x, digits=4, ...)
{

if(x$method=="Wald")
 {methvar<-"Wald"}
if(x$method=="ADD1")
 {methvar<-"Add-1"}
if(x$method=="ADD2")
 {methvar<-"Add-2"}

# The p-value of the maximum test

if(x$dist=="MVN")
{
cat("Union intersection test using",methvar," variance estimator:","\n")
cat("P-value of the maximum test:","\n")
print(round(x$pval, digits=digits))
cat(" ","\n")
}
else{
cat("Local p-values using",methvar," variance estimator:","\n")
cat(" ","\n")
}


# A table of testresults

testresult <- cbind(x$estimate, x$teststat, x$p.val.adj)

colnames(testresult) <- c("estimate","testat","p.val.adj") 

print(round(testresult, digits=digits))


cat(" ","\n")
cat("where proportions are the probability of the event", x$success, "\n")
cat(" ","\n")

invisible(x)

}

