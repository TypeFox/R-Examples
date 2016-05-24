"summary.binomRDci"<-function(object,...)
{

aargs<-list(...)
if(is.null(aargs$digits))
 {digits<-4}
else
 {digits<-aargs$digits}


cat("Summary statistics: \n")

summarystat<-rbind("number of successes"=object$x,
"number of trials"=object$n,
"estimated success probability"=round(object$x/object$n, digits=digits))

colnames(summarystat)<-object$names

print(summarystat, digits=digits)


cat("\n Contrast matrix: \n")

print(object$cmat, digits=digits)


if(object$dist=="MVN")
{
cat("\n The estimated correlation matrix of the contrasts is: \n")

print(object$corrmat, digits=digits)

}

cat("\n")

print(object, digits=digits)

invisible(object)
}