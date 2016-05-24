`summary.binomORci` <-
function(object,...)
{

aargs<-list(...)
if(is.null(aargs$digits))
 {digits<-4}
else
 {digits<-aargs$digits}

success<-object$success

cat("Summary statistics: \n")

summarystat<-rbind(object$x,
object$n,
round(object$x/object$n, digits=digits))

rownames(summarystat)<-c(paste("number of",success, collapse=" "),
"number of trials",
paste("estimated probability of",success, collapes=""))

colnames(summarystat)<-object$names

print(summarystat, digits=digits)


cat("\n Contrast matrix on the level of the logit link: \n")

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

