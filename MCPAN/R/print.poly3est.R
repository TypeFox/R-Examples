'print.poly3est' <-
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

sample.estimate<-rbind(x$Y, x$n, x$nadj, x$estimate)
rownames(sample.estimate)<-c("x","n", "adjusted n", "adjusted estimate")

cat("Raw and poly-", x$k,"-adjusted sample estimates:","\n", sep="")
print(round(sample.estimate, digits=digits))
cat(" ","\n")

invisible(x)

}

