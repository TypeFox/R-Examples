`print.poly3ci` <-
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

cat("Sample estimates, using poly-",x$k,"-adjustment","\n")
print(round(sample.estimate, digits=digits))
cat(" ","\n")

# The contrast matrix:

cat("Contrast matrix:","\n")
print(x$cmat)
cat(" ","\n")

# A table of confidence intervals


dist<-attr(x$quantile, which="dist")

if(dist=="MVN"){
cat("Simultaneous",x$conf.level*100,"percent confidence intervals using",methvar,"variance estimators:","\n")
}
else{
cat("Local",x$conf.level*100,"percent confidence intervals using",methvar,"variance estimators:","\n")
}

conf.int <- cbind(x$estimate, x$conf.int)

colnames(conf.int)[1] <- "estimate"

print(round(conf.int, digits=digits))

cat(" ","\n")

invisible(x)

}

