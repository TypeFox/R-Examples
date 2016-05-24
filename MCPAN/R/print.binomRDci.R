`print.binomRDci` <-
function(x, digits=4, ...)
{

if(x$method=="Wald")
 {methvar<-"Wald"}
if(x$method=="ADD1")
 {methvar<-"Add-1"}
if(x$method=="ADD2")
 {methvar<-"Add-2"}


# A table of confidence intervals


dist<-attr(x$quantile, which="dist")

if(dist=="MVN"){
cat("Simultaneous", round(x$conf.level*100,3),"percent",methvar,"-confidence intervals","\n",
"for the difference of proportions (RD)", "\n")
}
else{
cat("Local", round(x$conf.level*100,3),"percent",methvar,"-confidence intervals","\n",
"for the difference of proportions (RD)", "\n")
}

conf.int <- cbind(x$estimate, x$conf.int)

colnames(conf.int)[1] <- "estimate"

print(round(conf.int, digits=digits))

cat(" ","\n")
cat("where proportions are the probabilities to observe", x$success, "\n")
cat(" ","\n")

invisible(x)
}

