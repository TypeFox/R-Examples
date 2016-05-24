`print.binomRRci` <-
function(x, digits=4, ...)
{

METHOD<-x$method
# A table of confidence intervals


dist<-attr(x$quantile, which="dist")

if(dist=="MVN"){
cat("Simultaneous", round(x$conf.level*100,3),"percent-confidence intervals","\n",
"for the ratio of proportions (RR), based on a crude normal approximation", "\n")
}
else{
cat("Local", round(x$conf.level*100,3),"percent-confidence intervals","\n",
"for the ratio of proportions (RR), based on a crude normal approximation", "\n")
}

conf.int <- cbind(x$estimate, x$conf.int)

colnames(conf.int)[1] <- "estimate"

print(round(conf.int, digits=digits))

cat(" ","\n")
cat("where proportions are the probabilities to observe", x$success, "\n")
cat(" ","\n")

invisible(x)
}

