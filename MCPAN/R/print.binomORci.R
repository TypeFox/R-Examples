`print.binomORci` <-
function(x,...)
{

# A table of confidence intervals

aargs<-list(...)

if(is.null(aargs$digits))
 {digits<-4}
else
 {digits<-aargs$digits}

dist<-attr(x$quantile, which="dist")

if(dist=="MVN"){
cat("Simultaneous", round(x$conf.level*100,3),"percent-confidence intervals","\n",
"for the odds ratio (OR)", "\n")
}
else{
cat("Local", round(x$conf.level*100,3),"percent-confidence intervals","\n",
"for the odds ratio (OR)", "\n")
}

conf.int <- cbind(x$estimate, x$conf.int)

colnames(conf.int)[1] <- "estimate"

print(round(conf.int, digits=digits))

ORdef<-paste("p(",x$success,")/(1-p(",x$success,"))", sep="")

cat(" ","\n")
cat("where the odds is defined:", ORdef, "\n")
cat(" ","\n")

invisible(x)
}

