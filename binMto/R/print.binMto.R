"print.binMto" <-
function(x, digits=4, ...)
{
pargs<-list(...)

met<-x$method
conf <- x$conf.level
alt <- x$alternative
adj <- x$adj

if(x$adj=="Dunnett"){adjI<-"Dunnett-adjusted"}
if(x$adj=="Bonf"){adjI<-"Bonferroni-adjusted"}
if(x$adj=="Unadj"){adjI<-"Unadjusted"}


cat(" ", "\n")
cat(round(conf, digits=digits)*100, "-% confidence intervals for risk difference", "\n")
cat(adjI, "\n")
cat("used quantile: ", round(x$quantile, digits=digits), "\n")
cat(" ","\n")
pargs$x<-round(x$conf.int, digits=digits)
do.call("print", pargs)
cat(" ", "\n")
invisible(x)
}

