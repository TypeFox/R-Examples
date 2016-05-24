plot.bayes2=function(x,marginal=0,...)
if(marginal==0)image(as.numeric(dimnames(x$prob)[[1]]),
as.numeric(dimnames(x$prob)[[2]]),x$prob,
col=gray(1-(0:32)/32),...) else
if(marginal==1) barplot(apply(x$prob,1,sum),...) else
barplot(apply(x$prob,2,sum),...)