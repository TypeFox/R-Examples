
print.AccSampPlan = function(x,...){
print.default(x,...)
}

plot.AccSampPlan = function(x, y=NULL, ...){
    one.fig <- prod(par("mfcol")) == 1
plot(x$p, x$OC, type="l", 
ylab="Probability of Acceptance", 
xlab="Fraction Nonconforming p")

if(!is.null(x$ASN)){
if(one.fig){dev.new()}
plot(x$p, x$ASN, type="l", 
ylab="Uncurtailed average sample size", 
xlab="Fraction Nonconforming p")
}
if(!is.null(x$n)){
if(one.fig){dev.new()}
plot(x$p, x$n, type="l", 
ylab="Uncurtailed sample size", 
xlab="Fraction Nonconforming p")
}

if(one.fig){dev.new()}
plot(x$p, x$AOQ, type="l", 
ylab="AOQ", 
xlab="Fraction Nonconforming p")
title(paste("AOQL = ", formatC(max(x$AOQ))))

if(one.fig){dev.new()}
plot(x$p, x$ATI, type="l", 
ylab="ATI", 
xlab="Fraction Nonconforming p")
}
