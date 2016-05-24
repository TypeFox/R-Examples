corr.nn4pn<-function(lam, PN.cor){
x=rnorm(100000,0,1)
y=rnorm(100000,0,1)
xpois=qpois(pnorm(x),lam)
c=cor(xpois[order(xpois)],y[order(y)])/cor(x[order(x)],y[order(y)]) 
corrected=PN.cor/c 
return(corrected)
}
