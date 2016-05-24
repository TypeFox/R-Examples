corr.nn4pp<-function(lambda1, lambda2, PP.cor){
u = runif(100000, 0, 1)
lambda=c(lambda1,lambda2)
maxcor=cor(qpois(u, lambda1), qpois(u, lambda2))
mincor=cor(qpois(u, lambda1), qpois(1-u, lambda2))
a=-maxcor*mincor/(maxcor+mincor)
b=log((maxcor+a)/a, exp(1))
c=-a
corrected=log((PP.cor+a)/a, exp(1))/b
corrected=ifelse ((corrected>1 | corrected<(-1)),NA, corrected)
return(corrected)
}
