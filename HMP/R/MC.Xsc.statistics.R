MC.Xsc.statistics <-
function(Nrs, MC, fit, pi0=NULL, type="ha", siglev=0.05) {
if(missing(Nrs) || missing(MC) || missing(fit))
stop("Nrs, MC and/or fit missing.")

for(n in Nrs){
if(any(n!=n[1])){
warning("Unequal number of reads across samples.")
break
}
}

Xsc <- rep(0, MC)
for(i in 1:MC)
Xsc[i] <- Xsc.statistics.Hnull.Ha(Nrs, fit, type, pi0)

q.alpha <- qchisq(p=(1-siglev), df=length(fit$pi)-1, ncp=0, lower.tail=TRUE)
Xsc_pval <- sum((Xsc[Xsc != "NaN"] > q.alpha)/(sum(Xsc != "NaN")))

return(Xsc_pval)
}
