MC.ZT.statistics <-
function(Nrs, MC, fit, type="ha", siglev=0.05) {
if(missing(Nrs) || missing(MC) || missing(fit))
stop("Nrs, MC and/or fit missing.")

for(n in Nrs){
if(any(n!=n[1])){
warning("Unequal number of reads across samples.")
break
}
}

ZTt <- matrix(0, MC, 2)
for(i in 1:MC)
ZTt[i,] <- ZT.statistics.Hnull.Ha(Nrs, fit, type)

qAlpha <- qchisq(p=(1-siglev), df=length(fit$pi)-1, ncp=0, lower.tail=TRUE)

z <- ZTt[,1]
t <- ZTt[,2]
zpval <- sum((z[z != "NaN"] > qAlpha)/(sum(z != "NaN")))
tpval <- sum((t[t != "NaN"] > qAlpha)/(sum(t != "NaN")))

return(cbind(zpval, tpval))
}
