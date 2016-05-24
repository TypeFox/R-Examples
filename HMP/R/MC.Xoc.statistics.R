MC.Xoc.statistics <-
function(Nrs, MC, group.alphap, type="ha", siglev=0.05) {
if(missing(Nrs) || missing(group.alphap) || missing(MC))
stop("Nrs, group.alphap, and/or MC missing.")

for(n in Nrs){
if(any(n!=n[1])){
warning("Unequal number of reads across samples.")
break
}
}

numGroups <- length(Nrs)
Xoc <- rep(0, MC)
for(i in 1:MC)
Xoc[i] <- Xoc.statistics.Hnull.Ha(Nrs, group.alphap, numGroups, type)

q.alpha <- qchisq(p=(1-siglev), df=(numGroups-1), ncp=0, lower.tail=TRUE)
Xoc_pval <- sum((Xoc[Xoc != "NaN"] > q.alpha)/(sum(Xoc != "NaN")))

return(Xoc_pval)
}
