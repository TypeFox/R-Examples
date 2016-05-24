MC.Xdc.statistics <-
function(Nrs, MC, alphap, type="ha", siglev=0.05, est="mom") {
if(missing(Nrs) || missing(alphap) || missing(MC))
stop("Nrs, alphap and/or MC missing.")

if(tolower(type) == "hnull"){
K <- length(alphap)
}else if(tolower(type) == "ha"){
K <- ncol(alphap)
}else{
stop(sprintf("Type '%s' not found. Type must be 'ha' for power or 'hnull' for size.\n", as.character(type)))
}

for(n in Nrs){
if(any(n!=n[1])){
warning("Unequal number of reads across samples.")
break
}
}

numGroups <- length(Nrs)
Xdc <- rep(0, MC)
for(i in 1:MC)
Xdc[i] <- Xdc.statistics.Hnull.Ha(alphap, Nrs, numGroups, type, est)

q.alpha <- qchisq(p=(1-siglev), df=(numGroups-1)*K, ncp=0, lower.tail=TRUE)
Xdc_pval <- sum((Xdc[Xdc != "NaN"] > q.alpha)/(sum(Xdc != "NaN")))

return(Xdc_pval)
}
