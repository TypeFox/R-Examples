MC.Xmc.statistics <-
function(Nrs, MC, pi0, group.pi, group.theta, type="ha", siglev=0.05) {
if(missing(Nrs) || missing(MC) || missing(group.theta) || missing(pi0))
stop("Nrs, MC, pi0 and/or group.theta missing.")
if(missing(group.pi) && tolower(type) == "ha")
stop("group.pi missing.")

for(n in Nrs){
if(any(n!=n[1])){
warning("Unequal number of reads across samples.")
break
}
}

numGroups <- length(group.theta)
group.parameter <- vector("list", numGroups)
K <- length(pi0)

if(tolower(type) == "hnull"){
for (i in 1:numGroups)
group.parameter[[i]] <- c(pi0, group.theta[i], Nrs[[i]])
}else if(tolower(type) == "ha"){
for (i in 1:numGroups) 
group.parameter[[i]] <- c(group.pi[i,], group.theta[i], Nrs[[i]])
}else{
stop(sprintf("Type '%s' not found. Type must be 'ha' for power or 'hnull' for size.\n", as.character(type)))
}

Xmc <- rep(0, MC)
for(i in 1:MC)
Xmc[i] <- Xmc.statistics.Hnull.Ha(K, pi0, group.parameter)

q.alpha <- qchisq(p=(1-siglev), df=length(group.theta)*(K-1), ncp=0, lower.tail=TRUE)
Xmc_pval <- sum((Xmc[Xmc != "NaN"] > q.alpha)/(sum(Xmc != "NaN")))

return(Xmc_pval)
}
