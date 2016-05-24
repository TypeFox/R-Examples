MC.Xmcupo.statistics <-
function(Nrs, MC, pi0, group.pi, group.theta, type="ha", siglev=0.05) {
if(missing(Nrs) || missing(MC) || missing(group.theta))
stop("Nrs, MC, and/or group.theta missing.")
if(missing(group.pi) && tolower(type) == "ha")
stop("group.pi missing.")
if(missing(pi0) && tolower(type) == "hnull")
stop("pi0 missing.")

for(n in Nrs){
if(any(n!=n[1])){
warning("Unequal number of reads across samples.")
break
}
}

numGroups <- length(group.theta)
group.parameter <- vector("list", numGroups)

if(tolower(type) == "hnull"){
K <- length(pi0)
for(i in 1:numGroups)
group.parameter[[i]] <- c(pi0, group.theta[i], Nrs[[i]])
}else if(tolower(type) == "ha"){
K <- ncol(group.pi)
for(i in 1:numGroups) 
group.parameter[[i]] <- c(group.pi[i,], group.theta[i], Nrs[[i]])
}else{
stop(sprintf("Type '%s' not found. Type must be 'ha' for power or 'hnull' for size.\n", as.character(type)))
}

Xmcupo <- rep(0, MC)
for(i in 1:MC)
Xmcupo[i] <- Xmcupo.statistics.Hnull.Ha(K, group.parameter)

q.alpha <- qchisq(p=(1-siglev), df=length(group.theta)*(K-1), ncp=0, lower.tail=TRUE)
Xmcupo_pval <- sum((Xmcupo[Xmcupo != "NaN"] > q.alpha)/(sum(Xmcupo != "NaN")))

return(Xmcupo_pval)
}
