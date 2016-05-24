getPossPrices <- function(v,t,alpha,beta,kn,Kpa,type) {
K <- length(kn)
kn <- c(0,kn)
rslt <- vector("list",K)
for(k in 1:K) {
	a <- alpha[[k]](t)
	b <- beta[[k]](t)
	xlo <- kn[k]
	xhi <- kn[k+1]
	rslt[[k]] <- turnPts(a,b,v,Kpa,xlo,xhi,type)
}
if(type=="sip") {
    newres <- sort(unique(c(unlist(rslt),kn)))
} else {
  q <- length(v)
  newres <- vector("list",q)
  for(j in 1:q) {
    newres[[j]] <- lapply(1:K,function(i,x,j){x[[i]][[j]]},x=rslt,j=j)
  }
  newres <- lapply(newres,function(x,knots){unique(c(unlist(x),knots))},
                                            knots=kn)
  }
newres
}
