penplot <- function(base.parms, vbeta, variation="none", base.dist="Weibull", frailty.dist="gamma", depend=1, agemin=20){

x <- agemin:80
t0 <- x-agemin

if(variation=="secondgene") xbeta <- c( vbeta %*% t(expand.grid(c(0,1),c(0,1), c(0,1))) )
else xbeta <- c( vbeta %*% t(expand.grid(c(0,1),c(0,1)) ) )


s1 <- cumhaz(base.dist, t0, base.parms)*exp(xbeta[1]) # female non-carr
s2 <- cumhaz(base.dist, t0, base.parms)*exp(xbeta[2]) # male non-carr
s3 <- cumhaz(base.dist, t0, base.parms)*exp(xbeta[3]) # female carr
s4 <- cumhaz(base.dist, t0, base.parms)*exp(xbeta[4]) # male carr
s <- list(s1,s2,s3,s4)
if(variation=="secondgene"){
	# second gene carriers
s[[5]] <- cumhaz(base.dist, t0, base.parms)*exp(xbeta[5]) # female non-carr
s[[6]] <- cumhaz(base.dist, t0, base.parms)*exp(xbeta[6]) # male non-carr
s[[7]] <- cumhaz(base.dist, t0, base.parms)*exp(xbeta[7]) # female carr
s[[8]] <- cumhaz(base.dist, t0, base.parms)*exp(xbeta[8]) # male carr	
}

if(variation=="none" | variation=="secondgene") pen <- lapply(s, function(s) 1-exp(-s))
else if(variation=="frailty") pen <- lapply(s, function(dist,s,k) 1-laplace(dist,s,k), dist=frailty.dist, k=depend)
else stop("Unrecognized variation")
 
if(variation=="secondgene"){
par(mfrow=c(1,2))
	plot(x,pen[[5]], ylab="Penetrance", xlab="age at onset", ylim=c(0,1), type="l", lty=2, col="red", main=paste("Penetrances for second gene carriers \n", base.dist, "baseline"))
  lines(x, pen[[6]], lty=2, col="blue")
  lines(x, pen[[7]], lty=1, col="red")
  lines(x, pen[[8]], lty=1, col="blue")
  legend("topleft", c("male carrier", "female carrier", "male noncarrier", "female noncarrier"), bty="n", lty=c(1,1,2,2), col=c("blue","red","blue","red"))

	plot(x,pen[[1]], ylab="Penetrance", xlab="age at onset", ylim=c(0,1), type="l", lty=2, col="red", main=paste("Penetrances for second gene non-carriers \n", base.dist, "baseline"))
  lines(x, pen[[2]], lty=2, col="blue")
  lines(x, pen[[3]], lty=1, col="red")
  lines(x, pen[[4]], lty=1, col="blue")
  legend("topleft", c("male carrier", "female carrier", "male noncarrier", "female noncarrier"), bty="n", lty=c(1,1,2,2), col=c("blue","red","blue","red"))

}
else if(variation=="frailty"){
  plot(x,pen[[1]], ylab="Penetrance", xlab="age at onset", ylim=c(0,1), type="l", lty=2, col="red", main=paste("Penetrance curves \n", base.dist, "baseline and ", frailty.dist,"frailty"))
  lines(x, pen[[2]], lty=2, col="blue")
  lines(x, pen[[3]], lty=1, col="red")
  lines(x, pen[[4]], lty=1, col="blue")
  legend("topleft", c("male carrier", "female carrier", "male noncarrier", "female noncarrier"), bty="n", lty=c(1,1,2,2), col=c("blue","red","blue","red"))
}
else {
  plot(x,pen[[1]], ylab="Penetrance", xlab="age at onset", ylim=c(0,1), type="l", lty=2, col="red", main=paste("Penetrance curves \n", base.dist, "baseline"))
  lines(x, pen[[2]], lty=2, col="blue")
  lines(x, pen[[3]], lty=1, col="red")
  lines(x, pen[[4]], lty=1, col="blue")
  legend("topleft", c("male carrier", "female carrier", "male noncarrier", "female noncarrier"), bty="n", lty=c(1,1,2,2), col=c("blue","red","blue","red"))
}

}