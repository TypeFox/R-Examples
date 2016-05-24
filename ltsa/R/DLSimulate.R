`DLSimulate` <-
function(n,r,useC=TRUE, rand.gen=rnorm, ...)
{
    nr <- length(r)
    a <- rand.gen(n, ...)
    if (min(n,nr) <= 0) stop("input error")
    if (nr < n)
	r <- c(r, rep(0, n-nr))
    r <- r[1:n]
    EPS <- .Machine$double.eps # 1+EPS==1, machine epsilon
    z<-numeric(n)
    if (useC){
 	out<-.C("DLSim", z=as.double(z), as.double(a), as.integer(n), as.double(r), 
             as.double(EPS), fault = as.integer(1), PACKAGE="ltsa" )
    	fault<-out$fault
    	if (fault == 1) 
            stop("error: sequence not p.d.")
    	z<-out$z
    }
    else {
    sigmasqk <- r[1]
    error <- a[1]*sqrt(sigmasqk)
    phi <- r[2]/r[1]
    sigmasqk <- r[1] * (1 - phi^2)
    error <- a[2]*sqrt(sigmasqk)
    z[2] <- error + phi * z[1]
    sigmasqkm1<-sigmasqk
    for(k in 2:(n - 1)) {
            if (sigmasqkm1 < EPS) stop("r is not a p.d. sequence")
            phikk <- (r[k + 1] - phi %*% rev(r[2:k]))/sigmasqkm1
            sigmasqk <- sigmasqkm1 * (1 - phikk^2)
            phinew <- phi - phikk * rev(phi)
            phi <- c(phinew, phikk)
            sigmasqkm1 <- sigmasqk
            error <- a[k + 1]*sqrt(sigmasqk)
            z[k + 1] <- error + crossprod(phi, rev(z[1:k]))
            }
    }
 z
 }
