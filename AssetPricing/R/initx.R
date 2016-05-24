initx <- function (v,type) {
#
# Note that objects dS, gpr and alpha are assigned in the
# environment of initx.
#
# Set up objective function according to type:
    if(type=="sip") {
	sip <- function(x, vi, gpr) {
		q      <- length(vi)
		vq     <- vi[q]
		q      <- min(q,length(gpr))
		rv     <- c(vi[q:1],0)[-1]
		if(is.list(dS)) {
			S      <- lapply(dS[1:q],function(f,x){f(x,0)},x=x)
			dSdx   <- sapply(S,function(f){
					attr(f,"gradient")[,"x"]
				  })
			d2Sdx2 <- sapply(S,function(f){
					attr(f,"hessian")[,"x","x"]
		                  })
			S      <- unlist(S)
		} else {
			S      <- dS(x, 0, 1:q)
			dSdx   <- attr(S, "gradient")[, "x"]
			d2Sdx2 <- attr(S, "hessian")[, "x", "x"]
		}
		kappa  <- gpr[1:q]
		kappa[q] <- kappa[q] + alpha*(1-sum(kappa))
		fval     <- x * (sum((1:q) * dSdx * kappa)) +
                            sum((rv * dSdx + (1:q)*S)*kappa) - 
                            vq*sum(dSdx*kappa)
		jacobian <- x * (sum((1:q) * d2Sdx2*kappa)) +
                            sum((rv * d2Sdx2 + 2*(1:q)*dSdx)*kappa) - 
                            vq*sum(d2Sdx2*kappa)
		list(fval = fval, jacobian = jacobian)
    	}
    } else if(type=="dip") {
	dip <- function(x,vqmj,vq,j) {
		S <- if(is.list(dS)) dS[[j]](x,0) else dS(x,0,j)
		dSdx <- attr(S,"gradient")[,"x"]
		d2Sdx2 <- attr(S,"hessian")[,"x","x"]
		fval <- (j*x + vqmj - vq)*dSdx + j*S
		jacobian <- (x*j + vqmj + vq)*d2Sdx2 + 2*j*dSdx
		list(fval=fval,jacobian=jacobian)
	}
    } else {
	stop(paste("Type",type,"not recognized.\n"))
    }

qmax <- length(v)
jmax <- length(gpr)

# Solve according to type:
    if(type=="sip") {
	x <- numeric(qmax)
	for (i in 1:qmax) {
        	x[i] <- .newt(sip, v[1], vi = v[1:i], gpr=gpr)
	}
    } else {
	N <- jmax*(qmax - jmax/2 + 0.5)
	x <- numeric(N)
	for(q in 1:qmax) {
		vq   <- v[q]
		jtop <- min(q,jmax)
		for(j in 1:jtop) {
			vqmj <- if(j<q) v[q-j] else 0
			i <- qj2i(q,j,qmax)
			x[i] <- .newt(dip,v[1],vqmj=vqmj,vq=vq,j=j)
		}
	}
    }
    x
}
