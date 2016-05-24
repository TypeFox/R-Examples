vupdate <- function(x,tt,type) {
# Note that objects lambda, dS, gpr and alpha are assigned in the
# environment of vupdate.
N <- length(x)
if(type=="sip") {
	qmax <- N
	v    <- numeric(qmax+1)
	for(q in 1:qmax) {
		jtop <- min(q,length(gpr))
		rv   <- (v[q:1])[1:jtop]
		K    <- gpr[1:jtop]
		K[jtop] <- K[jtop] + alpha*(1-sum(K))
		if(is.list(dS)) {
			S      <- lapply(dS[1:jtop],function(f,x,t){
					f(x,t)},x=x[q],t=tt)
			dSdx   <- sapply(S,function(f){
				  attr(f,"gradient")[,"x"]
			  	  })
			S      <- unlist(S)
		} else {
			S      <- dS(x[q], tt, 1:jtop)
			dSdx   <- attr(S, "gradient")[, "x"]
		}
		num <- (x[q]*sum((1:jtop)*dSdx*K) + sum(rv*dSdx*K) +
                         sum((1:jtop)*S*K))
		den <- sum(dSdx*K)
		v[q+1] <- num/den
	}
	v <- v[-1]
} else if(type=="dip") {
	jmax <- length(gpr)
	qmax <- N/jmax + (jmax-1)/2
	x1 <- x[1:qmax]
	S    <- if(is.list(dS)) dS[1](x1,tt) else dS(x1,tt,1)
	dSdx <- attr(S,"gradient")[,"x"]
	v <- cumsum(x1+S/dSdx)
} else stop(paste("Type",type,"unrecognized.\n"))

vdot <- lambda(tt)*(-v + cev(x,tt,v,type))

list(v=v,vdot=vdot)
}
