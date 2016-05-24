scrG <- function(tt,x,parms,verbInt,tmax,info,...) {
# Note that the objects dS, gpr, alpha, and type are
# assigned in the environment of scrG.
#
# This is the value of xdot = {script G}(x,t) (equation (6) of paper.

vvdot <- vupdate(x,tt,type)
v     <- vvdot$v
vdot  <- vvdot$vdot
qmax  <- length(v)
if(type=="sip") {
	xdot <- numeric(qmax)
	for(q in 1:qmax) {
		vq      <- v[q]
		vdotq   <- vdot[q]
		jtop    <- min(q,length(gpr))
		jv      <- 1:jtop
		rv      <- (c(v[q:1],0)[-1])[1:jtop]
		rvdot   <- (c(vdot[q:1],0)[-1])[1:jtop]
		K       <- gpr[1:jtop]
		K[jtop]    <- K[jtop] + alpha*(1-sum(K))
		if(is.list(dS)) {
			S       <- lapply(dS[jv],function(f,xx,tt){f(xx,tt)},
                                          xx=x[q],tt=tt)
			dSdx    <- sapply(S,function(f){
				          attr(f,"gradient")[,"x"]
		 	                     })
			dSdt    <- sapply(S,function(f){
				          attr(f,"gradient")[,"t"]
		 	                     })
			d2Sdx2  <- sapply(S,function(f){
				          attr(f,"hessian")[,"x","x"]
	                                     })
			d2Sdxdt <- sapply(S,function(f){
				          attr(f,"hessian")[,"x","t"]
	                                     })
			S       <- unlist(S)
		} else {
			S       <- dS(x[q], tt, jv)
			dSdx    <- attr(S, "gradient")[, "x"]
			dSdt    <- attr(S, "gradient")[, "t"]
			d2Sdx2  <- attr(S, "hessian")[, "x", "x"]
			d2Sdxdt <- attr(S, "hessian")[, "x", "t"]
		}
		num     <- -(x[q]*sum(jv*d2Sdxdt*K) + sum((rv-vq)*d2Sdxdt*K) +
	                     sum((rvdot-vdotq)*dSdx*K) + sum(jv*dSdt*K))
		den     <- x[q]*sum(jv*d2Sdx2*K) + sum((rv-vq)*d2Sdx2*K) +
	                   2*sum(jv*dSdx*K)
		xdot[q] <- num/den
	}
} else if(type=="dip") {
	N    <- length(x)
	xdot <- numeric(N)
	jmax <- length(gpr)
	for(q in 1:qmax) {
		jtop <- min(q,jmax)
		for(j in 1:jtop) {
			i       <- qj2i(q,j,qmax)
			S       <- if(is.list(dS)) {
					dS[[j]](x[i],tt)
				   }  else dS(x[i],tt,j)
			dSdx    <- attr(S,"gradient")[,"x"]
			dSdt    <- attr(S,"gradient")[,"t"]
			d2Sdx2  <- attr(S,"hessian")[,"x","x"]
			d2Sdxdt <- attr(S,"hessian")[,"x","t"]
			vdqmj   <- if(j<q) vdot[q-j] else 0
			num     <- j*(S*d2Sdxdt - dSdx*dSdt) - 
		       	           (vdqmj - vdot[q])*dSdx^2
			den     <- j*(2*dSdx^2 - S*d2Sdx2)
			xdot[i] <- num/den
		}
	}
} else stop(paste("Type",type,"unrecognized.\n"))
if(verbInt > 0) progRep(info,verbInt,tt,tmax)
list(xdot=xdot,v=v,vdot=vdot)
}
