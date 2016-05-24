scrF <- function(tt,v,parms,verbInt,tmax,info,...) {
#
# The argument "parms" is a dummy, required by ode().  The "..."
# argument is not used.  The objects lambda, type, gpr and stabilize are
# (always) assigned in the environment of scrF.  When scrF is called
# by xsolve.disc(), the vector "x" of *possible* prices is assigned
# in the environment of scrF.  When scrF is called by vsolve()
# the pricing policy "x" is assigned in the environment of scrF.
# When scrF is called by xsolve.pwl() the lists "alpha" and "beta" of
# coefficient functions and knot vector "kn" for the piecewise linear
# representation of S(x,t) are assigned in the environment of scrF.
#
# The value of this function is a list whose entries are:
#
# (a) When this function is called by vsolve (whence the prices
# are given): just vdot = ``script F''(v,t), the (unstated :-( ))
# vectorized version of equation (2) of Banerjee and Turner (2012)
# (b) When this function is called by xsolve.disc (discrete prices)
# or by xsolve.pwl() (piecewise linear price sensitivity function):
# vdot, x = the vector of optimal prices and vdlit=vdot (again).
# The latter two list entries are the "global values that are required
# at each point".  See the description of the argument "func" to ode()
# in the help for ode().
E <- parent.env(environment())
if(is.null(E$x)) {
# Here scrF is being called by xsolve.pwl() and the price elasticity
# functions are piecewise linear (rather than discrete).  We want
# to maximize over the possible price values, but first we need to
# determine what the possible prices are.  Note that "x" is a list
# here, with the q-th entry of the list being the vector of possible
# optimal prices for stock size "q".
	qmax <- length(v)
	jmax <- length(gpr)
	x    <- vector("list",qmax)
	for(q in 1:qmax) {
		jtop <- min(q,jmax)
		Kpa  <- gpr[1:jtop]
		if(jtop < jmax)
			Kpa[jtop] <- Kpa[jtop] +
                                     environment(cev)$alpha*(1-sum(Kpa))
		x[[q]] <- getPossPrices(v[1:q],tt,E$alpha,E$beta,
                                        E$kn,Kpa,type=type)
	}
} else if(inherits(x,"flap")) {
# Here scrF is being called by vsolve() --- pricing policy is given
# so x supplies the actual prices. No optimization to be done,
# so just return the derivatives of the expected values.
	xx <- sapply(x,function(f,t){f(t)},t=tt)
        if(verbInt > 0) progRep(info,verbInt,tt,tmax)
	return(list(vdot=lambda(tt)*(-v + cev(xx,tt,v,type))))
}

# At this point either scrF was called by xsolve.pwl() and the
# vector of possible prices x has been constructed, or scrF is being
# called by xsolve.disc() and the set of possible (discrete) prices
# was provided in the vector x which was assigned in the environment
# of this function.  In either case the vector of possible prices
# is now available and we can maximize over this vector.
R <- try(cev(x,tt,v,type,maximize=TRUE))

# The object R consists of the expected values at the optimum prices
# (both discrete and pwl settings).  It has an attrbute consisting
# of the actual optimum prices; Return the vdot values *at* the
# optimum prices, and the optimum prices.
vdot <- lambda(tt)*(R - v)
xopt <- attr(R,"xopt")
if(stabilize) assign("xback",xopt,envir=environment(cev))
if(verbInt > 0) progRep(info,verbInt,tt,tmax)
list(vdot=vdot,x=xopt,vdlit=vdot) # "lit" for literal.
}
