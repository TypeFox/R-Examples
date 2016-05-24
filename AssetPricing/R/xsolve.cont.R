xsolve.cont <- function(S,lambda,gprob,tmax,qmax,nout,type,
                        alpha,salval,method,verbInt) {
#
# Function xsolve.cont to solve numerically the system of d.e.'s for
# the optimum price x_q(t) of an item of stock when there are q
# items remaining at time t.  Of necessity the price x_q(t) must be
# ***continuous*** for this function to make sense.  Note that time
# is residual time, i.e.  the time remaining until the expiry date,
# and hence decreases as the expiry date approaches. However we
# solve forward in time, starting from the expiry date, i.e. t=0.
#
# Note that
#           - S is an *expression* or *call*, or a *list* of such,
#             specifying the time-varying price sensitivity.
#           - lambda is the intensity of the Poisson arrival process
#           - gprob specifies the probability function for the size
#             of the arriving group of customers, either as a function
#             or as a vector of probabilities.
#           - tmax is the length of the time interval over which
#             we are solving the differential equation
#           - qmax is the initial number of items for sale, tmax
#             time units from the expiry date
#           - nout is the number of time points at which the values of
#             the solution are required. These points will be equispaced
#             on [0,tmax]
#           - type is the specification of the model type, either ``sip''
#             (singly indexed prices) or ``dip'' (doubly indexed prices).
#           - salval is the ``salvage value'' of an item of stock when
#             the expiry date or deadline is reached; used to determine
#             the initial conditions.
#
# Note that ``scrG'' means ``script G function''; the (vector) d.e.
# that we are solving is
#                        ``x.dot = {script G}(x,t)''
#

# Make sure the group size probabilities are OK.
    gpr <- if(is.function(gprob)) gprob(1:qmax) else gprob
    if(!is.numeric(gpr)) stop("Group size probabilities are not numeric.\n")
    if(any(gpr<0) | sum(gpr) > 1)
    	stop("Group size probabilities are not probabilities!\n")
    jmax <- max(which(gpr > sqrt(.Machine$double.eps)))
    jmax <- min(jmax,qmax)
    gpr  <- gpr[1:jmax]
    if(is.null(alpha)) {
       if(jmax > 1) stop(paste("When the maximum group size is great than 1,\n",
                               "\"alpha\" must be specified.\n"))
       alpha <- 1
    }

# If jmax = 1 we might as well set type equal to "sip" --- since
# indexing according to group size is "degenerate" in this case.
    if(jmax==1) type <- "sip"

# Make sure tmax is specified.
    if(is.null(tmax)) stop("Argument \"tmax\" was not specified.\n")

# Differentiate each price sensitivity function.  If there
# is only one such, raise it to the power "n" and differentiate
# the result, making "n" a function argument.
    if(is.list(S)) {
    	if(length(S) < jmax)
    		stop(paste("Length of \"S\" as a list must be at least ",
                               "\"jmax\" = ",jmax,".\n",sep=""))
    	dS <- list()
    	for(i in 1:jmax) {
    		xxx  <- deriv3(S[[i]],c("x","t"),function.arg=c("x","t"))
                environment(xxx) <- new.env()
    		pars <- attr(S[[i]],"parvec")
    		for(nm in names(pars)) {
    			assign(nm,pars[nm],envir=environment(xxx))
    		}
    		dS[[i]] <- xxx
    	}
    } else {
    	pars <- attr(S,"parvec")
    	S <- substitute(a^b,list(a=S[[1]],b=quote(n)))
    	dS <- deriv3(S,c("x","t"),function.arg=c("x","t","n"))
        environment(dS) <- new.env()
    	for(nm in names(pars)) {
    		assign(nm,pars[nm],envir=environment(dS))
    	}
    }

# Renew the environments of the functions into which objects
# are assigned to prevent old remnants hanging around and thereby
# instigating spurious results.
environment(vupdate) <- new.env()
environment(scrG) <- new.env()
environment(initx) <- new.env()
environment(cev) <- new.env()

#
    assign("dS",dS,envir=environment(vupdate))
    assign("dS",dS,envir=environment(scrG))
    assign("dS",dS,envir=environment(initx))
    assign("dS",dS,envir=environment(cev))
#
    assign("gpr",gpr,envir=environment(vupdate))
    assign("gpr",gpr,envir=environment(scrG))
    assign("gpr",gpr,envir=environment(initx))
    assign("gpr",gpr,envir=environment(cev))
#
    assign("alpha",alpha,envir=environment(vupdate))
    assign("alpha",alpha,envir=environment(scrG))
    assign("alpha",alpha,envir=environment(initx))
    assign("alpha",alpha,envir=environment(cev))
#
    assign("lambda",lambda,envir=environment(vupdate))
#
    assign("type",type,envir=environment(scrG))

# Do some setting up/initializing:
    tvec  <- seq(0,tmax,length=nout)
    v     <- (1:qmax)*salval
    x     <- initx(v,type)
    info  <- new.env()
    info$st.first <- info$st.last <- Sys.time()

# Solve the differential equation.
odeRslt <- ode(x,tvec,scrG,parms=NULL,method=method,verbInt=verbInt,
               tmax=tmax,info=info)
putAway(odeRslt,type,jmax,qmax,soltype="cont",x=NULL,prices=NULL)
}
