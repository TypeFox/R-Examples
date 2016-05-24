xsolve.disc <- function(S,lambda,gprob,tmax,qmax,prices,nout,type,
                        alpha,salval,epsilon,method,verbInt) {
#
# Function xsolve.disc to solve numerically the system of d.e.'s
# for the value v_q(t) of a stock of q items at time t, using the
# method of Runge-Kutta, in the setting in which prices vary over a
# ***discrete*** set.  The optimal prices are then determined by
# maximizing over this discrete set.  (Note that time is thought
# of as ***decreasing*** toward the ``departure time'' of 0.)
# But we solve ``forward in time'', starting from 0.
#
# We are solving the vector system of differential equations
#
#	vdot(t) = {script F}(v,t)
#
# In what follows, ``script F'' is denoted by ``scrF()''.
# In addition to "v" and "t", scrF() has an additional auxilliary
# argument "op".  This is a logical scalar which defaults to FALSE;
# if "op" is TRUE then the value returned by scrF() has an attribute
# "xopt" giving the vector of values of the optimal prices (chosen
# from amongst the finite set of possible prices.
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

# If there is only one price sensitivity function, create a new
# function equal to the original function raised to the power "n", with
# "n" a function argument.

    if(is.list(S)) {
    	if(length(S) < jmax)
    		stop(paste("Length of \"S\" as a list must be at least ",
                               "\"jmax\" = ",jmax,".\n",sep=""))
    	if(!all(sapply(S,is.function)))
    		stop("At least one entry of \"S\" is NOT a function.\n")
    	} else if(is.function(S)) {
    	oldS <- S
    	S    <- with(list(oldS=S),function(x,t,n) {oldS(x,t)^n})
    } else {
    	stop("Argument \"S\" must be either a function or a list of functions.\n")
    }

# Renew the environments of the functions into which objects
# are assigned to prevent old remnants hanging around and thereby
# instigating spurious results.
environment(scrF) <- new.env()
environment(cev) <- new.env()

#
    assign("dS",S,envir=environment(cev))  # "dS" to make notation compatible with
                                           # the "smooth" case.
    assign("gpr",gpr,envir=environment(cev))
    assign("alpha",alpha,envir=environment(cev))
    assign("epsilon",epsilon,envir=environment(cev))
#
    assign("x",prices,envir=environment(scrF))
    assign("lambda",lambda,envir=environment(scrF))
    assign("type",type,envir=environment(scrF))
    assign("gpr",gpr,envir=environment(scrF))
    assign("stabilize",epsilon>0,envir=environment(scrF))

# Do some setting up/initializing.
    tvec  <- seq(0,tmax,length=nout)
    v     <- (1:qmax)*salval
    info  <- new.env()
    info$st.first <- info$st.last <- Sys.time()

# Solve the differental equation:
    odeRslt <- ode(v,tvec,scrF,parms=NULL,method=method,verbInt=verbInt,
                   tmax=tmax,info=info)
    putAway(odeRslt,type,jmax,qmax,soltype="disc",x=NULL,prices=prices)
}
