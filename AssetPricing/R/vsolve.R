vsolve <- function(S,lambda,gprob=NULL,tmax=NULL,x,nout=300,
                   alpha=NULL,salval=0,method="lsoda",verbInt=0) {
#
# Function vsolve to solve numerically the system of d.e.'s for the
# value v_q(t) of a stock of q items at time t given a ***general***
# (not necessarily optimal; either discrete or continuous) pricing
# policy ``x'', specified as an object of class ``flap''.  Uses the
# method of Runge-Kutta.  (Note that time is ***residual*** time,
# decreasing toward the ``departure time'' of 0.)  But we solve
# ``forward in time'', starting from 0.
#

# Check on the given pricing policy.
    if(!inherits(x,"flap")) 
        stop("Argument x must be an object of class \"flap\".\n")
    qmax <- attr(x,"qmax")
    type <- if(inherits(x,"di.flap")) "dip" else "sip"

# If S is a piecewise linear price sensitivity function, obtain
# tmax, if NULL, as the corresponding attribute of that function
# otherwise make sure that it is smaller than that attribute.
    funtype <- attr(S,"funtype")
    if(!is.null(funtype) && funtype=="pwl") {
        if(is.null(tmax)) {
            tmax <- attr(S,"tmax")
        } else if(tmax > attr(S,"tmax")) {
                stop(paste("Argument \"tmax\" is greater than the \"tmax\" attribute\n",
                           "of the pwl price sensitivity function specified as\n",
                           "argument \"S\".\n"))
        }
    }

# Check that the functions in x were built on an interval at least
# as large as [0,tmax].  (Since they are step functions, their
# values are not really meaningful for arguments outside of the
# tlim attribute of x.)
    if(is.null(tmax)) {
        tmax <- attr(x,"tlim")[2]
    } else if(tmax > attr(x,"tlim")[2])
        stop("Argument \"x\" has attribute tlim[2] < tmax.\n")

# Make sure the group size probabilities are OK.
    if(is.null(gprob)) gprob <- 1
    gpr <- if(is.function(gprob)) gprob(1:qmax) else gprob
    if(!is.numeric(gpr)) stop("Group size probabilities are not numeric.\n")
    if(any(gpr<0) | sum(gpr) > 1)
        stop("Group size probabilities are not probabilities!\n")

# Find the value of jmax as determined by gprob.  If there is
# doubly indexed pricing, check that this value of jmax is less than 
# or equal to the jmax attribute of x.
    jmax <- max(which(gpr > sqrt(.Machine$double.eps)))
    jmax <- min(jmax,qmax)
    if(type=="dip" & jmax > attr(x,"jmax"))
        stop(paste("Pricing is group-size dependent and the maximum\n",
                   "customer group size is larger than the\n",
                   "\"jmax\" attribute of \"x\".\n"))

# Clip the group size probability vector to be of length jmax.
    gpr  <- gpr[1:jmax]

# Check up on alpha.
    if(is.null(alpha)) {
        if(jmax > 1) {
            stop(paste("Argument \"alpha\" must be specified if there is\n",
                       "non-zero probability of a group of size greater than 1.\n"))
        }
        alpha <- 1
    }

    N  <- length(x)
    nc <- if(type=="dip") jmax*(qmax - (jmax-1)/2) else qmax
    if(N != nc)
        stop("Length of \"x\" incommensurate with \"qmax\" and \"jmax\".\n")

# If S is a list make sure that it is of the right length.  Then
# check that all entries are either expressions (smooth functions
# for continuous prices) or functions (for discrete prices).  If S
# consists of only one function or expression, raise it to the power
# "n" making "n" a function argument.  In the expressions case
# differentiate the expression(s).
    if(is.list(S)) {
        if(length(S) != jmax)
            stop(paste("Length of \"S\" as a list must equal ",
                               "\"jmax\" = ",jmax,".\n",sep=""))
        if(all(sapply(S,is.expression))) {
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
        } else if(!all(sapply(S,is.function))) {
            stop("At least one entry of \"S\" is neither an expression nor a function.\n")
        }
    } else if(is.expression(S)) {
# Note that here S is an expression but even so the (counter-intuitive)
# syntax ``S[[1]]'' (a) makes ``sense'' and (b) is needed.
        pars <- attr(S,"parvec")
        S <- substitute(a^b,list(a=S[[1]],b=quote(n)))
        dS <- deriv3(S,c("x","t"),function.arg=c("x","t","n"))
        environment(dS) <- new.env()
        for(nm in names(pars)) {
            assign(nm,pars[nm],envir=environment(dS))
        }
    } else if(is.function(S)) {
        dS <- with(list(S=S),function(x,t,n) {S(x,t)^n})
    } else {
        stop(paste("Argument \"S\" must be either an expression or a list of such,\n",
                   "or a function or a list of such.\n"))
    }
    if(is.numeric(lambda)) {
        if(length(lambda) != 1 || lambda <=0)
            stop("When \"lambda\" is numeric it must be a positive scalar.\n")
        lambda <- with(list(lambda=lambda),function(t){rep(lambda,length(t))})
    }


# Renew the environment of cev() and scrF() to prevent old remnants hanging
# around and thereby instigating spurious results.
environment(cev) <- new.env()
environment(scrF) <- new.env()

#
    assign("dS",dS,envir=environment(cev))
    assign("gpr",gpr,envir=environment(cev))
    assign("alpha",alpha,envir=environment(cev))
#
    assign("x",x,envir=environment(scrF))
    assign("lambda",lambda,envir=environment(scrF))
    assign("type",type,envir=environment(scrF))

# Do some setting up/initializing:
    tvec  <- seq(0,tmax,length=nout)
    v     <- (1:qmax)*salval
    info  <- new.env()
    info$st.first <- info$st.last <- Sys.time()

# Solve the differential equation.
odeRslt <- ode(v,tvec,scrF,parms=NULL,method=method,
               verbInt=verbInt,tmax=tmax,info=info)

# The functions in the object x are ``non-parametric'' functions;
# they could have been defined over a larger interval than [0,tmax]
# which is currently being investigated.  If so we need to reset
# the tlim and ylim values based on the current value of tmax.
#
    if(tmax < attr(x,"tlim")[2]) {
        ttt <- seq(0,tmax,length=nout)
        foo <- function(f,tt) {
                return(range(f(tt)))
        }
        tstor <- lapply(x,foo,tt=ttt)
        attr(x,'tlim') <- c(0,tmax)
        attr(x,'ylim') <- range(unlist(tstor))
    }
    comment(x) <- c(comment(x),"Prices not necessarily optimal.")
    putAway(odeRslt,type,jmax,qmax,soltype="vsolve",x=x,prices=NULL)
}
