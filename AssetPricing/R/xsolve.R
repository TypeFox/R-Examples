xsolve <- function(S,lambda,gprob=1,tmax=NULL,qmax,prices=NULL,nout=300,
                   type="sip",alpha=NULL,salval=0,epsilon=NULL,
                   method="lsoda",verbInt=0) {
#
# Dispatch the problem to one of xsolve.cont(), xsolve.disc(),
# xsolve.pwl(), depending on the nature of the price sensitivity
# function "S".
#
    soltype <- findSolType(S,prices)

    if(is.numeric(lambda)) {
        if(length(lambda) != 1 || lambda <=0)
            stop("When \"lambda\" is numeric it must be a positive scalar.\n")
        lambda <- with(list(lambda=lambda),function(t){rep(lambda,length(t))})
    }

    if(is.null(epsilon)) {
        epsilon <- switch(EXPR=soltype,cont=NULL,
                                       disc = (.Machine$double.eps)^0.25,
                                       pwl  = (.Machine$double.eps)^0.5)
    }

    switch(EXPR = soltype,
	cont = xsolve.cont(S,lambda,gprob,tmax,qmax,nout,type,
                           alpha,salval,method=method,verbInt=verbInt),
	disc = xsolve.disc(S,lambda,gprob,tmax,qmax,prices,nout,type,
                           alpha,salval,epsilon,method=method,verbInt=verbInt),
	pwl  = xsolve.pwl(S,lambda,gprob,tmax,qmax,nout,type,
                           alpha,salval,epsilon,method=method,verbInt=verbInt)
    )
}
