buildS <- function(alpha,beta,kn,tmax) {
    K <- length(kn)
# Check that the resulting S function is continuous in x for all
# t in c(0,tmax).
    if(K > 1) {
	OK <- rep(TRUE,K-1)
        for(k in 1:(K-1)) {
	    foo <- function(t) {
		(alpha[[k]](t) - alpha[[k+1]](t) +
                (beta[[k]](t) -  beta[[k+1]](t))*kn[k])^2
            }
            chk <- optimize(foo,c(0,tmax),maximum=TRUE,tol=1e-12)
            if(chk$objective > .Machine$double.eps) OK[k] <- FALSE
        }
	if(any(!OK)) {
		NOK <- !OK
		nbad <- sum(NOK)
		kbad <- paste((1:(K-1))[NOK],collapse=" ")
		stop(paste("S is discontinuous in x at",
                           ngettext(nbad,"knot","each of the knots"),
			   kbad,"for some t.\n"))
	}
    }

# Check that S(x,t) is *non-increasing* in x for all t,
    OK <- rep(TRUE,K)
    for(k in 1:K) {
        chk <- optimize(beta[[k]],c(0,tmax),maximum=TRUE,tol=1e-12)
        if(chk$objective > .Machine$double.eps) OK[k] <- FALSE
    }
    if(any(!OK)) {
        NOK <- !OK
        nbad <- sum(NOK)
        kbad <- paste((1:(K-1))[NOK],collapse=" ")
        stop(paste("The \"beta\"",ngettext(nbad,"function","functions"),
            "numbered",kbad,ngettext(nbad,"is","are"),
            "positive for some values of \"t\" \n",
            "whence S(x,t) would fail to be non-increasing in \"x\".\n"))
    }

# Check that S(x,t) is *non-negative* for all t.
    foo <- function(t) {
        alpha[[K]](t) + beta[[K]](t)*kn[K]
    }
    chk <- optimize(foo,c(0,tmax),tol=1e-12)
    if(chk$objective < -(.Machine$double.eps))
        stop(paste("S(x,t) fails to be always non-negative and hence\n",
                   "does not define a probability.\n"))

# Check that S(0,t) = 1 for all t,
    foo <- function(t) {
        (alpha[[1]](t) - 1)^2
    }
    chk <- optimize(foo,c(0,tmax),maximum=TRUE,tol=1e-12)
    if(chk$objective > .Machine$double.eps)
        stop("S(0,t) is not equal to 1 for some \"t\".\n")

# OK, we're good to go.
    S <- function(x,t) {
	K <- length(kn)
	eps <- sqrt(.Machine$double.eps)
        if(any(x < -eps | x > kn[K]+eps))
		stop("At least one price value out of range.\n")
        #if(any(t < -eps | t > tmax+eps))
	#	stop("At least one time value out of range.\n")
        k <- cut(x,c(0,kn),include.lowest=TRUE,labels=1:K)
        k <- as.numeric(levels(k)[k])
        a <- lapply(k,function(i,alpha,t){alpha[[i]](t)},alpha=alpha,t=t)
        b <- lapply(k,function(i,beta,t){beta[[i]](t)},beta=beta,t=t)
        a <- matrix(unlist(a),nrow=length(x),byrow=TRUE)
        b <- matrix(unlist(b),nrow=length(x),byrow=TRUE)
        m <- a + b*x
        if(any(dim(m)==1)) as.vector(m) else m
    }
    environment(S) <- new.env()
    assign("alpha",alpha,envir=environment(S))
    assign("beta",beta,envir=environment(S))
    assign("kn",kn,envir=environment(S))
    assign("tmax",tmax,envir=environment(S))
    attr(S,"tmax") <- tmax
    attr(S,"funtype") <- "pwl"
# Why? Can't remember ....
    parent.env(environment(S)) <- globalenv()
    return(S)
}
