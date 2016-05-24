# depends on getpars and nobs
setMethod("logLik",signature(object="depmix"),
	#function(object,method="lystig") { 
	function(object,method=c("fb","lystig","classification"),na.allow=TRUE) { 
	    #4/5/2012: set to fb as this is now in C
	    method <- match.arg(method)
		if(method=="fb") ll <- fb(init=object@init,A=object@trDens,B=object@dens,ntimes=object@ntimes,homogeneous=object@homogeneous,na.allow=na.allow)$logLike
		if(method=="lystig") ll <- lystig(init=object@init,A=object@trDens,B=object@dens,ntimes=object@ntimes,homogeneous=object@homogeneous)$logLike
		if(method=="classification") {
		    ns <- nstates(object)
		    ntimes <- ntimes(object)
		    vstate <- viterbi(object)[,1]
		    B <- object@dens
		    if(na.allow) B[is.na(B)] <- 1
		    ll <- sum(log((apply(B,c(1,3),prod))[cbind(1:sum(ntimes),vstate)]))
		}
	    attr(ll, "df") <- freepars(object)
		attr(ll, "nobs") <- nobs(object)
		class(ll) <- "logLik"
		ll
	}
)

setMethod("logLik",signature(object="depmix.fitted.classLik"),
    function(object,method=c("classification","fb","lystig"),na.allow=TRUE) {
        method <- match.arg(method)
        callNextMethod(object=object,method=method,na.allow=na.allow)
     }
)

# depends on getpars and nobs
setMethod("logLik",signature(object="mix"),
	#function(object,method="lystig") { 
	function(object,method=c("fb","lystig","classification"),na.allow=TRUE) {
	    method <- match.arg(method)
		if(method=="fb") ll <- fb(init=object@init,A=matrix(0,1,1),B=object@dens,ntimes=object@ntimes,homogeneous=TRUE)$logLike
		if(method=="lystig") ll <- lystig(init=object@init,A=matrix(0,1,1),B=object@dens,ntimes=object@ntimes,homogeneous=TRUE)$logLike
		if(method=="classification") {
		    ntimes <- ntimes(object)
		    gamma <- fb(init=object@init,A=matrix(0,1,1),B=object@dens,ntimes=ntimes,homogeneous=TRUE)$gamma
		    vstate <- t(apply(gamma,1,ind.max))
		    B <- object@dens
		    if(na.allow) B[is.na(B)] <- 1
		    ll <- sum(log((apply(B,c(1,3),prod))[cbind(1:sum(ntimes(object)),vstate)]))
		}
		attr(ll, "df") <- freepars(object)
		attr(ll, "nobs") <- nobs(object)
		class(ll) <- "logLik"
		ll
	}
)

setMethod("logLik",signature(object="mix.fitted.classLik"),
    function(object,method=c("classification","fb","lystig"),na.allow=TRUE) {
        method <- match.arg(method)
        callNextMethod(object=object,method=method,na.allow=na.allow)
     }
)
