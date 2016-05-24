setClass("POISSONresponse",contains="GLMresponse")

# method 'fit'
# use: in EM (M step)
# returns: (fitted) response with (new) estimates of parameters

# methods 'logDens' & dens
# use: instead of density slot in rModel
# returns: matrix with log(p(y|x,parameters))
setMethod("logDens","POISSONresponse",
	function(object) {
		dpois(x=object@y,lambda=predict(object),log=TRUE)
	}
)

setMethod("dens","POISSONresponse",
	function(object,log=FALSE) {
		dpois(x=object@y,lambda=predict(object),log=log)
	}
)

setMethod("simulate",signature(object="POISSONresponse"),
  function(object,nsim=1,seed=NULL,times) {
    if(!is.null(seed)) set.seed(seed)
    if(missing(times)) {
      # draw in one go
      lambda <- predict(object)
    } else {
      lambda <- predict(object)[times,]
    }
    nt <- nrow(lambda)
    response <- rpois(nt*nsim,lambda=lambda)
    #if(nsim > 1) response <- matrix(response,ncol=nsim)
    response <- as.matrix(response)
    return(response)
  }
)
