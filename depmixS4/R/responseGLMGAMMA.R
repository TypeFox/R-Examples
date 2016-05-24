setClass("GAMMAresponse",contains="GLMresponse")

# method 'fit'
# use: in EM (M step)
# returns: (fitted) response with (new) estimates of parameters

# methods 'logDens' & dens
# use: instead of density slot in rModel
# returns: matrix with log(p(y|x,parameters))
setMethod("logDens","GAMMAresponse",
	function(object) {
		dgamma(x=object@y,shape=predict(object),log=TRUE)
	}
)

setMethod("dens","GAMMAresponse",
	function(object,log=FALSE) {
		dgamma(x=object@y,shape=predict(object),log=log)
	}
)

setMethod("simulate",signature(object="GAMMAresponse"),
  function(object,nsim=1,seed=NULL,times) {
    if(!is.null(seed)) set.seed(seed)
    if(missing(times)) {
      # draw in one go
      shape <- predict(object)
    } else {
      shape <- predict(object)[times,]
    }
    nt <- nrow(shape)
    response <- rgamma(nt*nsim,shape=shape)
    # if(nsim > 1) response <- matrix(response,ncol=nsim)
    response <- as.matrix(response)
    return(response)
  }
)
