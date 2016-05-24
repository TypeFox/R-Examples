"SSVSquantreg" <-
  function(formula, data=NULL, tau=0.5, include=NULL, burnin=1000, mcmc = 10000, thin=1,
           verbose = 0, seed = sample(1:1000000,1), pi0a0=1, pi0b0=1, 
           ...) {
    
    ## checks
    check.offset(list(...))
    if (pi0a0<=0 || pi0b0<=0){
	stop("Parameters pi0a0 and pi0b0 must be positive. \nPlease respecify and call again.\n")
	}
    cl <- match.call()
    if (tau<=0 || tau>=1){
	stop("tau must be in (0,1). \nPlease respecify and call again.\n")
	}

    
    ## seeds
    seeds <- form.seeds(seed) 
    lecuyer <- seeds[[1]]
    seed.array <- seeds[[2]]
    lecuyer.stream <- seeds[[3]]

    ## form response and model matrices
    holder <- parse.formula(formula, data=data)
    Y <- holder[[1]]
    X <- holder[[2]]
    xnames <- holder[[3]]    
    K <- ncol(X)  # number of covariates
    q <- length(include) #number of covariates that are pre-specified to appear in the model

    if (!is.null(include)){
	if(is.character(include)){
	  if(!all(include%in%xnames)){
		include.positions<-NA
		}
          else{
		include.positions<-match(include, xnames)
		}
	  }
	else{
	  if(max(include)>length(xnames)){
		include.positions<-NA
		}
	  else{  
	  	include.positions<-include
	  	}
	  }
	if (any(is.na(include.positions))){
	  stop("One or more covariates to be included are not present in the design matrix\n or one or more indices are out of range. Please respecify and call again.\n")
	  }

    ## Bring those covariates that are pre-specified to appear in the model to the first positions in the X matrix
	X <- cbind(X[,include.positions], X[,-include.positions])
	xnames <- c(xnames[include.positions], xnames[-include.positions])
	}
        
    ## define holder for posterior sample
    sample <- matrix(data=0, mcmc/thin, 2*K)
    posterior <- NULL 
    
    ## call C++ code to draw sample
    auto.Scythe.call(output.object="posterior", cc.fun.name="SSVSquantreg", 
                     sample.nonconst=sample, tau=as.double(tau), Y=Y, X=X, q=as.integer(q), 
                     burnin=as.integer(burnin), mcmc=as.integer(mcmc), thin=as.integer(thin),
                     lecuyer=as.integer(lecuyer), 
                     seedarray=as.integer(seed.array),
                     lecuyerstream=as.integer(lecuyer.stream),
                     verbose=as.integer(verbose),
                     pi0a0 = as.double(pi0a0), pi0b0=as.double(pi0b0),
		     package="MCMCpack")
 

output <- form.mcmc.object(posterior,                                                names=rep(xnames, times=2),
                               title="SSVSquantreg Posterior Sample",
                               y=Y, call=cl
                               )
    
gammasample<-output[,1:K]
attr(gammasample, "tau")<-tau
attr(gammasample, "xnames")<-xnames
attr(gammasample, "class")<-"qrssvs"
betasample<-output[,-(1:K)]
attr(betasample,"tau")<-tau
attr(gammasample, "call") <- attr(betasample, "call") <- cl
return(list(gamma=gammasample,beta=betasample))
}




