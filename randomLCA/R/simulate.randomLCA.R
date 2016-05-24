simulate.randomLCA <-
    function(object, nsim = 1, seed = as.integer(runif(1, 0, .Machine$integer.max)),
            ...)
{
    if(!exists(".Random.seed", envir = .GlobalEnv, inherits=FALSE))
        runif(1)		     # initialize the RNG if necessary
    RNGstate <- get(".Random.seed", envir = .GlobalEnv, inherits=FALSE)
    set.seed(seed)

# create missing data array
	ismissing <- ifelse(is.na(object$patterns[rep(1:dim(object$patterns)[1],times=object$freq),]),NA,1)
		
	value <- lapply(1:nsim, function(x)  {
		cumclassp <- cumsum(object$classp)
		u <- runif(object$nobs)
		rclass <- rep(1,object$nobs)
		if (object$nclass>1) {
			for (i in 1:(object$nclass-1)) {
				currp <- cumclassp[i]
				rclass <- ifelse(u>currp,i+1,rclass)
			}
		}
#		browser()
		# turn this into one per outcome
		alloutcomep <- as.vector(apply(matrix(rclass,ncol=1),1,function(x) {
			object$outcomep[x,]
		}))
		if (object$random) {
	# convert outcomep to appropriate scale
			if (object$probit) alloutcomex <- qnorm(alloutcomep)
			else alloutcomex <- log(alloutcomep/(1-alloutcomep))
	# create lambda randoms
			rlambda <- rep(rnorm(object$nobs),each=dim(object$patterns)[2])		
			if (!object$byclass) nrepeats <- dim(object$patterns)[2]/length(object$lambdacoef)
      else nrepeats <- dim(object$patterns)[2]/dim(object$lambdacoef)[2]
	
	if (object$level2) {
				if (object$byclass) alltau <- 
							as.vector(apply(matrix(rclass,ncol=1),1,function(x) {
								object$taucoef[x]
							}))
				else alltau <- rep(object$taucoef,object$nobs)
#				browser()
				rlambda <-  rlambda+rep(rnorm(object$nobs*dim(object$patterns)[2]/object$level2size),
						each=object$level2size)*rep(alltau,each=dim(object$patterns)[2])
			}
# browser()
			if (object$byclass) {
				alllambda <- as.vector(apply(matrix(rclass,ncol=1),1,function(x) {
								rep(object$lambdacoef[x,],nrepeats)
							}))
				alloutcomex <-  alloutcomex+rlambda*alllambda
			}
			else alloutcomex <- alloutcomex+
				rlambda*rep(object$lambdacoef,times=nrepeats*object$nobs)
			
			if (object$probit) alloutcomep <- pnorm(alloutcomex)
			else alloutcomep <- exp(alloutcomex)/(1+exp(alloutcomex))
		}
		# determine binomial random
#		browser()
		simrandom <- runif(length(alloutcomep))
#		print(simrandom[1:20])
		simoutcome <- t(matrix(ifelse(simrandom<alloutcomep,1,0),ncol=object$nobs))
	
	# now need to apply the missing values
	
		simoutcome <- as.data.frame(simoutcome*ismissing)
		names(simoutcome) <- NULL
		names(simoutcome) <- names(object$patterns)
#		print(simoutcome[88,])
		simoutcome
	})
    attr(value, "seed") <- seed
    class(value) <- "simulate.randomLCA"
    assign(".Random.seed", RNGstate, envir = .GlobalEnv, inherits=FALSE)
    value
}
