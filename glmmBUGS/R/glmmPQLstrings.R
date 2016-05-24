`glmmPQLstrings` <-
function(effects, covariates, observations, data=NULL, 
  family=c("bernoulli", "binomial", "poisson", "gaussian"), ...) {
  # call glmmPQL using strings rather than formulas

  # check to see if there are any factors in the model
  anyfactors = apply(data[,unlist(covariates), drop=FALSE], 2, is.factor)
  if(any(anyfactors))
    warning("Winbugs doesn't allow factors, \n convert them to indicator variables with the \n class.ind function in the nnet package")
  if(is.vector(family))
    family = family[1]
  

   
  if(family=="binomial" & (length(observations)>1) ) {
    theresp = as.matrix(data[,observations[1:2] ])
    theresp[,2] = theresp[,2] - theresp[,1]
    dimnames(theresp) = list(NULL, c("ones", "zeros"))
    observations = observations[-2]
  } else {
    theresp = data[,observations[1]] 
  }
  theformula = theresp ~ 1
  if(length(covariates)) {
    theformula = update(theformula, as.formula(paste("~ . +", 
      paste(unlist(covariates), collapse="+"), sep="")))
  }
  # add the offset
  if(length(observations) > 1 )
    theformula = update(theformula, as.formula(paste("~  . + offset(", observations[2], ")", sep="")))  
    
  therandom = as.formula(paste("~1|", paste(effects, collapse="/"), sep=""))

  if(is.character(family))
    if(family=="bernoulli")
        family = binomial

#return(list(theformula, therandom, data, theresp))

thepql = MASS::glmmPQL(theformula, random = therandom, family=family,
		data=data, verbose=FALSE,...)

thepql$effects = effects
thepql$covariates = covariates
thepql$observations = observations

return(thepql)

}

