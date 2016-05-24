#####################
#warnings for GLM SA
#####################
warnings <- function(formula,     #formula: assume treatment is 1st term on rhs
                     resp.family,  #family for GLM of model for response
                     trt.family,  #family for GLM of model for treatment
                     theta,       #Pr(U=1) for binomial model
                     grid.dim,    #final dimensions of output grid
                     standardize,  #Logical: should variables be standardized?
                     nsim,        #number of simulated Us to average over per cell in grid
                     zero.loc,		#location of zero along line y=x, as fraction in [0,1], or "full" if full range is desired
                     verbose,
                     buffer, 		  #restriction to range of coef on U to ensure stability around the edges
                     zetay.range,  	#custom range for zeta^y, e.g.(0,10), zero.loc will be overridden.
                     zetaz.range,  	#custom range for zeta^z, e.g.(-2,2), zero.loc will be overridden.
                     weights,     #some user-specified vector or "ATE", "ATT", or "ATC" for GLM.sens to create weights.
                     data) {
  
  if(is.null(data)) stop(paste("Either a matrix or data.frame object must be specified in data field."))
  
  #check that data is a data frame
  if(identical(class(data),"matrix")) {
    if(verbose) warning("coerced matrix to data frame")
    data = data.frame(data)
  }
  else if(!identical(class(data),"data.frame")) {
    stop(paste("Data is not a data.frame object"))
  } 
  
  #Code for listwise deletion
  postdata=na.omit(data)
  nobs.deleted = dim(data)[1] - dim(postdata)[1]
  if ((verbose) & nobs.deleted>0) warning(nobs.deleted, " observations were deleted listwise.\n") 
  
  #extract variables from formula
  form.vars <- parse.formula(formula, data)
  Y = form.vars$resp
  Z = form.vars$trt
  X = form.vars$covars 
   
  #check trt.family
  if(is.binary(Z)){ #binary treatment
    
    if(identical(class(trt.family),"family")) {
      if(identical(trt.family$family,"binomial")) {
        if(identical(trt.family$link,"probit")) {
          if(verbose) cat("Binomial family with probit link function is assumed in the treatment model.\n")
        }else if(identical(trt.family$link,"logit")) {
          warning("GLM.sens is not compatible with logistic link. Binomial family with probit link function is assumed in the treatment model.")
          trt.family = binomial(link="probit")
        }else{
          stop(paste("binomial(link=\"probit\") is the only available option for the binary treatment."))
        }
      }else{
        stop(paste("binomial(link=\"probit\") is the only available option for the binary treatment."))
      }
    }
    
    if(identical(class(trt.family),"function")) {
      if(identical(trt.family,binomial)) {
        if(verbose) cat("Binomial family with probit link function is assumed in the treatment model.\n")
        trt.family = binomial(link="probit")        
      }else{
        stop(paste("binomial(link=\"probit\") is the only available option with a binary treatment."))
      }
    }
    
    if(identical(class(trt.family),"character")) {
      if(identical(trt.family,"probit")||identical(trt.family,"binomial")||identical(trt.family,"binary")){
        if(verbose) cat("Binomial family with probit link function is assumed in the treatment model.\n")
        trt.family = binomial(link="probit")
      }else if(identical(trt.family,"logit")||identical(trt.family,"logistic")){
        warning("GLM.sens is not compatible with logistic link. Binomial family with probit link function is assumed in the treatment model.")
        trt.family = binomial(link="probit")        
      }else{
        stop(paste("binomial(link=\"probit\") is the only available option with a binary treatment."))
      }
    }
  }else{ #continuous treatment
    
    if(identical(class(trt.family),"family")) {
      if(identical(trt.family$family,"gaussian")) {
        if(identical(trt.family$link,"identity")) {
          if(verbose) cat("Gaussian family with identity link function is assumed in the treatment model.\n")
          trt.family = gaussian
        }else{
          stop(paste("Gaussian family with identity link function is the only available option with a continuous treatment."))
        }
      }else{
        stop(paste("Gaussian family with identity link function is the only available option with a continuous treatment."))
      }
    }
    
    if(identical(class(trt.family),"function")) {
      if(identical(trt.family,gaussian)) {
        if(verbose) cat("Gaussian family with identity link function is assumed in the treatment model.\n")
      }else{
        stop(paste("Gaussian family is the only available option with a continuous treatment."))
      }
    }
    
    if(identical(class(trt.family),"character")) {
      if(identical(trt.family,"gaussian")||identical(trt.family,"Gaussian")||identical(trt.family,"normal")||identical(trt.family,"identity")||identical(trt.family,"continuous")) {
        if(verbose) cat("Gaussian family with identity link function is assumed in the treatment model.\n")
        trt.family = gaussian     
      }else{
        stop(paste("Gaussian family is the only available option with a continuous treatment."))
      }
    }
  }
  
  
  #check resp.family
  if(is.binary(Y)){ #binary outcome
    stop(paste("An outcome variable needs to be continuous."))
  }else{ #continuous outcome
    
    if(identical(class(resp.family),"family")) {
      if(identical(resp.family$family,"gaussian")) {
        if(identical(resp.family$link,"identity")) {
          if(verbose) cat("Gaussian family with identity link function is assumed in the response model.\n")
          resp.family = gaussian
        }else{
          stop(paste("Gaussian family with identity link function is the only available option in the response model."))
        }
      }else{
        stop(paste("Gaussian family with identity link function is the only available option in the response model."))
      }
    }
    
    if(identical(class(resp.family),"function")) {
      if(identical(resp.family,gaussian)) {
        if(verbose) cat("Gaussian family with identity link function is assumed in the response model.\n")        
      }else{
        stop(paste("Gaussian family with identity link function is the only available option in the response model."))
      }
    }
    
    if(identical(class(resp.family),"character")) {
      if(identical(resp.family,"normal")||identical(resp.family,"continuous")||identical(resp.family,"gaussian")) {
        if(verbose) cat("Gaussian family with identity link function is assumed in the response model.\n")        
        resp.family = gaussian
      }else{
        stop(paste("Gaussian family with identity link function is the only available option in the response model."))
      }
    }
  }
  
  
  #check whether the dimentions of grid are at least 2.
  if(!is.null(grid.dim) && (length(grid.dim) != 2)) {
    stop(paste("Error: grid dimenstions must a vector of length 2"))
  }
  
  #If single value entered for range with 1 cell grid, expand to length 2 to avoid errors
  if(!is.null(zetay.range) && length(zetay.range) == 1 && (grid.dim == 1 | all.equal(grid.dim, c(1,1)))) 
    zetay.range = c(zetay.range, zetay.range)
  if(!is.null(zetaz.range) && length(zetaz.range) == 1 && (grid.dim == 1 | all.equal(grid.dim, c(1,1)))) 
    zetaz.range = c(zetaz.range, zetaz.range)
  
  #If single value entered for range with >1 cell grid, return an error
  if(!is.null(zetay.range) && (length(zetay.range) != 2)) {
    stop(paste("Error: zeta.y range must a vector of length 2"))
  }
  if(!is.null(zetaz.range) && (length(zetaz.range) != 2)) {
    stop(paste("Error: zeta.z range must a vector of length 2"))
  }
  
  return(list(formula=formula,
              resp.family=resp.family,
              trt.family=trt.family,
              theta=theta,
              grid.dim=grid.dim,
              standardize=standardize,
              nsim=nsim,
              zero.loc=zero.loc,
              verbose=verbose,
              buffer=buffer,
              weights=weights,
              zetay.range=zetay.range,
              zetaz.range=zetaz.range,
              data=postdata))
}            


#################
#warnings for BART SA
#################

warningsBART <- function(formula,     #formula: assume treatment is 1st term on rhs
                         grid.dim,    #final dimensions of output grid
                         verbose,
                         zetay.range,  	#custom range for zeta^y, e.g.(0,10), zero.loc will be overridden.
                         zetaz.range,  	#custom range for zeta^z, e.g.(-2,2), zero.loc will be overridden.
                         weights,     #some user-specified vector or "ATE", "ATT", or "ATC" for GLM.sens to create weights.
                         data) {
  
  if(is.null(data)) stop(paste("Either a matrix or data.frame object must be specified in data field."))
  
  #check that data is a data frame
  if(identical(class(data),"matrix")) {
    if(verbose) warning("coerced matrix to data frame")
    data = data.frame(data)
  }
  else if(!identical(class(data),"data.frame")) {
    stop(paste("Data is not a data.frame object"))
  } 
  
  #Code for listwise deletion
  postdata=na.omit(data)
  nobs.deleted = dim(data)[1] - dim(postdata)[1]
  if ((verbose) & nobs.deleted>0) warning(nobs.deleted, " observations were deleted listwise.\n") 
  
  #extract variables from formula
  form.vars <- parse.formula(formula, data)
  Y = form.vars$resp
  Z = form.vars$trt
  X = form.vars$covars 
  
  #check whether the dimentions of grid are at least 2.
  if(!is.null(grid.dim) && (length(grid.dim) != 2)) {
    stop(paste("Error: grid dimenstions must a vector of length 2"))
  }
  
  #If single value entered for range with 1 cell grid, expand to length 2 to avoid errors
  if(!is.null(zetay.range) && length(zetay.range) == 1 && (grid.dim == 1 | all.equal(grid.dim, c(1,1)))) 
    zetay.range = c(zetay.range, zetay.range)
  if(!is.null(zetaz.range) && length(zetaz.range) == 1 && (grid.dim == 1 | all.equal(grid.dim, c(1,1)))) 
    zetaz.range = c(zetaz.range, zetaz.range)
  
  #If single value entered for range with >1 cell grid, return an error
  if(!is.null(zetay.range) && (length(zetay.range) != 2)) {
    stop(paste("Error: zeta.y range must a vector of length 2"))
  }
  if(!is.null(zetaz.range) && (length(zetaz.range) != 2)) {
    stop(paste("Error: zeta.z range must a vector of length 2"))
  }
  
  return(list(formula=formula,
              grid.dim=grid.dim,
              zetay.range=zetay.range,
              zetaz.range=zetaz.range,
              data=postdata))
}            
