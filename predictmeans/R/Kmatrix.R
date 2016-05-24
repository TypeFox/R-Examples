Kmatrix <- function(model, modelterm, covariate=NULL, prtnum=FALSE)
{
   if (inherits(model, "mer") || inherits(model, "merMod")) {
    if(!lme4::isLMM(model) && !lme4::isGLMM(model))
      stop("Can't handle a nonlinear mixed model")
    thecall <- slot(model, "call")
    contrasts <- attr(model.matrix(model), "contrasts")
  }else if (inherits(model, "lme")) {
    thecall <- model$call
    contrasts <- model$contrasts
  }else if (inherits(model, "gls")) {
    thecall <- model$call
    contrasts <- model$contrasts
  }else if (inherits(model, "lm")) {  
    thecall <- model$call
    contrasts <- attr(model.matrix(model), "contrasts")
  }else stop(paste("Can't handle an model of class", class(model)[1]))

  cov.reduce <- function(x, name) mean(x)
  
  fac.reduce <- function(coefs, lev) apply(coefs, 2, mean) 
  
  # Get model formula and index of response
  Terms <- terms(model)
  yname <- as.character(attr(Terms, "variables"))[[2]]
  Terms <- delete.response(Terms)
  
  # get the pure formula w/o extra stuff
  formrhs <- formula(Terms)
  
  # All the variables in the model
  nm <- all.vars(formrhs)
  nm <- nm[nm!="pi"]
  
  # Figure out if any are coerced to factor or ordered
  anm <- all.names(formrhs)    
  coerced <- anm[1 + grep("factor|as.factor|ordered|as.ordered", anm)]
  
  # Obtain a simplified formula -- needed to recover the data in the model
  form <- as.formula(paste("~", paste(nm, collapse = "+")))
  envir <- attr(Terms, ".Environment")
  X <- model.frame(form, eval(thecall$data, envir=envir), 
                  subset = eval(thecall$subset, enclos=envir),
                  na.action = na.omit, drop.unused.levels = TRUE)
  preddf <- X
  baselevs <- xlev <- matdat <- list()  
  all.var.names <- names(X)

  for (xname in all.var.names) {
    obj <- X[[xname]]
    if (is.factor(obj)) {            
      xlev[[xname]] <- levels(obj)
      baselevs[[xname]] <- levels(obj)
    }
    else if (is.matrix(obj)) {
      # Matrices -- reduce columns thereof, but don't add to baselevs
      matdat[[xname]] <- apply(obj, 2, cov.reduce, xname)
    }
    else {
      # single numeric pred but coerced to a factor - use unique values
      if (length(grep(xname, coerced)) > 0)             
        baselevs[[xname]] <- sort(unique(obj))
      
      # Ordinary covariates - summarize if not in 'at' arg
      else baselevs[[xname]] <- cov.reduce(obj, xname)       
    }
  }
 
  covlevname <- setdiff(names(baselevs), c(names(xlev), coerced))

  if ((!is.null(covariate) && !covariate%in%c("NULL", "")) && is.numeric(covariate)) baselevs[covlevname] <- as.list(covariate)
  if ((!is.null(covariate) && !covariate%in%c("NULL", "")) && is.character(covariate) && covariate%in%covlevname) baselevs[[covariate]] <- seq(min(X[[covariate]]), max(X[[covariate]]), length=50)

  if (all(length(covlevname)!=0, prtnum)) {
    cat("\n", "The predicted means are estimated at \n\n")
    print(round( unlist(baselevs[covlevname]), 4))
    cat("\n")
  }

  # OK. Now make a grid of the factor levels of interest, along w/ covariate "at" values
  grid <- do.call("expand.grid", baselevs)

  # add any matrices
  for (nm in names(matdat))
    grid[[nm]] <- matrix(rep(matdat[[nm]], each=nrow(grid)), nrow=nrow(grid))

  # Now make a new dataset with just the factor combs and covariate values we want for prediction
  # WARNING -- This will overwrite X, so get anything you need from X BEFORE we get here
   
  m <- model.frame(Terms, grid, na.action = na.pass, xlev = xlev)
  X <- model.matrix(Terms, m, contrasts.arg = contrasts)

  # All factors (excluding covariates)
  # version 1.10 - no longer excluding covariates
  allFacs <- all.var.names
  
  ### Array of indexes for rows of X, organized by dimensions
  row.indexes <- array(seq_len(nrow(X)), sapply(baselevs, length))  

  # convert a string to a formula
  form <- as.formula(paste("~",modelterm))

  # These are the variables involved; and the label to use in the results
  facs <- all.vars(form)
  if ((!is.null(covariate) && !covariate%in%c("NULL", "")) && all(is.character(covariate), !covariate%in%facs)) facs <- c(facs, covariate)   
     
  if (any(sapply(facs, function(nm) length(grep(nm, allFacs)) == 0)))
    stop(paste("Unknown factor(s) in specification:", paste(form, collapse=" ")))
  
  # create the grid of factor combinations
  levs <- list()
  for (f in facs) levs[[f]] <- baselevs[[f]]
  
  combs <- do.call("expand.grid", levs)
  
  fctnames <- do.call("expand.grid", levs[rev(names(levs))])
  fctnames <- fctnames[, rev(names(fctnames)), drop=FALSE]
  rnK <- do.call("paste", c(fctnames, sep=":")) 
  
  K <- plyr::alply(row.indexes, match(facs, names(baselevs)), function(idx) {
    fac.reduce(X[idx, , drop=FALSE], "")
  })
  K <- as.matrix(as.data.frame(K))
  dimnames(K)[[2]] <- do.call("paste", c(combs, sep=":"))      
  K <-t(K)
  K <- K[rnK, , drop=FALSE]
  
  return(list(K=K, fctnames=fctnames, response=yname, preddf=preddf))
}
