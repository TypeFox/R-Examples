stability <- function(x, y, EV, threshold = 0.75, B = 100, fraction = 0.5,
                      model.selector = lasso.firstq,
                      args.model.selector = NULL,
                      parallel = FALSE, ncores = getOption("mc.cores", 2L),
                      verbose = FALSE)
{
  ## Purpose:
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Lukas Meier, Date: 28 Mar 2013, 12:49

  ## error checking

  if(threshold > 1 | threshold < 0.5)
    stop("threshold has to be in (0.5, 1)")
  ##-   if(any(EV >= q))
  ##-     stop("q should be larger than EV")
  ##- 
  ##-   if(length(q) > 1) ## q has to be a scalar (EV can be a vector)
  ##-     stop("q must be a scalar")

  ##

  n <- nrow(x)
  p <- ncol(x)
  
  col.nam <- colnames(x)

  q <- ceiling(sqrt(EV * p * (2 * threshold - 1)))
  
  ##thresholds <- 0.5 * (1 + q^2 / (p * EV)) ## vector of thresholds

  ##if(any(thresholds > 1))
  ##  warning("Some thresholds larger than 1. Decrease q or increase EV.")

  ## Matrix of selected models:
  ## rows = subsamples
  ## cols = predictors
  sel.mat <- matrix(FALSE, nrow = B, ncol = p) 

  sel.n <- floor(fraction * n)

  oneSample <- function(...){
    sel <- sample(1:n, sel.n, replace = FALSE)
    ## Current sub-sampled data
    x.sel <- x[sel,]
    y.sel <- y[sel]

    ## Get selected model
    sel.model <- do.call(model.selector, c(list(x = x.sel, y = y.sel, q = q),
                                           args.model.selector))
    out <- logical(ncol(x))
    out[sel.model] <- TRUE
    out
  }

  if(parallel){
    if(verbose)
      cat("...starting parallelization of bootstrap samples\n")
    sel.mat <- matrix(unlist(mclapply(1:B, oneSample, mc.cores = ncores)),
                      nrow = B, byrow = TRUE)
  }
  else{
    ## Subsampling
    for(b in 1:B){
      if(verbose)
        cat("...Subsample", b, "\n")
      sel.mat[b, ] <- oneSample()
      ##-     sel <- sample(1:n, sel.n, replace = FALSE)
      ##- 
      ##-     ## Current sub-sampled data
      ##-     x.sel <- x[sel,]
      ##-     y.sel <- y[sel]
      ##- 
      ##-     ## Get selected model
      ##-     sel.model <- do.call(model.selector, c(list(x = x.sel, y = y.sel, q = q),
      ##-                                            args.model.selector))
      ##-     sel.mat[b, sel.model] <- TRUE
    }
  }

  ## Get selection frequencies
  freq <- colMeans(sel.mat)
  names(freq) <- col.nam
        
  out <- list(); ##length(out) <- length(EV)
  
  ##for(i in 1:length(EV)){
  sel.current        <- which(freq >= threshold)
  names(sel.current) <- col.nam[sel.current]
    
  if(length(sel.current) == 0)
    sel.current <- NULL ## for safety reasons...
    
  out <- sel.current
  ##}

  out <- list(selected   = sel.current,
              EV         = EV,
              threshold  = threshold,
              freq       = freq,
              q          = q,
              method     = "stability",
              call       = match.call())
  
  class(out) <- "hdi"
  return(out)
}
