hdi <- function(x, y, method = "multi.split", B = NULL,
                fraction = 0.5,
                model.selector = NULL,
                EV = NULL, threshold = 0.75,         ## stability args
                gamma = seq(0.05, 0.99, by = 0.01),  ## multi-split args
                classical.fit = NULL,                ## multi-split args
                args.model.selector = NULL,          ## list
                args.classical.fit = NULL,           ## list
                verbose = FALSE, ...)
{
  ## Purpose:
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Lukas Meier, Date: 25 Mar 2013, 16:52

  ## Error checking and argument parsing

  ## ...

  #####################
  ##                 ##
  ## Multi-splitting ##
  ##                 ##
  #####################
  
  if(method %in% c("multi.split", "multi-split")){ ## for compatibility reasons
    ## argument checking
    if(is.null(B))
      B <- 100

    if(is.null(model.selector))
      model.selector <- lasso.cv
    
    if(is.null(classical.fit))
      classical.fit <- lm.pval

    ## perform multi-splitting
    out <- multi.split(x = x, y = y, B = B, fraction = fraction,
                       model.selector = model.selector,
                       classical.fit = classical.fit,
                       gamma = gamma,
                       args.model.selector = args.model.selector,
                       args.classical.fit = args.classical.fit,
                       verbose = verbose,...)
  }
  
  #########################
  ##                     ##
  ## Stability selection ##
  ##                     ##
  #########################

  else if(method == "stability"){
    ## argument checking
    if(is.null(B))
      B <- 100

    if(is.null(EV))
        stop("EV not defined")
    
    ##if(is.null(threshold))
    ##  stop("threshold not defined")
    
    if(is.null(model.selector))
      model.selector <- lasso.firstq

    ## perform stability-selection
    out <- stability(x = x, y = y, B = B, fraction = fraction,
                     threshold = threshold, model.selector = model.selector,
                     EV = EV,
                     args.model.selector = args.model.selector,
                     verbose = verbose,...)
  }
  else{
    stop("Method not (yet) defined")
  }
  ## Output parsing
  ## ...
  
  ## The following overwrites 
  out$call <- match.call()
  
  class(out) <- "hdi"
  out
}



