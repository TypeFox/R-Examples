## #############################################################
## A function for checking the balance in covariate distribution
## between treatment groups if class of object is 'data.frame'
## #############################################################
ps.balance.data.frame <- function(object,
                                  sel           = NULL,
                                  treat         = NULL,
                                  stratum.index = NULL,
                                  match.index   = NULL,
                                  method        = "classical",
                                  cat.levels    = 2,
                                  alpha         = 5,
                                  equal         = TRUE) 
{
  data <- object

  ## ##############
  ## find selection
  if(is.null(sel)){
    sel <- data
  }else{
    sel <- find.sel(data, sel)
  }
  
  ## ##########
  ## find treat
  if (is.null(treat)){
    stop("Argument 'treat' is needed.")
  }else{
    if (is.character(treat) | is.numeric(treat)){
      A <- find.treat(data=data,
                      treat=treat)
      treat      <- A[[1]]
      name.treat <- A[[2]]
    }else{
      stop("Argument 'treat' has to be either numeric or a string.") 
    }    
  }

  ## #################################
  ## find stratum.index or match.index
  if (is.null(stratum.index) & is.null(match.index)){

    stop("Either argument 'stratum.index' or 'match.index' is needed.")

  }else{

    if (!is.null(stratum.index) & !is.null(match.index)){     

      warning("Both arguments 'stratum.index' and 'match.index' are given: 'match.index' is used.")

    }

    if (!is.null(match.index)){

      stratum.index      <- NULL
      name.stratum.index <- NULL
      
      if (is.character(match.index) | is.numeric(match.index)){

        A <- find.strata(data    = data,
                         strata  = match.index,
                         match   = TRUE)
        
        match.index        <- A[[1]]
        name.match.index   <- A[[2]]

        match.strata <- c(rep(1,nrow(data)),
                          rep(2,nrow(data[match.index!=0,])))  
        data.new <- rbind(data, data[match.index!=0,])
        data.new$match.strata <- match.strata
          
        sel.match        <- as.data.frame(sel[match.index!=0,])
        names(sel.match) <- names(sel)
        sel <- rbind(sel, sel.match)
        treat1 <- c(treat, treat[match.index!=0])
        
      }else{
        stop("Argument 'match.index' has to be either numeric or a string.")
      }
     
    }else{ ## is.null(match)

      name.match.index <- NULL

      if (is.character(stratum.index) | is.numeric(stratum.index)){

        A <- find.strata(data   = data,
                         strata = stratum.index,
                         match  = FALSE)
        
        stratum.index        <- A[[1]]
        name.stratum.index   <- A[[2]]
        levels.stratum.index <- A[[3]]
        
      }else{
        stop("Argument 'stratum.index' has to be either numeric or a string.")
      }
    } 
  }

  if (is.null(stratum.index)){
    if (any(names(sel) == name.treat) |
        any(names(sel) == name.match.index)){
      warning("Argmuent 'sel' contains argument 'treat' or 'match.index'.")
    }
  }else{
    if (any(names(sel) == name.treat) |
        any(names(sel) == name.stratum.index)){
      warning("Argmuent 'sel' contains argument 'treat' or 'stratum.index'.")
    }
  }

  
  ## ###########
  ## Check alpha
  if (!is.numeric(alpha)){
    stop("Argument 'alpha' must be numeric.")
  }else{
    if (alpha < 1 | alpha > 100){
      stop("Argument 'alpha' must be in [1,100].")
    }
  }

  ## ################
  ## selection method
  if (is.null(stratum.index)){
    
    treat <- treat1
    index <- data.new$match.strata
    match.T <- TRUE

  }else{
    
    index <- stratum.index
    match.T <- FALSE
  }
    

  if (is.character(method) | is.numeric(method)){

    bal.test <- find.balance.method(sel        = sel,
                                    treat      = treat,
                                    index      = index,
                                    method     = method,
                                    cat.levels = cat.levels,
                                    match.T    = match.T,
                                    alpha      = alpha,
                                    equal      = equal)
  }
  
  object <- list(data                 = data,
                 treat                = treat,
                 name.treat           = name.treat,
                 stratum.index        = stratum.index,
                 name.stratum.index   = name.stratum.index,
                 #levels.stratum.index = levels.stratum.index,
                 match.index          = match.index,
                 name.match.index     = name.match.index,
                 bal.test             = bal.test)

  class(object) <- c("bal.data.frame",
                     class(object)[class(object)!="bal.data.frame"])
  

  return(object)
  
}
