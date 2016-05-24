## #############################################################
## A function for checking the balance in covariate distribution
## between treatment groups if class of object is
## 'stratified.data.frame'
## #############################################################
ps.balance.stratified.data.frame <- function(object,
                                             sel           = NULL,
                                             treat         = NULL,
                                             stratum.index = NULL,
                                             match.index   = NULL,
                                             method        = "classical",
                                             cat.levels    = 2,
                                             alpha         = 5,
                                             equal         = TRUE)
{
  data <- object$data

  ## ########
  ## find sel
  if(is.null(sel)){
    sel <- data
  }else{
    sel <- find.sel(data = data,
                    sel  = sel)
  }
  
  ## ##########
  ## find treat
  if (is.null(treat)){
    stop("Argument 'treat' is needed.")
  }else{
    if (is.character(treat) | is.numeric(treat)){
      A <- find.treat(data  = data,
                      treat = treat)
      treat      <- A[[1]]
      name.treat <- A[[2]]
  }else{
      stop("Argument 'treat' must be either numeric or a string.")
    }
  }
  
  ## ##################
  ## find stratum.index  
  if (!is.character(stratum.index) & !is.numeric(stratum.index) & !is.null(stratum.index)){
    stop("Argument 'stratum.index' has to be either numeric or a string.")
  }else{
    if (is.null(stratum.index)){
      stratum.index        <- object$stratum.index
      name.stratum.index   <- object$name.stratum.index
      levels.stratum.index <- levels(as.factor(stratum.index))
    }else{
      A <- find.strata(data   = data,
                       strata = stratum.index,
                       match  = FALSE)    
      stratum.index        <- A[[1]]
      name.stratum.index   <- A[[2]]
      levels.stratum.index <- A[[3]]
    }
  }
  
  if (any(names(sel) == name.stratum.index) |
      any(names(sel) == name.treat)){
    warning("Argmuent 'sel' contains argument 'stratum.index' or 'treat'.")}


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
    if (is.character(method) | is.numeric(method)){
      bal.test <- find.balance.method(sel        = sel,
                                      treat      = treat,
                                      index      = stratum.index,
                                      method     = method,
                                      cat.levels = cat.levels,
                                      match.T    = FALSE,
                                      alpha      = alpha,
                                      equal      = equal)      
    }

  object$bal.test             <- bal.test
  object$treat                <- treat
  object$name.treat           <- name.treat
  object$stratum.index        <- stratum.index
  object$name.stratum.index   <- name.stratum.index
  #object$levels.stratum.index <- levels.stratum.index

  class(object) <- c("bal.stratified.data.frame",
                     class(object)[class(object)!="bal.stratified.data.frame"])
  
  return(object)

}
