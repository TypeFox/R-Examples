## ###################################################################
## A function for checking the balance in covariate distribution
## between treatment groups if class of object is 'matched.data.frame'
## ###################################################################
ps.balance.matched.data.frame <- function(object,
                                          sel           = NULL,
                                          treat         = NULL,
                                          stratum.index = NULL,
                                          match.index   = NULL,
                                          method        = "classical",
                                          cat.levels    = 2,
                                          alpha         = 5,
                                          equal         = TRUE) 
{
  
  ## #################
  ## Check match.index
  if (is.null(match.index)){
    d1 <- object$data
    d2 <- object$data.matched
 
    d1[,"match.strata"] <- rep(1,nrow(d1))  ## unmatched
    d2[,"match.strata"] <- rep(2,nrow(d2))  ## matched
  
    data <- (rbind(d1,d2))
    data$match.strata <- as.factor(data$match.strata)

    match.index <- data$match.strata
    name.match.index <- "match.strata"
    
  }else{
    if (is.character(match.index) | is.numeric(match.index)){
      if (is.character(match.index)){      

        A <- find.strata(data   = object$data,
                         strata = match.index,
                         match  = TRUE)

        d1 <- object$data             ## A[[1]] contains match index
        d2 <- object$data[A[[1]]>0,]  ## A[[1]]>0 : matched obs

        d1[,"match.strata"] <- rep(1,nrow(d1))  ## unmatched
        d2[,"match.strata"] <- rep(2,nrow(d2))  ## matched
        
        data <- (rbind(d1,d2))
        data$match.strata <- as.factor(data$match.strata)

        match.index <- data$match.strata
        name.match.index <- A[[2]]
        
      }else{

        A <- find.strata(data   = object$data,
                         strata = match.index,
                         match  = TRUE)    

        d1 <- object$data
        d2 <- object$data[A[[1]]>0,]  ## A[[1]]>0 : matched obs

        d1[,"match.strata"] <- rep(1,nrow(d1))  ## unmatched
        d2[,"match.strata"] <- rep(2,nrow(d2))  ## matched
        
        data <- (rbind(d1,d2))
        data$match.strata <- as.factor(data$match.strata)

        match.index <- data$match.strata
        name.match.index <- A[[2]]
        
      }
    }else{
      stop("Argument 'match.index' must be numeric or a string.")
    } 
  }


  ## ########
  ## find sel
  if(is.null(sel)){
    sel <- data[,-which(c(names(data)==name.match.index))]
  }else{
    sel <- find.sel(data, sel)
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
        stop("Argument 'treat' has to be either numeric or a string.")
      }
  }

  if (any(names(sel)==name.treat))
    warning("Argmuent 'sel' contains argument 'treat'.")

  
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
                                    index      = match.index,
                                    method     = method,
                                    cat.levels = cat.levels,
                                    match.T    = TRUE,
                                    alpha      = alpha,
                                    equal      = equal)
  }
  
  object$data.matched      <- d2[,-which(c(names(d2)=="match.strata"))]
  object$bal.test          <- bal.test
  object$treat             <- treat
  object$name.treat        <- name.treat
  object$match.index       <- object$match.index
  object$name.match.index  <- object$name.match.index

  class(object) <- c("bal.matched.data.frame",
                     class(object)[class(object)!="bal.matched.data.frame"])
  
  
  return(object)
  
}
