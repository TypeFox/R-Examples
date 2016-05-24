## ###################################################
## Function to match data if object is of class pscore
## ###################################################

ps.match.pscore <- function(object,
                            object.control     = NULL,  
                            matched.by         = NULL,
                            control.matched.by = matched.by,
                            who.treated        = 1,
                            treat              = NULL,
                            name.match.index   = "match.index",
                            ratio              = 1,
                            caliper            = "logit",
                            x                  = 0.2,
                            givenTmatchingC    = TRUE,
                            bestmatch.first    = TRUE,
                            setseed            = FALSE,
                            combine.output     = TRUE)
{ 
  ## ############
  ## Extract data
  data <- object$data       

  ## ######################
  ## Check name.match.index
  if(any(names(data) == name.match.index))
    stop(paste("Argument 'name.match.index' =",
               name.match.index,
               " already exists in data.", sep=""))
  
  ## ################
  ## Check matched.by
  if ( is.null(matched.by) ){
    match.vec  <- object$pscore
    matched.by <- object$name.pscore
  }else{
    if (is.character(matched.by) | is.numeric(matched.by)){     
      A <- find.sel(data     = data,
                    sel      = matched.by,
                    sel.name = "matched.by")
      match.vec  <- A[,1]
      matched.by <- names(A)[1]
    }else{
      stop("Argument 'matched.by' must be either numeric or a string.")
    }
  }
 
  ## #######################
  ## Extract values of treat
  if ( is.null(treat) ){
    name.treat   <- object$name.treat
    treat        <- data[, name.treat]
    treat.values <- levels(as.factor(treat))
  }else{
    if (is.character(treat) | is.numeric(treat)){
      A <- find.treat(data  = data,
                      treat = treat)
      treat      <- A[[1]]
      name.treat <- A[[2]]
      treat.values <- levels(as.factor(treat))    }else{
      stop("Argument 'treat' has to be either numeric or a string.")
    }
  }
  
  if (any(treat.values == who.treated)){
    tvect <- data[,name.treat] == treat.values[treat.values == who.treated]
    ## TRUE if treated
    cvect <- data[,name.treat] == treat.values[treat.values != who.treated]
    ## TRUE if control
  }else{
    stop("Who was treated? Define argument 'who.treated'.")
  }

  ## #########################################
  ## Separate data regarding treated/untreated
  data1 <- data[tvect,]
  data2 <- data[cvect,]

 
  ## ######################
  ## Call match function
  match <- ps.matchcaliper(vect1 = data1[, matched.by],
                           vect2 = data2[, matched.by],
                           ratio,
                           caliper,
                           x,
                           givenTmatchingC,
                           bestmatch.first,
                           setseed )

  ## #############
  ## Manage output

  ## create new column where match.info is included
  data[,name.match.index] <- rep(NA,nrow(data))

  tvect[is.na(tvect)] <- cvect[is.na(cvect)] <- TRUE  

  data[tvect == TRUE, name.match.index] <- match$pairvect1
  data[cvect == TRUE, name.match.index] <- match$pairvect2

  match.index <- data[, name.match.index]
                      
  match.parameters <- list(caliper            = match$caliper,
                           ratio              = match$ratio,
                           who.treated        = who.treated,
                           givenTmatchingC    = match$givenTmatchingC,
                           bestmatch.first    = match$bestmatch.first)
  
  object$data                    <- data
  object$data[,name.match.index] <- match.index
  object$data.matched            <- data[data[, name.match.index] != 0, ]
  object$name.match.index        <- name.match.index
  object$match.index             <- match.index
  object$match.parameters        <- match.parameters
  object$matched.by              <- matched.by
  
  object$treat                   <- treat       ## if name.treat is
                                                ## specified,
  object$name.treat              <- name.treat  ## those arguments
                                                ## must be modified
                                                
  class(object) <- c("matched.pscore",
                     class(object)[class(object)!="matched.pscore"])
  
  
  return(object)

}

