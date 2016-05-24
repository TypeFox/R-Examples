## ########################################
## Function to match data from a data.frame
## ########################################
##
## It is possible to give one or two df's: if only one df is given, it
## must include at least one obs of each treatment group

ps.match.data.frame <- function(object,
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
  ## ################
  ## Check matched.by
  if ( is.null(matched.by) ){
    stop("Argument 'matched.by' is needed.")
  }else{
    if (is.character(matched.by) | is.numeric(matched.by)){
      A <- find.sel(data     = object,       
                    sel      = matched.by,  
                    sel.name = "matched.by")
      matched.by <- names(A)[1]
    }else{
      paste("Argument 'matched.by' must be either numeric or a string.")
    }
  }

 
  ## ####################
  ## Check object.control
  if ( is.null(object.control) ){

    ## ###########
    ## Check treat
    if (is.null(treat)){
      stop("Argument 'treat' is needed.") 
    }else{

      if (is.character(treat) | is.numeric(treat)){
          A <- find.treat(data  = object,
                          treat = treat)
          treat      <- A[[1]]
          name.treat <- A[[2]]
        }else{
          stop("Argument 'treat' has to be either numeric or a string.")
        }
      
      if (any(names(object) == name.treat)){
        treat.values <- levels(as.factor(object[, name.treat]))
      
        if (length(treat.values) == 1){
          stop(paste("Argument 'object' includes only observations with ",
                     name.treat, " = ", treat.values,".", sep=""))
        }else{

          ## #################
          ## Check who.treated
          if (any(treat.values == who.treated)){
            
            tvect <- object[, name.treat] == treat.values[treat.values == who.treated] 
            cvect <- object[, name.treat] == treat.values[treat.values != who.treated]

            object.treated <- object[tvect,]
            object.control <- object[cvect,]            

          }else{           
            stop("Who was treated? Define argument 'who.treated'.")
          }
        }
      }

      ## #######################
      ## Check name.match.index
      if(any(names(object) == name.match.index)){
        stop(paste("Argument 'name.match.index'=",
                   name.match.index,
                   " already exists in object.", sep=""))
      }
    } 
  }else{ ## not missing object.control
    
    object.treated <- object
    
    if( !is.data.frame(object.control) )
      stop("Argument 'object.control' has to be of class 'data.frame'.")

    
    ## ###########
    ## Check treat
    if ( is.null(treat) ){
      stop("Argument 'treat' is needed.") 
    }else{
      
      if ( is.character(treat) | is.numeric(treat) ){

        A1 <- find.treat.match(data      = object, 
                               treat     = treat,
                               obj.name  = "treated")
        treat1      <- A1[[1]]
        name.treat1 <- A1[[2]]
        
        A2 <- find.treat.match(data      = object.control, 
                               treat     = treat,
                               obj.name  = "control")
        treat2      <- A2[[1]]
        name.treat2 <- A2[[2]]

        if ( name.treat1 != name.treat2 ){
          stop("The name of argument 'treat' differs in input objects.")

        }else{
          if (length(levels(as.factor(c(treat1,treat2)))) != 2)
            stop("Argument 'treat' has less or more than two values.")

          treat <- c(treat1,treat2)
          name.treat <- name.treat1
        }
          
      }else{
        stop("Argument 'treat' has to be either numeric or a string.")
      }
    }
  
      
    

    ## #######################
    ## Check name.match.index
    if(any(names(object.treated) == name.match.index) ||
       any(names(object.control) == name.match.index))
      stop(paste("Argument 'name.match.index'=",
                 name.match.index,
                 " already exists in data.", sep=""))
    
    ## ########################
    ## Check control.matched.by
    if (is.character(control.matched.by) || is.numeric(control.matched.by)){
      A <- find.sel(data     = object.control,       ## nur fuer Fehlermeldung, falls 
                    sel      = control.matched.by,   ## control.matched.by nicht in Daten
                    sel.name = "control.matched.by")
      control.matched.by <- names(A)[1]
    }else{
      paste("Argument 'control.matched.by' must be either numeric or a string.")
    }
  } 

  ## ######################
  ## Call matching function
  match <- ps.matchcaliper(vect1=object.treated[,matched.by],
                           vect2=object.control[,control.matched.by],
                           ratio,
                           caliper,
                           x,
                           givenTmatchingC,
                           bestmatch.first,
                           setseed)


  ## #############
  ## Manage output
  object.treated[,name.match.index] <- match$pairvect1
  object.control[,name.match.index] <- match$pairvect2

  match.parameters <- list(caliper            = match$caliper,
                           ratio              = match$ratio,
                           who.treated        = who.treated,
                           givenTmatchingC    = match$givenTmatchingC,
                           bestmatch.first    = match$bestmatch.first)

  output <- list()

  
  if(combine.output){
      
    output$data              <- rbind(object.treated, object.control)
    output$data.matched      <- output$data[output$data[, name.match.index] != 0, ]
    output$match.index       <- output$data[, name.match.index]
    output$name.match.index  <- name.match.index
    output$match.parameters  <- match.parameters
    output$matched.by        <- matched.by
    output$name.treat        <- name.treat
    
    if (matched.by != control.matched.by){
      output$control.matched.by <- control.matched.by
    }
    
    class(output) <- c("matched.data.frame")
    
    
  }else{ ## combine.output = FALSE
    
    output$data             <- list(object.treated, object.control)
    output$data.matched     <- list(object.treated[object.treated[, name.match.index] != 0, ],
                                    object.control[object.control[, name.match.index] != 0, ])
    output$match.index      <- list(object.treated[, name.match.index],
                                    object.control[, name.match.index])
    output$name.match.index <- name.match.index
    output$match.parameters <- match.parameters
    output$matched.by       <- matched.by
    output$name.treat       <- name.treat
    
    if (matched.by != control.matched.by){
      output$control.matched.by <- control.matched.by
    }
    
    class(output) <- c("matched.data.frames")
    
  }
  
  return(output)
  
}
