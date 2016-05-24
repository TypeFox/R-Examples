
dist.plot.matched.data.frame <- function(object,
                                         sel           = NULL,
                                         treat         = NULL,
                                         stratum.index = NULL,
                                         match.index   = NULL,
                                         plot.type     = 1,
                                         compare       = FALSE,
                                         ...) 
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


  ## ###############################################
  ## check plot.type & calculate value to be plotted
  if (!is.numeric(plot.type)){

    stop("Argument 'plot.type' must be numeric.")

  }else{

    if (plot.type == 1){

      dist.plot.bars(sel        = sel,
                     treat      = treat,
                     name.treat = name.treat,
                     index      = match.index,
                     name.index = name.match.index,
                     compare    = compare,
                     match.T    = TRUE,
                     ...)
     
    }else{
      if(plot.type == 2){
           
        dist.plot.hist(sel        = sel,
                       treat      = treat,
                       name.treat = name.treat,
                       index      = match.index,
                       name.index = name.match.index,
                       compare    = compare,
                       match.T    = TRUE,
                       ...)
      
      }else
      stop("Argument 'plottype' must be either 1 or 2.")
    }
  } 
}


#dist.plot.matched.data.frame <- function(object,
#                                         sel           = NULL,
#                                         treat         = NULL,
#                                         stratum.index = NULL,
#                                         match.index   = NULL,
#                                         plot.type     = 1,
#                                         cat.levels    = 10,
#                                         plot.levels   = 5,
#                                         compare       = FALSE,
#                                         label.match   = NULL,
#                                         label.stratum = "Stratum",
#                                         with.legend  = TRUE,                
#                                         legend.title = NULL,
#                                         legend.cex   = 0.9,              
#                                         myoma        = c(2,2,2,2),
#                                         mymar        = c(2,4,1,2),
#                                         width        = 1,
#                                         xlim         = NULL,
#                                         ylim         = NULL,
#                                         col          = NULL,
#                                         las          = 1,
#                                         font.main    = 2,
#                                         font         = 1,
#                                         main         = NULL,
#                                         main.cex     = 1.2,
#                                         sub.cex      = 0.9,
#                                         bar.cex      = 0.8,
#                                         ...) 
