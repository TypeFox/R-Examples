
dist.plot.data.frame <- function(object,
                                 sel           = NULL,
                                 treat         = NULL,
                                 stratum.index = NULL,
                                 match.index   = NULL,
                                 plot.type     = 1,
                                 compare       = FALSE,
                                 ...) 
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


  ## ##############
  ## define match.T
  if (is.null(stratum.index)){
    
    treat <- treat1
    index <- data.new$match.strata
    name.index <- name.match.index
    match.T <- TRUE

  }else{
    
    index <- stratum.index
    name.index <- name.stratum.index
    match.T <- FALSE
  }

    
  ## ###############################################
  ## check plot.type & calculate value to be plotted
  if (!is.numeric(plot.type)){

    stop("Argument 'plot.type' must be numeric.")

  }else{

    if (plot.type == 1){

      dist.plot.bars(sel        = sel,
                     treat      = treat,
                     name.treat = name.treat,
                     index      = as.factor(index),
                     name.index = name.index,
                     compare    = compare,
                     match.T    = match.T,
                     ...)
     
    }else{
      if(plot.type == 2){
           
        dist.plot.hist(sel        = sel,
                       treat      = treat,
                       name.treat = name.treat,
                       index      = as.factor(index),
                       name.index = name.index,
                       compare    = compare,
                       match.T    = match.T,
                       ...)
      
      }else
      stop("Argument 'plottype' must be either 1 or 2.")
    }
  } 

  
}

#dist.plot.data.frame <- function(object,
#                                 sel           = NULL,
#                                 treat         = NULL,
#                                 stratum.index = NULL,
#                                 match.index   = NULL,
#                                 plot.type     = 1,
#                                 cat.levels    = 10,
#                                 plot.levels   = 5,
#                                 compare       = FALSE,
#                                 label.match   = NULL,
#                                 label.stratum = "Stratum",
#                                 with.legend  = TRUE,                
#                                 legend.title = NULL,
#                                 legend.cex   = 0.9,              
#                                 myoma        = c(2,2,2,2),
#                                 mymar        = c(2,4,1,2),
#                                 width        = 1,
#                                 xlim         = NULL,
#                                 ylim         = NULL,
#                                 col          = NULL,
#                                 las          = 1,
#                                 font.main    = 2,
#                                 font         = 1,
#                                 main         = NULL,
#                                 main.cex     = 1.2,
#                                 sub.cex      = 0.9,
#                                 bar.cex      = 0.8,
#                                 ...) 
