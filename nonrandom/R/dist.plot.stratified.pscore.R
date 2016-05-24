
dist.plot.stratified.pscore <- function(object,
                                        sel           = NULL,
                                        treat         = NULL,
                                        stratum.index = NULL,
                                        match.index   = NULL,
                                        plot.type     = 1,
                                        compare       = FALSE,
                                        ...)
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
  if (!is.character(treat) & !is.numeric(treat) & !is.null(treat)){
    stop("Argument 'treat' must be either numeric or a string.")
  }else{
    if (is.null(treat)){   
      name.treat <- object$name.treat
      treat      <- object$treat
    }else{
      A <- find.treat(data  = data,
                      treat = treat)
      treat      <- A[[1]]
      name.treat <- A[[2]]
    }
  }

  
  ## ##################
  ## find stratum.index  
  if (!is.character(stratum.index) &
      !is.numeric(stratum.index) &
      !is.null(stratum.index)){
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


  ## ###############################################
  ## check plot.type & calculate value to be plotted
  if (!is.numeric(plot.type)){    
    stop("Argument 'plot.type' must be numeric.")
  }else{
    if (plot.type == 1){
      dist.plot.bars(sel        = sel,
                     treat      = treat,
                     name.treat = name.treat,
                     index      = stratum.index,
                     name.index = name.stratum.index,
                     compare    = compare,
                     match.T    = FALSE,
                     ...)
    }else{
      if(plot.type == 2)
        dist.plot.hist(sel        = sel,
                       treat      = treat,
                       name.treat = name.treat,
                       index      = stratum.index,
                       name.index = name.stratum.index,
                       compare    = compare,
                       match.T    = FALSE,
                       ...)
      else
        stop("Argument 'plottype' must be either 1 or 2.")
    }
  }
}


#dist.plot.stratified.pscore <- function(object,
#                                        sel           = NULL,
#                                        treat         = NULL,
#                                        stratum.index = NULL,
#                                        #match.index   = NULL,
#                                        plot.type     = 1,
#                                        #cat.levels    = 10,
#                                        #plot.levels   = 5,
#                                        compare       = FALSE,
#                                        #label.match   = NULL,
#                                        #label.stratum = "Stratum",
#                                        #with.legend  = TRUE,                
#                                        #legend.title = NULL,
#                                        #legend.cex   = 0.9,              
#                                        #myoma        = c(2,2,2,2),
#                                        #mymar        = c(2,4,1,2),
#                                        #width        = 1,
#                                        #xlim         = NULL,
#                                        #ylim         = NULL,
#                                        #col          = NULL,
#                                        #las          = 1,
#                                        #font.main    = 2,
#                                        #font         = 1,
#                                        #main         = NULL,
#                                        #main.cex     = 1.2,
#                                        #sub.cex      = 0.9,
#                                        #bar.cex      = 0.8,
#                                        ...)
