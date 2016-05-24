## plot.type=1
dist.plot.bars <- function(sel,
                           treat,
                           name.treat,
                           name.index,
                           index,
                           compare,
                           match.T,
                           cat.levels  = 2,
                           label.match = NULL,
                           ...)
{
  if (match.T){
    if (is.null(label.match)){
      levels(index) <- c("Original", "Matched")
    }else{
      levels(index) <- label.match
    }
  }
  
  ## #########################################
  ## (1) define function to calculate means
  
  ## (1.1) Mean for continuous variables per treatment group
  func1 <- function(x) {
    tapply(as.numeric(sel[,x]), list(treat), mean)
  }
  
  ## (1.2) Mean for continuous variables per treatment group and per
  ## stratum  
  func2 <- function(x) {
    tapply(as.numeric(sel[,x]), list(treat, index), mean)
  }
  
  ## (1.3) Distribution (means) of categorical variables per treatment
  ## group
  func3 <- function(x) {
    t <- table(sel[,x], treat)
    for(i in 1:length(levels(as.factor(treat)))) {
      t[,i] <- t[,i]/sum(t[,i])
    }
    return(t)
  }
  
  ## (1.4) Distribution (means) of categorical variables per treatment
  ## group and per stratum
  func4 <- function(x) {
    t <- table(sel[,x], treat, index)
    
    for(i in 1:length(levels(as.factor(index)))) {
      for(j in 1:length(levels(as.factor(treat)))) 
        t[,j,i] <- t[,j,i]/sum(t[,j,i])
    }    
    return(t)
  }
  
  
  ## #######################################################
  ## (2) distingiush between categorical/continuous variables
  
  cat.index <- apply(sel,2, function(x) nlevels(as.factor(x))<=cat.levels)
  
  var.noncat <- names(sel)[1:length(sel)][cat.index == FALSE]
  var.cat    <- names(sel)[1:length(sel)][cat.index == TRUE]

  
  ## ###################
  ## (3) calculate means
  list.noncat <- list()
  list.cat <- list()

  if (!match.T){ ## stratification 
  
    if(compare == FALSE){ ## w/o unstratified data

      if(length(var.noncat) > 0)
        list.noncat <- lapply(var.noncat,func2) ## per stratum, continuous
      
      if(length(var.cat) > 0)
        list.cat <- lapply(var.cat,func4) ## per stratum, categorical
      
    }else{ ## with unstratified data
      
      if(length(var.noncat) > 0)
        list.noncat <- list(lapply(var.noncat,func1), ## continuous
                            lapply(var.noncat,func2)) ## per stratum, continuous
      
      if(length(var.cat) > 0)
        list.cat <- list(lapply(var.cat,func3),  ## categorical
                         lapply(var.cat,func4)) ## per stratum, categorical
    }
    
  }else{ ## matching
      
    if(length(var.noncat) > 0)
      list.noncat <- lapply(var.noncat,func2) ## continuous (unmatched, matched)
    
    if(length(var.cat) > 0)
      list.cat <- lapply(var.cat,func4) ## categorical (unmatched, matched)
    
  }

  res <- list(sel        = sel,
              treat      = treat,
              index      = index,
              var.noncat = var.noncat,
              var.cat    = var.cat,
              mean       = list.noncat,
              frequency  = list.cat)
  
  dist.plot.bar.plot(res,
                     name.treat,
                     name.index,
                     compare,
                     match.T,
                     ...)
  
}


