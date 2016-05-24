##generic
dictab <- function(cand.set, modnames = NULL, sort = TRUE, ...) {
  ##format list according to model class
  cand.set <- formatCands(cand.set)
  UseMethod("dictab", cand.set)
}


##default
dictab.default <- function(cand.set, modnames = NULL, sort = TRUE, ...) {
  stop("\nFunction not yet defined for this object class\n")
}



##bugs
dictab.AICbugs <- function(cand.set, modnames = NULL, sort = TRUE, ...) {
  
  ##check if named list if modnames are not supplied
  if(is.null(modnames)) {
    if(is.null(names(cand.set))) {
      modnames <- paste("Mod", 1:length(cand.set), sep = "")
      warning("\nModel names have been supplied automatically in the table\n")
    } else {
      modnames <- names(cand.set)
    }
  }
    
  
  Results <- data.frame(Modnames = modnames)                    #assign model names to first column
  Results$pD <- unlist(lapply(cand.set, DIC, return.pD = TRUE))     #extract number of parameters
  Results$DIC <- unlist(lapply(cand.set, DIC, return.pD = FALSE))  #extract DIC                                      #
  Results$Delta_DIC <- Results$DIC - min(Results$DIC)            #compute delta DIC
  Results$ModelLik <- exp(-0.5*Results$Delta_DIC)                #compute model likelihood required to compute Akaike weights
  Results$DICWt <- Results$ModelLik/sum(Results$ModelLik)        #compute Akaike weights
  Results$Deviance <- unlist(lapply(X = cand.set, FUN = function(i) i$mean$deviance))
  
  ##check if some models are redundant
  if(length(unique(Results$DIC)) != length(cand.set)) warning("\nCheck model structure carefully as some models may be redundant\n")
    
  if(sort)  {
    Results <- Results[rev(order(Results[, 6])),] 	  #if sort=TRUE, models are ranked based on Akaike weights
    Results$Cum.Wt <- cumsum(Results[, 6])                        #display cumulative sum of Akaike weights
  } else {Results$Cum.Wt <- NULL}
  
  class(Results) <- c("dictab", "data.frame")
  return(Results)
}




##rjags
dictab.AICrjags <- function(cand.set, modnames = NULL, sort = TRUE, ...) {

    ##check if named list if modnames are not supplied
    if(is.null(modnames)) {
      if(is.null(names(cand.set))) {
        modnames <- paste("Mod", 1:length(cand.set), sep = "")
        warning("\nModel names have been supplied automatically in the table\n")
      }
      modnames <- names(cand.set)
    }
    
    Results <- data.frame(Modnames = modnames)                    #assign model names to first column
    Results$pD <- unlist(lapply(cand.set, DIC, return.pD = TRUE))     #extract number of parameters
    Results$DIC <- unlist(lapply(cand.set, DIC, return.pD = FALSE))  #extract DIC                                      #
    Results$Delta_DIC <- Results$DIC - min(Results$DIC)            #compute delta DIC
    Results$ModelLik <- exp(-0.5*Results$Delta_DIC)                #compute model likelihood required to compute Akaike weights
    Results$DICWt <- Results$ModelLik/sum(Results$ModelLik)        #compute Akaike weights
    Results$Deviance <- unlist(lapply(X = cand.set, FUN = function(i) i$mean$deviance))
    
    ##check if some models are redundant
    if(length(unique(Results$DIC)) != length(cand.set)) warning("\nCheck model structure carefully as some models may be redundant\n")
    
    if(sort)  {
      Results <- Results[rev(order(Results[, 6])),] 	  #if sort=TRUE, models are ranked based on Akaike weights
      Results$Cum.Wt <- cumsum(Results[, 6])                        #display cumulative sum of Akaike weights
    } else {Results$Cum.Wt <- NULL}
   
    class(Results) <- c("dictab", "data.frame")
    return(Results)
  }


print.dictab <-
  function(x, digits = 2, deviance = TRUE, ...) {
    cat("\nModel selection based on", colnames(x)[3], ":\n")
    cat("\n")

    #check if Cum.Wt should be printed
    if(any(names(x) == "Cum.Wt")) {
      nice.tab <- cbind(x[, 2], x[, 3], x[, 4], x[, 6], x[, "Cum.Wt"], x[, 7])
      colnames(nice.tab) <- c(colnames(x)[c(2, 3, 4, 6)], "Cum.Wt", colnames(x)[7])
      rownames(nice.tab) <- x[, 1]
    } else {
      nice.tab <- cbind(x[, 2], x[, 3], x[, 4], x[, 6], x[, 7])
      colnames(nice.tab) <- c(colnames(x)[c(2, 3, 4, 6, 7)])
      rownames(nice.tab) <- x[, 1]
    }

    ##if deviance==FALSE
    if(identical(deviance, FALSE)) {
      names.cols <- colnames(nice.tab)
      sel.dev <- which(attr(regexpr(pattern = "Deviance", text = names.cols), "match.length") > 1)
      nice.tab <- nice.tab[, -sel.dev]
    }
    
    print(round(nice.tab, digits = digits)) #select rounding off with digits argument
    cat("\n")
  }
