##create generic aictab
aictab <- function(cand.set, modnames = NULL, second.ord = TRUE, nobs = NULL, sort = TRUE, ...) {
  ##format list according to model class
  cand.set <- formatCands(cand.set)
  UseMethod("aictab", cand.set)
}



##default to indicate when object class not supported
aictab.default <- function(cand.set, modnames = NULL, second.ord = TRUE, nobs = NULL, sort = TRUE, ...) {
    stop("\nFunction not yet defined for this object class\n")
  }



##aov
aictab.AICaov.lm <-
  function(cand.set, modnames = NULL, second.ord = TRUE, nobs = NULL, sort = TRUE, ...){

    ##check if named list if modnames are not supplied
    if(is.null(modnames)) {
      if(is.null(names(cand.set))) {
        modnames <- paste("Mod", 1:length(cand.set), sep = "")
        warning("\nModel names have been supplied automatically in the table\n")
      } else {
        modnames <- names(cand.set)
      }
    }
      

    ##add check to see whether response variable is the same for all models
    check.resp <- lapply(X = cand.set, FUN = function(b) formula(b)[2])
    if(length(unique(check.resp)) > 1) stop("\nYou must use the same response variable for all models\n")

    Results <- NULL
    Results <- data.frame(Modnames = modnames)                    #assign model names to first column
    Results$K <- unlist(lapply(X = cand.set, FUN = AICc, return.K = TRUE,
                               second.ord = second.ord, nobs = nobs))     #extract number of parameters
    Results$AICc <- unlist(lapply(X = cand.set, FUN = AICc, return.K = FALSE, 
                                  second.ord = second.ord, nobs = nobs))  #extract AICc                                      #
    Results$Delta_AICc <- Results$AICc - min(Results$AICc)            #compute delta AICc
    Results$ModelLik <- exp(-0.5*Results$Delta_AICc)                #compute model likelihood required to compute Akaike weights
    Results$AICcWt <- Results$ModelLik/sum(Results$ModelLik)        #compute Akaike weights

    ##check if some models are redundant
    if(length(unique(Results$AICc)) != length(cand.set)) warning("\nCheck model structure carefully as some models may be redundant\n")

    
    ##extract LL
    Results$LL <- unlist(lapply(X = cand.set, FUN = function(i) logLik(i)[1]))

    ##rename correctly to AIC
    if(second.ord == FALSE) {
      colnames(Results)[1:6] <- c("Modnames", "K", "AIC", "Delta_AIC", "ModelLik", "AICWt")
    }

    
    if(sort)  {
      Results <- Results[rev(order(Results[, 6])),] 	  #if sort=TRUE, models are ranked based on Akaike weights
      Results$Cum.Wt <- cumsum(Results[, 6])                        #display cumulative sum of Akaike weights
    } else {Results$Cum.Wt <- NULL}

    
    class(Results) <- c("aictab", "data.frame")
    return(Results)
  }



##betareg
aictab.AICbetareg <-
  function(cand.set, modnames = NULL, second.ord = TRUE, nobs = NULL, sort = TRUE, ...){

    ##check if named list if modnames are not supplied
    if(is.null(modnames)) {
      if(is.null(names(cand.set))) {
        modnames <- paste("Mod", 1:length(cand.set), sep = "")
        warning("\nModel names have been supplied automatically in the table\n")
      } else {
        modnames <- names(cand.set)
      }
    }
      

    ##add check to see whether response variable is the same for all models
    check.resp <- lapply(X = cand.set, FUN = function(b) formula(b)[2])
    if(length(unique(check.resp)) > 1) stop("\nYou must use the same response variable for all models\n")

    Results <- NULL
    Results <- data.frame(Modnames = modnames)                    #assign model names to first column
    Results$K <- unlist(lapply(X = cand.set, FUN = AICc, return.K = TRUE,
                               second.ord = second.ord, nobs = nobs))     #extract number of parameters
    Results$AICc <- unlist(lapply(X = cand.set, FUN = AICc, return.K = FALSE, 
                                  second.ord = second.ord, nobs = nobs))  #extract AICc                                      #
    Results$Delta_AICc <- Results$AICc - min(Results$AICc)            #compute delta AICc
    Results$ModelLik <- exp(-0.5*Results$Delta_AICc)                #compute model likelihood required to compute Akaike weights
    Results$AICcWt <- Results$ModelLik/sum(Results$ModelLik)        #compute Akaike weights

    ##check if some models are redundant
    if(length(unique(Results$AICc)) != length(cand.set)) warning("\nCheck model structure carefully as some models may be redundant\n")
         

    ##extract LL
    Results$LL <- unlist(lapply(X = cand.set, FUN = function(i) logLik(i)[1]))

    
    ##rename correctly to AIC
    if(second.ord == FALSE) {
      colnames(Results)[1:6] <- c("Modnames", "K", "AIC", "Delta_AIC", "ModelLik", "AICWt")
    }  

                
    if(sort)  {
      Results <- Results[rev(order(Results[, 6])),] 	  #if sort=TRUE, models are ranked based on Akaike weights
      Results$Cum.Wt <- cumsum(Results[, 6])                        #display cumulative sum of Akaike weights
    } else {Results$Cum.Wt <- NULL}

    
    class(Results) <- c("aictab", "data.frame")
    return(Results)
  }



##clm
aictab.AICsclm.clm <-
  function(cand.set, modnames = NULL, second.ord = TRUE, nobs = NULL, sort = TRUE, ...){
    
    ##check if named list if modnames are not supplied
    if(is.null(modnames)) {
      if(is.null(names(cand.set))) {
        modnames <- paste("Mod", 1:length(cand.set), sep = "")
        warning("\nModel names have been supplied automatically in the table\n")
      } else {
        modnames <- names(cand.set)
      }
    }
        
    ##add check to see whether response variable is the same for all models
    check.resp <- lapply(X = cand.set, FUN = function(b) formula(b)[2])
    if(length(unique(check.resp)) > 1) stop("\nYou must use the same response variable for all models\n")
       
    Results <- data.frame(Modnames = modnames)                    #assign model names to first column
    Results$K <- unlist(lapply(cand.set, AICc, return.K = TRUE, second.ord = second.ord, nobs = nobs))     #extract number of parameters
    Results$AICc <- unlist(lapply(cand.set, AICc, return.K = FALSE, second.ord = second.ord, nobs = nobs))  #extract AICc                                      #
    Results$Delta_AICc <- Results$AICc - min(Results$AICc)            #compute delta AICc
    Results$ModelLik <- exp(-0.5*Results$Delta_AICc)                #compute model likelihood required to compute Akaike weights
    Results$AICcWt <- Results$ModelLik/sum(Results$ModelLik)        #compute Akaike weights

   ##check if some models are redundant
    if(length(unique(Results$AICc)) != length(cand.set)) warning("\nCheck model structure carefully as some models may be redundant\n")
    
    Results$LL <- unlist(lapply(X= cand.set, FUN = function(i) logLik(i)[1]))      

    ##rename correctly to AIC
    if(second.ord == FALSE) {
      colnames(Results)[1:6] <- c("Modnames", "K", "AIC", "Delta_AIC", "ModelLik", "AICWt")
    }

    if(sort)  {
      Results <- Results[rev(order(Results[, 6])),] 	  #if sort=TRUE, models are ranked based on Akaike weights
      Results$Cum.Wt <- cumsum(Results[, 6])                        #display cumulative sum of Akaike weights
    } else {Results$Cum.Wt <- NULL}
   
    
    class(Results) <- c("aictab", "data.frame")
    return(Results)
  }



##clmm
aictab.AICclmm <-
  function(cand.set, modnames = NULL, second.ord = TRUE, nobs = NULL, sort = TRUE, ...){  

    ##check if named list if modnames are not supplied
    if(is.null(modnames)) {
      if(is.null(names(cand.set))) {
        modnames <- paste("Mod", 1:length(cand.set), sep = "")
        warning("\nModel names have been supplied automatically in the table\n")
      } else {
        modnames <- names(cand.set)
      }
    }
    
    
    ##add check to see whether response variable is the same for all models
    check.resp <- lapply(X = cand.set, FUN = function(b) formula(b)[2])
    if(length(unique(check.resp)) > 1) stop("\nYou must use the same response variable for all models\n")
       
    Results <- data.frame(Modnames = modnames)                    #assign model names to first column
    Results$K <- unlist(lapply(cand.set, AICc, return.K = TRUE, second.ord = second.ord, nobs = nobs))     #extract number of parameters
    Results$AICc <- unlist(lapply(cand.set, AICc, return.K = FALSE, second.ord = second.ord, nobs = nobs))  #extract AICc                                      #
    Results$Delta_AICc <- Results$AICc - min(Results$AICc)            #compute delta AICc
    Results$ModelLik <- exp(-0.5*Results$Delta_AICc)                #compute model likelihood required to compute Akaike weights
    Results$AICcWt <- Results$ModelLik/sum(Results$ModelLik)        #compute Akaike weights

   ##check if some models are redundant
    if(length(unique(Results$AICc)) != length(cand.set)) warning("\nCheck model structure carefully as some models may be redundant\n")
    
    Results$LL <- unlist(lapply(X= cand.set, FUN = function(i) logLik(i)[1]))      

    ##rename correctly to AIC
    if(second.ord == FALSE) {
      colnames(Results)[1:6] <- c("Modnames", "K", "AIC", "Delta_AIC", "ModelLik", "AICWt")
    }

    if(sort)  {
      Results <- Results[rev(order(Results[, 6])),] 	  #if sort=TRUE, models are ranked based on Akaike weights
      Results$Cum.Wt <- cumsum(Results[, 6])                        #display cumulative sum of Akaike weights
    } else {Results$Cum.Wt <- NULL}
   
    
    class(Results) <- c("aictab", "data.frame")
    return(Results)
  }



##coxme
aictab.AICcoxme <- function(cand.set, modnames = NULL, second.ord = TRUE, nobs = NULL, sort = TRUE, ...){  

  ##check if named list if modnames are not supplied
    if(is.null(modnames)) {
      if(is.null(names(cand.set))) {
        modnames <- paste("Mod", 1:length(cand.set), sep = "")
        warning("\nModel names have been supplied automatically in the table\n")
      } else {
        modnames <- names(cand.set)
      }
    }

    
  ##add check to see whether response variable is the same for all models
  check.resp <- lapply(X = cand.set, FUN = function(b) formula(b)$fixed[2])
  if(length(unique(check.resp)) > 1) stop("\nYou must use the same response variable for all models\n")
    
  ##arrange in table
  Results <- data.frame(Modnames=modnames)                    #assign model names to first column
  Results$K <- unlist(lapply(cand.set, AICc, return.K = TRUE, second.ord = second.ord, nobs = nobs))     #extract number of parameters
  Results$AICc <- unlist(lapply(cand.set, AICc, return.K = FALSE, second.ord = second.ord, nobs = nobs))  #extract AICc                    
  Results$Delta_AICc <- Results$AICc - min(Results$AICc)            #compute delta AICc
  Results$ModelLik <- exp(-0.5*Results$Delta_AICc)                #compute model likelihood required to compute Akaike weights
  Results$AICcWt <- Results$ModelLik/sum(Results$ModelLik)        #compute Akaike weights
    
  ##check if some models are redundant
  if(length(unique(Results$AICc)) != length(cand.set)) warning("\nCheck model structure carefully as some models may be redundant\n")
  
  ##add partial log-likelihood column
  Results$LL <- unlist(lapply(X = cand.set, FUN=function(i) extractLL(i)[1]))      
  

  ##rename correctly to AIC
  if(second.ord==FALSE) {
    colnames(Results)[1:6] <- c("Modnames", "K", "AIC", "Delta_AIC", "ModelLik", "AICWt")
  }

  if(sort)  {
    Results <- Results[rev(order(Results[, 6])),] 	  #if sort=TRUE, models are ranked based on Akaike weights
    Results$Cum.Wt <- cumsum(Results[, 6])                        #display cumulative sum of Akaike weights
  } else {Results$Cum.Wt <- NULL}
    
  class(Results) <- c("aictab", "data.frame")
  return(Results)
}



##coxph and clogit
aictab.AICcoxph <- function(cand.set, modnames = NULL, second.ord = TRUE, nobs = NULL, sort = TRUE, ...){  

  ##check if named list if modnames are not supplied
    if(is.null(modnames)) {
      if(is.null(names(cand.set))) {
        modnames <- paste("Mod", 1:length(cand.set), sep = "")
        warning("\nModel names have been supplied automatically in the table\n")
      } else {
        modnames <- names(cand.set)
      }
    }

    
  ##add check to see whether response variable is the same for all models
  check.resp <- lapply(X = cand.set, FUN = function(b) formula(b)[2])
  if(length(unique(check.resp)) > 1) stop("\nYou must use the same response variable for all models\n")
    
  ##arrange in table
  Results <- data.frame(Modnames=modnames)                    #assign model names to first column
  Results$K <- unlist(lapply(cand.set, AICc, return.K = TRUE, second.ord = second.ord, nobs = nobs))     #extract number of parameters
  Results$AICc <- unlist(lapply(cand.set, AICc, return.K = FALSE, second.ord = second.ord, nobs = nobs))  #extract AICc                    
  Results$Delta_AICc <- Results$AICc - min(Results$AICc)            #compute delta AICc
  Results$ModelLik <- exp(-0.5*Results$Delta_AICc)                #compute model likelihood required to compute Akaike weights
  Results$AICcWt <- Results$ModelLik/sum(Results$ModelLik)        #compute Akaike weights
    
  ##check if some models are redundant
  if(length(unique(Results$AICc)) != length(cand.set)) warning("\nCheck model structure carefully as some models may be redundant\n")
  
  ##add partial log-likelihood column
  Results$LL <- unlist(lapply(X = cand.set, FUN = function(i) extractLL(i)[1]))      
  

  ##rename correctly to AIC
  if(second.ord == FALSE) {
    colnames(Results)[1:6] <- c("Modnames", "K", "AIC", "Delta_AIC", "ModelLik", "AICWt")
  }

  if(sort)  {
    Results <- Results[rev(order(Results[, 6])),] 	  #if sort=TRUE, models are ranked based on Akaike weights
    Results$Cum.Wt <- cumsum(Results[, 6])                        #display cumulative sum of Akaike weights
  } else {Results$Cum.Wt <- NULL}
    
  class(Results) <- c("aictab", "data.frame")
  return(Results)
}



##fitdist (from fitdistrplus)
aictab.AICfitdist <-
  function(cand.set, modnames = NULL, second.ord = TRUE, nobs = NULL, sort = TRUE, ...){

    ##check if named list if modnames are not supplied
    if(is.null(modnames)) {
      if(is.null(names(cand.set))) {
        modnames <- paste("Mod", 1:length(cand.set), sep = "")
        warning("\nModel names have been supplied automatically in the table\n")
      } else {
        modnames <- names(cand.set)
      }
    }
      

    ##add check to see whether response variable is the same for all models
    #check.resp <- lapply(X = cand.set, FUN = function(b) formula(b)[2])
    #if(length(unique(check.resp)) > 1) stop("\nYou must use the same response variable for all models\n")

    Results <- NULL
    Results <- data.frame(Modnames = modnames)                    #assign model names to first column
    Results$K <- unlist(lapply(X = cand.set, FUN = AICc, return.K = TRUE,
                               second.ord = second.ord, nobs = nobs))     #extract number of parameters
    Results$AICc <- unlist(lapply(X = cand.set, FUN = AICc, return.K = FALSE, 
                                  second.ord = second.ord, nobs = nobs))  #extract AICc                                      #
    Results$Delta_AICc <- Results$AICc - min(Results$AICc)            #compute delta AICc
    Results$ModelLik <- exp(-0.5*Results$Delta_AICc)                #compute model likelihood required to compute Akaike weights
    Results$AICcWt <- Results$ModelLik/sum(Results$ModelLik)        #compute Akaike weights

    ##check if some models are redundant
    if(length(unique(Results$AICc)) != length(cand.set)) warning("\nCheck model structure carefully as some models may be redundant\n")

    
    ##extract LL
    Results$LL <- unlist(lapply(X = cand.set, FUN = function(i) logLik(i)[1]))

    ##rename correctly to AIC
    if(second.ord == FALSE) {
      colnames(Results)[1:6] <- c("Modnames", "K", "AIC", "Delta_AIC", "ModelLik", "AICWt")
    }

    
    if(sort)  {
      Results <- Results[rev(order(Results[, 6])),] 	  #if sort=TRUE, models are ranked based on Akaike weights
      Results$Cum.Wt <- cumsum(Results[, 6])                        #display cumulative sum of Akaike weights
    } else {Results$Cum.Wt <- NULL}

    
    class(Results) <- c("aictab", "data.frame")
    return(Results)
  }



##fitdistr (from MASS)
aictab.AICfitdistr <-
  function(cand.set, modnames = NULL, second.ord = TRUE, nobs = NULL, sort = TRUE, ...){

    ##check if named list if modnames are not supplied
    if(is.null(modnames)) {
      if(is.null(names(cand.set))) {
        modnames <- paste("Mod", 1:length(cand.set), sep = "")
        warning("\nModel names have been supplied automatically in the table\n")
      } else {
        modnames <- names(cand.set)
      }
    }
      

    ##add check to see whether response variable is the same for all models
    #check.resp <- lapply(X = cand.set, FUN = function(b) formula(b)[2])
    #if(length(unique(check.resp)) > 1) stop("\nYou must use the same response variable for all models\n")

    Results <- NULL
    Results <- data.frame(Modnames = modnames)                    #assign model names to first column
    Results$K <- unlist(lapply(X = cand.set, FUN = AICc, return.K = TRUE,
                               second.ord = second.ord, nobs = nobs))     #extract number of parameters
    Results$AICc <- unlist(lapply(X = cand.set, FUN = AICc, return.K = FALSE, 
                                  second.ord = second.ord, nobs = nobs))  #extract AICc                                      #
    Results$Delta_AICc <- Results$AICc - min(Results$AICc)            #compute delta AICc
    Results$ModelLik <- exp(-0.5*Results$Delta_AICc)                #compute model likelihood required to compute Akaike weights
    Results$AICcWt <- Results$ModelLik/sum(Results$ModelLik)        #compute Akaike weights

    ##check if some models are redundant
    if(length(unique(Results$AICc)) != length(cand.set)) warning("\nCheck model structure carefully as some models may be redundant\n")

    
    ##extract LL
    Results$LL <- unlist(lapply(X = cand.set, FUN = function(i) logLik(i)[1]))

    ##rename correctly to AIC
    if(second.ord == FALSE) {
      colnames(Results)[1:6] <- c("Modnames", "K", "AIC", "Delta_AIC", "ModelLik", "AICWt")
    }

    
    if(sort)  {
      Results <- Results[rev(order(Results[, 6])),] 	  #if sort=TRUE, models are ranked based on Akaike weights
      Results$Cum.Wt <- cumsum(Results[, 6])                        #display cumulative sum of Akaike weights
    } else {Results$Cum.Wt <- NULL}

    
    class(Results) <- c("aictab", "data.frame")
    return(Results)
  }



##glm
aictab.AICglm.lm <-
  function(cand.set, modnames = NULL, second.ord = TRUE, nobs = NULL, sort = TRUE, c.hat = 1, ...){
    
    ##check if named list if modnames are not supplied
    if(is.null(modnames)) {
      if(is.null(names(cand.set))) {
        modnames <- paste("Mod", 1:length(cand.set), sep = "")
        warning("\nModel names have been supplied automatically in the table\n")
      } else {
        modnames <- names(cand.set)
      }
    }
    

    ##add check to see whether response variable is the same for all models
    check.resp <- lapply(X = cand.set, FUN = function(b) formula(b)[2])
    if(length(unique(check.resp)) > 1) stop("\nYou must use the same response variable for all models\n")

    Results <- NULL
    Results <- data.frame(Modnames = modnames)                    #assign model names to first column
    Results$K <- unlist(lapply(X = cand.set, FUN = AICc, return.K = TRUE, second.ord = second.ord,
                               nobs = nobs, c.hat = c.hat))     #extract number of parameters
    Results$AICc <- unlist(lapply(X = cand.set, FUN = AICc, return.K = FALSE, second.ord = second.ord,
                                  nobs = nobs, c.hat = c.hat))  #extract AICc                                      
    Results$Delta_AICc <- Results$AICc - min(Results$AICc)            #compute delta AICc
    Results$ModelLik <- exp(-0.5*Results$Delta_AICc)                #compute model likelihood required to compute Akaike weights
    Results$AICcWt <- Results$ModelLik/sum(Results$ModelLik)        #compute Akaike weights

    ##check if some models are redundant
    if(length(unique(Results$AICc)) != length(cand.set)) warning("\nCheck model structure carefully as some models may be redundant\n")
         
    ##check if AICc and c.hat = 1
    if(second.ord == TRUE && c.hat == 1) {
      Results$LL <- unlist(lapply(X = cand.set, FUN=function(i) logLik(i)[1]))
    }
    
    ##rename correctly to QAICc and add column for c-hat
    if(second.ord == TRUE && c.hat > 1) {
      colnames(Results) <- c("Modnames", "K", "QAICc", "Delta_QAICc", "ModelLik", "QAICcWt")
      LL <- unlist(lapply(X=cand.set, FUN = function(i) logLik(i)[1]))
      Results$Quasi.LL <- LL/c.hat
      Results$c_hat <- c.hat
    }      

    ##rename correctly to AIC
    if(second.ord == FALSE && c.hat == 1) {
      colnames(Results)<-c("Modnames", "K", "AIC", "Delta_AIC", "ModelLik", "AICWt")
      Results$LL <- unlist(lapply(X=cand.set, FUN=function(i) logLik(i)[1]))      
    }  

    ##rename correctly to QAIC and add column for c-hat
    if(second.ord == FALSE && c.hat > 1) {
      colnames(Results)<-c("Modnames", "K", "QAIC", "Delta_QAIC", "ModelLik", "QAICWt")
      LL <- unlist(lapply(X=cand.set, FUN=function(i) logLik(i)[1]))
      Results$Quasi.LL <- LL/c.hat
      Results$c_hat<-c.hat
    }

    
    if(sort)  {
      Results <- Results[rev(order(Results[, 6])),] 	  #if sort=TRUE, models are ranked based on Akaike weights
      Results$Cum.Wt <- cumsum(Results[, 6])                        #display cumulative sum of Akaike weights
    } else {Results$Cum.Wt <- NULL}

    
    class(Results) <- c("aictab", "data.frame")
    return(Results)
  }



##gls
aictab.AICgls <-
  function(cand.set, modnames = NULL, second.ord = TRUE, nobs = NULL, sort = TRUE, ...){  #specify whether table should be sorted or not by delta AICc

    ##check if named list if modnames are not supplied
    if(is.null(modnames)) {
      if(is.null(names(cand.set))) {
        modnames <- paste("Mod", 1:length(cand.set), sep = "")
        warning("\nModel names have been supplied automatically in the table\n")
      } else {
        modnames <- names(cand.set)
      }
    }
    
    ##add check to see whether response variable is the same for all models
    check.resp <- lapply(X = cand.set, FUN = function(b) formula(b)[2])
    if(length(unique(check.resp)) > 1) stop("\nYou must use the same response variable for all models\n")
       
    Results <- NULL
    ##check if models were fit with same method (REML or ML)
    check_ML <- unlist(lapply(cand.set, FUN = function(i) i$method))
    
    if (any(check_ML != "ML")) {
      warning("\nModel selection for fixed effects is only appropriate with method='ML':", "\n",
              "REML (default) should only be used to select random effects for a constant set of fixed effects\n")
    }
    
    check.method <- unique(check_ML)
    
    if(identical(check.method, c("ML", "REML"))) {
      stop("\nYou should not have models fit with REML and ML in the same candidate model set\n")
    }
    
    Results <- data.frame(Modnames = modnames)                    #assign model names to first column
    Results$K <- unlist(lapply(cand.set, AICc, return.K = TRUE, second.ord = second.ord, nobs = nobs))     #extract number of parameters
    Results$AICc <- unlist(lapply(cand.set, AICc, return.K = FALSE, second.ord = second.ord, nobs = nobs))  #extract AICc                                      #
    Results$Delta_AICc <- Results$AICc - min(Results$AICc)            #compute delta AICc
    Results$ModelLik <- exp(-0.5*Results$Delta_AICc)                #compute model likelihood required to compute Akaike weights
    Results$AICcWt <- Results$ModelLik/sum(Results$ModelLik)        #compute Akaike weights

    ##check if some models are redundant
    if(length(unique(Results$AICc)) != length(cand.set)) warning("\nCheck model structure carefully as some models may be redundant\n")
    
    #check if ML or REML used and add column accordingly
    if(identical(check.method, "ML")) {
      Results$LL <- unlist(lapply(X = cand.set, FUN=function(i) logLik(i)[1]))      
    }

    if(identical(check.method, "REML")) {
      Results$Res.LL <- unlist(lapply(X = cand.set, FUN = function(i) logLik(i)[1]))      
    }

    
    #rename correctly to AIC
    if(second.ord == FALSE) {
      colnames(Results)[1:6] <- c("Modnames", "K", "AIC", "Delta_AIC", "ModelLik", "AICWt")
    }

    if(sort)  {
      Results <- Results[rev(order(Results[, 6])),] 	  #if sort=TRUE, models are ranked based on Akaike weights
      Results$Cum.Wt <- cumsum(Results[, 6])                        #display cumulative sum of Akaike weights
    } else {Results$Cum.Wt <- NULL}
   
    
    class(Results) <- c("aictab", "data.frame")
    return(Results)
  }



##gnls
aictab.AICgnls.gls <-
  function(cand.set, modnames = NULL, second.ord = TRUE, nobs = NULL, sort = TRUE, ...){  #specify whether table should be sorted or not by delta AICc

    ##check if named list if modnames are not supplied
    if(is.null(modnames)) {
      if(is.null(names(cand.set))) {
        modnames <- paste("Mod", 1:length(cand.set), sep = "")
        warning("\nModel names have been supplied automatically in the table\n")
      } else {
        modnames <- names(cand.set)
      }
    }
    
    ##add check to see whether response variable is the same for all models
    check.resp <- lapply(X = cand.set, FUN = function(b) formula(b)[2])
    if(length(unique(check.resp)) > 1) stop("\nYou must use the same response variable for all models\n")
       
    Results <- data.frame(Modnames=modnames)                    #assign model names to first column
    Results$K <- unlist(lapply(cand.set, AICc, return.K = TRUE, second.ord = second.ord, nobs = nobs))     #extract number of parameters
    Results$AICc <- unlist(lapply(cand.set, AICc, return.K = FALSE, second.ord = second.ord, nobs = nobs))  #extract AICc                                      #
    Results$Delta_AICc <- Results$AICc - min(Results$AICc)            #compute delta AICc
    Results$ModelLik <- exp(-0.5*Results$Delta_AICc)                #compute model likelihood required to compute Akaike weights
    Results$AICcWt <- Results$ModelLik/sum(Results$ModelLik)        #compute Akaike weights
    Results$LL <- unlist(lapply(X = cand.set, FUN = function(i) logLik(i)[1]))
    
    ##check if some models are redundant
    if(length(unique(Results$AICc)) != length(cand.set)) warning("\nCheck model structure carefully as some models may be redundant\n")
    
    ##rename correctly to AIC
    if(second.ord == FALSE) {
      colnames(Results)[1:6] <- c("Modnames", "K", "AIC", "Delta_AIC", "ModelLik", "AICWt")
    }

    if(sort)  {
      Results <- Results[rev(order(Results[, 6])),] 	  #if sort=TRUE, models are ranked based on Akaike weights
      Results$Cum.Wt <- cumsum(Results[, 6])                        #display cumulative sum of Akaike weights
    } else {Results$Cum.Wt <- NULL}
   
    
    class(Results) <- c("aictab", "data.frame")
    return(Results)
  }



##hurdle
aictab.AIChurdle <-
  function(cand.set, modnames = NULL, second.ord = TRUE, nobs = NULL, sort = TRUE, ...){  #specify whether table should be sorted or not by delta AICc

    ##check if named list if modnames are not supplied
    if(is.null(modnames)) {
      if(is.null(names(cand.set))) {
        modnames <- paste("Mod", 1:length(cand.set), sep = "")
        warning("\nModel names have been supplied automatically in the table\n")
      } else {
        modnames <- names(cand.set)
      }
    }
    
    #if(c.hat != 1) stop("\nThis function does not support overdispersion in \'zeroinfl\' models\n")
    
    ##add check to see whether response variable is the same for all models
    check.resp <- lapply(X = cand.set, FUN = function(b) formula(b)[2])                                        
    if(length(unique(check.resp)) > 1) stop("\nYou must use the same data set for all models\n")
       
    Results <- data.frame(Modnames = modnames)                    #assign model names to first column
    Results$K <- unlist(lapply(cand.set, AICc, return.K = TRUE, second.ord = second.ord, nobs = nobs))     #extract number of parameters
    Results$AICc <- unlist(lapply(cand.set, AICc, return.K = FALSE, second.ord = second.ord, nobs = nobs))  #extract AICc                                      #
    Results$Delta_AICc <- Results$AICc - min(Results$AICc)            #compute delta AICc
    Results$ModelLik <- exp(-0.5*Results$Delta_AICc)                #compute model likelihood required to compute Akaike weights
    Results$AICcWt <- Results$ModelLik/sum(Results$ModelLik)        #compute Akaike weights

   ##check if some models are redundant
    if(length(unique(Results$AICc)) != length(cand.set)) warning("\nCheck model structure carefully as some models may be redundant\n")
    
    Results$LL <- unlist(lapply(X= cand.set, FUN = function(i) logLik(i)))      

    ##rename correctly to AIC
    if(second.ord == FALSE) {
      colnames(Results)[1:6] <- c("Modnames", "K", "AIC", "Delta_AIC", "ModelLik", "AICWt")
    }

    if(sort)  {
      Results <- Results[rev(order(Results[, 6])),] 	  #if sort=TRUE, models are ranked based on Akaike weights
      Results$Cum.Wt <- cumsum(Results[, 6])                        #display cumulative sum of Akaike weights
    } else {Results$Cum.Wt <- NULL}
   
    
    class(Results) <- c("aictab", "data.frame")
    return(Results)
  }



##lavaan
aictab.AIClavaan <-
  function(cand.set, modnames = NULL, second.ord = TRUE, nobs = NULL, sort = TRUE, ...){

    ##check if named list if modnames are not supplied
    if(is.null(modnames)) {
      if(is.null(names(cand.set))) {
        modnames <- paste("Mod", 1:length(cand.set), sep = "")
        warning("\nModel names have been supplied automatically in the table\n")
      } else {
        modnames <- names(cand.set)
      }
    }
      

    ##add check to see whether observed variables are the same for all models
    check.obs <- unlist(lapply(X = cand.set, FUN = function(b) b@Data@ov.names[[1]]))
    ##frequency of each observed variable
    freq.obs <- table(check.obs)
    if(length(unique(freq.obs)) > 1) stop("\nModels with different sets of observed variables are not directly comparable\n")

    Results <- NULL
    Results <- data.frame(Modnames = modnames)                    #assign model names to first column
    Results$K <- unlist(lapply(X = cand.set, FUN = AICc, return.K = TRUE,
                               second.ord = second.ord, nobs = nobs))     #extract number of parameters
    Results$AICc <- unlist(lapply(X = cand.set, FUN = AICc, return.K = FALSE, 
                                  second.ord = second.ord, nobs = nobs))  #extract AICc                                      #
    Results$Delta_AICc <- Results$AICc - min(Results$AICc)            #compute delta AICc
    Results$ModelLik <- exp(-0.5*Results$Delta_AICc)                #compute model likelihood required to compute Akaike weights
    Results$AICcWt <- Results$ModelLik/sum(Results$ModelLik)        #compute Akaike weights

    ##check if some models are redundant
    if(length(unique(Results$AICc)) != length(cand.set)) warning("\nCheck model structure carefully as some models may be redundant\n")

    
    ##extract LL
    Results$LL <- unlist(lapply(X = cand.set, FUN = function(i) logLik(i)[1]))

    ##rename correctly to AIC
    if(second.ord == FALSE) {
      colnames(Results)[1:6] <- c("Modnames", "K", "AIC", "Delta_AIC", "ModelLik", "AICWt")
    }

    
    if(sort)  {
      Results <- Results[rev(order(Results[, 6])),] 	  #if sort=TRUE, models are ranked based on Akaike weights
      Results$Cum.Wt <- cumsum(Results[, 6])                        #display cumulative sum of Akaike weights
    } else {Results$Cum.Wt <- NULL}

    
    class(Results) <- c("aictab", "data.frame")
    return(Results)
  }



##lm
aictab.AIClm <-
  function(cand.set, modnames = NULL, second.ord = TRUE, nobs = NULL, sort = TRUE, ...){

    ##check if named list if modnames are not supplied
    if(is.null(modnames)) {
      if(is.null(names(cand.set))) {
        modnames <- paste("Mod", 1:length(cand.set), sep = "")
        warning("\nModel names have been supplied automatically in the table\n")
      } else {
        modnames <- names(cand.set)
      }
    }
      

    ##add check to see whether response variable is the same for all models
    check.resp <- lapply(X = cand.set, FUN = function(b) formula(b)[2])
    if(length(unique(check.resp)) > 1) stop("\nYou must use the same response variable for all models\n")

    Results <- NULL
    Results <- data.frame(Modnames = modnames)                    #assign model names to first column
    Results$K <- unlist(lapply(X = cand.set, FUN = AICc, return.K = TRUE,
                               second.ord = second.ord, nobs = nobs))     #extract number of parameters
    Results$AICc <- unlist(lapply(X = cand.set, FUN = AICc, return.K = FALSE, 
                                  second.ord = second.ord, nobs = nobs))  #extract AICc                                      #
    Results$Delta_AICc <- Results$AICc - min(Results$AICc)            #compute delta AICc
    Results$ModelLik <- exp(-0.5*Results$Delta_AICc)                #compute model likelihood required to compute Akaike weights
    Results$AICcWt <- Results$ModelLik/sum(Results$ModelLik)        #compute Akaike weights

    ##check if some models are redundant
    if(length(unique(Results$AICc)) != length(cand.set)) warning("\nCheck model structure carefully as some models may be redundant\n")
         
    ##extract LL
    Results$LL <- unlist(lapply(X = cand.set, FUN = function(i) logLik(i)[1]))

    ##rename correctly to AIC
    if(second.ord == FALSE) {
      colnames(Results)[1:6] <- c("Modnames", "K", "AIC", "Delta_AIC", "ModelLik", "AICWt")
    }  


    if(sort)  {
      Results <- Results[rev(order(Results[, 6])),] 	  #if sort=TRUE, models are ranked based on Akaike weights
      Results$Cum.Wt <- cumsum(Results[, 6])                        #display cumulative sum of Akaike weights
    } else {Results$Cum.Wt <- NULL}

    
    class(Results) <- c("aictab", "data.frame")
    return(Results)
  }



##lme
aictab.AIClme <-
  function(cand.set, modnames = NULL, second.ord = TRUE, nobs = NULL, sort = TRUE, ...){  #specify whether table should be sorted or not by delta AICc

    ##check if named list if modnames are not supplied
    if(is.null(modnames)) {
      if(is.null(names(cand.set))) {
        modnames <- paste("Mod", 1:length(cand.set), sep = "")
        warning("\nModel names have been supplied automatically in the table\n")
      } else {
        modnames <- names(cand.set)
      }
    }
    
    ##add check to see whether response variable is the same for all models
    check.resp <- lapply(X = cand.set, FUN = function(b) formula(b)[2])
    if(length(unique(check.resp)) > 1) stop("\nYou must use the same response variable for all models\n")
       
    Results <- NULL
#check if models were fit with same method (REML or ML)
    check_ML <- unlist(lapply(cand.set, FUN = function(i) i$method))
    
    if (any(check_ML != "ML")) {
      warning("\nModel selection for fixed effects is only appropriate with method='ML':", "\n",
                    "REML (default) should only be used to select random effects for a constant set of fixed effects\n")
    }
    
    check.method <- unique(check_ML)
    
    if(identical(check.method, c("ML", "REML"))) {
      stop("\nYou should not have models fit with REML and ML in the same candidate model set\n")
    }
    
    Results <- data.frame(Modnames = modnames)                    #assign model names to first column
    Results$K <- unlist(lapply(cand.set, AICc, return.K = TRUE, second.ord = second.ord, nobs = nobs))     #extract number of parameters
    Results$AICc <- unlist(lapply(cand.set, AICc, return.K = FALSE, second.ord = second.ord, nobs = nobs))  #extract AICc                          
    Results$Delta_AICc <- Results$AICc-min(Results$AICc)            #compute delta AICc
    Results$ModelLik <- exp(-0.5*Results$Delta_AICc)                #compute model likelihood required to compute Akaike weights
    Results$AICcWt <- Results$ModelLik/sum(Results$ModelLik)        #compute Akaike weights

   ##check if some models are redundant
    if(length(unique(Results$AICc)) != length(cand.set)) warning("\nCheck model structure carefully as some models may be redundant\n")
    
    #check if ML or REML used and add column accordingly
    if(identical(check.method, "ML")) {
      Results$LL <- unlist(lapply(X= cand.set, FUN = function(i) logLik(i)[1]))      
    }

    if(identical(check.method, "REML")) {
      Results$Res.LL <- unlist(lapply(X = cand.set, FUN = function(i) logLik(i)[1]))      
    }

    
    #rename correctly to AIC
    if(second.ord == FALSE) {
      colnames(Results)[1:6] <- c("Modnames", "K", "AIC", "Delta_AIC", "ModelLik", "AICWt")
    }

    if(sort)  {
      Results <- Results[rev(order(Results[, 6])),] 	  #if sort=TRUE, models are ranked based on Akaike weights
      Results$Cum.Wt <- cumsum(Results[, 6])                        #display cumulative sum of Akaike weights
    } else {Results$Cum.Wt <- NULL}
   
    
    class(Results) <- c("aictab", "data.frame")
    return(Results)
  }



##lmekin
aictab.AIClmekin <-
  function(cand.set, modnames = NULL, second.ord = TRUE, nobs = NULL, sort = TRUE, ...){  #specify whether table should be sorted or not by delta AICc

    ##check if named list if modnames are not supplied
    if(is.null(modnames)) {
      if(is.null(names(cand.set))) {
        modnames <- paste("Mod", 1:length(cand.set), sep = "")
        warning("\nModel names have been supplied automatically in the table\n")
      } else {
        modnames <- names(cand.set)
      }
    }
    
    ##add check to see whether response variable is the same for all models
    check.resp <- lapply(X = cand.set, FUN = function(b) formula(b)[2])
    if(length(unique(check.resp)) > 1) stop("\nYou must use the same response variable for all models\n")
       
    Results <- NULL
#check if models were fit with same method (REML or ML)
    check_ML <- unlist(lapply(cand.set, FUN = function(i) i$method))
    
    if (any(check_ML != "ML")) {
      warning("\nModel selection for fixed effects is only appropriate with method='ML':", "\n",
                    "REML (default) should only be used to select random effects for a constant set of fixed effects\n")
    }
    
    check.method <- unique(check_ML)
    
    if(identical(check.method, c("ML", "REML"))) {
      stop("\nYou should not have models fit with REML and ML in the same candidate model set\n")
    }
    
    Results <- data.frame(Modnames = modnames)                    #assign model names to first column
    Results$K <- unlist(lapply(cand.set, AICc, return.K = TRUE, second.ord = second.ord, nobs = nobs))     #extract number of parameters
    Results$AICc <- unlist(lapply(cand.set, AICc, return.K = FALSE, second.ord = second.ord, nobs = nobs))  #extract AICc                                      #
    Results$Delta_AICc <- Results$AICc-min(Results$AICc)            #compute delta AICc
    Results$ModelLik <- exp(-0.5*Results$Delta_AICc)                #compute model likelihood required to compute Akaike weights
    Results$AICcWt <- Results$ModelLik/sum(Results$ModelLik)        #compute Akaike weights

   ##check if some models are redundant
    if(length(unique(Results$AICc)) != length(cand.set)) warning("\nCheck model structure carefully as some models may be redundant\n")
    
    #check if ML or REML used and add column accordingly
    if(identical(check.method, "ML")) {
      Results$LL <- unlist(lapply(X= cand.set, FUN = function(i) extractLL(i)[1]))      
    }

    if(identical(check.method, "REML")) {
      Results$Res.LL <- unlist(lapply(X = cand.set, FUN = function(i) extractLL(i)[1]))      
    }

    
    #rename correctly to AIC
    if(second.ord == FALSE) {
      colnames(Results)[1:6] <- c("Modnames", "K", "AIC", "Delta_AIC", "ModelLik", "AICWt")
    }

    if(sort)  {
      Results <- Results[rev(order(Results[, 6])),] 	  #if sort=TRUE, models are ranked based on Akaike weights
      Results$Cum.Wt <- cumsum(Results[, 6])                        #display cumulative sum of Akaike weights
    } else {Results$Cum.Wt <- NULL}
   
    
    class(Results) <- c("aictab", "data.frame")
    return(Results)
  }



##maxlike
aictab.AICmaxlikeFit.list <-
  function(cand.set, modnames = NULL, second.ord = TRUE, nobs = NULL, sort = TRUE, c.hat = 1, ...){  #specify whether table should be sorted or not by delta AICc

    ##check if named list if modnames are not supplied
    if(is.null(modnames)) {
      if(is.null(names(cand.set))) {
        modnames <- paste("Mod", 1:length(cand.set), sep = "")
        warning("\nModel names have been supplied automatically in the table\n")
      } else {
        modnames <- names(cand.set)
      }
    }
    
    if(c.hat != 1) stop("\nThis function does not support overdispersion in \'maxlikeFit\' models\n")
    
    ##add check to see whether response variable is the same for all models
    #check.resp <- lapply(X = cand.set, FUN = function(b) nrow(b$points.retained))
    #if(length(unique(check.resp)) > 1) stop("\nYou must use the same data set for all models\n")
       
    Results <- data.frame(Modnames = modnames)                    #assign model names to first column
    Results$K <- unlist(lapply(cand.set, AICc, return.K = TRUE, second.ord = second.ord, nobs = nobs))     #extract number of parameters
    Results$AICc <- unlist(lapply(cand.set, AICc, return.K = FALSE, second.ord = second.ord, nobs = nobs))  #extract AICc                          
    Results$Delta_AICc <- Results$AICc - min(Results$AICc)            #compute delta AICc
    Results$ModelLik <- exp(-0.5*Results$Delta_AICc)                #compute model likelihood required to compute Akaike weights
    Results$AICcWt <- Results$ModelLik/sum(Results$ModelLik)        #compute Akaike weights

   ##check if some models are redundant
    if(length(unique(Results$AICc)) != length(cand.set)) warning("\nCheck model structure carefully as some models may be redundant\n")
    
    Results$LL <- unlist(lapply(X= cand.set, FUN = function(i) logLik(i)))      

    ##rename correctly to AIC
    if(second.ord == FALSE) {
      colnames(Results)[1:6] <- c("Modnames", "K", "AIC", "Delta_AIC", "ModelLik", "AICWt")
    }

    if(sort)  {
      Results <- Results[rev(order(Results[, 6])),] 	  #if sort=TRUE, models are ranked based on Akaike weights
      Results$Cum.Wt <- cumsum(Results[, 6])                        #display cumulative sum of Akaike weights
    } else {Results$Cum.Wt <- NULL}
   
    
    class(Results) <- c("aictab", "data.frame")
    return(Results)
  }



##mer - lme4 version < 1
aictab.AICmer <-
  function(cand.set, modnames = NULL, second.ord = TRUE, nobs = NULL, sort = TRUE, ...){  #specify whether table should be sorted or not by delta AICc
    ##check if named list if modnames are not supplied
    if (is.null(modnames)) {
      if (is.null(names(cand.set))) {
        modnames <- paste("Mod", 1:length(cand.set), sep = "")
        warning("\nModel names have been supplied automatically in the table\n")
      } else {
        modnames <- names(cand.set)
      }
    }
    
    ##add check to see whether response variable is the same for all models
    check.resp <- lapply(X = cand.set, FUN = function(b) formula(b)[2])
    if(length(unique(check.resp)) > 1) stop("\nYou must use the same response variable for all models\n")
       
    Results <- NULL
#check if models were fit with same method (REML or ML)
    check_bin <- unlist(lapply(cand.set, FUN = function(i) i@dims["REML"]))
    check_ML <- ifelse(check_bin == 1, "REML", "ML")
    
    if (any(check_ML != "ML")) {
      warning("\nModel selection for fixed effects is only appropriate with ML estimation:", "\n",
                    "REML (default) should only be used to select random effects for a constant set of fixed effects\n")
    }
    
    check.method <- unique(check_ML)
    
    if(length(check.method) > 1) {
      stop("\nYou should not have models fit with REML and ML in the same candidate model set\n")
    }
    
    Results <- data.frame(Modnames = modnames)                    #assign model names to first column
    Results$K <- unlist(lapply(cand.set, AICc, return.K = TRUE, second.ord = second.ord, nobs = nobs))     #extract number of parameters
    Results$AICc <- unlist(lapply(cand.set, AICc, return.K = FALSE, second.ord = second.ord, nobs = nobs))  #extract AICc                                      #
    Results$Delta_AICc <- Results$AICc - min(Results$AICc)            #compute delta AICc
    Results$ModelLik <- exp(-0.5*Results$Delta_AICc)                #compute model likelihood required to compute Akaike weights
    Results$AICcWt <- Results$ModelLik/sum(Results$ModelLik)        #compute Akaike weights

    ##check if some models are redundant
    if(length(unique(Results$AICc)) != length(cand.set)) warning("\nCheck model structure carefully as some models may be redundant\n") 
    
    #check if ML or REML used and add column accordingly
    if(identical(check.method, "ML")) {
      Results$LL <- unlist(lapply(X = cand.set, FUN=function(i) logLik(i)[1]))      
    }

    if(identical(check.method, "REML")) {
      Results$Res.LL <- unlist(lapply(X = cand.set, FUN=function(i) logLik(i)[1]))      
    }

    
    #rename correctly to AIC
    if(second.ord == FALSE) {
      colnames(Results)[1:6] <- c("Modnames", "K", "AIC", "Delta_AIC", "ModelLik", "AICWt")
    }

    if(sort)  {
      Results <- Results[rev(order(Results[, 6])),] 	  #if sort=TRUE, models are ranked based on Akaike weights
      Results$Cum.Wt <- cumsum(Results[, 6])                        #display cumulative sum of Akaike weights
    } else {Results$Cum.Wt <- NULL}
   
    
    class(Results) <- c("aictab", "data.frame")
    return(Results)
  }



##lmerMod
aictab.AIClmerMod <-
  function(cand.set, modnames = NULL, second.ord = TRUE, nobs = NULL, sort = TRUE, ...){  #specify whether table should be sorted or not by delta AICc

    ##check if named list if modnames are not supplied
    if(is.null(modnames)) {
      if(is.null(names(cand.set))) {
        modnames <- paste("Mod", 1:length(cand.set), sep = "")
        warning("\nModel names have been supplied automatically in the table\n")
      } else {
        modnames <- names(cand.set)
      }
    }
    
    ##check for subclass of object
    sub.class <- lapply(X = cand.set, FUN = class)

    ##add check to see whether response variable is the same for all models
    check.resp <- lapply(X = cand.set, FUN = function(b) formula(b)[2])
    if(length(unique(check.resp)) > 1) stop("\nYou must use the same response variable for all models\n")

    Results <- NULL
    ##check if models were fit with same method (REML or ML)
    check_REML <- unlist(lapply(cand.set, FUN = function(i) isREML(i)))
    check_ML <- ifelse(check_REML, "REML", "ML")
    
    if (any(check_REML)) {
      warning("\nModel selection for fixed effects is only appropriate with ML estimation:", "\n",
                    "REML (default) should only be used to select random effects for a constant set of fixed effects\n")
    }

    check.method <- unique(check_ML)
    if(length(check.method) > 1) {
      stop("\nYou should not have models fit with REML and ML in the same candidate model set\n")
    }
    
    Results <- data.frame(Modnames = modnames)                    #assign model names to first column
    Results$K <- unlist(lapply(cand.set, AICc, return.K = TRUE, second.ord = second.ord, nobs = nobs))     #extract number of parameters
    Results$AICc <- unlist(lapply(cand.set, AICc, return.K = FALSE, second.ord = second.ord, nobs = nobs))  #extract AICc                          
    Results$Delta_AICc <- Results$AICc - min(Results$AICc)            #compute delta AICc
    Results$ModelLik <- exp(-0.5*Results$Delta_AICc)                #compute model likelihood required to compute Akaike weights
    Results$AICcWt <- Results$ModelLik/sum(Results$ModelLik)        #compute Akaike weights

    ##check if some models are redundant
    if(length(unique(Results$AICc)) != length(cand.set)) warning("\nCheck model structure carefully as some models may be redundant\n") 
    
    #check if ML or REML used and add column accordingly
    if(identical(check.method, "ML")) {
      Results$LL <- unlist(lapply(X = cand.set, FUN = function(i) logLik(i)[1]))      
    }

    if(identical(check.method, "REML")) {
      Results$Res.LL <- unlist(lapply(X = cand.set, FUN = function(i) logLik(i)[1]))      
    }

    
    #rename correctly to AIC
    if(second.ord == FALSE) {
      colnames(Results)[1:6] <- c("Modnames", "K", "AIC", "Delta_AIC", "ModelLik", "AICWt")
    }

    if(sort)  {
      Results <- Results[rev(order(Results[, 6])),] 	  #if sort=TRUE, models are ranked based on Akaike weights
      Results$Cum.Wt <- cumsum(Results[, 6])                        #display cumulative sum of Akaike weights
    } else {Results$Cum.Wt <- NULL}
   
    
    class(Results) <- c("aictab", "data.frame")
    return(Results)
  }



##glmerMod
aictab.AICglmerMod <-
  function(cand.set, modnames = NULL, second.ord = TRUE, nobs = NULL, sort = TRUE, ...){  #specify whether table should be sorted or not by delta AICc

    ##check if named list if modnames are not supplied
    if(is.null(modnames)) {
      if(is.null(names(cand.set))) {
        modnames <- paste("Mod", 1:length(cand.set), sep = "")
        warning("\nModel names have been supplied automatically in the table\n")
      } else {
        modnames <- names(cand.set)
      }
    }
    
        
    ##add check to see whether response variable is the same for all models
    check.resp <- lapply(X = cand.set, FUN = function(b) formula(b)[2])
    if(length(unique(check.resp)) > 1) stop("\nYou must use the same response variable for all models\n")

    
    Results <- data.frame(Modnames = modnames)                    #assign model names to first column
    Results$K <- unlist(lapply(cand.set, AICc, return.K = TRUE, second.ord = second.ord, nobs = nobs))     #extract number of parameters
    Results$AICc <- unlist(lapply(cand.set, AICc, return.K = FALSE, second.ord = second.ord, nobs = nobs))  #extract AICc                          
    Results$Delta_AICc <- Results$AICc - min(Results$AICc)            #compute delta AICc
    Results$ModelLik <- exp(-0.5*Results$Delta_AICc)                #compute model likelihood required to compute Akaike weights
    Results$AICcWt <- Results$ModelLik/sum(Results$ModelLik)        #compute Akaike weights

    ##check if some models are redundant
    if(length(unique(Results$AICc)) != length(cand.set)) warning("\nCheck model structure carefully as some models may be redundant\n") 
    
    ##add LL
    Results$LL <- unlist(lapply(X = cand.set, FUN = function(i) logLik(i)[1]))      
    
    
    ##rename correctly to AIC
    if(second.ord == FALSE) {
      colnames(Results)[1:6] <- c("Modnames", "K", "AIC", "Delta_AIC", "ModelLik", "AICWt")
    }

    if(sort)  {
      Results <- Results[rev(order(Results[, 6])),] 	  #if sort=TRUE, models are ranked based on Akaike weights
      Results$Cum.Wt <- cumsum(Results[, 6])                        #display cumulative sum of Akaike weights
    } else {Results$Cum.Wt <- NULL}
   
    
    class(Results) <- c("aictab", "data.frame")
    return(Results)
  }



##nlme
aictab.AICnlme.lme <-
  function(cand.set, modnames = NULL, second.ord = TRUE, nobs = NULL, sort = TRUE, ...){  #specify whether table should be sorted or not by delta AICc

    ##check if named list if modnames are not supplied
    if(is.null(modnames)) {
      if(is.null(names(cand.set))) {
        modnames <- paste("Mod", 1:length(cand.set), sep = "")
        warning("\nModel names have been supplied automatically in the table\n")
      } else {
        modnames <- names(cand.set)
      }
    }
    
    ##add check to see whether response variable is the same for all models
    check.resp <- lapply(X = cand.set, FUN = function(b) formula(b)[2])
    if(length(unique(check.resp)) > 1) stop("\nYou must use the same response variable for all models\n")
       
    Results <- NULL
#check if models were fit with same method (REML or ML)
    check_ML <- unlist(lapply(cand.set, FUN = function(i) i$method))
    
    if (any(check_ML != "ML")) {
      warning("\nModel selection for fixed effects is only appropriate with method='ML':", "\n",
                    "REML (default) should only be used to select random effects for a constant set of fixed effects\n")
    }
    
    check.method <- unique(check_ML)
    
    if(identical(check.method, c("ML", "REML"))) {
      stop("\nYou should not have models fit with REML and ML in the same candidate model set\n")
    }
    
    Results <- data.frame(Modnames = modnames)                    #assign model names to first column
    Results$K <- unlist(lapply(cand.set, AICc, return.K = TRUE, second.ord = second.ord, nobs = nobs))     #extract number of parameters
    Results$AICc <- unlist(lapply(cand.set, AICc, return.K = FALSE, second.ord = second.ord, nobs = nobs))  #extract AICc                          
    Results$Delta_AICc <- Results$AICc-min(Results$AICc)            #compute delta AICc
    Results$ModelLik <- exp(-0.5*Results$Delta_AICc)                #compute model likelihood required to compute Akaike weights
    Results$AICcWt <- Results$ModelLik/sum(Results$ModelLik)        #compute Akaike weights

   ##check if some models are redundant
    if(length(unique(Results$AICc)) != length(cand.set)) warning("\nCheck model structure carefully as some models may be redundant\n")
    
    #check if ML or REML used and add column accordingly
    if(identical(check.method, "ML")) {
      Results$LL <- unlist(lapply(X= cand.set, FUN = function(i) logLik(i)[1]))      
    }

    if(identical(check.method, "REML")) {
      Results$Res.LL <- unlist(lapply(X = cand.set, FUN = function(i) logLik(i)[1]))      
    }

    
    #rename correctly to AIC
    if(second.ord == FALSE) {
      colnames(Results)[1:6] <- c("Modnames", "K", "AIC", "Delta_AIC", "ModelLik", "AICWt")
    }

    if(sort)  {
      Results <- Results[rev(order(Results[, 6])),] 	  #if sort=TRUE, models are ranked based on Akaike weights
      Results$Cum.Wt <- cumsum(Results[, 6])                        #display cumulative sum of Akaike weights
    } else {Results$Cum.Wt <- NULL}
   
    
    class(Results) <- c("aictab", "data.frame")
    return(Results)
  }



##nlmerMod
aictab.AICnlmerMod <-
  function(cand.set, modnames = NULL, second.ord = TRUE, nobs = NULL, sort = TRUE, ...){  #specify whether table should be sorted or not by delta AICc

    ##check if named list if modnames are not supplied
    if(is.null(modnames)) {
      if(is.null(names(cand.set))) {
        modnames <- paste("Mod", 1:length(cand.set), sep = "")
        warning("\nModel names have been supplied automatically in the table\n")
      } else {
        modnames <- names(cand.set)
      }
    }
    
    ##add check to see whether response variable is the same for all models
    check.resp <- lapply(X = cand.set, FUN = function(b) unlist(strsplit(x = as.character(formula(b)), split = "~"))[2])
    if(length(unique(check.resp)) > 1) stop("\nYou must use the same response variable for all models\n")

    Results <- data.frame(Modnames = modnames)                    #assign model names to first column
    Results$K <- unlist(lapply(cand.set, AICc, return.K = TRUE, second.ord = second.ord, nobs = nobs))     #extract number of parameters
    Results$AICc <- unlist(lapply(cand.set, AICc, return.K = FALSE, second.ord = second.ord, nobs = nobs))  #extract AICc                          
    Results$Delta_AICc <- Results$AICc - min(Results$AICc)            #compute delta AICc
    Results$ModelLik <- exp(-0.5*Results$Delta_AICc)                #compute model likelihood required to compute Akaike weights
    Results$AICcWt <- Results$ModelLik/sum(Results$ModelLik)        #compute Akaike weights

    ##check if some models are redundant
    if(length(unique(Results$AICc)) != length(cand.set)) warning("\nCheck model structure carefully as some models may be redundant\n") 

    ##add LL
    Results$LL <- unlist(lapply(X = cand.set, FUN = function(i) logLik(i)[1]))      
    
       
    #rename correctly to AIC
    if(second.ord == FALSE) {
      colnames(Results)[1:6] <- c("Modnames", "K", "AIC", "Delta_AIC", "ModelLik", "AICWt")
    }

    if(sort)  {
      Results <- Results[rev(order(Results[, 6])),] 	  #if sort=TRUE, models are ranked based on Akaike weights
      Results$Cum.Wt <- cumsum(Results[, 6])                        #display cumulative sum of Akaike weights
    } else {Results$Cum.Wt <- NULL}
   
    
    class(Results) <- c("aictab", "data.frame")
    return(Results)
  }



##multinom
aictab.AICmultinom.nnet <-
  function(cand.set, modnames = NULL, second.ord = TRUE, nobs = NULL, sort = TRUE, c.hat = 1, ...){  #specify whether table should be sorted or not by delta AICc

    ##check if named list if modnames are not supplied
    if(is.null(modnames)) {
      if(is.null(names(cand.set))) {
        modnames <- paste("Mod", 1:length(cand.set), sep = "")
        warning("\nModel names have been supplied automatically in the table\n")
      } else {
        modnames <- names(cand.set)
      }
    }
    
    ##add check to see whether response variable is the same for all models
    check.resp <- lapply(X = cand.set, FUN = function(b) formula(b)[2])
    if(length(unique(check.resp)) > 1) stop("\nYou must use the same response variable for all models\n")
    
    Results <- data.frame(Modnames = modnames)                    #assign model names to first column
    Results$K <- unlist(lapply(X = cand.set, FUN = AICc, return.K = TRUE, 
                               second.ord = second.ord, nobs = nobs, c.hat = c.hat))     #extract number of parameters
    Results$AICc <- unlist(lapply(X = cand.set, FUN = AICc, return.K = FALSE, 
                                  second.ord = second.ord, nobs = nobs, c.hat = c.hat))  #extract AICc                                      #
    Results$Delta_AICc <- Results$AICc - min(Results$AICc)            #compute delta AICc
    Results$ModelLik <- exp(-0.5*Results$Delta_AICc)                #compute model likelihood required to compute Akaike weights
    Results$AICcWt <- Results$ModelLik/sum(Results$ModelLik)        #compute Akaike weights

    ##check if some models are redundant
    if(length(unique(Results$AICc)) != length(cand.set)) warning("\nCheck model structure carefully as some models may be redundant\n")
    
    ##check if LL computed
    if(second.ord==TRUE && c.hat == 1) {
      Results$LL <- unlist(lapply(X = cand.set, FUN = function(i) logLik(i)[1]))
    }
    
    
    ##rename correctly to QAICc and add column for c-hat
    if(second.ord == TRUE && c.hat > 1) {
      colnames(Results) <- c("Modnames", "K", "QAICc", "Delta_QAICc", "ModelLik", "QAICcWt")
      LL <- unlist(lapply(X = cand.set, FUN = function(i) logLik(i)[1]))
      Results$Quasi.LL <- LL/c.hat
      Results$c_hat <- c.hat
    }      
    

    ##rename correctly to AIC
    if(second.ord == FALSE && c.hat == 1) {
      colnames(Results) <- c("Modnames", "K", "AIC", "Delta_AIC", "ModelLik", "AICWt")
      Results$LL <- unlist(lapply(X = cand.set, FUN = function(i) logLik(i)[1]))      
    }  

    ##rename correctly to QAIC and add column for c-hat
    if(second.ord == FALSE && c.hat > 1) {
      colnames(Results) <- c("Modnames", "K", "QAIC", "Delta_QAIC", "ModelLik", "QAICWt")
      LL <- unlist(lapply(X = cand.set, FUN = function(i) logLik(i)[1]))
      Results$Quasi.LL <- LL/c.hat
      Results$c_hat <- c.hat
    }      


    if(sort)  {
      Results <- Results[rev(order(Results[, 6])),] 	  #if sort=TRUE, models are ranked based on Akaike weights
      Results$Cum.Wt <- cumsum(Results[, 6])                        #display cumulative sum of Akaike weights
    } else {Results$Cum.Wt <- NULL}

    
    class(Results) <- c("aictab", "data.frame")
    return(Results)
  }



##lme
aictab.AIClme <-
  function(cand.set, modnames = NULL, second.ord = TRUE, nobs = NULL, sort = TRUE, ...){  #specify whether table should be sorted or not by delta AICc

    ##check if named list if modnames are not supplied
    if(is.null(modnames)) {
      if(is.null(names(cand.set))) {
        modnames <- paste("Mod", 1:length(cand.set), sep = "")
        warning("\nModel names have been supplied automatically in the table\n")
      } else {
        modnames <- names(cand.set)
      }
    }
    
    ##add check to see whether response variable is the same for all models
    check.resp <- lapply(X = cand.set, FUN = function(b) formula(b)[2])
    if(length(unique(check.resp)) > 1) stop("\nYou must use the same response variable for all models\n")
       
    Results <- NULL
    ##check if models were fit with same method (REML or ML)
    check_ML <- unlist(lapply(cand.set, FUN = function(i) i$method))
    
    if (any(check_ML != "ML")) {
      warning("\nModel selection for fixed effects is only appropriate with method=ML:", "\n",
                    "REML (default) should only be used to select random effects for a constant set of fixed effects\n")
    }
    
    check.method <- unique(check_ML)
    
    if(identical(check.method, c("ML", "REML"))) {
      stop("\nYou should not have models fit with REML and ML in the same candidate model set\n")
    }
    
    Results <- data.frame(Modnames = modnames)                    #assign model names to first column
    Results$K <- unlist(lapply(cand.set, AICc, return.K = TRUE,
                               second.ord = second.ord, nobs = nobs))     #extract number of parameters
    Results$AICc <- unlist(lapply(cand.set, AICc, return.K = FALSE,
                                  second.ord = second.ord, nobs = nobs))  #extract AICc                                      #
    Results$Delta_AICc <- Results$AICc - min(Results$AICc)            #compute delta AICc
    Results$ModelLik <- exp(-0.5*Results$Delta_AICc)                #compute model likelihood required to compute Akaike weights
    Results$AICcWt <- Results$ModelLik/sum(Results$ModelLik)        #compute Akaike weights

    ##check if some models are redundant
    if(length(unique(Results$AICc)) != length(cand.set)) warning("\nCheck model structure carefully as some models may be redundant\n")
    
    #check if ML or REML used and add column accordingly
    if(identical(check.method, "ML")) {
      Results$LL <- unlist(lapply(X= cand.set, FUN = function(i) logLik(i)[1]))      
    }

    if(identical(check.method, "REML")) {
      Results$Res.LL <- unlist(lapply(X = cand.set, FUN = function(i) logLik(i)[1]))      
    }

    
    #rename correctly to AIC
    if(second.ord == FALSE) {
      colnames(Results)[1:6] <- c("Modnames", "K", "AIC", "Delta_AIC", "ModelLik", "AICWt")
    }

    if(sort)  {
      Results <- Results[rev(order(Results[, 6])),] 	  #if sort=TRUE, models are ranked based on Akaike weights
      Results$Cum.Wt <- cumsum(Results[, 6])                        #display cumulative sum of Akaike weights
    } else {Results$Cum.Wt <- NULL}
   
    
    class(Results) <- c("aictab", "data.frame")
    return(Results)
  }



##nls
aictab.AICnls <-
  function(cand.set, modnames = NULL, second.ord = TRUE, nobs = NULL, sort = TRUE, ...){  #specify whether table should be sorted or not by delta AICc

    ##check if named list if modnames are not supplied
    if(is.null(modnames)) {
      if(is.null(names(cand.set))) {
        modnames <- paste("Mod", 1:length(cand.set), sep = "")
        warning("\nModel names have been supplied automatically in the table\n")
      } else {
        modnames <- names(cand.set)
      }
    }
    
    ##add check to see whether response variable is the same for all models
    check.resp <- lapply(X = cand.set, FUN = function(b) formula(b)[2])
    if(length(unique(check.resp)) > 1) stop("\nYou must use the same response variable for all models\n")
       
    Results <- data.frame(Modnames=modnames)                    #assign model names to first column
    Results$K <- unlist(lapply(cand.set, AICc, return.K = TRUE, second.ord = second.ord, nobs = nobs))     #extract number of parameters
    Results$AICc <- unlist(lapply(cand.set, AICc, return.K = FALSE, second.ord = second.ord, nobs = nobs))  #extract AICc                                      #
    Results$Delta_AICc <- Results$AICc - min(Results$AICc)            #compute delta AICc
    Results$ModelLik <- exp(-0.5*Results$Delta_AICc)                #compute model likelihood required to compute Akaike weights
    Results$AICcWt <- Results$ModelLik/sum(Results$ModelLik)        #compute Akaike weights
    Results$LL <- unlist(lapply(X = cand.set, FUN = function(i) logLik(i)[1]))
    
    ##check if some models are redundant
    if(length(unique(Results$AICc)) != length(cand.set)) warning("\nCheck model structure carefully as some models may be redundant\n")
    
    ##rename correctly to AIC
    if(second.ord == FALSE) {
      colnames(Results)[1:6] <- c("Modnames", "K", "AIC", "Delta_AIC", "ModelLik", "AICWt")
    }

    if(sort)  {
      Results <- Results[rev(order(Results[, 6])),] 	  #if sort=TRUE, models are ranked based on Akaike weights
      Results$Cum.Wt <- cumsum(Results[, 6])                        #display cumulative sum of Akaike weights
    } else {Results$Cum.Wt <- NULL}
   
    
    class(Results) <- c("aictab", "data.frame")
    return(Results)
  }



##polr
aictab.AICpolr <-
  function(cand.set, modnames = NULL, second.ord = TRUE, nobs = NULL, sort = TRUE, ...){  #specify whether table should be sorted or not by delta AICc

    ##check if named list if modnames are not supplied
    if(is.null(modnames)) {
      if(is.null(names(cand.set))) {
        modnames <- paste("Mod", 1:length(cand.set), sep = "")
        warning("\nModel names have been supplied automatically in the table\n")
      } else {
        modnames <- names(cand.set)
      }
    }
    
    ##add check to see whether response variable is the same for all models
    check.resp <- lapply(X = cand.set, FUN = function(b) formula(b)[2])
    if(length(unique(check.resp)) > 1) stop("\nYou must use the same response variable for all models\n")
    
    Results <- NULL
    Results <- data.frame(Modnames=modnames)                    #assign model names to first column
    Results$K <- unlist(lapply(X = cand.set, FUN = AICc, return.K = TRUE, 
                               second.ord = second.ord, nobs = nobs))     #extract number of parameters
    Results$AICc <- unlist(lapply(X = cand.set, FUN = AICc, return.K = FALSE, 
                                  second.ord = second.ord, nobs = nobs))  #extract AICc                                      #
    Results$Delta_AICc <- Results$AICc - min(Results$AICc)            #compute delta AICc
    Results$ModelLik <- exp(-0.5*Results$Delta_AICc)                #compute model likelihood required to compute Akaike weights
    Results$AICcWt <- Results$ModelLik/sum(Results$ModelLik)        #compute Akaike weights
    Results$LL <- unlist(lapply(X = cand.set, FUN = function(i) logLik(i)[1]))      

    ##check if some models are redundant
    if(length(unique(Results$AICc)) != length(cand.set)) warning("\nCheck model structure carefully as some models may be redundant\n")
    
    #rename correctly to AIC
    if(second.ord==FALSE) {
      colnames(Results)[1:6]<-c("Modnames", "K", "AIC", "Delta_AIC", "ModelLik",
                           "AICWt")
    }


    if(sort)  {
      Results <- Results[rev(order(Results[, 6])),] 	  #if sort=TRUE, models are ranked based on Akaike weights
      Results$Cum.Wt <- cumsum(Results[, 6])                        #display cumulative sum of Akaike weights
    } else {Results$Cum.Wt <- NULL}


    
    class(Results) <- c("aictab", "data.frame")
    return(Results)
  }



##rlm
aictab.AICrlm.lm <-
  function(cand.set, modnames = NULL, second.ord = TRUE, nobs = NULL, sort = TRUE, ...){  #specify whether table should be sorted or not by delta AICc

    ##check if named list if modnames are not supplied
    if(is.null(modnames)) {
      if(is.null(names(cand.set))) {
        modnames <- paste("Mod", 1:length(cand.set), sep = "")
        warning("\nModel names have been supplied automatically in the table\n")
      } else {
        modnames <- names(cand.set)
      }
    }
    
    ##add check to see whether response variable is the same for all models
    check.resp <- lapply(X = cand.set, FUN = function(b) formula(b)[2])
    if(length(unique(check.resp)) > 1) stop("\nYou must use the same response variable for all models\n")
       
    Results <- data.frame(Modnames = modnames)                    #assign model names to first column
    Results$K <- unlist(lapply(cand.set, AICc, return.K = TRUE,
                               second.ord = second.ord, nobs = nobs))     #extract number of parameters
    Results$AICc <- unlist(lapply(cand.set, AICc, return.K = FALSE,
                                  second.ord = second.ord, nobs = nobs))  #extract AICc                                      
    Results$Delta_AICc <- Results$AICc - min(Results$AICc)            #compute delta AICc
    Results$ModelLik <- exp(-0.5*Results$Delta_AICc)                #compute model likelihood required to compute Akaike weights
    Results$AICcWt <- Results$ModelLik/sum(Results$ModelLik)        #compute Akaike weights

   ##check if some models are redundant
    if(length(unique(Results$AICc)) != length(cand.set)) warning("\nCheck model structure carefully as some models may be redundant\n")
    
    Results$LL <- unlist(lapply(X= cand.set, FUN = function(i) logLik(i)[1]))      

    ##rename correctly to AIC
    if(second.ord == FALSE) {
      colnames(Results)[1:6] <- c("Modnames", "K", "AIC", "Delta_AIC", "ModelLik", "AICWt")
    }

    if(sort)  {
      Results <- Results[rev(order(Results[, 6])),] 	  #if sort=TRUE, models are ranked based on Akaike weights
      Results$Cum.Wt <- cumsum(Results[, 6])                        #display cumulative sum of Akaike weights
    } else {Results$Cum.Wt <- NULL}
   
    
    class(Results) <- c("aictab", "data.frame")
    return(Results)
  }



##survreg
aictab.AICsurvreg <-
  function(cand.set, modnames = NULL, second.ord = TRUE, nobs = NULL, sort = TRUE, ...){

    ##check if named list if modnames are not supplied
    if(is.null(modnames)) {
      if(is.null(names(cand.set))) {
        modnames <- paste("Mod", 1:length(cand.set), sep = "")
        warning("\nModel names have been supplied automatically in the table\n")
      } else {
        modnames <- names(cand.set)
      }
    }
      

    ##add check to see whether response variable is the same for all models
    check.resp <- lapply(X = cand.set, FUN = function(b) formula(b)[2])
    if(length(unique(check.resp)) > 1) stop("\nYou must use the same response variable for all models\n")

    Results <- NULL
    Results <- data.frame(Modnames = modnames)                    #assign model names to first column
    Results$K <- unlist(lapply(X = cand.set, FUN = AICc, return.K = TRUE,
                               second.ord = second.ord, nobs = nobs))     #extract number of parameters
    Results$AICc <- unlist(lapply(X = cand.set, FUN = AICc, return.K = FALSE, 
                                  second.ord = second.ord, nobs = nobs))  #extract AICc                                      #
    Results$Delta_AICc <- Results$AICc - min(Results$AICc)            #compute delta AICc
    Results$ModelLik <- exp(-0.5*Results$Delta_AICc)                #compute model likelihood required to compute Akaike weights
    Results$AICcWt <- Results$ModelLik/sum(Results$ModelLik)        #compute Akaike weights

    ##check if some models are redundant
    if(length(unique(Results$AICc)) != length(cand.set)) warning("\nCheck model structure carefully as some models may be redundant\n")
         
    ##extract LL
    Results$LL <- unlist(lapply(X = cand.set, FUN = function(i) logLik(i)[1]))
        
    ##rename correctly to AIC
    if(second.ord == FALSE) {
      colnames(Results)[1:6] <- c("Modnames", "K", "AIC", "Delta_AIC", "ModelLik", "AICWt")
    }  

      
    if(sort)  {
      Results <- Results[rev(order(Results[, 6])),] 	  #if sort=TRUE, models are ranked based on Akaike weights
      Results$Cum.Wt <- cumsum(Results[, 6])                        #display cumulative sum of Akaike weights
    } else {Results$Cum.Wt <- NULL}

    
    class(Results) <- c("aictab", "data.frame")
    return(Results)
  }



##occu
aictab.AICunmarkedFitOccu <- function(cand.set, modnames = NULL, second.ord = TRUE, nobs = NULL, sort = TRUE, c.hat = 1, ...){  #specify whether table should be sorted or not by delta AICc

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
  Results$K <- unlist(lapply(cand.set, AICc, return.K = TRUE,
                             second.ord = second.ord, nobs = nobs, c.hat = c.hat))     #extract number of parameters
  Results$AICc <- unlist(lapply(cand.set, AICc, return.K = FALSE,
                                second.ord = second.ord, nobs = nobs, c.hat = c.hat))  #extract AICc                                      #
  Results$Delta_AICc <- Results$AICc-min(Results$AICc)            #compute delta AICc
  Results$ModelLik <- exp(-0.5*Results$Delta_AICc)                #compute model likelihood required to compute Akaike weights
  Results$AICcWt <- Results$ModelLik/sum(Results$ModelLik)        #compute Akaike weights

  ##check if some models are redundant
  if(length(unique(Results$AICc)) != length(cand.set)) warning("\nCheck model structure carefully as some models may be redundant\n")


  ##check if AICc and c.hat = 1
  if(second.ord == TRUE && c.hat == 1) {
    Results$LL <- unlist(lapply(X = cand.set, FUN = function(i) extractLL(i)))      
  }
  
  ##rename correctly to QAICc and add column for c-hat
  if(second.ord == TRUE && c.hat > 1) {
    colnames(Results) <- c("Modnames", "K", "QAICc", "Delta_QAICc", "ModelLik", "QAICcWt")
    LL <- unlist(lapply(X = cand.set, FUN = function(i) extractLL(i))) 
    Results$Quasi.LL <- LL/c.hat
    Results$c_hat <- c.hat
  }      

  ##rename correctly to AIC
  if(second.ord == FALSE && c.hat==1) {
    colnames(Results)<-c("Modnames", "K", "AIC", "Delta_AIC", "ModelLik", "AICWt")
    Results$LL <- unlist(lapply(X = cand.set, FUN = function(i) extractLL(i))) 
  }  
  
  ##rename correctly to QAIC and add column for c-hat
  if(second.ord == FALSE && c.hat > 1) {
    colnames(Results)<-c("Modnames", "K", "QAIC", "Delta_QAIC", "ModelLik", "QAICWt")
    LL <- unlist(lapply(X = cand.set, FUN = function(i) extractLL(i))) 
    Results$Quasi.LL <- LL/c.hat
    Results$c_hat<-c.hat
  }      


  if(sort)  {
    Results <- Results[rev(order(Results[, 6])),] 	  #if sort=TRUE, models are ranked based on Akaike weights
    Results$Cum.Wt <- cumsum(Results[, 6])                        #display cumulative sum of Akaike weights
  } else {Results$Cum.Wt <- NULL}
  
  class(Results) <- c("aictab", "data.frame")
  return(Results)
}


##colext
aictab.AICunmarkedFitColExt <- function(cand.set, modnames = NULL, second.ord = TRUE, nobs = NULL, sort = TRUE, c.hat = 1, ...){  #specify whether table should be sorted or not by delta AICc

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
  Results$K <- unlist(lapply(cand.set, AICc, return.K = TRUE,
                             second.ord = second.ord, nobs = nobs, c.hat = c.hat))     #extract number of parameters
  Results$AICc <- unlist(lapply(cand.set, AICc, return.K = FALSE, 
                                second.ord = second.ord, nobs = nobs, c.hat = c.hat))  #extract AICc                                      #
  Results$Delta_AICc <- Results$AICc-min(Results$AICc)            #compute delta AICc
  Results$ModelLik <- exp(-0.5*Results$Delta_AICc)                #compute model likelihood required to compute Akaike weights
  Results$AICcWt <- Results$ModelLik/sum(Results$ModelLik)        #compute Akaike weights

  ##check if some models are redundant
  if(length(unique(Results$AICc)) != length(cand.set)) warning("\nCheck model structure carefully as some models may be redundant\n")


  ##check if AICc and c.hat = 1
  if(second.ord == TRUE && c.hat == 1) {
    Results$LL <- unlist(lapply(X = cand.set, FUN = function(i) extractLL(i)))      
  }
  
  ##rename correctly to QAICc and add column for c-hat
  if(second.ord == TRUE && c.hat > 1) {
    colnames(Results) <- c("Modnames", "K", "QAICc", "Delta_QAICc", "ModelLik", "QAICcWt")
    LL <- unlist(lapply(X = cand.set, FUN = function(i) extractLL(i))) 
    Results$Quasi.LL <- LL/c.hat
    Results$c_hat <- c.hat
  }      

  ##rename correctly to AIC
  if(second.ord == FALSE && c.hat==1) {
    colnames(Results)<-c("Modnames", "K", "AIC", "Delta_AIC", "ModelLik", "AICWt")
    Results$LL <- unlist(lapply(X = cand.set, FUN = function(i) extractLL(i))) 
  }  
  
  ##rename correctly to QAIC and add column for c-hat
  if(second.ord == FALSE && c.hat > 1) {
    colnames(Results)<-c("Modnames", "K", "QAIC", "Delta_QAIC", "ModelLik", "QAICWt")
    LL <- unlist(lapply(X = cand.set, FUN = function(i) extractLL(i))) 
    Results$Quasi.LL <- LL/c.hat
    Results$c_hat<-c.hat
  }      


  if(sort)  {
    Results <- Results[rev(order(Results[, 6])),] 	  #if sort=TRUE, models are ranked based on Akaike weights
    Results$Cum.Wt <- cumsum(Results[, 6])                        #display cumulative sum of Akaike weights
  } else {Results$Cum.Wt <- NULL}
  
  class(Results) <- c("aictab", "data.frame")
  return(Results)
}


##occuRN
aictab.AICunmarkedFitOccuRN <- function(cand.set, modnames = NULL, second.ord = TRUE, nobs = NULL, sort = TRUE, c.hat = 1, ...){  #specify whether table should be sorted or not by delta AICc

  ##check if named list if modnames are not supplied
  if(is.null(modnames)) {
    if(is.null(names(cand.set))) {
      modnames <- paste("Mod", 1:length(cand.set), sep = "")
      warning("\nModel names have been supplied automatically in the table\n")
    } else {
      modnames <- names(cand.set)
    }
  }
    

  ##add check for use of c-hat
  #if(c.hat > 1) stop("\nThe correction for overdispersion is not appropriate with Royle-Nichols heterogeneity models\n")

  Results <- data.frame(Modnames = modnames)                    #assign model names to first column
  Results$K <- unlist(lapply(cand.set, AICc, return.K = TRUE,
                             second.ord = second.ord, nobs = nobs, c.hat = c.hat))     #extract number of parameters
  Results$AICc <- unlist(lapply(cand.set, AICc, return.K = FALSE,
                                second.ord = second.ord, nobs = nobs, c.hat = c.hat))  #extract AICc                                      
  Results$Delta_AICc <- Results$AICc-min(Results$AICc)            #compute delta AICc
  Results$ModelLik <- exp(-0.5*Results$Delta_AICc)                #compute model likelihood required to compute Akaike weights
  Results$AICcWt <- Results$ModelLik/sum(Results$ModelLik)        #compute Akaike weights

  ##check if some models are redundant
  if(length(unique(Results$AICc)) != length(cand.set)) warning("\nCheck model structure carefully as some models may be redundant\n")

  ##check if AICc and c.hat = 1
  if(second.ord == TRUE && c.hat == 1) {
    Results$LL <- unlist(lapply(X = cand.set, FUN = function(i) extractLL(i)))      
  }
  
  ##rename correctly to QAICc and add column for c-hat
  if(second.ord == TRUE && c.hat > 1) {
    colnames(Results) <- c("Modnames", "K", "QAICc", "Delta_QAICc", "ModelLik", "QAICcWt")
    LL <- unlist(lapply(X = cand.set, FUN = function(i) extractLL(i))) 
    Results$Quasi.LL <- LL/c.hat
    Results$c_hat <- c.hat
  }      

  ##rename correctly to AIC
  if(second.ord == FALSE && c.hat==1) {
    colnames(Results)<-c("Modnames", "K", "AIC", "Delta_AIC", "ModelLik", "AICWt")
    Results$LL <- unlist(lapply(X = cand.set, FUN = function(i) extractLL(i))) 
  }  
  
  ##rename correctly to QAIC and add column for c-hat
  if(second.ord == FALSE && c.hat > 1) {
    colnames(Results)<-c("Modnames", "K", "QAIC", "Delta_QAIC", "ModelLik", "QAICWt")
    LL <- unlist(lapply(X = cand.set, FUN = function(i) extractLL(i))) 
    Results$Quasi.LL <- LL/c.hat
    Results$c_hat<-c.hat
  }      

  
  if(sort)  {
    Results <- Results[rev(order(Results[, 6])),] 	  #if sort=TRUE, models are ranked based on Akaike weights
    Results$Cum.Wt <- cumsum(Results[, 6])                        #display cumulative sum of Akaike weights
  } else {Results$Cum.Wt <- NULL}
  
  class(Results) <- c("aictab", "data.frame")
  return(Results)
}



##pcount
aictab.AICunmarkedFitPCount <- function(cand.set, modnames = NULL, second.ord = TRUE, nobs = NULL, sort = TRUE, c.hat = 1, ...){  #specify whether table should be sorted or not by delta AICc

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
  Results$K <- unlist(lapply(cand.set, AICc, return.K = TRUE,
                             second.ord = second.ord, nobs = nobs, c.hat = c.hat))     #extract number of parameters
  Results$AICc <- unlist(lapply(cand.set, AICc, return.K = FALSE,
                                second.ord = second.ord, nobs = nobs, c.hat = c.hat))  #extract AICc                                      #
  Results$Delta_AICc <- Results$AICc - min(Results$AICc)            #compute delta AICc
  Results$ModelLik <- exp(-0.5*Results$Delta_AICc)                #compute model likelihood required to compute Akaike weights
  Results$AICcWt <- Results$ModelLik/sum(Results$ModelLik)        #compute Akaike weights

  ##check if some models are redundant
  if(length(unique(Results$AICc)) != length(cand.set)) warning("\nCheck model structure carefully as some models may be redundant\n")


  ##check if AICc and c.hat = 1
  if(second.ord == TRUE && c.hat == 1) {
    Results$LL <- unlist(lapply(X = cand.set, FUN = function(i) extractLL(i)))      
  }
  
  ##rename correctly to QAICc and add column for c-hat
  if(second.ord == TRUE && c.hat > 1) {
    colnames(Results) <- c("Modnames", "K", "QAICc", "Delta_QAICc", "ModelLik", "QAICcWt")
    LL <- unlist(lapply(X = cand.set, FUN = function(i) extractLL(i))) 
    Results$Quasi.LL <- LL/c.hat
    Results$c_hat <- c.hat
  }      

  ##rename correctly to AIC
  if(second.ord == FALSE && c.hat == 1) {
    colnames(Results)<-c("Modnames", "K", "AIC", "Delta_AIC", "ModelLik", "AICWt")
    Results$LL <- unlist(lapply(X = cand.set, FUN = function(i) extractLL(i))) 
  }  
  
  ##rename correctly to QAIC and add column for c-hat
  if(second.ord == FALSE && c.hat > 1) {
    colnames(Results)<-c("Modnames", "K", "QAIC", "Delta_QAIC", "ModelLik", "QAICWt")
    LL <- unlist(lapply(X = cand.set, FUN = function(i) extractLL(i))) 
    Results$Quasi.LL <- LL/c.hat
    Results$c_hat<-c.hat
  }      


  if(sort)  {
    Results <- Results[rev(order(Results[, 6])),] 	  #if sort=TRUE, models are ranked based on Akaike weights
    Results$Cum.Wt <- cumsum(Results[, 6])                        #display cumulative sum of Akaike weights
  } else {Results$Cum.Wt <- NULL}
  
  class(Results) <- c("aictab", "data.frame")
  return(Results)
}



##same function as that for objects created by pcount( )
aictab.AICunmarkedFitPCO <- function(cand.set, modnames = NULL, second.ord = TRUE, nobs = NULL, sort = TRUE, c.hat = 1, ...){  #specify whether table should be sorted or not by delta AICc

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
  Results$K <- unlist(lapply(cand.set, AICc, return.K = TRUE,
                             second.ord = second.ord, nobs = nobs, c.hat = c.hat))     #extract number of parameters
  Results$AICc <- unlist(lapply(cand.set, AICc, return.K = FALSE,
                                second.ord = second.ord, nobs = nobs, c.hat = c.hat))  #extract AICc                                      
  Results$Delta_AICc <- Results$AICc - min(Results$AICc)            #compute delta AICc
  Results$ModelLik <- exp(-0.5*Results$Delta_AICc)                #compute model likelihood required to compute Akaike weights
  Results$AICcWt <- Results$ModelLik/sum(Results$ModelLik)        #compute Akaike weights

  ##check if some models are redundant
  if(length(unique(Results$AICc)) != length(cand.set)) warning("\nCheck model structure carefully as some models may be redundant\n")


  ##check if AICc and c.hat = 1
  if(second.ord == TRUE && c.hat == 1) {
    Results$LL <- unlist(lapply(X = cand.set, FUN = function(i) extractLL(i)))      
  }
  
  ##rename correctly to QAICc and add column for c-hat
  if(second.ord == TRUE && c.hat > 1) {
    colnames(Results) <- c("Modnames", "K", "QAICc", "Delta_QAICc", "ModelLik", "QAICcWt")
    LL <- unlist(lapply(X = cand.set, FUN = function(i) extractLL(i))) 
    Results$Quasi.LL <- LL/c.hat
    Results$c_hat <- c.hat
  }      

  ##rename correctly to AIC
  if(second.ord == FALSE && c.hat == 1) {
    colnames(Results)<-c("Modnames", "K", "AIC", "Delta_AIC", "ModelLik", "AICWt")
    Results$LL <- unlist(lapply(X = cand.set, FUN = function(i) extractLL(i))) 
  }  
  
  ##rename correctly to QAIC and add column for c-hat
  if(second.ord == FALSE && c.hat > 1) {
    colnames(Results)<-c("Modnames", "K", "QAIC", "Delta_QAIC", "ModelLik", "QAICWt")
    LL <- unlist(lapply(X = cand.set, FUN = function(i) extractLL(i))) 
    Results$Quasi.LL <- LL/c.hat
    Results$c_hat<-c.hat
  }      


  if(sort)  {
    Results <- Results[rev(order(Results[, 6])),] 	  #if sort=TRUE, models are ranked based on Akaike weights
    Results$Cum.Wt <- cumsum(Results[, 6])                        #display cumulative sum of Akaike weights
  } else {Results$Cum.Wt <- NULL}
  
  class(Results) <- c("aictab", "data.frame")
  return(Results)
}


##distsamp
aictab.AICunmarkedFitDS <- function(cand.set, modnames = NULL, second.ord = TRUE, nobs = NULL, sort = TRUE, c.hat = 1, ...){  #specify whether table should be sorted or not by delta AICc

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
  Results$K <- unlist(lapply(cand.set, AICc, return.K = TRUE,
                             second.ord = second.ord, nobs = nobs, c.hat = c.hat))     #extract number of parameters
  Results$AICc <- unlist(lapply(cand.set, AICc, return.K = FALSE,
                                second.ord = second.ord, nobs = nobs, c.hat = c.hat))  #extract AICc                                      #
  Results$Delta_AICc <- Results$AICc - min(Results$AICc)            #compute delta AICc
  Results$ModelLik <- exp(-0.5*Results$Delta_AICc)                #compute model likelihood required to compute Akaike weights
  Results$AICcWt <- Results$ModelLik/sum(Results$ModelLik)        #compute Akaike weights

  ##check if some models are redundant
  if(length(unique(Results$AICc)) != length(cand.set)) warning("\nCheck model structure carefully as some models may be redundant\n")

  ##check if AICc and c.hat = 1
  if(second.ord == TRUE && c.hat == 1) {
    Results$LL <- unlist(lapply(X = cand.set, FUN = function(i) extractLL(i)))      
  }
  
  ##rename correctly to QAICc and add column for c-hat
  if(second.ord == TRUE && c.hat > 1) {
    colnames(Results) <- c("Modnames", "K", "QAICc", "Delta_QAICc", "ModelLik", "QAICcWt")
    LL <- unlist(lapply(X = cand.set, FUN = function(i) extractLL(i))) 
    Results$Quasi.LL <- LL/c.hat
    Results$c_hat <- c.hat
  }      

  ##rename correctly to AIC
  if(second.ord == FALSE && c.hat == 1) {
    colnames(Results)<-c("Modnames", "K", "AIC", "Delta_AIC", "ModelLik", "AICWt")
    Results$LL <- unlist(lapply(X = cand.set, FUN = function(i) extractLL(i))) 
  }  
  
  ##rename correctly to QAIC and add column for c-hat
  if(second.ord == FALSE && c.hat > 1) {
    colnames(Results)<-c("Modnames", "K", "QAIC", "Delta_QAIC", "ModelLik", "QAICWt")
    LL <- unlist(lapply(X = cand.set, FUN = function(i) extractLL(i))) 
    Results$Quasi.LL <- LL/c.hat
    Results$c_hat<-c.hat
  }
  
  if(sort)  {
    Results <- Results[rev(order(Results[, 6])),] 	  #if sort=TRUE, models are ranked based on Akaike weights
    Results$Cum.Wt <- cumsum(Results[, 6])                        #display cumulative sum of Akaike weights
  } else {Results$Cum.Wt <- NULL}
  
  class(Results) <- c("aictab", "data.frame")
  return(Results)
}



##gdistsamp
aictab.AICunmarkedFitGDS <- function(cand.set, modnames = NULL, second.ord = TRUE, nobs = NULL, sort = TRUE, c.hat = 1, ...){  #specify whether table should be sorted or not by delta AICc

  ##check if named list if modnames are not supplied
  if(is.null(modnames)) {
    if(is.null(names(cand.set))) {
      modnames <- paste("Mod", 1:length(cand.set), sep = "")
      warning("\nModel names have been supplied automatically in the table\n")
    } else {
      modnames <- names(cand.set)
    }
  }
    

  ##add check for use of c-hat
  #if(c.hat > 1) stop("\nThe correction for overdispersion is not appropriate for distance sampling models\n")

  Results <- data.frame(Modnames = modnames)                    #assign model names to first column
  Results$K <- unlist(lapply(cand.set, AICc, return.K = TRUE,
                             second.ord = second.ord, nobs = nobs, c.hat = c.hat))     #extract number of parameters
  Results$AICc <- unlist(lapply(cand.set, AICc, return.K = FALSE,
                                second.ord = second.ord, nobs = nobs, c.hat = c.hat))  #extract AICc
  Results$Delta_AICc <- Results$AICc - min(Results$AICc)            #compute delta AICc
  Results$ModelLik <- exp(-0.5*Results$Delta_AICc)                #compute model likelihood required to compute Akaike weights
  Results$AICcWt <- Results$ModelLik/sum(Results$ModelLik)        #compute Akaike weights

  ##check if some models are redundant
  if(length(unique(Results$AICc)) != length(cand.set)) warning("\nCheck model structure carefully as some models may be redundant\n")

  ##check if AICc and c.hat = 1
  if(second.ord == TRUE && c.hat == 1) {
    Results$LL <- unlist(lapply(X = cand.set, FUN = function(i) extractLL(i)))      
  }
  
  ##rename correctly to QAICc and add column for c-hat
  if(second.ord == TRUE && c.hat > 1) {
    colnames(Results) <- c("Modnames", "K", "QAICc", "Delta_QAICc", "ModelLik", "QAICcWt")
    LL <- unlist(lapply(X = cand.set, FUN = function(i) extractLL(i))) 
    Results$Quasi.LL <- LL/c.hat
    Results$c_hat <- c.hat
  }      

  ##rename correctly to AIC
  if(second.ord == FALSE && c.hat == 1) {
    colnames(Results)<-c("Modnames", "K", "AIC", "Delta_AIC", "ModelLik", "AICWt")
    Results$LL <- unlist(lapply(X = cand.set, FUN = function(i) extractLL(i))) 
  }  
  
  ##rename correctly to QAIC and add column for c-hat
  if(second.ord == FALSE && c.hat > 1) {
    colnames(Results)<-c("Modnames", "K", "QAIC", "Delta_QAIC", "ModelLik", "QAICWt")
    LL <- unlist(lapply(X = cand.set, FUN = function(i) extractLL(i))) 
    Results$Quasi.LL <- LL/c.hat
    Results$c_hat<-c.hat
  }      

  
  if(sort)  {
    Results <- Results[rev(order(Results[, 6])),] 	  #if sort=TRUE, models are ranked based on Akaike weights
    Results$Cum.Wt <- cumsum(Results[, 6])                        #display cumulative sum of Akaike weights
  } else {Results$Cum.Wt <- NULL}
  
  class(Results) <- c("aictab", "data.frame")
  return(Results)
}



##occuFP
aictab.AICunmarkedFitOccuFP <- function(cand.set, modnames = NULL, second.ord = TRUE, nobs = NULL, sort = TRUE, c.hat = 1, ...){  #specify whether table should be sorted or not by delta AICc

  ##check if named list if modnames are not supplied
  if(is.null(modnames)) {
    if(is.null(names(cand.set))) {
      modnames <- paste("Mod", 1:length(cand.set), sep = "")
      warning("\nModel names have been supplied automatically in the table\n")
    } else {
      modnames <- names(cand.set)
    }
  }
    
  ##add check for use of c-hat
  if(c.hat > 1) stop("\nThe correction for overdispersion is not yet implemented for false-positive occupancy models\n")

  Results <- data.frame(Modnames = modnames)                    #assign model names to first column
  Results$K <- unlist(lapply(cand.set, AICc, return.K = TRUE,
                             second.ord = second.ord, nobs = nobs, c.hat = c.hat))     #extract number of parameters
  Results$AICc <- unlist(lapply(cand.set, AICc, return.K = FALSE,
                                second.ord = second.ord, nobs = nobs, c.hat = c.hat))  #extract AICc
  Results$Delta_AICc <- Results$AICc - min(Results$AICc)            #compute delta AICc
  Results$ModelLik <- exp(-0.5*Results$Delta_AICc)                #compute model likelihood required to compute Akaike weights
  Results$AICcWt <- Results$ModelLik/sum(Results$ModelLik)        #compute Akaike weights

  ##check if some models are redundant
  if(length(unique(Results$AICc)) != length(cand.set)) warning("\nCheck model structure carefully as some models may be redundant\n")

  Results$LL <- unlist(lapply(X = cand.set, FUN = function(i) extractLL(i)))      
    
  ##rename correctly to AIC
  if(second.ord == FALSE) {
    colnames(Results)[1:6] <- c("Modnames", "K", "AIC", "Delta_AIC", "ModelLik", "AICWt")
  }  
  
  if(sort)  {
    Results <- Results[rev(order(Results[, 6])),] 	  #if sort=TRUE, models are ranked based on Akaike weights
    Results$Cum.Wt <- cumsum(Results[, 6])                        #display cumulative sum of Akaike weights
  } else {Results$Cum.Wt <- NULL}
  
  class(Results) <- c("aictab", "data.frame")
  return(Results)
}


##multinomPois
aictab.AICunmarkedFitMPois <- function(cand.set, modnames = NULL, second.ord = TRUE, nobs = NULL, sort = TRUE, c.hat = 1, ...){  #specify whether table should be sorted or not by delta AICc

  ##check if named list if modnames are not supplied
  if(is.null(modnames)) {
    if(is.null(names(cand.set))) {
      modnames <- paste("Mod", 1:length(cand.set), sep = "")
      warning("\nModel names have been supplied automatically in the table\n")
    } else {
      modnames <- names(cand.set)
    }
  }
    
  ##add check for use of c-hat
  #if(c.hat > 1) stop("\nThe correction for overdispersion is not yet implemented for multinomial Poisson models\n")

  Results <- data.frame(Modnames = modnames)                    #assign model names to first column
  Results$K <- unlist(lapply(cand.set, AICc, return.K = TRUE,
                             second.ord = second.ord, nobs = nobs, c.hat = c.hat))     #extract number of parameters
  Results$AICc <- unlist(lapply(cand.set, AICc, return.K = FALSE,
                                second.ord = second.ord, nobs = nobs, c.hat = c.hat))  #extract AICc                                      #
  Results$Delta_AICc <- Results$AICc - min(Results$AICc)            #compute delta AICc
  Results$ModelLik <- exp(-0.5*Results$Delta_AICc)                #compute model likelihood required to compute Akaike weights
  Results$AICcWt <- Results$ModelLik/sum(Results$ModelLik)        #compute Akaike weights

  ##check if some models are redundant
  if(length(unique(Results$AICc)) != length(cand.set)) warning("\nCheck model structure carefully as some models may be redundant\n")


  ##check if AICc and c.hat = 1
  if(second.ord == TRUE && c.hat == 1) {
    Results$LL <- unlist(lapply(X = cand.set, FUN = function(i) extractLL(i)))      
  }
  
  ##rename correctly to QAICc and add column for c-hat
  if(second.ord == TRUE && c.hat > 1) {
    colnames(Results) <- c("Modnames", "K", "QAICc", "Delta_QAICc", "ModelLik", "QAICcWt")
    LL <- unlist(lapply(X = cand.set, FUN = function(i) extractLL(i))) 
    Results$Quasi.LL <- LL/c.hat
    Results$c_hat <- c.hat
  }      

  ##rename correctly to AIC
  if(second.ord == FALSE && c.hat == 1) {
    colnames(Results)<-c("Modnames", "K", "AIC", "Delta_AIC", "ModelLik", "AICWt")
    Results$LL <- unlist(lapply(X = cand.set, FUN = function(i) extractLL(i))) 
  }  
  
  ##rename correctly to QAIC and add column for c-hat
  if(second.ord == FALSE && c.hat > 1) {
    colnames(Results)<-c("Modnames", "K", "QAIC", "Delta_QAIC", "ModelLik", "QAICWt")
    LL <- unlist(lapply(X = cand.set, FUN = function(i) extractLL(i))) 
    Results$Quasi.LL <- LL/c.hat
    Results$c_hat<-c.hat
  }      

  
  if(sort)  {
    Results <- Results[rev(order(Results[, 6])),] 	  #if sort=TRUE, models are ranked based on Akaike weights
    Results$Cum.Wt <- cumsum(Results[, 6])                        #display cumulative sum of Akaike weights
  } else {Results$Cum.Wt <- NULL}
  
  class(Results) <- c("aictab", "data.frame")
  return(Results)
}


##gmultmix
aictab.AICunmarkedFitGMM <- function(cand.set, modnames = NULL, second.ord = TRUE, nobs = NULL, sort = TRUE, c.hat = 1, ...){  #specify whether table should be sorted or not by delta AICc

  ##check if named list if modnames are not supplied
  if(is.null(modnames)) {
    if(is.null(names(cand.set))) {
      modnames <- paste("Mod", 1:length(cand.set), sep = "")
      warning("\nModel names have been supplied automatically in the table\n")
    } else {
      modnames <- names(cand.set)
    }
  }
    

  ##add check for use of c-hat
  #if(c.hat > 1) stop("\nThe correction for overdispersion is not yet implemented for generalized multinomial mixture models\n")

  Results <- data.frame(Modnames = modnames)                    #assign model names to first column
  Results$K <- unlist(lapply(cand.set, AICc, return.K = TRUE,
                             second.ord = second.ord, nobs = nobs, c.hat = c.hat))     #extract number of parameters
  Results$AICc <- unlist(lapply(cand.set, AICc, return.K = FALSE,
                                second.ord = second.ord, nobs = nobs, c.hat = c.hat))  #extract AICc
  Results$Delta_AICc <- Results$AICc - min(Results$AICc)            #compute delta AICc
  Results$ModelLik <- exp(-0.5*Results$Delta_AICc)                #compute model likelihood required to compute Akaike weights
  Results$AICcWt <- Results$ModelLik/sum(Results$ModelLik)        #compute Akaike weights

  ##check if some models are redundant
  if(length(unique(Results$AICc)) != length(cand.set)) warning("\nCheck model structure carefully as some models may be redundant\n")


  ##check if AICc and c.hat = 1
  if(second.ord == TRUE && c.hat == 1) {
    Results$LL <- unlist(lapply(X = cand.set, FUN = function(i) extractLL(i)))      
  }
  
  ##rename correctly to QAICc and add column for c-hat
  if(second.ord == TRUE && c.hat > 1) {
    colnames(Results) <- c("Modnames", "K", "QAICc", "Delta_QAICc", "ModelLik", "QAICcWt")
    LL <- unlist(lapply(X = cand.set, FUN = function(i) extractLL(i))) 
    Results$Quasi.LL <- LL/c.hat
    Results$c_hat <- c.hat
  }      

  ##rename correctly to AIC
  if(second.ord == FALSE && c.hat == 1) {
    colnames(Results) <- c("Modnames", "K", "AIC", "Delta_AIC", "ModelLik", "AICWt")
    Results$LL <- unlist(lapply(X = cand.set, FUN = function(i) extractLL(i))) 
  }  
  
  ##rename correctly to QAIC and add column for c-hat
  if(second.ord == FALSE && c.hat > 1) {
    colnames(Results) <- c("Modnames", "K", "QAIC", "Delta_QAIC", "ModelLik", "QAICWt")
    LL <- unlist(lapply(X = cand.set, FUN = function(i) extractLL(i))) 
    Results$Quasi.LL <- LL/c.hat
    Results$c_hat <- c.hat
  }      


  if(sort)  {
    Results <- Results[rev(order(Results[, 6])),] 	  #if sort=TRUE, models are ranked based on Akaike weights
    Results$Cum.Wt <- cumsum(Results[, 6])                        #display cumulative sum of Akaike weights
  } else {Results$Cum.Wt <- NULL}
  
  class(Results) <- c("aictab", "data.frame")
  return(Results)
}


##gpcount
aictab.AICunmarkedFitGPC <- function(cand.set, modnames = NULL, second.ord = TRUE, nobs = NULL, sort = TRUE, c.hat = 1, ...){  #specify whether table should be sorted or not by delta AICc

  ##check if named list if modnames are not supplied
  if(is.null(modnames)) {
    if(is.null(names(cand.set))) {
      modnames <- paste("Mod", 1:length(cand.set), sep = "")
      warning("\nModel names have been supplied automatically in the table\n")
    } else {
      modnames <- names(cand.set)
    }
  }
    

  ##add check for use of c-hat
  #if(c.hat > 1) stop("\nThe correction for overdispersion is not yet implemented for generalized binomial mixture models\n")

  Results <- data.frame(Modnames = modnames)                    #assign model names to first column
  Results$K <- unlist(lapply(cand.set, AICc, return.K = TRUE,
                             second.ord = second.ord, nobs = nobs, c.hat = c.hat))     #extract number of parameters
  Results$AICc <- unlist(lapply(cand.set, AICc, return.K = FALSE,
                                second.ord = second.ord, nobs = nobs, c.hat = c.hat))  #extract AICc
  Results$Delta_AICc <- Results$AICc - min(Results$AICc)            #compute delta AICc
  Results$ModelLik <- exp(-0.5*Results$Delta_AICc)                #compute model likelihood required to compute Akaike weights
  Results$AICcWt <- Results$ModelLik/sum(Results$ModelLik)        #compute Akaike weights

  ##check if some models are redundant
  if(length(unique(Results$AICc)) != length(cand.set)) warning("\nCheck model structure carefully as some models may be redundant\n")


  ##check if AICc and c.hat = 1
  if(second.ord == TRUE && c.hat == 1) {
    Results$LL <- unlist(lapply(X = cand.set, FUN = function(i) extractLL(i)))      
  }
  
  ##rename correctly to QAICc and add column for c-hat
  if(second.ord == TRUE && c.hat > 1) {
    colnames(Results) <- c("Modnames", "K", "QAICc", "Delta_QAICc", "ModelLik", "QAICcWt")
    LL <- unlist(lapply(X = cand.set, FUN = function(i) extractLL(i))) 
    Results$Quasi.LL <- LL/c.hat
    Results$c_hat <- c.hat
  }      

  ##rename correctly to AIC
  if(second.ord == FALSE && c.hat == 1) {
    colnames(Results) <- c("Modnames", "K", "AIC", "Delta_AIC", "ModelLik", "AICWt")
    Results$LL <- unlist(lapply(X = cand.set, FUN = function(i) extractLL(i))) 
  }  
  
  ##rename correctly to QAIC and add column for c-hat
  if(second.ord == FALSE && c.hat > 1) {
    colnames(Results) <- c("Modnames", "K", "QAIC", "Delta_QAIC", "ModelLik", "QAICWt")
    LL <- unlist(lapply(X = cand.set, FUN = function(i) extractLL(i))) 
    Results$Quasi.LL <- LL/c.hat
    Results$c_hat <- c.hat
  }      

  
  if(sort)  {
    Results <- Results[rev(order(Results[, 6])),] 	  #if sort=TRUE, models are ranked based on Akaike weights
    Results$Cum.Wt <- cumsum(Results[, 6])                        #display cumulative sum of Akaike weights
  } else {Results$Cum.Wt <- NULL}
  
  class(Results) <- c("aictab", "data.frame")
  return(Results)
}



##vglm
aictab.AICvglm <-
  function(cand.set, modnames = NULL, second.ord = TRUE, nobs = NULL, sort = TRUE, c.hat = 1, ...){  #specify whether table should be sorted or not by delta AICc
    ##check if named list if modnames are not supplied
    if(is.null(modnames)) {
      if(is.null(names(cand.set))) {
        modnames <- paste("Mod", 1:length(cand.set), sep = "")
        warning("\nModel names have been supplied automatically in the table\n")
      } else {
        modnames <- names(cand.set)
      }
    }
    

    ##add check to see whether response variable is the same for all models
    check.resp <- lapply(X = cand.set, FUN = function(b) formula(b)[2])
    if(length(unique(check.resp)) > 1) stop("\nYou must use the same response variable for all models\n")

    ##changed to AICcmodavg:::AICc.vglm to avoid conflicts with VGAM
    Results <- data.frame(Modnames=modnames)                    #assign model names to first column
    Results$K <- unlist(lapply(X = cand.set, FUN = AICc, return.K = TRUE, 
                               second.ord = second.ord, nobs = nobs, c.hat = c.hat))     #extract number of parameters
    Results$AICc <- unlist(lapply(X = cand.set, FUN = AICc, return.K = FALSE, 
                                  second.ord = second.ord, nobs = nobs, c.hat = c.hat))  #extract AICc
    Results$Delta_AICc <- Results$AICc - min(Results$AICc)            #compute delta AICc
    Results$ModelLik <- exp(-0.5*Results$Delta_AICc)                #compute model likelihood required to compute Akaike weights
    Results$AICcWt <- Results$ModelLik/sum(Results$ModelLik)        #compute Akaike weights

    ##check if some models are redundant
    if(length(unique(Results$AICc)) != length(cand.set)) warning("\nCheck model structure carefully as some models may be redundant\n")
         
    ##check if AICc and c.hat = 1
    if(second.ord==TRUE && c.hat == 1) {
      Results$LL <- unlist(lapply(X = cand.set, FUN = function(i) logLik(i)))
    }
    
    ##rename correctly to QAICc and add column for c-hat
    if(second.ord == TRUE && c.hat > 1) {
      colnames(Results) <- c("Modnames", "K", "QAICc", "Delta_QAICc", "ModelLik", "QAICcWt")
      LL <- unlist(lapply(X = cand.set, FUN = function(i) logLik(i)))
      Results$Quasi.LL <- LL/c.hat
      Results$c_hat <- c.hat
    }

    ##rename correctly to AIC
    if(second.ord == FALSE && c.hat == 1) {
      colnames(Results)<-c("Modnames", "K", "AIC", "Delta_AIC", "ModelLik", "AICWt")
      Results$LL <- unlist(lapply(X = cand.set, FUN = function(i) logLik(i)))      
    }  

    ##rename correctly to QAIC and add column for c-hat
    if(second.ord == FALSE && c.hat > 1) {
      colnames(Results)<-c("Modnames", "K", "QAIC", "Delta_QAIC", "ModelLik", "QAICWt")
      LL <- unlist(lapply(X = cand.set, FUN = function(i) logLik(i)))
      Results$Quasi.LL <- LL/c.hat
      Results$c_hat <- c.hat
    }      

    
    if(sort)  {
      Results <- Results[rev(order(Results[, 6])),] 	  #if sort=TRUE, models are ranked based on Akaike weights
      Results$Cum.Wt <- cumsum(Results[, 6])                        #display cumulative sum of Akaike weights
    } else {Results$Cum.Wt <- NULL}

    
    class(Results) <- c("aictab", "data.frame")
    return(Results)
  }



##zeroinfl
aictab.AICzeroinfl <-
  function(cand.set, modnames = NULL, second.ord = TRUE, nobs = NULL, sort = TRUE, ...){  #specify whether table should be sorted or not by delta AICc

    ##check if named list if modnames are not supplied
    if(is.null(modnames)) {
      if(is.null(names(cand.set))) {
        modnames <- paste("Mod", 1:length(cand.set), sep = "")
        warning("\nModel names have been supplied automatically in the table\n")
      } else {
        modnames <- names(cand.set)
      }
    }
    
    #if(c.hat != 1) stop("\nThis function does not support overdispersion in \'zeroinfl\' models\n")
    
    ##add check to see whether response variable is the same for all models
    check.resp <- lapply(X = cand.set, FUN = function(b) formula(b)[2])                                        
    if(length(unique(check.resp)) > 1) stop("\nYou must use the same data set for all models\n")
       
    Results <- data.frame(Modnames = modnames)                    #assign model names to first column
    Results$K <- unlist(lapply(cand.set, AICc, return.K = TRUE, second.ord = second.ord, nobs = nobs))     #extract number of parameters
    Results$AICc <- unlist(lapply(cand.set, AICc, return.K = FALSE, second.ord = second.ord, nobs = nobs))  #extract AICc                                      #
    Results$Delta_AICc <- Results$AICc - min(Results$AICc)            #compute delta AICc
    Results$ModelLik <- exp(-0.5*Results$Delta_AICc)                #compute model likelihood required to compute Akaike weights
    Results$AICcWt <- Results$ModelLik/sum(Results$ModelLik)        #compute Akaike weights

   ##check if some models are redundant
    if(length(unique(Results$AICc)) != length(cand.set)) warning("\nCheck model structure carefully as some models may be redundant\n")
    
    Results$LL <- unlist(lapply(X= cand.set, FUN = function(i) logLik(i)))      

    ##rename correctly to AIC
    if(second.ord == FALSE) {
      colnames(Results)[1:6] <- c("Modnames", "K", "AIC", "Delta_AIC", "ModelLik", "AICWt")
    }

    if(sort)  {
      Results <- Results[rev(order(Results[, 6])),] 	  #if sort=TRUE, models are ranked based on Akaike weights
      Results$Cum.Wt <- cumsum(Results[, 6])                        #display cumulative sum of Akaike weights
    } else {Results$Cum.Wt <- NULL}
   
    
    class(Results) <- c("aictab", "data.frame")
    return(Results)
  }



print.aictab <-
  function(x, digits = 2, LL = TRUE, ...) {
    cat("\nModel selection based on ", colnames(x)[3], ":\n", sep = "")
    if (any(names(x) == "c_hat")) {cat("(c-hat estimate = ", x$c_hat[1], ")\n", sep = "")}
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
    

    #if LL==FALSE
    if(identical(LL, FALSE)) {
      names.cols <- colnames(nice.tab)
      sel.LL <- which(attr(regexpr(pattern = "LL", text = names.cols), "match.length") > 1)
      nice.tab <- nice.tab[, -sel.LL]
    }
    
    print(round(nice.tab, digits = digits)) #select rounding off with digits argument
    cat("\n")
  }

