##generic
modavgEffect <- function(cand.set, modnames = NULL, newdata, second.ord = TRUE,
                         nobs = NULL, uncond.se = "revised", conf.level = 0.95,
                         ...){
  cand.set <- formatCands(cand.set)
  UseMethod("modavgEffect", cand.set)
}



##default
modavgEffect.default <- function(cand.set, modnames = NULL, newdata, second.ord = TRUE,
                                 nobs = NULL, uncond.se = "revised", conf.level = 0.95,
                                 ...){
  stop("\nFunction not yet defined for this object class\n")
}



##aov
modavgEffect.AICaov.lm <- function(cand.set, modnames = NULL, newdata, second.ord = TRUE,
                                   nobs = NULL, uncond.se = "revised", conf.level = 0.95,
                                   ...){
  ##newdata is data frame with exact structure of the original data frame (same variable names and type)

  ##check if named list if modnames are not supplied
  if(is.null(modnames)) {
    if(is.null(names(cand.set))) {
      modnames <- paste("Mod", 1:length(cand.set), sep = "")
      warning("\nModel names have been supplied automatically in the table\n")
    } else {
      modnames <- names(cand.set)
    }
  }
    
  ##determine number of observations in new data set
  nobserv <- nrow(newdata)

  if(nobserv > 2) stop("\nCurrent maximum number of groups compared is 2:\nmodify newdata argument accordingly\n")

  ##determine number of columns in new data set
  ncolumns <- ncol(newdata)
        
  ##if only 1 column, add an additional column to avoid problems in computation
  if(ncolumns == 1) newdata$blank.fake.column.NAs <- NA

  
  ##determine which column varies
  uniques <- apply(X = newdata, MARGIN = 2, FUN = unique)
  lengths <- lapply(X = uniques, FUN = length)
  varies <- sapply(X = lengths, FUN = function(i) i > 1)

    
  ##extract name of column
  if(sum(varies) == 1) {
    var.id <- names(varies)[which(varies == TRUE)]
    
    ##determine name of groups compared
    group1 <- as.character(newdata[,paste(var.id)][1])
    group2 <- as.character(newdata[,paste(var.id)][2])

  } else {
    ##warn that no single variable defines groups
    warning("\nGroups do not seem to be defined by a single variable.\n Function proceeding with generic group names\n")
    ##use generic names
    var.id <- "Groups"
    group1 <- "group 1"
    group2 <- "group 2"
  }

  
    
  ##number of models
  nmods <- length(modnames)

  ##compute fitted values
  fit <- matrix(data = unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata)$fit)),
                nrow = nmods, ncol = 2, byrow = TRUE)
  
  ##compute SE's on fitted values
  SE <- matrix(data = unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata)$se.fit)),
               nrow = nmods, ncol = 2, byrow = TRUE)
    
  ##difference between groups 
  diff <- fit[, 1] - fit[, 2]
    
  ##SE on difference
  SE.diff <- sqrt(SE[, 1]^2 + SE[, 2]^2)

  ##store AICc table
  AICctab <- aictab(cand.set = cand.set, modnames = modnames, second.ord = second.ord, nobs = nobs, sort = FALSE)

  ##create object to hold Model-averaged estimates and unconditional SE's
  Mod.avg.out <- matrix(NA, nrow = 1, ncol = 2)
  colnames(Mod.avg.out) <- c("Mod.avg.diff", "Uncond.SE")

  ##begin loop - AICc
  if(second.ord == TRUE){
                   
    ##create temporary data.frame to store fitted values and SE 
    AICctmp <- AICctab
    AICctmp$diff <- diff
    AICctmp$SE.diff <- SE.diff
      
    ##compute model averaged prediction and store in output matrix
    Mod.avg.out[, 1] <- sum(AICctmp$AICcWt*AICctmp$diff)
    ##compute unconditional SE and store in output matrix
      
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Mod.avg.out[, 2] <- sum(AICctmp$AICcWt*sqrt(AICctmp$SE.diff^2 + (AICctmp$diff - Mod.avg.out[, 1])^2))
    }

      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Mod.avg.out[, 2] <- sqrt(sum(AICctmp$AICcWt*(AICctmp$SE.diff^2 + (AICctmp$diff - Mod.avg.out[, 1])^2)))
      }
      ##store table
      AICc.out <- AICctmp
    }
  


    ##create temporary data.frame to store fitted values and SE - AIC
    if(second.ord==FALSE) {
      
      AICtmp <- AICctab
      AICtmp$diff <- diff
      AICtmp$SE.diff <- SE.diff

      ##compute model averaged prediction and store in output matrix
      Mod.avg.out[, 1] <- sum(AICtmp$AICWt*AICtmp$diff)

      ##compute unconditional SE and store in output matrix
      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Mod.avg.out[, 2] <- sum(AICtmp$AICWt*sqrt(AICtmp$SE.diff^2 + (AICtmp$diff - Mod.avg.out[, 1])^2))
      }

      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Mod.avg.out[, 2] <- sqrt(sum(AICtmp$AICWt*(AICtmp$SE.diff^2 + (AICtmp$diff - Mod.avg.out[, 1])^2)))
      }  
      ##store table
      AICc.out <- AICtmp
    }

    ##indicate scale of predictions
    type <- "response"

    zcrit <- qnorm(p=(1-conf.level)/2, lower.tail=FALSE)
    Lower.CL <- Mod.avg.out[, 1] - zcrit * Mod.avg.out[, 2]
    Upper.CL <- Mod.avg.out[, 1] + zcrit * Mod.avg.out[, 2]
    Mod.eff.list <- list("Group.variable" = var.id, "Group1" = group1,
                         "Group2" = group2, "Type" = type, "Mod.avg.table" = AICc.out, "Mod.avg.eff" = Mod.avg.out[,1], 
                          "Uncond.se" = Mod.avg.out[,2], "Conf.level" = conf.level, "Lower.CL" = Lower.CL,
                          "Upper.CL" = Upper.CL)
    class(Mod.eff.list) <- c("modavgEffect", "list")
    return(Mod.eff.list)
}



##glm
modavgEffect.AICglm.lm <- function(cand.set, modnames = NULL, newdata, second.ord = TRUE,
                                   nobs = NULL, uncond.se = "revised", conf.level = 0.95,
                                   type = "response", c.hat = 1, gamdisp = NULL,
                                   ...){
  
  ##check if named list if modnames are not supplied
  if(is.null(modnames)) {
    if(is.null(names(cand.set))) {
      modnames <- paste("Mod", 1:length(cand.set), sep = "")
      warning("\nModel names have been supplied automatically in the table\n")
    } else {
      modnames <- names(cand.set)
    }
  }

  
  ##check on newdata
  ##determine number of observations in new data set
  nobserv <- nrow(newdata)
  if(nobserv > 2) stop("\nCurrent maximum number of groups compared is 2:\nmodify newdata argument accordingly\n")

  ##determine number of columns in new data set
  ncolumns <- ncol(newdata)
  
  ##if only 1 column, add an additional column to avoid problems in computation
  if(ncolumns == 1) newdata$blank.fake.column.NAs <- NA
  
  ##determine which column varies
  uniques <- apply(X = newdata, MARGIN = 2, FUN = unique)
  lengths <- lapply(X = uniques, FUN = length)
  varies <- sapply(X = lengths, FUN = function(i) i > 1)

  ##extract name of column
  if(sum(varies) == 1) {
    var.id <- names(varies)[which(varies == TRUE)]
                
    ##determine name of groups compared
    group1 <- as.character(newdata[,paste(var.id)][1])
    group2 <- as.character(newdata[,paste(var.id)][2])

  } else {
    ##warn that no single variable defines groups
    warning("\nGroups do not seem to be defined by a single variable.\n Function proceeding with generic group names\n")
    ##use generic names
    var.id <- "Groups"
    group1 <- "group 1"
    group2 <- "group 2"
  }

  
  ##newdata is data frame with exact structure of the original data frame (same variable names and type)
  if(type == "terms") {stop("\nThe terms argument is not defined for this function\n")}

  ##check family of glm to avoid problems when requesting predictions with argument 'dispersion'
  fam.type <- unlist(lapply(cand.set, FUN=function(i) family(i)$family))
  fam.unique <- unique(fam.type)
  if(identical(fam.unique, "gaussian")) {
    dispersion <- NULL  #set to NULL if gaussian is used
  } else{dispersion <- c.hat}
  ##poisson and binomial defaults to 1 (no separate parameter for variance)

  ##for negative binomial - reset to NULL
  if(any(regexpr("Negative Binomial", fam.type) != -1)) {
    dispersion <- NULL
    ##check for mixture of negative binomial and other
    ##number of models with negative binomial
    negbin.num <- sum(regexpr("Negative Binomial", fam.type) != -1)
    if(negbin.num < length(fam.type)) {
      stop("Function does not support mixture of negative binomial with other distribution")
    }
  }
  
  
    
###################CHANGES####
##############################
  if(c.hat > 1) {dispersion <- c.hat }
  if(!is.null(gamdisp)) {dispersion <- gamdisp}
  if(c.hat > 1 && !is.null(gamdisp)) {stop("\nYou cannot specify values for both \'c.hat\' and \'gamdisp\'\n")}
  ##dispersion is the dispersion parameter - this influences the SE's (to specify dispersion parameter for either overdispersed Poisson or Gamma glm)
  ##type enables to specify either "response" (original scale = point estimate) or "link" (linear predictor)
  
  ##check if object is of "lm" or "glm" class
  ##extract classes
  mod.class <- unlist(lapply(X = cand.set, FUN = class))
  ##check if all are identical
  check.class <- unique(mod.class)

  ##check that link function is the same for all models if linear predictor is used
  if(identical(type, "link")) {
    check.link <- unlist(lapply(X = cand.set, FUN = function(i) i$family$link))
    unique.link <- unique(x = check.link)
    if(length(unique.link) > 1) {stop("\nIt is not appropriate to compute a model averaged beta estimate\n",
                                          "with different link functions\n")}
  }

 
  ##check if model uses gamma distribution
  gam1 <- unlist(lapply(cand.set, FUN = function(i) family(i)$family[1] == "Gamma")) #check for gamma regression models
  ##correct SE's for estimates of gamma regressions when gamdisp is specified
  if(any(gam1) == TRUE)  {
    ##check for specification of gamdisp argument
    if(is.null(gamdisp)) stop("\nYou must specify a gamma dispersion parameter with gamma generalized linear models\n")
  }
  
  ##number of models
  nmods <- length(modnames)

  ##compute fitted values
  fit <- matrix(data = unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata, type = type,
                                                     dispersion = dispersion)$fit)), nrow = nmods, ncol = 2, byrow = TRUE)
    
  ##compute SE's on fitted values
  SE <- matrix(data = unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata, type = type,
                                                    dispersion = dispersion)$se.fit)), nrow = nmods, ncol = 2, byrow = TRUE)
  
  ##difference between groups 
  diff <- fit[, 1] - fit[, 2]
  
  ##SE on difference
  SE.diff <- sqrt(SE[, 1]^2 + SE[, 2]^2)
  
  ##store AICc table
  AICctab <- aictab(cand.set = cand.set, modnames = modnames, second.ord = second.ord,
                    nobs = nobs, sort = FALSE, c.hat = c.hat)

  ##create object to hold Model-averaged estimates and unconditional SE's
  Mod.avg.out <- matrix(NA, nrow = 1, ncol = 2)
  colnames(Mod.avg.out) <- c("Mod.avg.diff", "Uncond.SE")

  ##begin loop - AICc
  if(second.ord == TRUE && c.hat == 1){
             
    ##create temporary data.frame to store fitted values and SE 
    AICctmp <- AICctab
    AICctmp$diff <- diff
    AICctmp$SE.diff <- SE.diff
    
    ##compute model averaged prediction and store in output matrix
    Mod.avg.out[, 1] <- sum(AICctmp$AICcWt*AICctmp$diff)
    ##compute unconditional SE and store in output matrix
      
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Mod.avg.out[, 2] <- sum(AICctmp$AICcWt*sqrt(AICctmp$SE.diff^2 + (AICctmp$diff - Mod.avg.out[, 1])^2))
    }

    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Mod.avg.out[, 2] <- sqrt(sum(AICctmp$AICcWt*(AICctmp$SE.diff^2 + (AICctmp$diff - Mod.avg.out[, 1])^2)))
    }
    ##store table
    AICc.out <- AICctmp
  }
  
  ##create temporary data.frame to store fitted values and SE - QAICc
  if(second.ord==TRUE && c.hat > 1) {
    
    QAICctmp <- AICctab
    QAICctmp$diff <- diff
    QAICctmp$SE.diff <- SE.diff
    
    ##compute model averaged prediction and store in output matrix
    Mod.avg.out[, 1] <- sum(QAICctmp$QAICcWt*QAICctmp$diff)
      
    ##compute unconditional SE and store in output matrix
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Mod.avg.out[, 2] <- sum(QAICctmp$QAICcWt*sqrt(QAICctmp$SE.diff^2 + (QAICctmp$diff - Mod.avg.out[, 1])^2))
    }
      
    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Mod.avg.out[, 2] <- sqrt(sum(QAICctmp$QAICcWt*(QAICctmp$SE.diff^2 + (QAICctmp$diff - Mod.avg.out[, 1])^2)))  
    }
    ##store table
    AICc.out <- QAICctmp
    
  }

  ##create temporary data.frame to store fitted values and SE - AIC
  if(second.ord == FALSE && c.hat == 1) {
    
    AICtmp <- AICctab
    AICtmp$diff <- diff
    AICtmp$SE.diff <- SE.diff

    ##compute model averaged prediction and store in output matrix
    Mod.avg.out[, 1] <- sum(AICtmp$AICWt*AICtmp$diff)
    
    ##compute unconditional SE and store in output matrix
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Mod.avg.out[, 2] <- sum(AICtmp$AICWt*sqrt(AICtmp$SE.diff^2 + (AICtmp$diff - Mod.avg.out[, 1])^2))
    }

    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Mod.avg.out[, 2] <- sqrt(sum(AICtmp$AICWt*(AICtmp$SE.diff^2 + (AICtmp$diff - Mod.avg.out[, 1])^2)))
    }  
    ##store table
    AICc.out <- AICtmp
  }
  
  ##create temporary data.frame to store fitted values and SE - QAIC
  if(second.ord == FALSE && c.hat > 1) {
          
    QAICtmp <- AICctab
    QAICtmp$diff <- diff
    QAICtmp$SE.diff <- SE.diff
      
    ##compute model averaged prediction and store in output matrix
    Mod.avg.out[, 1] <- sum(QAICtmp$QAICWt*QAICtmp$diff)
      
    ##compute unconditional SE and store in output matrix
    if(identical(uncond.se, "old")) {
      Mod.avg.out[, 2] <- sum(QAICtmp$QAICWt*sqrt(QAICtmp$SE.diff^2 + (QAICtmp$diff - Mod.avg.out[, 1])^2))
    }
      
    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Mod.avg.out[, 2] <- sqrt(sum(QAICtmp$QAICWt*(QAICtmp$SE.diff^2 + (QAICtmp$diff - Mod.avg.out[, 1])^2)))
    }  
    ##store table
    AICc.out <- QAICtmp
  }

  zcrit <- qnorm(p=(1-conf.level)/2, lower.tail=FALSE)
  Lower.CL <- Mod.avg.out[, 1] - zcrit * Mod.avg.out[, 2]
  Upper.CL <- Mod.avg.out[, 1] + zcrit * Mod.avg.out[, 2]
  Mod.eff.list <- list("Group.variable" = var.id, "Group1" = group1,
                       "Group2" = group2, "Type" = type, "Mod.avg.table" = AICc.out, "Mod.avg.eff" = Mod.avg.out[,1],
                       "Uncond.se" = Mod.avg.out[,2], "Conf.level" = conf.level, "Lower.CL" = Lower.CL,
                       "Upper.CL" = Upper.CL)
  class(Mod.eff.list) <- c("modavgEffect", "list")
  return(Mod.eff.list)
}



##gls
modavgEffect.AICgls <- function(cand.set, modnames = NULL, newdata, second.ord = TRUE,
                                nobs = NULL, uncond.se = "revised", conf.level = 0.95,
                                ...){
  ##newdata is data frame with exact structure of the original data frame (same variable names and type)

  ##check if named list if modnames are not supplied
  if(is.null(modnames)) {
    if(is.null(names(cand.set))) {
      modnames <- paste("Mod", 1:length(cand.set), sep = "")
      warning("\nModel names have been supplied automatically in the table\n")
    } else {
      modnames <- names(cand.set)
    }
  }
    
    ##determine number of observations in new data set
    nobserv <- nrow(newdata)

    if(nobserv > 2) stop("\nCurrent maximum number of groups compared is 2:\nmodify newdata argument accordingly\n")

    ##determine number of columns in new data set
    ncolumns <- ncol(newdata)
        
    ##if only 1 column, add an additional column to avoid problems in computation
    if(ncolumns == 1) newdata$blank.fake.column.NAs <- NA

  
  ##determine which column varies
  uniques <- apply(X = newdata, MARGIN = 2, FUN = unique)
  lengths <- lapply(X = uniques, FUN = length)
  varies <- sapply(X = lengths, FUN = function(i) i > 1)

    
  ##extract name of column
  if(sum(varies) == 1) {
    var.id <- names(varies)[which(varies == TRUE)]
                
    ##determine name of groups compared
    group1 <- as.character(newdata[,paste(var.id)][1])
    group2 <- as.character(newdata[,paste(var.id)][2])

  } else {
    ##warn that no single variable defines groups
    warning("\nGroups do not seem to be defined by a single variable.\n Function proceeding with generic group names\n")
    ##use generic names
    var.id <- "Groups"
    group1 <- "group 1"
    group2 <- "group 2"
  }

    
    ##number of models
    nmods <- length(modnames)

    ##compute fitted values
    fit <- matrix(data = unlist(lapply(X = cand.set, FUN = function(i)predictSE(i, se.fit = TRUE, newdata = newdata)$fit)),
                  nrow = nmods, ncol = 2, byrow = TRUE)
    
    ##compute SE's on fitted values
    SE <- matrix(data = unlist(lapply(X = cand.set, FUN = function(i)predictSE(i, se.fit = TRUE, newdata = newdata)$se.fit)),
                 nrow = nmods, ncol = 2, byrow = TRUE)
    
    ##difference between groups 
    diff <- fit[, 1] - fit[, 2]
    
    ##SE on difference
    SE.diff <- sqrt(SE[, 1]^2 + SE[, 2]^2)


    
    ##store AICc table
    AICctab <- aictab(cand.set = cand.set, modnames = modnames, second.ord = second.ord, nobs = nobs, sort = FALSE)

    #create object to hold Model-averaged estimates and unconditional SE's
    Mod.avg.out <- matrix(NA, nrow = 1, ncol = 2)
    colnames(Mod.avg.out) <- c("Mod.avg.diff", "Uncond.SE")

    ##begin loop - AICc
    if(second.ord == TRUE){
                   
      ##create temporary data.frame to store fitted values and SE 
      AICctmp <- AICctab
      AICctmp$diff <- diff
      AICctmp$SE.diff <- SE.diff
      
      ##compute model averaged prediction and store in output matrix
      Mod.avg.out[, 1] <- sum(AICctmp$AICcWt*AICctmp$diff)
      ##compute unconditional SE and store in output matrix
      
      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Mod.avg.out[, 2] <- sum(AICctmp$AICcWt*sqrt(AICctmp$SE.diff^2 + (AICctmp$diff - Mod.avg.out[, 1])^2))
      }

      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Mod.avg.out[, 2] <- sqrt(sum(AICctmp$AICcWt*(AICctmp$SE.diff^2 + (AICctmp$diff - Mod.avg.out[, 1])^2)))
      }
      ##store table
      AICc.out <- AICctmp
    }
  


    ##create temporary data.frame to store fitted values and SE - AIC
    if(second.ord==FALSE) {
      
      AICtmp <- AICctab
      AICtmp$diff <- diff
      AICtmp$SE.diff <- SE.diff

      ##compute model averaged prediction and store in output matrix
      Mod.avg.out[, 1] <- sum(AICtmp$AICWt*AICtmp$diff)

      ##compute unconditional SE and store in output matrix
      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Mod.avg.out[, 2] <- sum(AICtmp$AICWt*sqrt(AICtmp$SE.diff^2 + (AICtmp$diff - Mod.avg.out[, 1])^2))
      }

      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Mod.avg.out[, 2] <- sqrt(sum(AICtmp$AICWt*(AICtmp$SE.diff^2 + (AICtmp$diff - Mod.avg.out[, 1])^2)))
      }  
      ##store table
      AICc.out <- AICtmp
    }

    ##indicate scale of predictions
    type <- "response"

    zcrit <- qnorm(p=(1-conf.level)/2, lower.tail=FALSE)
    Lower.CL <- Mod.avg.out[, 1] - zcrit * Mod.avg.out[, 2]
    Upper.CL <- Mod.avg.out[, 1] + zcrit * Mod.avg.out[, 2]
    Mod.eff.list <- list("Group.variable" = var.id, "Group1" = group1,
                         "Group2" = group2, "Type" = type, "Mod.avg.table" = AICc.out, "Mod.avg.eff" = Mod.avg.out[,1], 
                          "Uncond.se" = Mod.avg.out[,2], "Conf.level" = conf.level, "Lower.CL" = Lower.CL,
                          "Upper.CL" = Upper.CL)
    class(Mod.eff.list) <- c("modavgEffect", "list")
    return(Mod.eff.list)
}



##lm
modavgEffect.AIClm <- function(cand.set, modnames = NULL, newdata, second.ord = TRUE,
                               nobs = NULL, uncond.se = "revised", conf.level = 0.95,
                               ...){
  ##newdata is data frame with exact structure of the original data frame (same variable names and type)

  ##check if named list if modnames are not supplied
  if(is.null(modnames)) {
    if(is.null(names(cand.set))) {
      modnames <- paste("Mod", 1:length(cand.set), sep = "")
      warning("\nModel names have been supplied automatically in the table\n")
    } else {
      modnames <- names(cand.set)
    }
  }
    
    ##determine number of observations in new data set
    nobserv <- nrow(newdata)

    if(nobserv > 2) stop("\nCurrent maximum number of groups compared is 2:\nmodify newdata argument accordingly\n")

    ##determine number of columns in new data set
    ncolumns <- ncol(newdata)
        
    ##if only 1 column, add an additional column to avoid problems in computation
    if(ncolumns == 1) newdata$blank.fake.column.NAs <- NA

  
  ##determine which column varies
  uniques <- apply(X = newdata, MARGIN = 2, FUN = unique)
  lengths <- lapply(X = uniques, FUN = length)
  varies <- sapply(X = lengths, FUN = function(i) i > 1)

    
  ##extract name of column
  if(sum(varies) == 1) {
    var.id <- names(varies)[which(varies == TRUE)]
                
    ##determine name of groups compared
    group1 <- as.character(newdata[,paste(var.id)][1])
    group2 <- as.character(newdata[,paste(var.id)][2])

  } else {
    ##warn that no single variable defines groups
    warning("\nGroups do not seem to be defined by a single variable.\n Function proceeding with generic group names\n")
    ##use generic names
    var.id <- "Groups"
    group1 <- "group 1"
    group2 <- "group 2"
  }

    
    ##number of models
    nmods <- length(modnames)

    ##compute fitted values
    fit <- matrix(data = unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata)$fit)),
                  nrow = nmods, ncol = 2, byrow = TRUE)
    
    ##compute SE's on fitted values
    SE <- matrix(data = unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata)$se.fit)),
                 nrow = nmods, ncol = 2, byrow = TRUE)
    
    ##difference between groups 
    diff <- fit[, 1] - fit[, 2]
    
    ##SE on difference
    SE.diff <- sqrt(SE[, 1]^2 + SE[, 2]^2)

    ##store AICc table
    AICctab <- aictab(cand.set = cand.set, modnames = modnames, second.ord = second.ord, nobs = nobs, sort = FALSE)

    #create object to hold Model-averaged estimates and unconditional SE's
    Mod.avg.out <- matrix(NA, nrow = 1, ncol = 2)
    colnames(Mod.avg.out) <- c("Mod.avg.diff", "Uncond.SE")

    ##begin loop - AICc
    if(second.ord == TRUE){
                   
      ##create temporary data.frame to store fitted values and SE 
      AICctmp <- AICctab
      AICctmp$diff <- diff
      AICctmp$SE.diff <- SE.diff
      
      ##compute model averaged prediction and store in output matrix
      Mod.avg.out[, 1] <- sum(AICctmp$AICcWt*AICctmp$diff)
      ##compute unconditional SE and store in output matrix
      
      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Mod.avg.out[, 2] <- sum(AICctmp$AICcWt*sqrt(AICctmp$SE.diff^2 + (AICctmp$diff - Mod.avg.out[, 1])^2))
      }

      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Mod.avg.out[, 2] <- sqrt(sum(AICctmp$AICcWt*(AICctmp$SE.diff^2 + (AICctmp$diff - Mod.avg.out[, 1])^2)))
      }
      ##store table
      AICc.out <- AICctmp
    }
  


    ##create temporary data.frame to store fitted values and SE - AIC
    if(second.ord==FALSE) {
      
      AICtmp <- AICctab
      AICtmp$diff <- diff
      AICtmp$SE.diff <- SE.diff

      ##compute model averaged prediction and store in output matrix
      Mod.avg.out[, 1] <- sum(AICtmp$AICWt*AICtmp$diff)

      ##compute unconditional SE and store in output matrix
      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Mod.avg.out[, 2] <- sum(AICtmp$AICWt*sqrt(AICtmp$SE.diff^2 + (AICtmp$diff - Mod.avg.out[, 1])^2))
      }

      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Mod.avg.out[, 2] <- sqrt(sum(AICtmp$AICWt*(AICtmp$SE.diff^2 + (AICtmp$diff - Mod.avg.out[, 1])^2)))
      }  
      ##store table
      AICc.out <- AICtmp
    }

    ##indicate scale of predictions
    type <- "response"

    zcrit <- qnorm(p=(1-conf.level)/2, lower.tail=FALSE)
    Lower.CL <- Mod.avg.out[, 1] - zcrit * Mod.avg.out[, 2]
    Upper.CL <- Mod.avg.out[, 1] + zcrit * Mod.avg.out[, 2]
    Mod.eff.list <- list("Group.variable" = var.id, "Group1" = group1,
                         "Group2" = group2, "Type" = type, "Mod.avg.table" = AICc.out, "Mod.avg.eff" = Mod.avg.out[,1], 
                          "Uncond.se" = Mod.avg.out[,2], "Conf.level" = conf.level, "Lower.CL" = Lower.CL,
                          "Upper.CL" = Upper.CL)
    class(Mod.eff.list) <- c("modavgEffect", "list")
    return(Mod.eff.list)
}



##lme
modavgEffect.AIClme <-
function(cand.set, modnames = NULL, newdata, second.ord = TRUE, nobs = NULL,
         uncond.se = "revised", conf.level = 0.95, ...){
  ##newdata is data frame with exact structure of the original data frame (same variable names and type)

  ##check if named list if modnames are not supplied
  if(is.null(modnames)) {
    if(is.null(names(cand.set))) {
      modnames <- paste("Mod", 1:length(cand.set), sep = "")
      warning("\nModel names have been supplied automatically in the table\n")
    } else {
      modnames <- names(cand.set)
    }
  }
    
  ##check on newdata
  ##determine number of observations in new data set
  nobserv <- nrow(newdata)
  if(nobserv > 2) stop("\nCurrent maximum number of groups compared is 2:\nmodify newdata argument accordingly\n")

  ##determine number of columns in new data set
  ncolumns <- ncol(newdata)

  ##if only 1 column, add an additional column to avoid problems in computation
  if(ncolumns == 1) newdata$blank.fake.column.NAs <- NA

  
  ##determine which column varies
  uniques <- apply(X = newdata, MARGIN = 2, FUN = unique)
  lengths <- lapply(X = uniques, FUN = length)
  varies <- sapply(X = lengths, FUN = function(i) i > 1)

    
  ##extract name of column
  if(sum(varies) == 1) {
    var.id <- names(varies)[which(varies == TRUE)]
                
    ##determine name of groups compared
    group1 <- as.character(newdata[,paste(var.id)][1])
    group2 <- as.character(newdata[,paste(var.id)][2])

  } else {
    ##warn that no single variable defines groups
    warning("\nGroups do not seem to be defined by a single variable.\n Function proceeding with generic group names\n")
    ##use generic names
    var.id <- "Groups"
    group1 <- "group 1"
    group2 <- "group 2"
  }
  
    
  ##number of models
  nmods <- length(modnames)

  ##compute fitted values
  fit <- matrix(data = unlist(lapply(X = cand.set, FUN = function(i)predictSE(i, se.fit = TRUE, newdata = newdata)$fit)),
                nrow = nmods, ncol = 2, byrow = TRUE)
    
  ##compute SE's on fitted values
  SE <- matrix(data = unlist(lapply(X = cand.set, FUN = function(i)predictSE(i, se.fit = TRUE, newdata = newdata)$se.fit)),
               nrow = nmods, ncol = 2, byrow = TRUE)
    
  ##difference between groups 
  diff <- fit[, 1] - fit[, 2]
    
  ##SE on difference
  SE.diff <- sqrt(SE[, 1]^2 + SE[, 2]^2)
  
  ##store AICc table
  AICctab <- aictab(cand.set = cand.set, modnames = modnames, second.ord = second.ord, nobs = nobs, sort = FALSE)
  
  
  ##create object to hold Model-averaged estimates and unconditional SE's
  Mod.avg.out <- matrix(NA, nrow = 1, ncol = 2)
  colnames(Mod.avg.out) <- c("Mod.avg.diff", "Uncond.SE")
  
  ##begin loop - AICc
  if(second.ord == TRUE){
    
    ##create temporary data.frame to store fitted values and SE 
    AICctmp <- AICctab
    AICctmp$diff <- diff
    AICctmp$SE.diff <- SE.diff
    
    ##compute model averaged prediction and store in output matrix
    Mod.avg.out[, 1] <- sum(AICctmp$AICcWt*AICctmp$diff)
    ##compute unconditional SE and store in output matrix
    
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Mod.avg.out[, 2] <- sum(AICctmp$AICcWt*sqrt(AICctmp$SE.diff^2 + (AICctmp$diff - Mod.avg.out[, 1])^2))
    }

    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Mod.avg.out[, 2] <- sqrt(sum(AICctmp$AICcWt*(AICctmp$SE.diff^2 + (AICctmp$diff - Mod.avg.out[, 1])^2)))
    }
    ##store table
    AICc.out <- AICctmp
  }
  
  ##create temporary data.frame to store fitted values and SE - AIC
  if(second.ord==FALSE) {
    
    AICtmp <- AICctab
    AICtmp$diff <- diff
    AICtmp$SE.diff <- SE.diff

    ##compute model averaged prediction and store in output matrix
    Mod.avg.out[, 1] <- sum(AICtmp$AICWt*AICtmp$diff)
    
    ##compute unconditional SE and store in output matrix
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Mod.avg.out[, 2] <- sum(AICtmp$AICWt*sqrt(AICtmp$SE.diff^2 + (AICtmp$diff - Mod.avg.out[, 1])^2))
    }

    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Mod.avg.out[, 2] <- sqrt(sum(AICtmp$AICWt*(AICtmp$SE.diff^2 + (AICtmp$diff - Mod.avg.out[, 1])^2)))
    }  
    ##store table
    AICc.out <- AICtmp
  }
  
  ##indicate scale of predictions
  type <- "response"
  
  zcrit <- qnorm(p=(1-conf.level)/2, lower.tail=FALSE)
  Lower.CL <- Mod.avg.out[, 1] - zcrit * Mod.avg.out[, 2]
  Upper.CL <- Mod.avg.out[, 1] + zcrit * Mod.avg.out[, 2]
  Mod.eff.list <- list("Group.variable" = var.id, "Group1" = group1,
                       "Group2" = group2, "Type" = type, "Mod.avg.table" = AICc.out, "Mod.avg.eff" = Mod.avg.out[,1],
                       "Uncond.se" = Mod.avg.out[,2], "Conf.level" = conf.level, "Lower.CL" = Lower.CL,
                       "Upper.CL" = Upper.CL)
  class(Mod.eff.list) <- c("modavgEffect", "list")
  return(Mod.eff.list)
}



##mer  - lme4 version < 1
modavgEffect.AICmer <- function(cand.set, modnames = NULL, newdata, second.ord = TRUE,
                                nobs = NULL, uncond.se = "revised", conf.level = 0.95,
                                type = "response", ...) {

  ##check if named list if modnames are not supplied
  if(is.null(modnames)) {
    if(is.null(names(cand.set))) {
      modnames <- paste("Mod", 1:length(cand.set), sep = "")
      warning("\nModel names have been supplied automatically in the table\n")
    } else {
      modnames <- names(cand.set)
    }
  }
    
  ##check on newdata
  ##determine number of observations in new data set
  nobserv <- nrow(newdata)
  if(nobserv > 2) stop("\nCurrent maximum number of groups compared is 2:\nmodify newdata argument accordingly\n")

  ##determine number of columns in new data set
  ncolumns <- ncol(newdata)
  
  ##if only 1 column, add an additional column to avoid problems in computation
  if(ncolumns == 1) newdata$blank.fake.column.NAs <- NA
  

  ##determine which column varies
  uniques <- apply(X = newdata, MARGIN = 2, FUN = unique)
  lengths <- lapply(X = uniques, FUN = length)
  varies <- sapply(X = lengths, FUN = function(i) i > 1)

    
  ##extract name of column
  if(sum(varies) == 1) {
    var.id <- names(varies)[which(varies == TRUE)]
                
    ##determine name of groups compared
    group1 <- as.character(newdata[,paste(var.id)][1])
    group2 <- as.character(newdata[,paste(var.id)][2])

  } else {
    ##warn that no single variable defines groups
    warning("\nGroups do not seem to be defined by a single variable.\n Function proceeding with generic group names\n")
    ##use generic names
    var.id <- "Groups"
    group1 <- "group 1"
    group2 <- "group 2"
  }

    
  ##extract classes
  mod.class <- unlist(lapply(X=cand.set, FUN=class))
  ##check if all are identical
  check.class <- unique(mod.class)

  ##check that link function is the same for all models if linear predictor is used
  if(identical(type, "link")) {
    link.list <- unlist(lapply(X = cand.set, FUN = function(i) fam.link.mer(i)$link))
    check.link <- unique(link.list)
    if(length(check.link) > 1) stop("\nIt is not appropriate to compute a model-averaged beta estimate\n",
                                        "from models using different link functions\n")
  }

 
  ##number of models
  nmods <- length(modnames)

    
  ##compute fitted values
  fit <- matrix(data = unlist(lapply(X = cand.set, FUN = function(i)predictSE(i, se.fit = TRUE, newdata = newdata, type = type)$fit)),
                nrow = nmods, ncol = 2, byrow = TRUE)
    
  ##compute SE's on fitted values
  SE <- matrix(data = unlist(lapply(X = cand.set, FUN = function(i)predictSE(i, se.fit = TRUE, newdata = newdata, type = type)$se.fit)),
               nrow = nmods, ncol = 2, byrow = TRUE)
  
  ##difference between groups 
  diff <- fit[, 1] - fit[, 2]
    
  ##SE on difference
  SE.diff <- sqrt(SE[, 1]^2 + SE[, 2]^2)
  
  
  ##store AICc table
  AICctab <- aictab(cand.set = cand.set, modnames = modnames, second.ord = second.ord,
                    nobs = nobs, sort = FALSE)

  ##create object to hold Model-averaged estimates and unconditional SE's
  Mod.avg.out <- matrix(NA, nrow = 1, ncol = 2)
  colnames(Mod.avg.out) <- c("Mod.avg.diff", "Uncond.SE")


  ##begin loop - AICc
  if(second.ord==TRUE){
    ##create temporary data.frame to store fitted values and SE 
    AICctmp <- AICctab
    AICctmp$diff <- diff
    AICctmp$SE.diff <- SE.diff
    
    ##compute model averaged prediction and store in output matrix
    Mod.avg.out[, 1] <- sum(AICctmp$AICcWt*AICctmp$diff)
    ##compute unconditional SE and store in output matrix
      
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Mod.avg.out[, 2] <- sum(AICctmp$AICcWt*sqrt(AICctmp$SE.diff^2 + (AICctmp$diff - Mod.avg.out[, 1])^2))
    }

    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Mod.avg.out[, 2] <- sqrt(sum(AICctmp$AICcWt*(AICctmp$SE.diff^2 + (AICctmp$diff - Mod.avg.out[, 1])^2)))
    }
    ##store table
    AICc.out <- AICctmp
  }
  

  ##create temporary data.frame to store fitted values and SE - AIC
  if(second.ord==FALSE) {
    AICtmp <- AICctab
    AICtmp$diff <- diff
    AICtmp$SE.diff <- SE.diff

    ##compute model averaged prediction and store in output matrix
    Mod.avg.out[, 1] <- sum(AICtmp$AICWt*AICtmp$diff)

    ##compute unconditional SE and store in output matrix
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Mod.avg.out[, 2] <- sum(AICtmp$AICWt*sqrt(AICtmp$SE.diff^2 + (AICtmp$diff - Mod.avg.out[, 1])^2))
    }

    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Mod.avg.out[, 2] <- sqrt(sum(AICtmp$AICWt*(AICtmp$SE.diff^2 + (AICtmp$diff - Mod.avg.out[, 1])^2)))
    }  
    ##store table
    AICc.out <- AICtmp
  }
  
  zcrit <- qnorm(p=(1-conf.level)/2, lower.tail=FALSE)
  Lower.CL <- Mod.avg.out[, 1] - zcrit * Mod.avg.out[, 2]
  Upper.CL <- Mod.avg.out[, 1] + zcrit * Mod.avg.out[, 2]
  Mod.eff.list <- list("Group.variable" = var.id, "Group1" = group1,
                        "Group2" = group2, "Type" = type, "Mod.avg.table" = AICc.out, "Mod.avg.eff" = Mod.avg.out[,1],
                        "Uncond.se" = Mod.avg.out[,2],  "Conf.level" = conf.level, "Lower.CL" = Lower.CL,
                        "Upper.CL" = Upper.CL)
  class(Mod.eff.list) <- c("modavgEffect", "list")
  return(Mod.eff.list)
}



##glmerMod
modavgEffect.AICglmerMod <- function(cand.set, modnames = NULL, newdata, second.ord = TRUE,
                                     nobs = NULL, uncond.se = "revised", conf.level = 0.95,
                                     type = "response", ...) {

  ##check if named list if modnames are not supplied
    if(is.null(modnames)) {
      if(is.null(names(cand.set))) {
        modnames <- paste("Mod", 1:length(cand.set), sep = "")
        warning("\nModel names have been supplied automatically in the table\n")
      }
      modnames <- names(cand.set)
    }
    
  ##newdata is data frame with exact structure of the original data frame (same variable names and type)
  ##check on newdata
  ##determine number of observations in new data set
  nobserv <- nrow(newdata)
  if(nobserv > 2) stop("\nCurrent maximum number of groups compared is 2:\nmodify newdata argument accordingly\n")

  ##determine number of columns in new data set
  ncolumns <- ncol(newdata)
  
  ##if only 1 column, add an additional column to avoid problems in computation
  if(ncolumns == 1) newdata$blank.fake.column.NAs <- NA

    
  ##determine which column varies
  uniques <- apply(X = newdata, MARGIN = 2, FUN = unique)
  lengths <- lapply(X = uniques, FUN = length)
  varies <- sapply(X = lengths, FUN = function(i) i > 1)

    
  ##extract name of column
  if(sum(varies) == 1) {
    var.id <- names(varies)[which(varies == TRUE)]
                
    ##determine name of groups compared
    group1 <- as.character(newdata[,paste(var.id)][1])
    group2 <- as.character(newdata[,paste(var.id)][2])

  } else {
    ##warn that no single variable defines groups
    warning("\nGroups do not seem to be defined by a single variable.\n Function proceeding with generic group names\n")
    ##use generic names
    var.id <- "Groups"
    group1 <- "group 1"
    group2 <- "group 2"
  }

    
    
  ##extract classes
  mod.class <- unlist(lapply(X=cand.set, FUN=class))
  ##check if all are identical
  check.class <- unique(mod.class)

  ##check that link function is the same for all models if linear predictor is used
  if(identical(type, "link")) {
    link.list <- unlist(lapply(X = cand.set, FUN = function(i) fam.link.mer(i)$link))
    check.link <- unique(link.list)
    if(length(check.link) > 1) stop("\nIt is not appropriate to compute a model-averaged beta estimate\n",
                                        "from models using different link functions\n")
  }

 
  ##number of models
  nmods <- length(modnames)

    
  ##compute fitted values
  fit <- matrix(data = unlist(lapply(X = cand.set, FUN = function(i)predictSE(i, se.fit = TRUE, newdata = newdata, type = type)$fit)),
                nrow = nmods, ncol = 2, byrow = TRUE)
    
  ##compute SE's on fitted values
  SE <- matrix(data = unlist(lapply(X = cand.set, FUN = function(i)predictSE(i, se.fit = TRUE, newdata = newdata, type = type)$se.fit)),
               nrow = nmods, ncol = 2, byrow = TRUE)
  
  ##difference between groups 
  diff <- fit[, 1] - fit[, 2]
    
  ##SE on difference
  SE.diff <- sqrt(SE[, 1]^2 + SE[, 2]^2)
  
  
  ##store AICc table
  AICctab <- aictab(cand.set = cand.set, modnames = modnames, second.ord = second.ord,
                    nobs = nobs, sort = FALSE)

  ##create object to hold Model-averaged estimates and unconditional SE's
  Mod.avg.out <- matrix(NA, nrow = 1, ncol = 2)
  colnames(Mod.avg.out) <- c("Mod.avg.diff", "Uncond.SE")


  ##begin loop - AICc
  if(second.ord==TRUE){
    ##create temporary data.frame to store fitted values and SE 
    AICctmp <- AICctab
    AICctmp$diff <- diff
    AICctmp$SE.diff <- SE.diff
    
    ##compute model averaged prediction and store in output matrix
    Mod.avg.out[, 1] <- sum(AICctmp$AICcWt*AICctmp$diff)
    ##compute unconditional SE and store in output matrix
      
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Mod.avg.out[, 2] <- sum(AICctmp$AICcWt*sqrt(AICctmp$SE.diff^2 + (AICctmp$diff - Mod.avg.out[, 1])^2))
    }

    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Mod.avg.out[, 2] <- sqrt(sum(AICctmp$AICcWt*(AICctmp$SE.diff^2 + (AICctmp$diff - Mod.avg.out[, 1])^2)))
    }
    ##store table
    AICc.out <- AICctmp
  }
  

  ##create temporary data.frame to store fitted values and SE - AIC
  if(second.ord==FALSE) {
    AICtmp <- AICctab
    AICtmp$diff <- diff
    AICtmp$SE.diff <- SE.diff

    ##compute model averaged prediction and store in output matrix
    Mod.avg.out[, 1] <- sum(AICtmp$AICWt*AICtmp$diff)

    ##compute unconditional SE and store in output matrix
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Mod.avg.out[, 2] <- sum(AICtmp$AICWt*sqrt(AICtmp$SE.diff^2 + (AICtmp$diff - Mod.avg.out[, 1])^2))
    }

    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Mod.avg.out[, 2] <- sqrt(sum(AICtmp$AICWt*(AICtmp$SE.diff^2 + (AICtmp$diff - Mod.avg.out[, 1])^2)))
    }  
    ##store table
    AICc.out <- AICtmp
  }
  
  zcrit <- qnorm(p=(1-conf.level)/2, lower.tail=FALSE)
  Lower.CL <- Mod.avg.out[, 1] - zcrit * Mod.avg.out[, 2]
  Upper.CL <- Mod.avg.out[, 1] + zcrit * Mod.avg.out[, 2]
  Mod.eff.list <- list("Group.variable" = var.id, "Group1" = group1,
                        "Group2" = group2, "Type" = type, "Mod.avg.table" = AICc.out, "Mod.avg.eff" = Mod.avg.out[,1],
                        "Uncond.se" = Mod.avg.out[,2],  "Conf.level" = conf.level, "Lower.CL" = Lower.CL,
                        "Upper.CL" = Upper.CL)
  class(Mod.eff.list) <- c("modavgEffect", "list")
  return(Mod.eff.list)
}



##lmerMod
modavgEffect.AIClmerMod <- function(cand.set, modnames = NULL, newdata, second.ord = TRUE,
                                    nobs = NULL, uncond.se = "revised", conf.level = 0.95,
                                    ...) {

  ##check if named list if modnames are not supplied
  if(is.null(modnames)) {
    if(is.null(names(cand.set))) {
      modnames <- paste("Mod", 1:length(cand.set), sep = "")
      warning("\nModel names have been supplied automatically in the table\n")
    }
    modnames <- names(cand.set)
  }
    
  ##newdata is data frame with exact structure of the original data frame (same variable names and type)
  ##check on newdata
  ##determine number of observations in new data set
  nobserv <- nrow(newdata)
  if(nobserv > 2) stop("\nCurrent maximum number of groups compared is 2:\nmodify newdata argument accordingly\n")

  ##determine number of columns in new data set
  ncolumns <- ncol(newdata)
  
  ##if only 1 column, add an additional column to avoid problems in computation
  if(ncolumns == 1) newdata$blank.fake.column.NAs <- NA
  

  ##determine which column varies
  uniques <- apply(X = newdata, MARGIN = 2, FUN = unique)
  lengths <- lapply(X = uniques, FUN = length)
  varies <- sapply(X = lengths, FUN = function(i) i > 1)

    
  ##extract name of column
  if(sum(varies) == 1) {
    var.id <- names(varies)[which(varies == TRUE)]
                
    ##determine name of groups compared
    group1 <- as.character(newdata[,paste(var.id)][1])
    group2 <- as.character(newdata[,paste(var.id)][2])

  } else {
    ##warn that no single variable defines groups
    warning("\nGroups do not seem to be defined by a single variable.\n Function proceeding with generic group names\n")
    ##use generic names
    var.id <- "Groups"
    group1 <- "group 1"
    group2 <- "group 2"
  }

    
  ##extract classes
  mod.class <- unlist(lapply(X=cand.set, FUN=class))
  ##check if all are identical
  check.class <- unique(mod.class)

  ##number of models
  nmods <- length(modnames)
  
    
  ##compute fitted values
  fit <- matrix(data = unlist(lapply(X = cand.set, FUN = function(i)predictSE(i, se.fit = TRUE, newdata = newdata)$fit)),
                nrow = nmods, ncol = 2, byrow = TRUE)
    
  ##compute SE's on fitted values
  SE <- matrix(data = unlist(lapply(X = cand.set, FUN = function(i)predictSE(i, se.fit = TRUE, newdata = newdata)$se.fit)),
               nrow = nmods, ncol = 2, byrow = TRUE)
  
  ##difference between groups 
  diff <- fit[, 1] - fit[, 2]
  
  ##SE on difference
  SE.diff <- sqrt(SE[, 1]^2 + SE[, 2]^2)
  
  ##store AICc table
  AICctab <- aictab(cand.set = cand.set, modnames = modnames, second.ord = second.ord,
                    nobs = nobs, sort = FALSE)

  ##create object to hold Model-averaged estimates and unconditional SE's
  Mod.avg.out <- matrix(NA, nrow = 1, ncol = 2)
  colnames(Mod.avg.out) <- c("Mod.avg.diff", "Uncond.SE")

  ##begin loop - AICc
  if(second.ord==TRUE){
    ##create temporary data.frame to store fitted values and SE 
    AICctmp <- AICctab
    AICctmp$diff <- diff
    AICctmp$SE.diff <- SE.diff
    
    ##compute model averaged prediction and store in output matrix
    Mod.avg.out[, 1] <- sum(AICctmp$AICcWt*AICctmp$diff)
    ##compute unconditional SE and store in output matrix
      
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Mod.avg.out[, 2] <- sum(AICctmp$AICcWt*sqrt(AICctmp$SE.diff^2 + (AICctmp$diff - Mod.avg.out[, 1])^2))
    }

    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Mod.avg.out[, 2] <- sqrt(sum(AICctmp$AICcWt*(AICctmp$SE.diff^2 + (AICctmp$diff - Mod.avg.out[, 1])^2)))
    }
    ##store table
    AICc.out <- AICctmp
  }
  
  
  ##create temporary data.frame to store fitted values and SE - AIC
  if(second.ord==FALSE) {
    AICtmp <- AICctab
    AICtmp$diff <- diff
    AICtmp$SE.diff <- SE.diff

    ##compute model averaged prediction and store in output matrix
    Mod.avg.out[, 1] <- sum(AICtmp$AICWt*AICtmp$diff)

    ##compute unconditional SE and store in output matrix
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Mod.avg.out[, 2] <- sum(AICtmp$AICWt*sqrt(AICtmp$SE.diff^2 + (AICtmp$diff - Mod.avg.out[, 1])^2))
    }

    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Mod.avg.out[, 2] <- sqrt(sum(AICtmp$AICWt*(AICtmp$SE.diff^2 + (AICtmp$diff - Mod.avg.out[, 1])^2)))
    }  
    ##store table
    AICc.out <- AICtmp
  }

  ##scale of predictions
  type <- "response"
    
  zcrit <- qnorm(p=(1-conf.level)/2, lower.tail=FALSE)
  Lower.CL <- Mod.avg.out[, 1] - zcrit * Mod.avg.out[, 2]
  Upper.CL <- Mod.avg.out[, 1] + zcrit * Mod.avg.out[, 2]
  Mod.eff.list <- list("Group.variable" = var.id, "Group1" = group1,
                        "Group2" = group2, "Type" = type, "Mod.avg.table" = AICc.out, "Mod.avg.eff" = Mod.avg.out[,1],
                        "Uncond.se" = Mod.avg.out[,2],  "Conf.level" = conf.level, "Lower.CL" = Lower.CL,
                        "Upper.CL" = Upper.CL)
  class(Mod.eff.list) <- c("modavgEffect", "list")
  return(Mod.eff.list)
}



##rlm
modavgEffect.AICrlm.lm <- function(cand.set, modnames = NULL, newdata, second.ord = TRUE,
                                   nobs = NULL, uncond.se = "revised", conf.level = 0.95,
                                   ...){
  ##newdata is data frame with exact structure of the original data frame (same variable names and type)

  ##check if named list if modnames are not supplied
  if(is.null(modnames)) {
    if(is.null(names(cand.set))) {
      modnames <- paste("Mod", 1:length(cand.set), sep = "")
      warning("\nModel names have been supplied automatically in the table\n")
    } else {
      modnames <- names(cand.set)
    }
  }

  ##check on newdata
  ##determine number of observations in new data set
  nobserv <- nrow(newdata)
  if(nobserv > 2) stop("\nCurrent maximum number of groups compared is 2:\nmodify newdata argument accordingly\n")

  ##determine number of columns in new data set
  ncolumns <- ncol(newdata)

  ##if only 1 column, add an additional column to avoid problems in computation
  if(ncolumns == 1) newdata$blank.fake.column.NAs <- NA

  
  ##determine which column varies
  uniques <- apply(X = newdata, MARGIN = 2, FUN = unique)
  lengths <- lapply(X = uniques, FUN = length)
  varies <- sapply(X = lengths, FUN = function(i) i > 1)

    
  ##extract name of column
  if(sum(varies) == 1) {
    var.id <- names(varies)[which(varies == TRUE)]
                
    ##determine name of groups compared
    group1 <- as.character(newdata[,paste(var.id)][1])
    group2 <- as.character(newdata[,paste(var.id)][2])

  } else {
    ##warn that no single variable defines groups
    warning("\nGroups do not seem to be defined by a single variable.\n Function proceeding with generic group names\n")
    ##use generic names
    var.id <- "Groups"
    group1 <- "group 1"
    group2 <- "group 2"
  }

  
  ##number of models
  nmods <- length(modnames)
    
  ##compute fitted values
  fit <- matrix(data = unlist(lapply(X = cand.set, FUN = function(i) predict(i, se.fit = TRUE, newdata = newdata)$fit)),
                nrow = nmods, ncol = 2, byrow = TRUE)
    
  ##compute SE's on fitted values
  SE <- matrix(data = unlist(lapply(X = cand.set, FUN = function(i) predict(i, se.fit = TRUE, newdata = newdata)$se.fit)),
               nrow = nmods, ncol = 2, byrow = TRUE)
    
  ##difference between groups 
  diff <- fit[, 1] - fit[, 2]
    
  ##SE on difference
  SE.diff <- sqrt(SE[, 1]^2 + SE[, 2]^2)
    
  ##store AICc table
  AICctab <- aictab(cand.set = cand.set, modnames = modnames, second.ord = second.ord, nobs = nobs, sort = FALSE)

  ##create object to hold Model-averaged estimates and unconditional SE's
  Mod.avg.out <- matrix(NA, nrow = 1, ncol = 2)
  colnames(Mod.avg.out) <- c("Mod.avg.diff", "Uncond.SE")
  
  ##begin loop - AICc
  if(second.ord == TRUE){
                   
    ##create temporary data.frame to store fitted values and SE 
    AICctmp <- AICctab
    AICctmp$diff <- diff
    AICctmp$SE.diff <- SE.diff
      
    ##compute model averaged prediction and store in output matrix
    Mod.avg.out[, 1] <- sum(AICctmp$AICcWt*AICctmp$diff)
    ##compute unconditional SE and store in output matrix
      
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Mod.avg.out[, 2] <- sum(AICctmp$AICcWt*sqrt(AICctmp$SE.diff^2 + (AICctmp$diff - Mod.avg.out[, 1])^2))
    }
    
    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Mod.avg.out[, 2] <- sqrt(sum(AICctmp$AICcWt*(AICctmp$SE.diff^2 + (AICctmp$diff - Mod.avg.out[, 1])^2)))
    }
    ##store table
    AICc.out <- AICctmp
  }
  
  ##create temporary data.frame to store fitted values and SE - AIC
  if(second.ord==FALSE) {
      
    AICtmp <- AICctab
    AICtmp$diff <- diff
    AICtmp$SE.diff <- SE.diff

    ##compute model averaged prediction and store in output matrix
    Mod.avg.out[, 1] <- sum(AICtmp$AICWt*AICtmp$diff)

    ##compute unconditional SE and store in output matrix
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Mod.avg.out[, 2] <- sum(AICtmp$AICWt*sqrt(AICtmp$SE.diff^2 + (AICtmp$diff - Mod.avg.out[, 1])^2))
    }
    
    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Mod.avg.out[, 2] <- sqrt(sum(AICtmp$AICWt*(AICtmp$SE.diff^2 + (AICtmp$diff - Mod.avg.out[, 1])^2)))
    }  
    ##store table
    AICc.out <- AICtmp
  }
  
  ##indicate scale of predictions
  type <- "response"

  zcrit <- qnorm(p=(1-conf.level)/2, lower.tail=FALSE)
  Lower.CL <- Mod.avg.out[, 1] - zcrit * Mod.avg.out[, 2]
  Upper.CL <- Mod.avg.out[, 1] + zcrit * Mod.avg.out[, 2]
  Mod.eff.list <- list("Group.variable" = var.id, "Group1" = group1,
                       "Group2" = group2, "Type" = type, "Mod.avg.table" = AICc.out, "Mod.avg.eff" = Mod.avg.out[,1],
                       "Uncond.se" = Mod.avg.out[,2], "Conf.level" = conf.level, "Lower.CL" = Lower.CL,
                       "Upper.CL" = Upper.CL)
  class(Mod.eff.list) <- c("modavgEffect", "list")
  return(Mod.eff.list)
}



##survreg
modavgEffect.AICsurvreg <- function(cand.set, modnames = NULL, newdata, second.ord = TRUE,
                                    nobs = NULL, uncond.se = "revised", conf.level = 0.95,
                                    type = "response", ...){
  ##newdata is data frame with exact structure of the original data frame (same variable names and type)

  ##check if named list if modnames are not supplied
  if(is.null(modnames)) {
    if(is.null(names(cand.set))) {
      modnames <- paste("Mod", 1:length(cand.set), sep = "")
      warning("\nModel names have been supplied automatically in the table\n")
    } else {
      modnames <- names(cand.set)
    }
  }

  ##check that distribution is the same for all models
  if(identical(type, "link")) {
    check.dist <- sapply(X = cand.set, FUN = function(i) i$dist)
    unique.dist <- unique(x = check.dist)
    if(length(unique.dist) > 1) stop("\nFunction does not support model-averaging effect size on link scale using different distributions\n")
  }
  
  
    ##determine number of observations in new data set
    nobserv <- nrow(newdata)

    if(nobserv > 2) stop("\nCurrent maximum number of groups compared is 2:\nmodify newdata argument accordingly\n")

    ##determine number of columns in new data set
    ncolumns <- ncol(newdata)
        
    ##if only 1 column, add an additional column to avoid problems in computation
    if(ncolumns == 1) newdata$blank.fake.column.NAs <- NA

  
  ##determine which column varies
  uniques <- apply(X = newdata, MARGIN = 2, FUN = unique)
  lengths <- lapply(X = uniques, FUN = length)
  varies <- sapply(X = lengths, FUN = function(i) i > 1)

    
  ##extract name of column
  if(sum(varies) == 1) {
    var.id <- names(varies)[which(varies == TRUE)]
                
    ##determine name of groups compared
    group1 <- as.character(newdata[,paste(var.id)][1])
    group2 <- as.character(newdata[,paste(var.id)][2])

  } else {
    ##warn that no single variable defines groups
    warning("\nGroups do not seem to be defined by a single variable.\n Function proceeding with generic group names\n")
    ##use generic names
    var.id <- "Groups"
    group1 <- "group 1"
    group2 <- "group 2"
  }
    
    
    ##number of models
    nmods <- length(modnames)

    ##compute fitted values
    fit <- matrix(data = unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata, type = type)$fit)),
                  nrow = nmods, ncol = 2, byrow = TRUE)
    
    ##compute SE's on fitted values
    SE <- matrix(data = unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata, type = type)$se.fit)),
                 nrow = nmods, ncol = 2, byrow = TRUE)
    
    ##difference between groups 
    diff <- fit[, 1] - fit[, 2]
    
    ##SE on difference
    SE.diff <- sqrt(SE[, 1]^2 + SE[, 2]^2)

    ##store AICc table
    AICctab <- aictab(cand.set = cand.set, modnames = modnames, second.ord = second.ord, nobs = nobs, sort = FALSE)

    #create object to hold Model-averaged estimates and unconditional SE's
    Mod.avg.out <- matrix(NA, nrow = 1, ncol = 2)
    colnames(Mod.avg.out) <- c("Mod.avg.diff", "Uncond.SE")

    ##begin loop - AICc
    if(second.ord == TRUE){
                   
      ##create temporary data.frame to store fitted values and SE 
      AICctmp <- AICctab
      AICctmp$diff <- diff
      AICctmp$SE.diff <- SE.diff
      
      ##compute model averaged prediction and store in output matrix
      Mod.avg.out[, 1] <- sum(AICctmp$AICcWt*AICctmp$diff)
      ##compute unconditional SE and store in output matrix
      
      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Mod.avg.out[, 2] <- sum(AICctmp$AICcWt*sqrt(AICctmp$SE.diff^2 + (AICctmp$diff - Mod.avg.out[, 1])^2))
      }

      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Mod.avg.out[, 2] <- sqrt(sum(AICctmp$AICcWt*(AICctmp$SE.diff^2 + (AICctmp$diff - Mod.avg.out[, 1])^2)))
      }
      ##store table
      AICc.out <- AICctmp
    }
  


    ##create temporary data.frame to store fitted values and SE - AIC
    if(second.ord==FALSE) {
      
      AICtmp <- AICctab
      AICtmp$diff <- diff
      AICtmp$SE.diff <- SE.diff

      ##compute model averaged prediction and store in output matrix
      Mod.avg.out[, 1] <- sum(AICtmp$AICWt*AICtmp$diff)

      ##compute unconditional SE and store in output matrix
      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Mod.avg.out[, 2] <- sum(AICtmp$AICWt*sqrt(AICtmp$SE.diff^2 + (AICtmp$diff - Mod.avg.out[, 1])^2))
      }

      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Mod.avg.out[, 2] <- sqrt(sum(AICtmp$AICWt*(AICtmp$SE.diff^2 + (AICtmp$diff - Mod.avg.out[, 1])^2)))
      }  
      ##store table
      AICc.out <- AICtmp
    }

    
    zcrit <- qnorm(p=(1-conf.level)/2, lower.tail=FALSE)
    Lower.CL <- Mod.avg.out[, 1] - zcrit * Mod.avg.out[, 2]
    Upper.CL <- Mod.avg.out[, 1] + zcrit * Mod.avg.out[, 2]
    Mod.eff.list <- list("Group.variable" = var.id, "Group1" = group1,
                         "Group2" = group2, "Type" = type, "Mod.avg.table" = AICc.out, "Mod.avg.eff" = Mod.avg.out[,1], 
                          "Uncond.se" = Mod.avg.out[,2], "Conf.level" = conf.level, "Lower.CL" = Lower.CL,
                          "Upper.CL" = Upper.CL)
    class(Mod.eff.list) <- c("modavgEffect", "list")
    return(Mod.eff.list)
}



##occu
modavgEffect.AICunmarkedFitOccu <- function(cand.set, modnames = NULL, newdata, second.ord = TRUE,
                                            nobs = NULL, uncond.se = "revised", conf.level = 0.95,
                                            type = "response", c.hat = 1, parm.type = NULL,
                                            ...) {

  ##check if named list if modnames are not supplied
  if(is.null(modnames)) {
    if(is.null(names(cand.set))) {
      modnames <- paste("Mod", 1:length(cand.set), sep = "")
      warning("\nModel names have been supplied automatically in the table\n")
    } else {
      modnames <- names(cand.set)
    }
  }

  ##rename values according to unmarked to extract from object
  ##psi
  if(identical(parm.type, "psi")) {
    parm.type1 <- "state"; parm.id <- "psi"
  }

  ##detect
  if(identical(parm.type, "detect")) {parm.type1 <- "det"; parm.id <- "p"}
     
  ##newdata is data frame with exact structure of the original data frame (same variable names and type)
  ##check on newdata
  ##determine number of observations in new data set
  nobserv <- nrow(newdata)
  if(nobserv > 2) stop("\nCurrent maximum number of groups compared is 2:\nmodify newdata argument accordingly\n")

  ##determine number of columns in new data set
  ncolumns <- ncol(newdata)

  ##if only 1 column, add an additional column to avoid problems in computation
  if(ncolumns == 1) newdata$blank.fake.column.NAs <- NA

  
  ##determine which column varies
  uniques <- apply(X = newdata, MARGIN = 2, FUN = unique)
  lengths <- lapply(X = uniques, FUN = length)
  varies <- sapply(X = lengths, FUN = function(i) i > 1)

    
  ##extract name of column
  if(sum(varies) == 1) {
    var.id <- names(varies)[which(varies == TRUE)]
                
    ##determine name of groups compared
    group1 <- as.character(newdata[,paste(var.id)][1])
    group2 <- as.character(newdata[,paste(var.id)][2])

  } else {
    ##warn that no single variable defines groups
    warning("\nGroups do not seem to be defined by a single variable.\n Function proceeding with generic group names\n")
    ##use generic names
    var.id <- "Groups"
    group1 <- "group 1"
    group2 <- "group 2"
  }

  
  ##number of models
  nmods <- length(modnames)


  ##compute predicted values
  ##point estimate
  if(identical(type, "response")) {
    ##extract fitted value for observation obs
    fit <- matrix(data = unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata,
                                                       type = parm.type1)$Predicted)),
                  nrow = nmods, ncol = 2, byrow = TRUE)

    ##extract SE for fitted value for observation obs
    SE <- matrix(data = unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata,
                                                      type = parm.type1)$SE)),
                 nrow = nmods, ncol = 2, byrow = TRUE)
  }

  
  ##link scale
  if(identical(type, "link")) {
    ##extract fitted value for observation obs
    fit <- matrix(data = unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata,
                                                       type = parm.type1, backTransform = FALSE)$Predicted)),
                  nrow = nmods, ncol = 2, byrow = TRUE)
    
    ##extract SE for fitted value for observation obs
    SE <- matrix(data = unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata,
                                                      type = parm.type1, backTransform = FALSE)$SE)),
                 nrow = nmods, ncol = 2, byrow = TRUE)
  }  

  
  ##difference between groups 
  diff <- fit[, 1] - fit[, 2]
    
  ##SE on difference
  SE.diff <- sqrt(SE[, 1]^2 + SE[, 2]^2)
  
  ##store AICc table
  AICctab <- aictab(cand.set = cand.set, modnames = modnames, second.ord = second.ord,
                    nobs = nobs, sort = FALSE, c.hat = c.hat)

  ##create object to hold Model-averaged estimates and unconditional SE's
  Mod.avg.out <- matrix(NA, nrow = 1, ncol = 2)
  colnames(Mod.avg.out) <- c("Mod.avg.diff", "Uncond.SE")

  ##begin loop - AICc
  if(second.ord == TRUE && c.hat == 1){
     
    ##create temporary data.frame to store fitted values and SE 
    AICctmp <- AICctab
    AICctmp$diff <- diff
    AICctmp$SE.diff <- SE.diff
    
    ##compute model averaged prediction and store in output matrix
    Mod.avg.out[, 1] <- sum(AICctmp$AICcWt*AICctmp$diff)
    ##compute unconditional SE and store in output matrix
      
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Mod.avg.out[, 2] <- sum(AICctmp$AICcWt*sqrt(AICctmp$SE.diff^2 + (AICctmp$diff - Mod.avg.out[, 1])^2))
    }

    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Mod.avg.out[, 2] <- sqrt(sum(AICctmp$AICcWt*(AICctmp$SE.diff^2 + (AICctmp$diff - Mod.avg.out[, 1])^2)))
    }
    ##store table
    AICc.out <- AICctmp
  }


  ##begin loop - QAICc
  if(second.ord == TRUE && c.hat > 1){
      ##create temporary data.frame to store fitted values and SE 
      QAICctmp <- AICctab
      QAICctmp$diff <- diff
      QAICctmp$SE.diff <- SE.diff * sqrt(c.hat)

      ##compute model averaged prediction and store in output matrix
      Mod.avg.out[, 1] <- sum(QAICctmp$QAICcWt*QAICctmp$diff)
      
      ##compute unconditional SE and store in output matrix
      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Mod.avg.out[, 2] <- sum(QAICctmp$QAICcWt*sqrt(QAICctmp$SE.diff^2 + (QAICctmp$diff - Mod.avg.out[, 1])^2))
      }
      
      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Mod.avg.out[, 2] <- sqrt(sum(QAICctmp$QAICcWt*(QAICctmp$SE.diff^2 + (QAICctmp$diff - Mod.avg.out[, 1])^2)))  
      }
      ##store table
      AICc.out <- QAICctmp
    
    }
  
  
  
  ##create temporary data.frame to store fitted values and SE - AIC
  if(second.ord == FALSE && c.hat == 1) {
    AICtmp <- AICctab
    AICtmp$diff <- diff
    AICtmp$SE.diff <- SE.diff

    ##compute model averaged prediction and store in output matrix
    Mod.avg.out[, 1] <- sum(AICtmp$AICWt*AICtmp$diff)

    ##compute unconditional SE and store in output matrix
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Mod.avg.out[, 2] <- sum(AICtmp$AICWt*sqrt(AICtmp$SE.diff^2 + (AICtmp$diff - Mod.avg.out[, 1])^2))
    }

    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Mod.avg.out[, 2] <- sqrt(sum(AICtmp$AICWt*(AICtmp$SE.diff^2 + (AICtmp$diff - Mod.avg.out[, 1])^2)))
    }  
    ##store table
    AICc.out <- AICtmp
  }

  
  ##begin loop - QAICc
  if(second.ord == FALSE && c.hat > 1){
      
    ##create temporary data.frame to store fitted values and SE 
    QAICtmp <- AICctab
    QAICtmp$diff <- diff
    QAICtmp$SE.diff <- SE.diff* sqrt(c.hat)
      
    ##compute model averaged prediction and store in output matrix
    Mod.avg.out[, 1] <- sum(QAICtmp$QAICWt*QAICtmp$diff)
      
    ##compute unconditional SE and store in output matrix
    if(identical(uncond.se, "old")) {
      Mod.avg.out[, 2] <- sum(QAICtmp$QAICWt*sqrt(QAICtmp$SE.diff^2 + (QAICtmp$diff - Mod.avg.out[, 1])^2))
    }
      
    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Mod.avg.out[, 2] <- sqrt(sum(QAICtmp$QAICWt*(QAICtmp$SE.diff^2 + (QAICtmp$diff - Mod.avg.out[, 1])^2)))
    }  
    ##store table
    AICc.out <- QAICtmp
  }

  Group.variable <- paste(parm.id, "(", var.id, ")")
  
  zcrit <- qnorm(p=(1-conf.level)/2, lower.tail=FALSE)
  Lower.CL <- Mod.avg.out[, 1] - zcrit * Mod.avg.out[, 2]
  Upper.CL <- Mod.avg.out[, 1] + zcrit * Mod.avg.out[, 2]
  Mod.eff.list <- list("Group.variable" = Group.variable, "Group1" = group1,
                       "Group2" = group2, "Type" = type, "Mod.avg.table" = AICc.out,
                       "Mod.avg.eff" = Mod.avg.out[,1], "Uncond.se" = Mod.avg.out[,2], "Conf.level" = conf.level,
                       "Lower.CL" = Lower.CL, "Upper.CL" = Upper.CL)
  class(Mod.eff.list) <- c("modavgEffect", "list")
  return(Mod.eff.list)  
}



##colext
modavgEffect.AICunmarkedFitColExt <- function(cand.set, modnames = NULL, newdata, second.ord = TRUE,
                                              nobs = NULL, uncond.se = "revised", conf.level = 0.95,
                                              type = "response", c.hat = 1, parm.type = NULL,
                                              ...) {

  ##check if named list if modnames are not supplied
  if(is.null(modnames)) {
    if(is.null(names(cand.set))) {
      modnames <- paste("Mod", 1:length(cand.set), sep = "")
      warning("\nModel names have been supplied automatically in the table\n")
    } else {
      modnames <- names(cand.set)
    }
  }

  ##rename values according to unmarked to extract from object
  ##psi
  if(identical(parm.type, "psi")) {
    parm.type1 <- "psi"; parm.id <- "psi"
  }

  ##gamma
  if(identical(parm.type, "gamma")) {
    parm.type1 <- "col"; parm.id <- "col"
  }

  ##epsilon
  if(identical(parm.type, "epsilon")) {
    parm.type1 <- "ext"; parm.id <- "ext"
  }

  ##detect
  if(identical(parm.type, "detect")) {parm.type1 <- "det"; parm.id <- "p"}
  
  ##newdata is data frame with exact structure of the original data frame (same variable names and type)
  ##check on newdata
  ##determine number of observations in new data set
  nobserv <- nrow(newdata)
  if(nobserv > 2) stop("\nCurrent maximum number of groups compared is 2:\nmodify newdata argument accordingly\n")

  ##determine number of columns in new data set
  ncolumns <- ncol(newdata)

  ##if only 1 column, add an additional column to avoid problems in computation
  if(ncolumns == 1) newdata$blank.fake.column.NAs <- NA

  
  ##determine which column varies
  uniques <- apply(X = newdata, MARGIN = 2, FUN = unique)
  lengths <- lapply(X = uniques, FUN = length)
  varies <- sapply(X = lengths, FUN = function(i) i > 1)

    
  ##extract name of column
  if(sum(varies) == 1) {
    var.id <- names(varies)[which(varies == TRUE)]
                
    ##determine name of groups compared
    group1 <- as.character(newdata[,paste(var.id)][1])
    group2 <- as.character(newdata[,paste(var.id)][2])

  } else {
    ##warn that no single variable defines groups
    warning("\nGroups do not seem to be defined by a single variable.\n Function proceeding with generic group names\n")
    ##use generic names
    var.id <- "Groups"
    group1 <- "group 1"
    group2 <- "group 2"
  }

  
  ##number of models
  nmods <- length(modnames)


  ##compute predicted values
  ##point estimate
  if(identical(type, "response")) {
    ##extract fitted value for observation obs
    fit <- matrix(data = unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata,
                                                       type = parm.type1)$Predicted)),
                  nrow = nmods, ncol = 2, byrow = TRUE)

    ##extract SE for fitted value for observation obs
    SE <- matrix(data = unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata,
                                                      type = parm.type1)$SE)),
                 nrow = nmods, ncol = 2, byrow = TRUE)
  }
  
  ##link scale
  if(identical(type, "link")) {
    ##extract fitted value for observation obs
    fit <- matrix(data = unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata,
                                                       type = parm.type1, backTransform = FALSE)$Predicted)),
                  nrow = nmods, ncol = 2, byrow = TRUE)

    ##extract SE for fitted value for observation obs
    SE <- matrix(data = unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata,
                                                      type = parm.type1, backTransform = FALSE)$SE)),
                 nrow = nmods, ncol = 2, byrow = TRUE)
  }  

  
  ##difference between groups 
  diff <- fit[, 1] - fit[, 2]
    
  ##SE on difference
  SE.diff <- sqrt(SE[, 1]^2 + SE[, 2]^2)
  
  ##store AICc table
  AICctab <- aictab(cand.set = cand.set, modnames = modnames, second.ord = second.ord,
                    nobs = nobs, sort = FALSE, c.hat = c.hat)

  ##create object to hold Model-averaged estimates and unconditional SE's
  Mod.avg.out <- matrix(NA, nrow = 1, ncol = 2)
  colnames(Mod.avg.out) <- c("Mod.avg.diff", "Uncond.SE")

  ##begin loop - AICc
  if(second.ord == TRUE && c.hat == 1){
     
    ##create temporary data.frame to store fitted values and SE 
    AICctmp <- AICctab
    AICctmp$diff <- diff
    AICctmp$SE.diff <- SE.diff
    
    ##compute model averaged prediction and store in output matrix
    Mod.avg.out[, 1] <- sum(AICctmp$AICcWt*AICctmp$diff)
    ##compute unconditional SE and store in output matrix
      
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Mod.avg.out[, 2] <- sum(AICctmp$AICcWt*sqrt(AICctmp$SE.diff^2 + (AICctmp$diff - Mod.avg.out[, 1])^2))
    }

    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Mod.avg.out[, 2] <- sqrt(sum(AICctmp$AICcWt*(AICctmp$SE.diff^2 + (AICctmp$diff - Mod.avg.out[, 1])^2)))
    }
    ##store table
    AICc.out <- AICctmp
  }


  ##begin loop - QAICc
  if(second.ord == TRUE && c.hat > 1){
      ##create temporary data.frame to store fitted values and SE 
      QAICctmp <- AICctab
      QAICctmp$diff <- diff
      QAICctmp$SE.diff <- SE.diff * sqrt(c.hat)

      ##compute model averaged prediction and store in output matrix
      Mod.avg.out[, 1] <- sum(QAICctmp$QAICcWt*QAICctmp$diff)
      
      ##compute unconditional SE and store in output matrix
      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Mod.avg.out[, 2] <- sum(QAICctmp$QAICcWt*sqrt(QAICctmp$SE.diff^2 + (QAICctmp$diff - Mod.avg.out[, 1])^2))
      }
      
      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Mod.avg.out[, 2] <- sqrt(sum(QAICctmp$QAICcWt*(QAICctmp$SE.diff^2 + (QAICctmp$diff - Mod.avg.out[, 1])^2)))  
      }
      ##store table
      AICc.out <- QAICctmp
    
    }
  
  
  
  ##create temporary data.frame to store fitted values and SE - AIC
  if(second.ord == FALSE && c.hat == 1) {
    AICtmp <- AICctab
    AICtmp$diff <- diff
    AICtmp$SE.diff <- SE.diff

    ##compute model averaged prediction and store in output matrix
    Mod.avg.out[, 1] <- sum(AICtmp$AICWt*AICtmp$diff)

    ##compute unconditional SE and store in output matrix
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Mod.avg.out[, 2] <- sum(AICtmp$AICWt*sqrt(AICtmp$SE.diff^2 + (AICtmp$diff - Mod.avg.out[, 1])^2))
    }

    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Mod.avg.out[, 2] <- sqrt(sum(AICtmp$AICWt*(AICtmp$SE.diff^2 + (AICtmp$diff - Mod.avg.out[, 1])^2)))
    }  
    ##store table
    AICc.out <- AICtmp
  }

  
  ##begin loop - QAICc
  if(second.ord == FALSE && c.hat > 1){
      
    ##create temporary data.frame to store fitted values and SE 
    QAICtmp <- AICctab
    QAICtmp$diff <- diff
    QAICtmp$SE.diff <- SE.diff* sqrt(c.hat)
      
    ##compute model averaged prediction and store in output matrix
    Mod.avg.out[, 1] <- sum(QAICtmp$QAICWt*QAICtmp$diff)
      
    ##compute unconditional SE and store in output matrix
    if(identical(uncond.se, "old")) {
      Mod.avg.out[, 2] <- sum(QAICtmp$QAICWt*sqrt(QAICtmp$SE.diff^2 + (QAICtmp$diff - Mod.avg.out[, 1])^2))
    }
      
    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Mod.avg.out[, 2] <- sqrt(sum(QAICtmp$QAICWt*(QAICtmp$SE.diff^2 + (QAICtmp$diff - Mod.avg.out[, 1])^2)))
    }  
    ##store table
    AICc.out <- QAICtmp
  }

  Group.variable <- paste(parm.id, "(", var.id, ")")
  
  zcrit <- qnorm(p=(1-conf.level)/2, lower.tail=FALSE)
  Lower.CL <- Mod.avg.out[, 1] - zcrit * Mod.avg.out[, 2]
  Upper.CL <- Mod.avg.out[, 1] + zcrit * Mod.avg.out[, 2]
  Mod.eff.list <- list("Group.variable" = Group.variable, "Group1" = group1,
                       "Group2" = group2, "Type" = type, "Mod.avg.table" = AICc.out,
                       "Mod.avg.eff" = Mod.avg.out[,1], "Uncond.se" = Mod.avg.out[,2], "Conf.level" = conf.level,
                       "Lower.CL" = Lower.CL, "Upper.CL" = Upper.CL)
  class(Mod.eff.list) <- c("modavgEffect", "list")
  return(Mod.eff.list)  
}



##occuRN
modavgEffect.AICunmarkedFitOccuRN <- function(cand.set, modnames = NULL, newdata, second.ord = TRUE,
                                              nobs = NULL, uncond.se = "revised", conf.level = 0.95,
                                              type = "response", c.hat = 1, parm.type = NULL,
                                              ...) {

  ##check if named list if modnames are not supplied
  if(is.null(modnames)) {
    if(is.null(names(cand.set))) {
      modnames <- paste("Mod", 1:length(cand.set), sep = "")
      warning("\nModel names have been supplied automatically in the table\n")
    } else {
      modnames <- names(cand.set)
    }
  }

  ##rename values according to unmarked to extract from object
  ##lambda
  if(identical(parm.type, "lambda")) {
    parm.type1 <- "state"; parm.id <- "lam"
  }

  ##detect
  if(identical(parm.type, "detect")) {parm.type1 <- "det"; parm.id <- "p"}
     
  ##newdata is data frame with exact structure of the original data frame (same variable names and type)
  ##check on newdata
  ##determine number of observations in new data set
  nobserv <- nrow(newdata)
  if(nobserv > 2) stop("\nCurrent maximum number of groups compared is 2:\nmodify newdata argument accordingly\n")

  ##determine number of columns in new data set
  ncolumns <- ncol(newdata)

  ##if only 1 column, add an additional column to avoid problems in computation
  if(ncolumns == 1) newdata$blank.fake.column.NAs <- NA
  

  ##determine which column varies
  uniques <- apply(X = newdata, MARGIN = 2, FUN = unique)
  lengths <- lapply(X = uniques, FUN = length)
  varies <- sapply(X = lengths, FUN = function(i) i > 1)

    
  ##extract name of column
  if(sum(varies) == 1) {
    var.id <- names(varies)[which(varies == TRUE)]
                
    ##determine name of groups compared
    group1 <- as.character(newdata[,paste(var.id)][1])
    group2 <- as.character(newdata[,paste(var.id)][2])

  } else {
    ##warn that no single variable defines groups
    warning("\nGroups do not seem to be defined by a single variable.\n Function proceeding with generic group names\n")
    ##use generic names
    var.id <- "Groups"
    group1 <- "group 1"
    group2 <- "group 2"
  }

  
  ##number of models
  nmods <- length(modnames)


  ##compute predicted values
  ##point estimate
  if(identical(type, "response")) {
    fit <- matrix(data = unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata,
                                                       type = parm.type1)$Predicted)),
                  nrow = nmods, ncol = 2, byrow = TRUE)

    ##extract SE for fitted value for observation obs
    SE <- matrix(data = unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata,
                                                      type = parm.type1)$SE)),
                 nrow = nmods, ncol = 2, byrow = TRUE)
  }
  
  ##link scale
  if(identical(type, "link")) {
    ##extract fitted value for observation obs
    fit <- matrix(data = unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata,
                                                       type = parm.type1, backTransform = FALSE)$Predicted)),
                  nrow = nmods, ncol = 2, byrow = TRUE)

    ##extract SE for fitted value for observation obs
    SE <- matrix(data = unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata,
                                                      type = parm.type1, backTransform = FALSE)$SE)),
                 nrow = nmods, ncol = 2, byrow = TRUE)
  }  

  
  ##difference between groups 
  diff <- fit[, 1] - fit[, 2]
    
  ##SE on difference
  SE.diff <- sqrt(SE[, 1]^2 + SE[, 2]^2)
  
  ##store AICc table
  AICctab <- aictab(cand.set = cand.set, modnames = modnames, second.ord = second.ord,
                    nobs = nobs, sort = FALSE, c.hat = c.hat)

  ##create object to hold Model-averaged estimates and unconditional SE's
  Mod.avg.out <- matrix(NA, nrow = 1, ncol = 2)
  colnames(Mod.avg.out) <- c("Mod.avg.diff", "Uncond.SE")

  ##begin loop - AICc
  if(second.ord == TRUE && c.hat == 1){
     
    ##create temporary data.frame to store fitted values and SE 
    AICctmp <- AICctab
    AICctmp$diff <- diff
    AICctmp$SE.diff <- SE.diff
    
    ##compute model averaged prediction and store in output matrix
    Mod.avg.out[, 1] <- sum(AICctmp$AICcWt*AICctmp$diff)
    ##compute unconditional SE and store in output matrix
      
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Mod.avg.out[, 2] <- sum(AICctmp$AICcWt*sqrt(AICctmp$SE.diff^2 + (AICctmp$diff - Mod.avg.out[, 1])^2))
    }

    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Mod.avg.out[, 2] <- sqrt(sum(AICctmp$AICcWt*(AICctmp$SE.diff^2 + (AICctmp$diff - Mod.avg.out[, 1])^2)))
    }
    ##store table
    AICc.out <- AICctmp
  }


  ##begin loop - QAICc
  if(second.ord == TRUE && c.hat > 1){
      ##create temporary data.frame to store fitted values and SE 
      QAICctmp <- AICctab
      QAICctmp$diff <- diff
      QAICctmp$SE.diff <- SE.diff * sqrt(c.hat)

      ##compute model averaged prediction and store in output matrix
      Mod.avg.out[, 1] <- sum(QAICctmp$QAICcWt*QAICctmp$diff)
      
      ##compute unconditional SE and store in output matrix
      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Mod.avg.out[, 2] <- sum(QAICctmp$QAICcWt*sqrt(QAICctmp$SE.diff^2 + (QAICctmp$diff - Mod.avg.out[, 1])^2))
      }
      
      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Mod.avg.out[, 2] <- sqrt(sum(QAICctmp$QAICcWt*(QAICctmp$SE.diff^2 + (QAICctmp$diff - Mod.avg.out[, 1])^2)))  
      }
      ##store table
      AICc.out <- QAICctmp
    
    }
  
  
  
  ##create temporary data.frame to store fitted values and SE - AIC
  if(second.ord == FALSE && c.hat == 1) {
    AICtmp <- AICctab
    AICtmp$diff <- diff
    AICtmp$SE.diff <- SE.diff

    ##compute model averaged prediction and store in output matrix
    Mod.avg.out[, 1] <- sum(AICtmp$AICWt*AICtmp$diff)

    ##compute unconditional SE and store in output matrix
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Mod.avg.out[, 2] <- sum(AICtmp$AICWt*sqrt(AICtmp$SE.diff^2 + (AICtmp$diff - Mod.avg.out[, 1])^2))
    }

    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Mod.avg.out[, 2] <- sqrt(sum(AICtmp$AICWt*(AICtmp$SE.diff^2 + (AICtmp$diff - Mod.avg.out[, 1])^2)))
    }  
    ##store table
    AICc.out <- AICtmp
  }

  
  ##begin loop - QAICc
  if(second.ord == FALSE && c.hat > 1){
      
    ##create temporary data.frame to store fitted values and SE 
    QAICtmp <- AICctab
    QAICtmp$diff <- diff
    QAICtmp$SE.diff <- SE.diff* sqrt(c.hat)
      
    ##compute model averaged prediction and store in output matrix
    Mod.avg.out[, 1] <- sum(QAICtmp$QAICWt*QAICtmp$diff)
      
    ##compute unconditional SE and store in output matrix
    if(identical(uncond.se, "old")) {
      Mod.avg.out[, 2] <- sum(QAICtmp$QAICWt*sqrt(QAICtmp$SE.diff^2 + (QAICtmp$diff - Mod.avg.out[, 1])^2))
    }
      
    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Mod.avg.out[, 2] <- sqrt(sum(QAICtmp$QAICWt*(QAICtmp$SE.diff^2 + (QAICtmp$diff - Mod.avg.out[, 1])^2)))
    }  
    ##store table
    AICc.out <- QAICtmp
  }

  Group.variable <- paste(parm.id, "(", var.id, ")")
  
  zcrit <- qnorm(p=(1-conf.level)/2, lower.tail=FALSE)
  Lower.CL <- Mod.avg.out[, 1] - zcrit * Mod.avg.out[, 2]
  Upper.CL <- Mod.avg.out[, 1] + zcrit * Mod.avg.out[, 2]
  Mod.eff.list <- list("Group.variable" = Group.variable, "Group1" = group1,
                       "Group2" = group2, "Type" = type, "Mod.avg.table" = AICc.out,
                       "Mod.avg.eff" = Mod.avg.out[,1], "Uncond.se" = Mod.avg.out[,2], "Conf.level" = conf.level,
                       "Lower.CL" = Lower.CL, "Upper.CL" = Upper.CL)
  class(Mod.eff.list) <- c("modavgEffect", "list")
  return(Mod.eff.list)  
}



##pcount
modavgEffect.AICunmarkedFitPCount <- function(cand.set, modnames = NULL, newdata, second.ord = TRUE,
                                              nobs = NULL, uncond.se = "revised", conf.level = 0.95,
                                              type = "response", c.hat = 1, parm.type = NULL,
                                              ...) {

  ##check if named list if modnames are not supplied
  if(is.null(modnames)) {
    if(is.null(names(cand.set))) {
      modnames <- paste("Mod", 1:length(cand.set), sep = "")
      warning("\nModel names have been supplied automatically in the table\n")
    } else {
      modnames <- names(cand.set)
    }
  }

  ##rename values according to unmarked to extract from object
  ##lambda
  if(identical(parm.type, "lambda")) {
    parm.type1 <- "state"; parm.id <- "lam"
    
    ##check mixture type for mixture models
    mixture.type <- sapply(X = cand.set, FUN = function(i) i@mixture)
    unique.mixture <- unique(mixture.type)
    if(length(unique.mixture) > 1) {
      if(any(unique.mixture == "ZIP")) stop("\nThis function does not yet support mixing ZIP with other distributions\n")
    } else {
      mixture.id <- unique(mixture.type)
      if(identical(unique.mixture, "ZIP")) {
        if(identical(type, "link")) stop("\nLink scale not yet supported for ZIP mixtures\n")
      }
    }
  }
   

  ##detect
  if(identical(parm.type, "detect")) {parm.type1 <- "det"; parm.id <- "p"}
     
  ##newdata is data frame with exact structure of the original data frame (same variable names and type)
  ##check on newdata
  ##determine number of observations in new data set
  nobserv <- nrow(newdata)
  if(nobserv > 2) stop("\nCurrent maximum number of groups compared is 2:\nmodify newdata argument accordingly\n")

  ##determine number of columns in new data set
  ncolumns <- ncol(newdata)

  ##if only 1 column, add an additional column to avoid problems in computation
  if(ncolumns == 1) newdata$blank.fake.column.NAs <- NA
  

  ##determine which column varies
  uniques <- apply(X = newdata, MARGIN = 2, FUN = unique)
  lengths <- lapply(X = uniques, FUN = length)
  varies <- sapply(X = lengths, FUN = function(i) i > 1)

    
  ##extract name of column
  if(sum(varies) == 1) {
    var.id <- names(varies)[which(varies == TRUE)]
                
    ##determine name of groups compared
    group1 <- as.character(newdata[,paste(var.id)][1])
    group2 <- as.character(newdata[,paste(var.id)][2])

  } else {
    ##warn that no single variable defines groups
    warning("\nGroups do not seem to be defined by a single variable.\n Function proceeding with generic group names\n")
    ##use generic names
    var.id <- "Groups"
    group1 <- "group 1"
    group2 <- "group 2"
  }

  
  
  ##number of models
  nmods <- length(modnames)


  ##compute predicted values
  ##point estimate
  if(identical(type, "response")) {
    ##extract fitted value for observation obs
    if(identical(parm.type, "lambda") && identical(mixture.id, "ZIP")) {
      fit <- matrix(data = unlist(lapply(X = cand.set, FUN = function(i)predictSE(i, se.fit = TRUE,
                                           newdata = newdata)$fit)),
                    nrow = nmods, ncol = 2, byrow = TRUE)
      
      SE <- matrix(data = unlist(lapply(X = cand.set, FUN = function(i)predictSE(i, se.fit = TRUE,
                                          newdata = newdata)$se.fit)),
                   nrow = nmods, ncol = 2, byrow = TRUE)
      
    } else {
      fit <- matrix(data = unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata,
                                                         type = parm.type1)$Predicted)),
                    nrow = nmods, ncol = 2, byrow = TRUE)

    ##extract SE for fitted value for observation obs
    SE <- matrix(data = unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata,
                                        type = parm.type1)$SE)),
                 nrow = nmods, ncol = 2, byrow = TRUE)
    }
  }
  
  ##link scale
  if(identical(type, "link")) {
    ##extract fitted value for observation obs
    fit <- matrix(data = unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata,
                                                       type = parm.type1, backTransform = FALSE)$Predicted)),
                  nrow = nmods, ncol = 2, byrow = TRUE)

    ##extract SE for fitted value for observation obs
    SE <- matrix(data = unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata,
                                                      type = parm.type1, backTransform = FALSE)$SE)),
                 nrow = nmods, ncol = 2, byrow = TRUE)
  }  

  
  ##difference between groups 
  diff <- fit[, 1] - fit[, 2]
    
  ##SE on difference
  SE.diff <- sqrt(SE[, 1]^2 + SE[, 2]^2)
  
  ##store AICc table
  AICctab <- aictab(cand.set = cand.set, modnames = modnames, second.ord = second.ord,
                    nobs = nobs, sort = FALSE, c.hat = c.hat)

  ##create object to hold Model-averaged estimates and unconditional SE's
  Mod.avg.out <- matrix(NA, nrow = 1, ncol = 2)
  colnames(Mod.avg.out) <- c("Mod.avg.diff", "Uncond.SE")

  ##begin loop - AICc
  if(second.ord == TRUE && c.hat == 1){
     
    ##create temporary data.frame to store fitted values and SE 
    AICctmp <- AICctab
    AICctmp$diff <- diff
    AICctmp$SE.diff <- SE.diff
    
    ##compute model averaged prediction and store in output matrix
    Mod.avg.out[, 1] <- sum(AICctmp$AICcWt*AICctmp$diff)
    ##compute unconditional SE and store in output matrix
      
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Mod.avg.out[, 2] <- sum(AICctmp$AICcWt*sqrt(AICctmp$SE.diff^2 + (AICctmp$diff - Mod.avg.out[, 1])^2))
    }

    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Mod.avg.out[, 2] <- sqrt(sum(AICctmp$AICcWt*(AICctmp$SE.diff^2 + (AICctmp$diff - Mod.avg.out[, 1])^2)))
    }
    ##store table
    AICc.out <- AICctmp
  }


  ##begin loop - QAICc
  if(second.ord == TRUE && c.hat > 1){
      ##create temporary data.frame to store fitted values and SE 
      QAICctmp <- AICctab
      QAICctmp$diff <- diff
      QAICctmp$SE.diff <- SE.diff * sqrt(c.hat)

      ##compute model averaged prediction and store in output matrix
      Mod.avg.out[, 1] <- sum(QAICctmp$QAICcWt*QAICctmp$diff)
      
      ##compute unconditional SE and store in output matrix
      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Mod.avg.out[, 2] <- sum(QAICctmp$QAICcWt*sqrt(QAICctmp$SE.diff^2 + (QAICctmp$diff - Mod.avg.out[, 1])^2))
      }
      
      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Mod.avg.out[, 2] <- sqrt(sum(QAICctmp$QAICcWt*(QAICctmp$SE.diff^2 + (QAICctmp$diff - Mod.avg.out[, 1])^2)))  
      }
      ##store table
      AICc.out <- QAICctmp
    
    }
  
  
  
  ##create temporary data.frame to store fitted values and SE - AIC
  if(second.ord == FALSE && c.hat == 1) {
    AICtmp <- AICctab
    AICtmp$diff <- diff
    AICtmp$SE.diff <- SE.diff

    ##compute model averaged prediction and store in output matrix
    Mod.avg.out[, 1] <- sum(AICtmp$AICWt*AICtmp$diff)

    ##compute unconditional SE and store in output matrix
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Mod.avg.out[, 2] <- sum(AICtmp$AICWt*sqrt(AICtmp$SE.diff^2 + (AICtmp$diff - Mod.avg.out[, 1])^2))
    }

    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Mod.avg.out[, 2] <- sqrt(sum(AICtmp$AICWt*(AICtmp$SE.diff^2 + (AICtmp$diff - Mod.avg.out[, 1])^2)))
    }  
    ##store table
    AICc.out <- AICtmp
  }

  
  ##begin loop - QAICc
  if(second.ord == FALSE && c.hat > 1){
      
    ##create temporary data.frame to store fitted values and SE 
    QAICtmp <- AICctab
    QAICtmp$diff <- diff
    QAICtmp$SE.diff <- SE.diff* sqrt(c.hat)
      
    ##compute model averaged prediction and store in output matrix
    Mod.avg.out[, 1] <- sum(QAICtmp$QAICWt*QAICtmp$diff)
      
    ##compute unconditional SE and store in output matrix
    if(identical(uncond.se, "old")) {
      Mod.avg.out[, 2] <- sum(QAICtmp$QAICWt*sqrt(QAICtmp$SE.diff^2 + (QAICtmp$diff - Mod.avg.out[, 1])^2))
    }
      
    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Mod.avg.out[, 2] <- sqrt(sum(QAICtmp$QAICWt*(QAICtmp$SE.diff^2 + (QAICtmp$diff - Mod.avg.out[, 1])^2)))
    }  
    ##store table
    AICc.out <- QAICtmp
  }

  Group.variable <- paste(parm.id, "(", var.id, ")")
  
  zcrit <- qnorm(p=(1-conf.level)/2, lower.tail=FALSE)
  Lower.CL <- Mod.avg.out[, 1] - zcrit * Mod.avg.out[, 2]
  Upper.CL <- Mod.avg.out[, 1] + zcrit * Mod.avg.out[, 2]
  Mod.eff.list <- list("Group.variable" = Group.variable, "Group1" = group1,
                       "Group2" = group2, "Type" = type, "Mod.avg.table" = AICc.out,
                       "Mod.avg.eff" = Mod.avg.out[,1], "Uncond.se" = Mod.avg.out[,2], "Conf.level" = conf.level,
                       "Lower.CL" = Lower.CL, "Upper.CL" = Upper.CL)
  class(Mod.eff.list) <- c("modavgEffect", "list")
  return(Mod.eff.list)  
}

##unmarkedFitPCO
modavgEffect.AICunmarkedFitPCO <- function(cand.set, modnames = NULL, newdata, second.ord = TRUE,
                                           nobs = NULL, uncond.se = "revised", conf.level = 0.95,
                                           type = "response", c.hat = 1, parm.type = NULL,
                                           ...) {
  
  ##check if named list if modnames are not supplied
  if(is.null(modnames)) {
    if(is.null(names(cand.set))) {
      modnames <- paste("Mod", 1:length(cand.set), sep = "")
      warning("\nModel names have been supplied automatically in the table\n")
    } else {
      modnames <- names(cand.set)
    }
  }

  ##rename values according to unmarked to extract from object
  ##gamma
  if(identical(parm.type, "gamma")) {
    parm.type1 <- "gamma"; parm.id <- "gam"
  }

  ##lambda
  if(identical(parm.type, "lambda")) {
    parm.type1 <- "lambda"; parm.id <- "lam"
    
    ##check mixture type for mixture models
    mixture.type <- sapply(X = cand.set, FUN = function(i) i@mixture)
    unique.mixture <- unique(mixture.type)
    if(length(unique.mixture) > 1) {
      if(any(unique.mixture == "ZIP")) stop("\nThis function does not yet support mixing ZIP with other distributions\n")
    } else {
      mixture.id <- unique(mixture.type)
      if(identical(unique.mixture, "ZIP")) {
        if(identical(type, "link")) stop("\nLink scale not yet supported for ZIP mixtures\n")
      }
    }
  }
  

  ##omega
  if(identical(parm.type, "omega")) {
    parm.type1 <- "omega"; parm.id <- "omega"
  }

  ##detect
  if(identical(parm.type, "detect")) {parm.type1 <- "det"; parm.id <- "p"}
  
    ##newdata is data frame with exact structure of the original data frame (same variable names and type)
  ##check on newdata
  ##determine number of observations in new data set
  nobserv <- nrow(newdata)
  if(nobserv > 2) stop("\nCurrent maximum number of groups compared is 2:\nmodify newdata argument accordingly\n")

  ##determine number of columns in new data set
  ncolumns <- ncol(newdata)

  ##if only 1 column, add an additional column to avoid problems in computation
  if(ncolumns == 1) newdata$blank.fake.column.NAs <- NA

  
  ##determine which column varies
  uniques <- apply(X = newdata, MARGIN = 2, FUN = unique)
  lengths <- lapply(X = uniques, FUN = length)
  varies <- sapply(X = lengths, FUN = function(i) i > 1)

    
  ##extract name of column
  if(sum(varies) == 1) {
    var.id <- names(varies)[which(varies == TRUE)]
                
    ##determine name of groups compared
    group1 <- as.character(newdata[,paste(var.id)][1])
    group2 <- as.character(newdata[,paste(var.id)][2])

  } else {
    ##warn that no single variable defines groups
    warning("\nGroups do not seem to be defined by a single variable.\n Function proceeding with generic group names\n")
    ##use generic names
    var.id <- "Groups"
    group1 <- "group 1"
    group2 <- "group 2"
  }

  
  
  ##number of models
  nmods <- length(modnames)


  ##compute predicted values
  ##point estimate
  if(identical(type, "response")) {
    ##extract fitted value for observation obs
    if(identical(parm.type, "lambda") && identical(mixture.id, "ZIP")) {
      fit <- matrix(data = unlist(lapply(X = cand.set, FUN = function(i)predictSE(i, se.fit = TRUE,
                                           newdata = newdata)$fit)),
                    nrow = nmods, ncol = 2, byrow = TRUE)
      
      SE <- matrix(data = unlist(lapply(X = cand.set, FUN = function(i)predictSE(i, se.fit = TRUE,
                                          newdata = newdata)$se.fit)),
                   nrow = nmods, ncol = 2, byrow = TRUE)
      
    } else {
      fit <- matrix(data = unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata,
                                                         type = parm.type1)$Predicted)),
                    nrow = nmods, ncol = 2, byrow = TRUE)

    ##extract SE for fitted value for observation obs
    SE <- matrix(data = unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata,
                                        type = parm.type1)$SE)),
                 nrow = nmods, ncol = 2, byrow = TRUE)
    }
  }
  
  ##link scale
  if(identical(type, "link")) {
    ##extract fitted value for observation obs
    fit <- matrix(data = unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata,
                                                       type = parm.type1, backTransform = FALSE)$Predicted)),
                  nrow = nmods, ncol = 2, byrow = TRUE)

    ##extract SE for fitted value for observation obs
    SE <- matrix(data = unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata,
                                                      type = parm.type1, backTransform = FALSE)$SE)),
                 nrow = nmods, ncol = 2, byrow = TRUE)
  }  

  
  ##difference between groups 
  diff <- fit[, 1] - fit[, 2]
    
  ##SE on difference
  SE.diff <- sqrt(SE[, 1]^2 + SE[, 2]^2)
  
  ##store AICc table
  AICctab <- aictab(cand.set = cand.set, modnames = modnames, second.ord = second.ord,
                    nobs = nobs, sort = FALSE, c.hat = c.hat)

  ##create object to hold Model-averaged estimates and unconditional SE's
  Mod.avg.out <- matrix(NA, nrow = 1, ncol = 2)
  colnames(Mod.avg.out) <- c("Mod.avg.diff", "Uncond.SE")

  ##begin loop - AICc
  if(second.ord == TRUE && c.hat == 1){
     
    ##create temporary data.frame to store fitted values and SE 
    AICctmp <- AICctab
    AICctmp$diff <- diff
    AICctmp$SE.diff <- SE.diff
    
    ##compute model averaged prediction and store in output matrix
    Mod.avg.out[, 1] <- sum(AICctmp$AICcWt*AICctmp$diff)
    ##compute unconditional SE and store in output matrix
      
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Mod.avg.out[, 2] <- sum(AICctmp$AICcWt*sqrt(AICctmp$SE.diff^2 + (AICctmp$diff - Mod.avg.out[, 1])^2))
    }

    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Mod.avg.out[, 2] <- sqrt(sum(AICctmp$AICcWt*(AICctmp$SE.diff^2 + (AICctmp$diff - Mod.avg.out[, 1])^2)))
    }
    ##store table
    AICc.out <- AICctmp
  }


  ##begin loop - QAICc
  if(second.ord == TRUE && c.hat > 1){
      ##create temporary data.frame to store fitted values and SE 
      QAICctmp <- AICctab
      QAICctmp$diff <- diff
      QAICctmp$SE.diff <- SE.diff * sqrt(c.hat)

      ##compute model averaged prediction and store in output matrix
      Mod.avg.out[, 1] <- sum(QAICctmp$QAICcWt*QAICctmp$diff)
      
      ##compute unconditional SE and store in output matrix
      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Mod.avg.out[, 2] <- sum(QAICctmp$QAICcWt*sqrt(QAICctmp$SE.diff^2 + (QAICctmp$diff - Mod.avg.out[, 1])^2))
      }
      
      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Mod.avg.out[, 2] <- sqrt(sum(QAICctmp$QAICcWt*(QAICctmp$SE.diff^2 + (QAICctmp$diff - Mod.avg.out[, 1])^2)))  
      }
      ##store table
      AICc.out <- QAICctmp
    
    }
  
  
  
  ##create temporary data.frame to store fitted values and SE - AIC
  if(second.ord == FALSE && c.hat == 1) {
    AICtmp <- AICctab
    AICtmp$diff <- diff
    AICtmp$SE.diff <- SE.diff

    ##compute model averaged prediction and store in output matrix
    Mod.avg.out[, 1] <- sum(AICtmp$AICWt*AICtmp$diff)

    ##compute unconditional SE and store in output matrix
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Mod.avg.out[, 2] <- sum(AICtmp$AICWt*sqrt(AICtmp$SE.diff^2 + (AICtmp$diff - Mod.avg.out[, 1])^2))
    }

    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Mod.avg.out[, 2] <- sqrt(sum(AICtmp$AICWt*(AICtmp$SE.diff^2 + (AICtmp$diff - Mod.avg.out[, 1])^2)))
    }  
    ##store table
    AICc.out <- AICtmp
  }

  
  ##begin loop - QAICc
  if(second.ord == FALSE && c.hat > 1){
      
    ##create temporary data.frame to store fitted values and SE 
    QAICtmp <- AICctab
    QAICtmp$diff <- diff
    QAICtmp$SE.diff <- SE.diff* sqrt(c.hat)
      
    ##compute model averaged prediction and store in output matrix
    Mod.avg.out[, 1] <- sum(QAICtmp$QAICWt*QAICtmp$diff)
      
    ##compute unconditional SE and store in output matrix
    if(identical(uncond.se, "old")) {
      Mod.avg.out[, 2] <- sum(QAICtmp$QAICWt*sqrt(QAICtmp$SE.diff^2 + (QAICtmp$diff - Mod.avg.out[, 1])^2))
    }
      
    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Mod.avg.out[, 2] <- sqrt(sum(QAICtmp$QAICWt*(QAICtmp$SE.diff^2 + (QAICtmp$diff - Mod.avg.out[, 1])^2)))
    }  
    ##store table
    AICc.out <- QAICtmp
  }

  Group.variable <- paste(parm.id, "(", var.id, ")")
  
  zcrit <- qnorm(p=(1-conf.level)/2, lower.tail=FALSE)
  Lower.CL <- Mod.avg.out[, 1] - zcrit * Mod.avg.out[, 2]
  Upper.CL <- Mod.avg.out[, 1] + zcrit * Mod.avg.out[, 2]
  Mod.eff.list <- list("Group.variable" = Group.variable, "Group1" = group1,
                       "Group2" = group2, "Type" = type, "Mod.avg.table" = AICc.out,
                       "Mod.avg.eff" = Mod.avg.out[,1], "Uncond.se" = Mod.avg.out[,2], "Conf.level" = conf.level,
                       "Lower.CL" = Lower.CL, "Upper.CL" = Upper.CL)
  class(Mod.eff.list) <- c("modavgEffect", "list")
  return(Mod.eff.list)  
}



##DS
modavgEffect.AICunmarkedFitDS <- function(cand.set, modnames = NULL, newdata, second.ord = TRUE,
                                          nobs = NULL, uncond.se = "revised", conf.level = 0.95,
                                          type = "response", c.hat = 1, parm.type = NULL,
                                          ...) {

  ##check if named list if modnames are not supplied
  if(is.null(modnames)) {
    if(is.null(names(cand.set))) {
      modnames <- paste("Mod", 1:length(cand.set), sep = "")
      warning("\nModel names have been supplied automatically in the table\n")
    } else {
      modnames <- names(cand.set)
    }
  }

  ##rename values according to unmarked to extract from object

  ##lambda
  if(identical(parm.type, "lambda")) {
    parm.type1 <- "state"; parm.id <- "lam"
    ##check mixture type for mixture models
  }

  ##detect
  if(identical(parm.type, "detect")) {
    stop("\nModel-averaging effect sizes on detection not yet supported for unmarkedFitDS class\n")
    ##parm.type1 <- "det"; parm.id <- "p"
  }
  ##newdata is data frame with exact structure of the original data frame (same variable names and type)
  ##check on newdata
  ##determine number of observations in new data set
  nobserv <- nrow(newdata)
  if(nobserv > 2) stop("\nCurrent maximum number of groups compared is 2:\nmodify newdata argument accordingly\n")

  ##determine number of columns in new data set
  ncolumns <- ncol(newdata)

  ##if only 1 column, add an additional column to avoid problems in computation
  if(ncolumns == 1) newdata$blank.fake.column.NAs <- NA
  

  ##determine which column varies
  uniques <- apply(X = newdata, MARGIN = 2, FUN = unique)
  lengths <- lapply(X = uniques, FUN = length)
  varies <- sapply(X = lengths, FUN = function(i) i > 1)

    
  ##extract name of column
  if(sum(varies) == 1) {
    var.id <- names(varies)[which(varies == TRUE)]
                
    ##determine name of groups compared
    group1 <- as.character(newdata[,paste(var.id)][1])
    group2 <- as.character(newdata[,paste(var.id)][2])

  } else {
    ##warn that no single variable defines groups
    warning("\nGroups do not seem to be defined by a single variable.\n Function proceeding with generic group names\n")
    ##use generic names
    var.id <- "Groups"
    group1 <- "group 1"
    group2 <- "group 2"
  }

  
  ##number of models
  nmods <- length(modnames)


  ##compute predicted values
  ##point estimate
  if(identical(type, "response")) {
    fit <- matrix(data = unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata,
                                                       type = parm.type1)$Predicted)),
                  nrow = nmods, ncol = 2, byrow = TRUE)

    ##extract SE for fitted value for observation obs
    SE <- matrix(data = unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata,
                                                      type = parm.type1)$SE)),
                 nrow = nmods, ncol = 2, byrow = TRUE)
  }

  
  ##link scale
  if(identical(type, "link")) {
    ##extract fitted value for observation obs
    fit <- matrix(data = unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata,
                                                       type = parm.type1, backTransform = FALSE)$Predicted)),
                  nrow = nmods, ncol = 2, byrow = TRUE)

    ##extract SE for fitted value for observation obs
    SE <- matrix(data = unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata,
                                                      type = parm.type1, backTransform = FALSE)$SE)),
                 nrow = nmods, ncol = 2, byrow = TRUE)
  }  

  
  ##difference between groups 
  diff <- fit[, 1] - fit[, 2]
    
  ##SE on difference
  SE.diff <- sqrt(SE[, 1]^2 + SE[, 2]^2)
  
  ##store AICc table
  AICctab <- aictab(cand.set = cand.set, modnames = modnames, second.ord = second.ord,
                    nobs = nobs, sort = FALSE, c.hat = c.hat)

  ##create object to hold Model-averaged estimates and unconditional SE's
  Mod.avg.out <- matrix(NA, nrow = 1, ncol = 2)
  colnames(Mod.avg.out) <- c("Mod.avg.diff", "Uncond.SE")

  ##begin loop - AICc
  if(second.ord == TRUE && c.hat == 1){
     
    ##create temporary data.frame to store fitted values and SE 
    AICctmp <- AICctab
    AICctmp$diff <- diff
    AICctmp$SE.diff <- SE.diff
    
    ##compute model averaged prediction and store in output matrix
    Mod.avg.out[, 1] <- sum(AICctmp$AICcWt*AICctmp$diff)
    ##compute unconditional SE and store in output matrix
      
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Mod.avg.out[, 2] <- sum(AICctmp$AICcWt*sqrt(AICctmp$SE.diff^2 + (AICctmp$diff - Mod.avg.out[, 1])^2))
    }

    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Mod.avg.out[, 2] <- sqrt(sum(AICctmp$AICcWt*(AICctmp$SE.diff^2 + (AICctmp$diff - Mod.avg.out[, 1])^2)))
    }
    ##store table
    AICc.out <- AICctmp
  }


  ##begin loop - QAICc
  if(second.ord == TRUE && c.hat > 1){
      ##create temporary data.frame to store fitted values and SE 
      QAICctmp <- AICctab
      QAICctmp$diff <- diff
      QAICctmp$SE.diff <- SE.diff * sqrt(c.hat)

      ##compute model averaged prediction and store in output matrix
      Mod.avg.out[, 1] <- sum(QAICctmp$QAICcWt*QAICctmp$diff)
      
      ##compute unconditional SE and store in output matrix
      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Mod.avg.out[, 2] <- sum(QAICctmp$QAICcWt*sqrt(QAICctmp$SE.diff^2 + (QAICctmp$diff - Mod.avg.out[, 1])^2))
      }
      
      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Mod.avg.out[, 2] <- sqrt(sum(QAICctmp$QAICcWt*(QAICctmp$SE.diff^2 + (QAICctmp$diff - Mod.avg.out[, 1])^2)))  
      }
      ##store table
      AICc.out <- QAICctmp
    
    }
  
  
  
  ##create temporary data.frame to store fitted values and SE - AIC
  if(second.ord == FALSE && c.hat == 1) {
    AICtmp <- AICctab
    AICtmp$diff <- diff
    AICtmp$SE.diff <- SE.diff

    ##compute model averaged prediction and store in output matrix
    Mod.avg.out[, 1] <- sum(AICtmp$AICWt*AICtmp$diff)

    ##compute unconditional SE and store in output matrix
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Mod.avg.out[, 2] <- sum(AICtmp$AICWt*sqrt(AICtmp$SE.diff^2 + (AICtmp$diff - Mod.avg.out[, 1])^2))
    }

    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Mod.avg.out[, 2] <- sqrt(sum(AICtmp$AICWt*(AICtmp$SE.diff^2 + (AICtmp$diff - Mod.avg.out[, 1])^2)))
    }  
    ##store table
    AICc.out <- AICtmp
  }

  
  ##begin loop - QAICc
  if(second.ord == FALSE && c.hat > 1){
      
    ##create temporary data.frame to store fitted values and SE 
    QAICtmp <- AICctab
    QAICtmp$diff <- diff
    QAICtmp$SE.diff <- SE.diff* sqrt(c.hat)
      
    ##compute model averaged prediction and store in output matrix
    Mod.avg.out[, 1] <- sum(QAICtmp$QAICWt*QAICtmp$diff)
      
    ##compute unconditional SE and store in output matrix
    if(identical(uncond.se, "old")) {
      Mod.avg.out[, 2] <- sum(QAICtmp$QAICWt*sqrt(QAICtmp$SE.diff^2 + (QAICtmp$diff - Mod.avg.out[, 1])^2))
    }
      
    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Mod.avg.out[, 2] <- sqrt(sum(QAICtmp$QAICWt*(QAICtmp$SE.diff^2 + (QAICtmp$diff - Mod.avg.out[, 1])^2)))
    }  
    ##store table
    AICc.out <- QAICtmp
  }

  Group.variable <- paste(parm.id, "(", var.id, ")")
  
  zcrit <- qnorm(p=(1-conf.level)/2, lower.tail=FALSE)
  Lower.CL <- Mod.avg.out[, 1] - zcrit * Mod.avg.out[, 2]
  Upper.CL <- Mod.avg.out[, 1] + zcrit * Mod.avg.out[, 2]
  Mod.eff.list <- list("Group.variable" = Group.variable, "Group1" = group1,
                       "Group2" = group2, "Type" = type, "Mod.avg.table" = AICc.out,
                       "Mod.avg.eff" = Mod.avg.out[,1], "Uncond.se" = Mod.avg.out[,2], "Conf.level" = conf.level,
                       "Lower.CL" = Lower.CL, "Upper.CL" = Upper.CL)
  class(Mod.eff.list) <- c("modavgEffect", "list")
  return(Mod.eff.list)  
}



##gdistsamp
modavgEffect.AICunmarkedFitGDS <- function(cand.set, modnames = NULL, newdata, second.ord = TRUE,
                                           nobs = NULL, uncond.se = "revised", conf.level = 0.95,
                                           type = "response", c.hat = 1, parm.type = NULL,
                                           ...) {
  
  ##check if named list if modnames are not supplied
  if(is.null(modnames)) {
    if(is.null(names(cand.set))) {
      modnames <- paste("Mod", 1:length(cand.set), sep = "")
      warning("\nModel names have been supplied automatically in the table\n")
    } else {
      modnames <- names(cand.set)
    }
  }

  ##rename values according to unmarked to extract from object

  ##lambda
  if(identical(parm.type, "lambda")) {
    parm.type1 <- "lambda"; parm.id <- "lam"
  }

  ##detect
  if(identical(parm.type, "detect")) {
    stop("\nModel-averaging effect sizes on detection not yet supported for unmarkedFitGDS class\n")
    ##parm.type1 <- "det"; parm.id <- "p"
  }
  ##availability
  if(identical(parm.type, "phi")) {parm.type1 <- "phi"; parm.id <- "phi"}
  
  ##newdata is data frame with exact structure of the original data frame (same variable names and type)
  ##check on newdata
  ##determine number of observations in new data set
  nobserv <- nrow(newdata)
  if(nobserv > 2) stop("\nCurrent maximum number of groups compared is 2:\nmodify newdata argument accordingly\n")

  ##determine number of columns in new data set
  ncolumns <- ncol(newdata)

  ##if only 1 column, add an additional column to avoid problems in computation
  if(ncolumns == 1) newdata$blank.fake.column.NAs <- NA
  

  ##determine which column varies
  uniques <- apply(X = newdata, MARGIN = 2, FUN = unique)
  lengths <- lapply(X = uniques, FUN = length)
  varies <- sapply(X = lengths, FUN = function(i) i > 1)

    
  ##extract name of column
  if(sum(varies) == 1) {
    var.id <- names(varies)[which(varies == TRUE)]
                
    ##determine name of groups compared
    group1 <- as.character(newdata[,paste(var.id)][1])
    group2 <- as.character(newdata[,paste(var.id)][2])

  } else {
    ##warn that no single variable defines groups
    warning("\nGroups do not seem to be defined by a single variable.\n Function proceeding with generic group names\n")
    ##use generic names
    var.id <- "Groups"
    group1 <- "group 1"
    group2 <- "group 2"
  }

  
  ##number of models
  nmods <- length(modnames)


  ##compute predicted values
  ##point estimate
  if(identical(type, "response")) {
    ##extract fitted value for observation obs
    fit <- matrix(data = unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata,
                                                       type = parm.type1)$Predicted)),
                  nrow = nmods, ncol = 2, byrow = TRUE)

    ##extract SE for fitted value for observation obs
    SE <- matrix(data = unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata,
                                        type = parm.type1)$SE)),
                 nrow = nmods, ncol = 2, byrow = TRUE)
  }
  
  
  ##link scale
  if(identical(type, "link")) {
    ##extract fitted value for observation obs
    fit <- matrix(data = unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata,
                                                       type = parm.type1, backTransform = FALSE)$Predicted)),
                  nrow = nmods, ncol = 2, byrow = TRUE)

    ##extract SE for fitted value for observation obs
    SE <- matrix(data = unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata,
                                                      type = parm.type1, backTransform = FALSE)$SE)),
                 nrow = nmods, ncol = 2, byrow = TRUE)
  }  

  
  ##difference between groups 
  diff <- fit[, 1] - fit[, 2]
    
  ##SE on difference
  SE.diff <- sqrt(SE[, 1]^2 + SE[, 2]^2)
  
  ##store AICc table
  AICctab <- aictab(cand.set = cand.set, modnames = modnames, second.ord = second.ord,
                    nobs = nobs, sort = FALSE, c.hat = c.hat)

  ##create object to hold Model-averaged estimates and unconditional SE's
  Mod.avg.out <- matrix(NA, nrow = 1, ncol = 2)
  colnames(Mod.avg.out) <- c("Mod.avg.diff", "Uncond.SE")

  ##begin loop - AICc
  if(second.ord == TRUE && c.hat == 1){
     
    ##create temporary data.frame to store fitted values and SE 
    AICctmp <- AICctab
    AICctmp$diff <- diff
    AICctmp$SE.diff <- SE.diff
    
    ##compute model averaged prediction and store in output matrix
    Mod.avg.out[, 1] <- sum(AICctmp$AICcWt*AICctmp$diff)
    ##compute unconditional SE and store in output matrix
      
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Mod.avg.out[, 2] <- sum(AICctmp$AICcWt*sqrt(AICctmp$SE.diff^2 + (AICctmp$diff - Mod.avg.out[, 1])^2))
    }

    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Mod.avg.out[, 2] <- sqrt(sum(AICctmp$AICcWt*(AICctmp$SE.diff^2 + (AICctmp$diff - Mod.avg.out[, 1])^2)))
    }
    ##store table
    AICc.out <- AICctmp
  }


  ##begin loop - QAICc
  if(second.ord == TRUE && c.hat > 1){
      ##create temporary data.frame to store fitted values and SE 
      QAICctmp <- AICctab
      QAICctmp$diff <- diff
      QAICctmp$SE.diff <- SE.diff * sqrt(c.hat)

      ##compute model averaged prediction and store in output matrix
      Mod.avg.out[, 1] <- sum(QAICctmp$QAICcWt*QAICctmp$diff)
      
      ##compute unconditional SE and store in output matrix
      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Mod.avg.out[, 2] <- sum(QAICctmp$QAICcWt*sqrt(QAICctmp$SE.diff^2 + (QAICctmp$diff - Mod.avg.out[, 1])^2))
      }
      
      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Mod.avg.out[, 2] <- sqrt(sum(QAICctmp$QAICcWt*(QAICctmp$SE.diff^2 + (QAICctmp$diff - Mod.avg.out[, 1])^2)))  
      }
      ##store table
      AICc.out <- QAICctmp
    
    }
  
  
  
  ##create temporary data.frame to store fitted values and SE - AIC
  if(second.ord == FALSE && c.hat == 1) {
    AICtmp <- AICctab
    AICtmp$diff <- diff
    AICtmp$SE.diff <- SE.diff

    ##compute model averaged prediction and store in output matrix
    Mod.avg.out[, 1] <- sum(AICtmp$AICWt*AICtmp$diff)

    ##compute unconditional SE and store in output matrix
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Mod.avg.out[, 2] <- sum(AICtmp$AICWt*sqrt(AICtmp$SE.diff^2 + (AICtmp$diff - Mod.avg.out[, 1])^2))
    }

    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Mod.avg.out[, 2] <- sqrt(sum(AICtmp$AICWt*(AICtmp$SE.diff^2 + (AICtmp$diff - Mod.avg.out[, 1])^2)))
    }  
    ##store table
    AICc.out <- AICtmp
  }

  
  ##begin loop - QAICc
  if(second.ord == FALSE && c.hat > 1){
      
    ##create temporary data.frame to store fitted values and SE 
    QAICtmp <- AICctab
    QAICtmp$diff <- diff
    QAICtmp$SE.diff <- SE.diff* sqrt(c.hat)
      
    ##compute model averaged prediction and store in output matrix
    Mod.avg.out[, 1] <- sum(QAICtmp$QAICWt*QAICtmp$diff)
      
    ##compute unconditional SE and store in output matrix
    if(identical(uncond.se, "old")) {
      Mod.avg.out[, 2] <- sum(QAICtmp$QAICWt*sqrt(QAICtmp$SE.diff^2 + (QAICtmp$diff - Mod.avg.out[, 1])^2))
    }
      
    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Mod.avg.out[, 2] <- sqrt(sum(QAICtmp$QAICWt*(QAICtmp$SE.diff^2 + (QAICtmp$diff - Mod.avg.out[, 1])^2)))
    }  
    ##store table
    AICc.out <- QAICtmp
  }

  Group.variable <- paste(parm.id, "(", var.id, ")")
  
  zcrit <- qnorm(p=(1-conf.level)/2, lower.tail=FALSE)
  Lower.CL <- Mod.avg.out[, 1] - zcrit * Mod.avg.out[, 2]
  Upper.CL <- Mod.avg.out[, 1] + zcrit * Mod.avg.out[, 2]
  Mod.eff.list <- list("Group.variable" = Group.variable, "Group1" = group1,
                       "Group2" = group2, "Type" = type, "Mod.avg.table" = AICc.out,
                       "Mod.avg.eff" = Mod.avg.out[,1], "Uncond.se" = Mod.avg.out[,2], "Conf.level" = conf.level,
                       "Lower.CL" = Lower.CL, "Upper.CL" = Upper.CL)
  class(Mod.eff.list) <- c("modavgEffect", "list")
  return(Mod.eff.list)  
}



##occuFP
modavgEffect.AICunmarkedFitOccuFP <- function(cand.set, modnames = NULL, newdata, second.ord = TRUE,
                                              nobs = NULL, uncond.se = "revised", conf.level = 0.95,
                                              type = "response", c.hat = 1, parm.type = NULL,
                                              ...) {

  ##check if named list if modnames are not supplied
  if(is.null(modnames)) {
    if(is.null(names(cand.set))) {
      modnames <- paste("Mod", 1:length(cand.set), sep = "")
      warning("\nModel names have been supplied automatically in the table\n")
    } else {
      modnames <- names(cand.set)
    }
  }

  ##rename values according to unmarked to extract from object
  ##psi
  if(identical(parm.type, "psi")) {
    parm.type1 <- "state"; parm.id <- "psi"
  }


  ##detect
  if(identical(parm.type, "detect")) {parm.type1 <- "det"; parm.id <- "p"}

  ##false positives
  if(identical(parm.type, "fp")) {parm.type1 <- "fp"; parm.id <- "fp"}
    
  ##newdata is data frame with exact structure of the original data frame (same variable names and type)
  ##check on newdata
  ##determine number of observations in new data set
  nobserv <- nrow(newdata)
  if(nobserv > 2) stop("\nCurrent maximum number of groups compared is 2:\nmodify newdata argument accordingly\n")

  ##determine number of columns in new data set
  ncolumns <- ncol(newdata)

  ##if only 1 column, add an additional column to avoid problems in computation
  if(ncolumns == 1) newdata$blank.fake.column.NAs <- NA
  

  ##determine which column varies
  uniques <- apply(X = newdata, MARGIN = 2, FUN = unique)
  lengths <- lapply(X = uniques, FUN = length)
  varies <- sapply(X = lengths, FUN = function(i) i > 1)

    
  ##extract name of column
  if(sum(varies) == 1) {
    var.id <- names(varies)[which(varies == TRUE)]
                
    ##determine name of groups compared
    group1 <- as.character(newdata[,paste(var.id)][1])
    group2 <- as.character(newdata[,paste(var.id)][2])

  } else {
    ##warn that no single variable defines groups
    warning("\nGroups do not seem to be defined by a single variable.\n Function proceeding with generic group names\n")
    ##use generic names
    var.id <- "Groups"
    group1 <- "group 1"
    group2 <- "group 2"
  }

  
  ##number of models
  nmods <- length(modnames)


  ##compute predicted values
  ##point estimate
  if(identical(type, "response")) {
    ##extract fitted value for observation obs
    fit <- matrix(data = unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata,
                                                       type = parm.type1)$Predicted)),
                  nrow = nmods, ncol = 2, byrow = TRUE)

    ##extract SE for fitted value for observation obs
    SE <- matrix(data = unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata,
                                                      type = parm.type1)$SE)),
                 nrow = nmods, ncol = 2, byrow = TRUE)
    }
  
  ##link scale
  if(identical(type, "link")) {
    ##extract fitted value for observation obs
    fit <- matrix(data = unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata,
                                                       type = parm.type1, backTransform = FALSE)$Predicted)),
                  nrow = nmods, ncol = 2, byrow = TRUE)

    ##extract SE for fitted value for observation obs
    SE <- matrix(data = unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata,
                                                      type = parm.type1, backTransform = FALSE)$SE)),
                 nrow = nmods, ncol = 2, byrow = TRUE)
  }  

  
  ##difference between groups 
  diff <- fit[, 1] - fit[, 2]
    
  ##SE on difference
  SE.diff <- sqrt(SE[, 1]^2 + SE[, 2]^2)
  
  ##store AICc table
  AICctab <- aictab(cand.set = cand.set, modnames = modnames, second.ord = second.ord,
                    nobs = nobs, sort = FALSE, c.hat = c.hat)

  ##create object to hold Model-averaged estimates and unconditional SE's
  Mod.avg.out <- matrix(NA, nrow = 1, ncol = 2)
  colnames(Mod.avg.out) <- c("Mod.avg.diff", "Uncond.SE")

  ##begin loop - AICc
  if(second.ord == TRUE && c.hat == 1){
     
    ##create temporary data.frame to store fitted values and SE 
    AICctmp <- AICctab
    AICctmp$diff <- diff
    AICctmp$SE.diff <- SE.diff
    
    ##compute model averaged prediction and store in output matrix
    Mod.avg.out[, 1] <- sum(AICctmp$AICcWt*AICctmp$diff)
    ##compute unconditional SE and store in output matrix
      
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Mod.avg.out[, 2] <- sum(AICctmp$AICcWt*sqrt(AICctmp$SE.diff^2 + (AICctmp$diff - Mod.avg.out[, 1])^2))
    }

    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Mod.avg.out[, 2] <- sqrt(sum(AICctmp$AICcWt*(AICctmp$SE.diff^2 + (AICctmp$diff - Mod.avg.out[, 1])^2)))
    }
    ##store table
    AICc.out <- AICctmp
  }


  ##begin loop - QAICc
  if(second.ord == TRUE && c.hat > 1){
      ##create temporary data.frame to store fitted values and SE 
      QAICctmp <- AICctab
      QAICctmp$diff <- diff
      QAICctmp$SE.diff <- SE.diff * sqrt(c.hat)

      ##compute model averaged prediction and store in output matrix
      Mod.avg.out[, 1] <- sum(QAICctmp$QAICcWt*QAICctmp$diff)
      
      ##compute unconditional SE and store in output matrix
      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Mod.avg.out[, 2] <- sum(QAICctmp$QAICcWt*sqrt(QAICctmp$SE.diff^2 + (QAICctmp$diff - Mod.avg.out[, 1])^2))
      }
      
      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Mod.avg.out[, 2] <- sqrt(sum(QAICctmp$QAICcWt*(QAICctmp$SE.diff^2 + (QAICctmp$diff - Mod.avg.out[, 1])^2)))  
      }
      ##store table
      AICc.out <- QAICctmp
    
    }
  
  
  
  ##create temporary data.frame to store fitted values and SE - AIC
  if(second.ord == FALSE && c.hat == 1) {
    AICtmp <- AICctab
    AICtmp$diff <- diff
    AICtmp$SE.diff <- SE.diff

    ##compute model averaged prediction and store in output matrix
    Mod.avg.out[, 1] <- sum(AICtmp$AICWt*AICtmp$diff)

    ##compute unconditional SE and store in output matrix
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Mod.avg.out[, 2] <- sum(AICtmp$AICWt*sqrt(AICtmp$SE.diff^2 + (AICtmp$diff - Mod.avg.out[, 1])^2))
    }

    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Mod.avg.out[, 2] <- sqrt(sum(AICtmp$AICWt*(AICtmp$SE.diff^2 + (AICtmp$diff - Mod.avg.out[, 1])^2)))
    }  
    ##store table
    AICc.out <- AICtmp
  }

  
  ##begin loop - QAICc
  if(second.ord == FALSE && c.hat > 1){
      
    ##create temporary data.frame to store fitted values and SE 
    QAICtmp <- AICctab
    QAICtmp$diff <- diff
    QAICtmp$SE.diff <- SE.diff* sqrt(c.hat)
      
    ##compute model averaged prediction and store in output matrix
    Mod.avg.out[, 1] <- sum(QAICtmp$QAICWt*QAICtmp$diff)
      
    ##compute unconditional SE and store in output matrix
    if(identical(uncond.se, "old")) {
      Mod.avg.out[, 2] <- sum(QAICtmp$QAICWt*sqrt(QAICtmp$SE.diff^2 + (QAICtmp$diff - Mod.avg.out[, 1])^2))
    }
      
    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Mod.avg.out[, 2] <- sqrt(sum(QAICtmp$QAICWt*(QAICtmp$SE.diff^2 + (QAICtmp$diff - Mod.avg.out[, 1])^2)))
    }  
    ##store table
    AICc.out <- QAICtmp
  }

  Group.variable <- paste(parm.id, "(", var.id, ")")
  
  zcrit <- qnorm(p=(1-conf.level)/2, lower.tail=FALSE)
  Lower.CL <- Mod.avg.out[, 1] - zcrit * Mod.avg.out[, 2]
  Upper.CL <- Mod.avg.out[, 1] + zcrit * Mod.avg.out[, 2]
  Mod.eff.list <- list("Group.variable" = Group.variable, "Group1" = group1,
                       "Group2" = group2, "Type" = type, "Mod.avg.table" = AICc.out,
                       "Mod.avg.eff" = Mod.avg.out[,1], "Uncond.se" = Mod.avg.out[,2], "Conf.level" = conf.level,
                       "Lower.CL" = Lower.CL, "Upper.CL" = Upper.CL)
  class(Mod.eff.list) <- c("modavgEffect", "list")
  return(Mod.eff.list)  
}



##multinomPois
modavgEffect.AICunmarkedFitMPois <- function(cand.set, modnames = NULL, newdata, second.ord = TRUE,
                                             nobs = NULL, uncond.se = "revised", conf.level = 0.95,
                                             type = "response", c.hat = 1, parm.type = NULL,
                                             ...) {

  ##check if named list if modnames are not supplied
  if(is.null(modnames)) {
    if(is.null(names(cand.set))) {
      modnames <- paste("Mod", 1:length(cand.set), sep = "")
      warning("\nModel names have been supplied automatically in the table\n")
    } else {
      modnames <- names(cand.set)
    }
  }

  ##rename values according to unmarked to extract from object
  ##lambda
  if(identical(parm.type, "lambda")) {
    parm.type1 <- "state"; parm.id <- "lam"
    ##set check to NULL for other models
    mixture.id <- NULL
  }
  

  ##detect
  if(identical(parm.type, "detect")) {parm.type1 <- "det"; parm.id <- "p"}
       
  ##newdata is data frame with exact structure of the original data frame (same variable names and type)
  ##check on newdata
  ##determine number of observations in new data set
  nobserv <- nrow(newdata)
  if(nobserv > 2) stop("\nCurrent maximum number of groups compared is 2:\nmodify newdata argument accordingly\n")

  ##determine number of columns in new data set
  ncolumns <- ncol(newdata)

  ##if only 1 column, add an additional column to avoid problems in computation
  if(ncolumns == 1) newdata$blank.fake.column.NAs <- NA
  

  ##determine which column varies
  uniques <- apply(X = newdata, MARGIN = 2, FUN = unique)
  lengths <- lapply(X = uniques, FUN = length)
  varies <- sapply(X = lengths, FUN = function(i) i > 1)

    
  ##extract name of column
  if(sum(varies) == 1) {
    var.id <- names(varies)[which(varies == TRUE)]
                
    ##determine name of groups compared
    group1 <- as.character(newdata[,paste(var.id)][1])
    group2 <- as.character(newdata[,paste(var.id)][2])

  } else {
    ##warn that no single variable defines groups
    warning("\nGroups do not seem to be defined by a single variable.\n Function proceeding with generic group names\n")
    ##use generic names
    var.id <- "Groups"
    group1 <- "group 1"
    group2 <- "group 2"
  }

  
  ##number of models
  nmods <- length(modnames)


  ##compute predicted values
  ##point estimate
  if(identical(type, "response")) {
    fit <- matrix(data = unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata,
                                                       type = parm.type1)$Predicted)),
                  nrow = nmods, ncol = 2, byrow = TRUE)
    
    ##extract SE for fitted value for observation obs
    SE <- matrix(data = unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata,
                                                      type = parm.type1)$SE)),
                 nrow = nmods, ncol = 2, byrow = TRUE)
  }
    
  ##link scale
  if(identical(type, "link")) {
    ##extract fitted value for observation obs
    fit <- matrix(data = unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata,
                                                       type = parm.type1, backTransform = FALSE)$Predicted)),
                  nrow = nmods, ncol = 2, byrow = TRUE)

    ##extract SE for fitted value for observation obs
    SE <- matrix(data = unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata,
                                                      type = parm.type1, backTransform = FALSE)$SE)),
                 nrow = nmods, ncol = 2, byrow = TRUE)
  }  

  
  ##difference between groups 
  diff <- fit[, 1] - fit[, 2]
    
  ##SE on difference
  SE.diff <- sqrt(SE[, 1]^2 + SE[, 2]^2)
  
  ##store AICc table
  AICctab <- aictab(cand.set = cand.set, modnames = modnames, second.ord = second.ord,
                    nobs = nobs, sort = FALSE, c.hat = c.hat)

  ##create object to hold Model-averaged estimates and unconditional SE's
  Mod.avg.out <- matrix(NA, nrow = 1, ncol = 2)
  colnames(Mod.avg.out) <- c("Mod.avg.diff", "Uncond.SE")

  ##begin loop - AICc
  if(second.ord == TRUE && c.hat == 1){
     
    ##create temporary data.frame to store fitted values and SE 
    AICctmp <- AICctab
    AICctmp$diff <- diff
    AICctmp$SE.diff <- SE.diff
    
    ##compute model averaged prediction and store in output matrix
    Mod.avg.out[, 1] <- sum(AICctmp$AICcWt*AICctmp$diff)
    ##compute unconditional SE and store in output matrix
      
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Mod.avg.out[, 2] <- sum(AICctmp$AICcWt*sqrt(AICctmp$SE.diff^2 + (AICctmp$diff - Mod.avg.out[, 1])^2))
    }

    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Mod.avg.out[, 2] <- sqrt(sum(AICctmp$AICcWt*(AICctmp$SE.diff^2 + (AICctmp$diff - Mod.avg.out[, 1])^2)))
    }
    ##store table
    AICc.out <- AICctmp
  }


  ##begin loop - QAICc
  if(second.ord == TRUE && c.hat > 1){
      ##create temporary data.frame to store fitted values and SE 
      QAICctmp <- AICctab
      QAICctmp$diff <- diff
      QAICctmp$SE.diff <- SE.diff * sqrt(c.hat)

      ##compute model averaged prediction and store in output matrix
      Mod.avg.out[, 1] <- sum(QAICctmp$QAICcWt*QAICctmp$diff)
      
      ##compute unconditional SE and store in output matrix
      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Mod.avg.out[, 2] <- sum(QAICctmp$QAICcWt*sqrt(QAICctmp$SE.diff^2 + (QAICctmp$diff - Mod.avg.out[, 1])^2))
      }
      
      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Mod.avg.out[, 2] <- sqrt(sum(QAICctmp$QAICcWt*(QAICctmp$SE.diff^2 + (QAICctmp$diff - Mod.avg.out[, 1])^2)))  
      }
      ##store table
      AICc.out <- QAICctmp
    
    }
  
  
  
  ##create temporary data.frame to store fitted values and SE - AIC
  if(second.ord == FALSE && c.hat == 1) {
    AICtmp <- AICctab
    AICtmp$diff <- diff
    AICtmp$SE.diff <- SE.diff

    ##compute model averaged prediction and store in output matrix
    Mod.avg.out[, 1] <- sum(AICtmp$AICWt*AICtmp$diff)

    ##compute unconditional SE and store in output matrix
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Mod.avg.out[, 2] <- sum(AICtmp$AICWt*sqrt(AICtmp$SE.diff^2 + (AICtmp$diff - Mod.avg.out[, 1])^2))
    }

    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Mod.avg.out[, 2] <- sqrt(sum(AICtmp$AICWt*(AICtmp$SE.diff^2 + (AICtmp$diff - Mod.avg.out[, 1])^2)))
    }  
    ##store table
    AICc.out <- AICtmp
  }

  
  ##begin loop - QAICc
  if(second.ord == FALSE && c.hat > 1){
      
    ##create temporary data.frame to store fitted values and SE 
    QAICtmp <- AICctab
    QAICtmp$diff <- diff
    QAICtmp$SE.diff <- SE.diff* sqrt(c.hat)
      
    ##compute model averaged prediction and store in output matrix
    Mod.avg.out[, 1] <- sum(QAICtmp$QAICWt*QAICtmp$diff)
      
    ##compute unconditional SE and store in output matrix
    if(identical(uncond.se, "old")) {
      Mod.avg.out[, 2] <- sum(QAICtmp$QAICWt*sqrt(QAICtmp$SE.diff^2 + (QAICtmp$diff - Mod.avg.out[, 1])^2))
    }
      
    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Mod.avg.out[, 2] <- sqrt(sum(QAICtmp$QAICWt*(QAICtmp$SE.diff^2 + (QAICtmp$diff - Mod.avg.out[, 1])^2)))
    }  
    ##store table
    AICc.out <- QAICtmp
  }

  Group.variable <- paste(parm.id, "(", var.id, ")")
  
  zcrit <- qnorm(p=(1-conf.level)/2, lower.tail=FALSE)
  Lower.CL <- Mod.avg.out[, 1] - zcrit * Mod.avg.out[, 2]
  Upper.CL <- Mod.avg.out[, 1] + zcrit * Mod.avg.out[, 2]
  Mod.eff.list <- list("Group.variable" = Group.variable, "Group1" = group1,
                       "Group2" = group2, "Type" = type, "Mod.avg.table" = AICc.out,
                       "Mod.avg.eff" = Mod.avg.out[,1], "Uncond.se" = Mod.avg.out[,2], "Conf.level" = conf.level,
                       "Lower.CL" = Lower.CL, "Upper.CL" = Upper.CL)
  class(Mod.eff.list) <- c("modavgEffect", "list")
  return(Mod.eff.list)  
}



##gmultmix
modavgEffect.AICunmarkedFitGMM <- function(cand.set, modnames = NULL, newdata, second.ord = TRUE,
                                           nobs = NULL, uncond.se = "revised", conf.level = 0.95,
                                           type = "response", c.hat = 1, parm.type = NULL,
                                           ...) {

  ##check if named list if modnames are not supplied
  if(is.null(modnames)) {
    if(is.null(names(cand.set))) {
      modnames <- paste("Mod", 1:length(cand.set), sep = "")
      warning("\nModel names have been supplied automatically in the table\n")
    } else {
      modnames <- names(cand.set)
    }
  }

  ##rename values according to unmarked to extract from object
  ##lambda
  if(identical(parm.type, "lambda")) {
    parm.type1 <- "lambda"; parm.id <- "lam"
  }

  ##detect
  if(identical(parm.type, "detect")) {parm.type1 <- "det"; parm.id <- "p"}
  
  ##availability
  if(identical(parm.type, "phi")) {parm.type1 <- "phi"; parm.id <- "phi"}
  
  ##newdata is data frame with exact structure of the original data frame (same variable names and type)
  ##check on newdata
  ##determine number of observations in new data set
  nobserv <- nrow(newdata)
  if(nobserv > 2) stop("\nCurrent maximum number of groups compared is 2:\nmodify newdata argument accordingly\n")

  ##determine number of columns in new data set
  ncolumns <- ncol(newdata)

  ##if only 1 column, add an additional column to avoid problems in computation
  if(ncolumns == 1) newdata$blank.fake.column.NAs <- NA

  
  ##determine which column varies
  uniques <- apply(X = newdata, MARGIN = 2, FUN = unique)
  lengths <- lapply(X = uniques, FUN = length)
  varies <- sapply(X = lengths, FUN = function(i) i > 1)

    
  ##extract name of column
  if(sum(varies) == 1) {
    var.id <- names(varies)[which(varies == TRUE)]
                
    ##determine name of groups compared
    group1 <- as.character(newdata[,paste(var.id)][1])
    group2 <- as.character(newdata[,paste(var.id)][2])

  } else {
    ##warn that no single variable defines groups
    warning("\nGroups do not seem to be defined by a single variable.\n Function proceeding with generic group names\n")
    ##use generic names
    var.id <- "Groups"
    group1 <- "group 1"
    group2 <- "group 2"
  }

  
  
  ##number of models
  nmods <- length(modnames)


  ##compute predicted values
  ##point estimate
  if(identical(type, "response")) {
    fit <- matrix(data = unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata,
                                                       type = parm.type1)$Predicted)),
                  nrow = nmods, ncol = 2, byrow = TRUE)
    
    ##extract SE for fitted value for observation obs
    SE <- matrix(data = unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata,
                                                      type = parm.type1)$SE)),
                 nrow = nmods, ncol = 2, byrow = TRUE)
  }
  
  ##link scale
  if(identical(type, "link")) {
    ##extract fitted value for observation obs
    fit <- matrix(data = unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata,
                                                       type = parm.type1, backTransform = FALSE)$Predicted)),
                  nrow = nmods, ncol = 2, byrow = TRUE)

    ##extract SE for fitted value for observation obs
    SE <- matrix(data = unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata,
                                                      type = parm.type1, backTransform = FALSE)$SE)),
                 nrow = nmods, ncol = 2, byrow = TRUE)
  }  
  
  
  ##difference between groups 
  diff <- fit[, 1] - fit[, 2]
    
  ##SE on difference
  SE.diff <- sqrt(SE[, 1]^2 + SE[, 2]^2)
  
  ##store AICc table
  AICctab <- aictab(cand.set = cand.set, modnames = modnames, second.ord = second.ord,
                    nobs = nobs, sort = FALSE, c.hat = c.hat)

  ##create object to hold Model-averaged estimates and unconditional SE's
  Mod.avg.out <- matrix(NA, nrow = 1, ncol = 2)
  colnames(Mod.avg.out) <- c("Mod.avg.diff", "Uncond.SE")

  ##begin loop - AICc
  if(second.ord == TRUE && c.hat == 1){
     
    ##create temporary data.frame to store fitted values and SE 
    AICctmp <- AICctab
    AICctmp$diff <- diff
    AICctmp$SE.diff <- SE.diff
    
    ##compute model averaged prediction and store in output matrix
    Mod.avg.out[, 1] <- sum(AICctmp$AICcWt*AICctmp$diff)
    ##compute unconditional SE and store in output matrix
      
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Mod.avg.out[, 2] <- sum(AICctmp$AICcWt*sqrt(AICctmp$SE.diff^2 + (AICctmp$diff - Mod.avg.out[, 1])^2))
    }

    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Mod.avg.out[, 2] <- sqrt(sum(AICctmp$AICcWt*(AICctmp$SE.diff^2 + (AICctmp$diff - Mod.avg.out[, 1])^2)))
    }
    ##store table
    AICc.out <- AICctmp
  }


  ##begin loop - QAICc
  if(second.ord == TRUE && c.hat > 1){
      ##create temporary data.frame to store fitted values and SE 
      QAICctmp <- AICctab
      QAICctmp$diff <- diff
      QAICctmp$SE.diff <- SE.diff * sqrt(c.hat)

      ##compute model averaged prediction and store in output matrix
      Mod.avg.out[, 1] <- sum(QAICctmp$QAICcWt*QAICctmp$diff)
      
      ##compute unconditional SE and store in output matrix
      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Mod.avg.out[, 2] <- sum(QAICctmp$QAICcWt*sqrt(QAICctmp$SE.diff^2 + (QAICctmp$diff - Mod.avg.out[, 1])^2))
      }
      
      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Mod.avg.out[, 2] <- sqrt(sum(QAICctmp$QAICcWt*(QAICctmp$SE.diff^2 + (QAICctmp$diff - Mod.avg.out[, 1])^2)))  
      }
      ##store table
      AICc.out <- QAICctmp
    
    }
  
  
  
  ##create temporary data.frame to store fitted values and SE - AIC
  if(second.ord == FALSE && c.hat == 1) {
    AICtmp <- AICctab
    AICtmp$diff <- diff
    AICtmp$SE.diff <- SE.diff

    ##compute model averaged prediction and store in output matrix
    Mod.avg.out[, 1] <- sum(AICtmp$AICWt*AICtmp$diff)

    ##compute unconditional SE and store in output matrix
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Mod.avg.out[, 2] <- sum(AICtmp$AICWt*sqrt(AICtmp$SE.diff^2 + (AICtmp$diff - Mod.avg.out[, 1])^2))
    }

    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Mod.avg.out[, 2] <- sqrt(sum(AICtmp$AICWt*(AICtmp$SE.diff^2 + (AICtmp$diff - Mod.avg.out[, 1])^2)))
    }  
    ##store table
    AICc.out <- AICtmp
  }

  
  ##begin loop - QAICc
  if(second.ord == FALSE && c.hat > 1){
      
    ##create temporary data.frame to store fitted values and SE 
    QAICtmp <- AICctab
    QAICtmp$diff <- diff
    QAICtmp$SE.diff <- SE.diff* sqrt(c.hat)
      
    ##compute model averaged prediction and store in output matrix
    Mod.avg.out[, 1] <- sum(QAICtmp$QAICWt*QAICtmp$diff)
      
    ##compute unconditional SE and store in output matrix
    if(identical(uncond.se, "old")) {
      Mod.avg.out[, 2] <- sum(QAICtmp$QAICWt*sqrt(QAICtmp$SE.diff^2 + (QAICtmp$diff - Mod.avg.out[, 1])^2))
    }
      
    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Mod.avg.out[, 2] <- sqrt(sum(QAICtmp$QAICWt*(QAICtmp$SE.diff^2 + (QAICtmp$diff - Mod.avg.out[, 1])^2)))
    }  
    ##store table
    AICc.out <- QAICtmp
  }

  Group.variable <- paste(parm.id, "(", var.id, ")")
  
  zcrit <- qnorm(p=(1-conf.level)/2, lower.tail=FALSE)
  Lower.CL <- Mod.avg.out[, 1] - zcrit * Mod.avg.out[, 2]
  Upper.CL <- Mod.avg.out[, 1] + zcrit * Mod.avg.out[, 2]
  Mod.eff.list <- list("Group.variable" = Group.variable, "Group1" = group1,
                       "Group2" = group2, "Type" = type, "Mod.avg.table" = AICc.out,
                       "Mod.avg.eff" = Mod.avg.out[,1], "Uncond.se" = Mod.avg.out[,2], "Conf.level" = conf.level,
                       "Lower.CL" = Lower.CL, "Upper.CL" = Upper.CL)
  class(Mod.eff.list) <- c("modavgEffect", "list")
  return(Mod.eff.list)  
}



##gpcount
modavgEffect.AICunmarkedFitGPC <- function(cand.set, modnames = NULL, newdata, second.ord = TRUE,
                                           nobs = NULL, uncond.se = "revised", conf.level = 0.95,
                                           type = "response", c.hat = 1, parm.type = NULL,
                                           ...) {

  ##check if named list if modnames are not supplied
  if(is.null(modnames)) {
    if(is.null(names(cand.set))) {
      modnames <- paste("Mod", 1:length(cand.set), sep = "")
      warning("\nModel names have been supplied automatically in the table\n")
    } else {
      modnames <- names(cand.set)
    }
  }

  ##rename values according to unmarked to extract from object

  ##lambda
  if(identical(parm.type, "lambda")) {
    parm.type1 <- "lambda"; parm.id <- "lam"
  }

  ##detect
  if(identical(parm.type, "detect")) {parm.type1 <- "det"; parm.id <- "p"}
  
  ##availability
  if(identical(parm.type, "phi")) {parm.type1 <- "phi"; parm.id <- "phi"}

  ##newdata is data frame with exact structure of the original data frame (same variable names and type)
  ##check on newdata
  ##determine number of observations in new data set
  nobserv <- nrow(newdata)
  if(nobserv > 2) stop("\nCurrent maximum number of groups compared is 2:\nmodify newdata argument accordingly\n")

  ##determine number of columns in new data set
  ncolumns <- ncol(newdata)

  ##if only 1 column, add an additional column to avoid problems in computation
  if(ncolumns == 1) newdata$blank.fake.column.NAs <- NA

  
  ##determine which column varies
  uniques <- apply(X = newdata, MARGIN = 2, FUN = unique)
  lengths <- lapply(X = uniques, FUN = length)
  varies <- sapply(X = lengths, FUN = function(i) i > 1)

    
  ##extract name of column
  if(sum(varies) == 1) {
    var.id <- names(varies)[which(varies == TRUE)]
                
    ##determine name of groups compared
    group1 <- as.character(newdata[,paste(var.id)][1])
    group2 <- as.character(newdata[,paste(var.id)][2])

  } else {
    ##warn that no single variable defines groups
    warning("\nGroups do not seem to be defined by a single variable.\n Function proceeding with generic group names\n")
    ##use generic names
    var.id <- "Groups"
    group1 <- "group 1"
    group2 <- "group 2"
  }

    
  ##number of models
  nmods <- length(modnames)


  ##compute predicted values
  ##point estimate
  if(identical(type, "response")) {
    fit <- matrix(data = unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata,
                                                       type = parm.type1)$Predicted)),
                  nrow = nmods, ncol = 2, byrow = TRUE)

    ##extract SE for fitted value for observation obs
    SE <- matrix(data = unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata,
                                                      type = parm.type1)$SE)),
                 nrow = nmods, ncol = 2, byrow = TRUE)
  }
  
  ##link scale
  if(identical(type, "link")) {
    ##extract fitted value for observation obs
    fit <- matrix(data = unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata,
                                                       type = parm.type1, backTransform = FALSE)$Predicted)),
                  nrow = nmods, ncol = 2, byrow = TRUE)

    ##extract SE for fitted value for observation obs
    SE <- matrix(data = unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata,
                                                      type = parm.type1, backTransform = FALSE)$SE)),
                 nrow = nmods, ncol = 2, byrow = TRUE)
  }  

  
  ##difference between groups 
  diff <- fit[, 1] - fit[, 2]
    
  ##SE on difference
  SE.diff <- sqrt(SE[, 1]^2 + SE[, 2]^2)
  
  ##store AICc table
  AICctab <- aictab(cand.set = cand.set, modnames = modnames, second.ord = second.ord,
                    nobs = nobs, sort = FALSE, c.hat = c.hat)

  ##create object to hold Model-averaged estimates and unconditional SE's
  Mod.avg.out <- matrix(NA, nrow = 1, ncol = 2)
  colnames(Mod.avg.out) <- c("Mod.avg.diff", "Uncond.SE")

  ##begin loop - AICc
  if(second.ord == TRUE && c.hat == 1){
     
    ##create temporary data.frame to store fitted values and SE 
    AICctmp <- AICctab
    AICctmp$diff <- diff
    AICctmp$SE.diff <- SE.diff
    
    ##compute model averaged prediction and store in output matrix
    Mod.avg.out[, 1] <- sum(AICctmp$AICcWt*AICctmp$diff)
    ##compute unconditional SE and store in output matrix
      
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Mod.avg.out[, 2] <- sum(AICctmp$AICcWt*sqrt(AICctmp$SE.diff^2 + (AICctmp$diff - Mod.avg.out[, 1])^2))
    }

    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Mod.avg.out[, 2] <- sqrt(sum(AICctmp$AICcWt*(AICctmp$SE.diff^2 + (AICctmp$diff - Mod.avg.out[, 1])^2)))
    }
    ##store table
    AICc.out <- AICctmp
  }


  ##begin loop - QAICc
  if(second.ord == TRUE && c.hat > 1){
      ##create temporary data.frame to store fitted values and SE 
      QAICctmp <- AICctab
      QAICctmp$diff <- diff
      QAICctmp$SE.diff <- SE.diff * sqrt(c.hat)

      ##compute model averaged prediction and store in output matrix
      Mod.avg.out[, 1] <- sum(QAICctmp$QAICcWt*QAICctmp$diff)
      
      ##compute unconditional SE and store in output matrix
      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Mod.avg.out[, 2] <- sum(QAICctmp$QAICcWt*sqrt(QAICctmp$SE.diff^2 + (QAICctmp$diff - Mod.avg.out[, 1])^2))
      }
      
      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Mod.avg.out[, 2] <- sqrt(sum(QAICctmp$QAICcWt*(QAICctmp$SE.diff^2 + (QAICctmp$diff - Mod.avg.out[, 1])^2)))  
      }
      ##store table
      AICc.out <- QAICctmp
    
    }
  
  
  
  ##create temporary data.frame to store fitted values and SE - AIC
  if(second.ord == FALSE && c.hat == 1) {
    AICtmp <- AICctab
    AICtmp$diff <- diff
    AICtmp$SE.diff <- SE.diff

    ##compute model averaged prediction and store in output matrix
    Mod.avg.out[, 1] <- sum(AICtmp$AICWt*AICtmp$diff)

    ##compute unconditional SE and store in output matrix
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Mod.avg.out[, 2] <- sum(AICtmp$AICWt*sqrt(AICtmp$SE.diff^2 + (AICtmp$diff - Mod.avg.out[, 1])^2))
    }

    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Mod.avg.out[, 2] <- sqrt(sum(AICtmp$AICWt*(AICtmp$SE.diff^2 + (AICtmp$diff - Mod.avg.out[, 1])^2)))
    }  
    ##store table
    AICc.out <- AICtmp
  }

  
  ##begin loop - QAICc
  if(second.ord == FALSE && c.hat > 1){
      
    ##create temporary data.frame to store fitted values and SE 
    QAICtmp <- AICctab
    QAICtmp$diff <- diff
    QAICtmp$SE.diff <- SE.diff* sqrt(c.hat)
      
    ##compute model averaged prediction and store in output matrix
    Mod.avg.out[, 1] <- sum(QAICtmp$QAICWt*QAICtmp$diff)
      
    ##compute unconditional SE and store in output matrix
    if(identical(uncond.se, "old")) {
      Mod.avg.out[, 2] <- sum(QAICtmp$QAICWt*sqrt(QAICtmp$SE.diff^2 + (QAICtmp$diff - Mod.avg.out[, 1])^2))
    }
      
    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Mod.avg.out[, 2] <- sqrt(sum(QAICtmp$QAICWt*(QAICtmp$SE.diff^2 + (QAICtmp$diff - Mod.avg.out[, 1])^2)))
    }  
    ##store table
    AICc.out <- QAICtmp
  }

  Group.variable <- paste(parm.id, "(", var.id, ")")
  
  zcrit <- qnorm(p=(1-conf.level)/2, lower.tail=FALSE)
  Lower.CL <- Mod.avg.out[, 1] - zcrit * Mod.avg.out[, 2]
  Upper.CL <- Mod.avg.out[, 1] + zcrit * Mod.avg.out[, 2]
  Mod.eff.list <- list("Group.variable" = Group.variable, "Group1" = group1,
                       "Group2" = group2, "Type" = type, "Mod.avg.table" = AICc.out,
                       "Mod.avg.eff" = Mod.avg.out[,1], "Uncond.se" = Mod.avg.out[,2], "Conf.level" = conf.level,
                       "Lower.CL" = Lower.CL, "Upper.CL" = Upper.CL)
  class(Mod.eff.list) <- c("modavgEffect", "list")
  return(Mod.eff.list)  
}



print.modavgEffect <- function(x, digits = 2, ...) {

  ##rework Group.variable labels
  old.type <- x$Group.variable
  stripped.type <- unlist(strsplit(old.type, split = "\\("))

  ic <- colnames(x$Mod.avg.table)[3]
  
  cat("\nModel-averaged effect size on the", x$Type, "scale based on entire model set:\n\n")
  
  ##extract elements
  if(length(stripped.type) == 1) {
    cat("\nMultimodel inference on \"", paste(x$Group.variable, x$Group1, sep = ""), " - ",
        paste(x$Group.variable, x$Group2, sep = ""), "\" based on ", ic, "\n", sep = "")
    
    ##if unmarkedFit model, then print differently
  } else {
    ##extract parameter name
    parm.type <- gsub("(^ +)|( +$)", "", stripped.type[1])
    
    ##extract Group.variable name
    var.id <- gsub("(^ +)|( +$)", "", unlist(strsplit(stripped.type[2], "\\)"))[1])

    cat("\nMultimodel inference on \"", paste(parm.type, "(", var.id, x$Group1, ")", sep = ""), " - ",
        paste(parm.type, "(", var.id, x$Group2, ")", sep = ""), "\" based on ", ic, "\n", sep = "")
  }

  
  cat("\n", ic, " table used to obtain model-averaged effect size:\n", sep = "")
  oldtab <- x$Mod.avg.table
  if (any(names(oldtab)=="c_hat")) {cat("\t(c-hat estimate = ", oldtab$c_hat[1], ")\n", sep = "")}
  cat("\n")
  if (any(names(oldtab)=="c_hat")) {
    nice.tab <- cbind(oldtab[,2], oldtab[,3], oldtab[,4], oldtab[,6],
                      oldtab[,9], oldtab[,10])
  } else {nice.tab <- cbind(oldtab[,2], oldtab[,3], oldtab[,4], oldtab[,6],
                            oldtab[,8], oldtab[,9])
        }

  colnames(nice.tab) <- c(colnames(oldtab)[c(2,3,4,6)], paste("Effect(", x$Group1, " - ", x$Group2, ")", sep = ""), "SE")
  rownames(nice.tab) <- oldtab[,1]
  print(round(nice.tab, digits=digits))
  cat("\nModel-averaged effect size:", eval(round(x$Mod.avg.eff, digits=digits)), "\n")
  cat("Unconditional SE:", eval(round(x$Uncond.se, digits=digits)), "\n")
  cat("",x$Conf.level * 100, "% Unconditional confidence interval: ", round(x$Lower.CL, digits=digits),
      ", ", round(x$Upper.CL, digits=digits), "\n\n", sep = "")
}

