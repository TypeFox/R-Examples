##generic
modavg <- function(cand.set, parm, modnames = NULL, second.ord = TRUE,
                   nobs = NULL, uncond.se = "revised", conf.level = 0.95, 
                   exclude = NULL, warn = TRUE, ...){
  cand.set <- formatCands(cand.set)
  UseMethod("modavg", cand.set)
}



##default
modavg.default <- function(cand.set, parm, modnames = NULL, second.ord = TRUE,
                           nobs = NULL,  uncond.se = "revised", conf.level = 0.95, 
                           exclude = NULL, warn = TRUE, ...){
  stop("\nFunction not yet defined for this object class\n")
}



##aov
modavg.AICaov.lm <-
  function(cand.set, parm, modnames = NULL, second.ord = TRUE,
           nobs = NULL,  uncond.se = "revised", conf.level = 0.95, 
           exclude = NULL, warn = TRUE, ...){

    ##check if named list if modnames are not supplied
    if(is.null(modnames)) {
      if(is.null(names(cand.set))) {
        modnames <- paste("Mod", 1:length(cand.set), sep = "")
        warning("\nModel names have been supplied automatically in the table\n")
      } else {
        modnames <- names(cand.set)
      }
    }


#####MODIFICATIONS BEGIN#######
    ##remove all leading and trailing white space and within parm
    parm <- gsub('[[:space:]]+', "", parm)
     
    ##reverse parm
    reversed.parm <- reverse.parm(parm)
    exclude <- reverse.exclude(exclude = exclude)
#####MODIFICATIONS END######

    ##extract model formula for each model in cand.set
    mod_formula <- lapply(cand.set, FUN = function(i) rownames(summary(i)$coefficients)) #extract model formula for each model in cand.set
    
    nmods <- length(cand.set)
  
    ##setup matrix to indicate presence of parm in the model
    include <- matrix(NA, nrow=nmods, ncol=1)
    ##add a check for multiple instances of same variable in given model (i.e., interactions)
    include.check <- matrix(NA, nrow=nmods, ncol=1)

    ##iterate over each formula in mod_formula list
    for (i in 1:nmods) {
      idents <- NULL
      idents.check <- NULL
      form <- mod_formula[[i]]

######################################################################################################
######################################################################################################
###MODIFICATIONS BEGIN
      ##iterate over each element of formula[[i]] in list
      if(is.null(reversed.parm)) {
        for (j in 1:length(form)) {
          idents[j] <- identical(parm, form[j])
          idents.check[j] <- ifelse(attr(regexpr(parm, form[j], fixed=TRUE), "match.length")== "-1" , 0, 1)  
        }
      } else {
        for (j in 1:length(form)) {
          idents[j] <- identical(parm, form[j]) | identical(reversed.parm, form[j])
          idents.check[j] <- ifelse(attr(regexpr(parm, form[j], fixed=TRUE), "match.length")=="-1" & attr(regexpr(reversed.parm, form[j],
                                                                  fixed=TRUE), "match.length")=="-1" , 0, 1)  
        }
      }
###MODIFICATIONS END
######################################################################################################
######################################################################################################


      include[i] <- ifelse(any(idents==1), 1, 0)
      include.check[i] <- ifelse(sum(idents.check)>1, "duplicates", "OK")
    }

#####################################################
    ##exclude == NULL; warn=TRUE:  warn that duplicates occur and stop
    if(is.null(exclude) && identical(warn, TRUE)) {
      if(any(include.check == "duplicates")) {
        stop("\nSome models possibly include more than one instance of the parameter of interest.\n",
             "This may be due to the presence of interaction/polynomial terms, or variables\n",
             "with similar names:\n",
             "\tsee \"?modavg\" for details on variable specification and \"exclude\" argument\n")
      }

    }

    ##exclude == NULL; warn=FALSE:  compute model-averaged beta estimate from models including variable of interest,
    ##assuming that the variable is not involved in interaction or higher order polynomial (x^2, x^3, etc...),
    ##warn that models were not excluded
    if(is.null(exclude) && identical(warn, FALSE)) {
      if(any(include.check == "duplicates")) {
        warning("Multiple instances of parameter of interest in given model is presumably\n",
                "not due to interaction or polynomial terms - these models will not be\n",
                "excluded from the computation of model-averaged estimate\n")
      }
      
    }



    ##warn if exclude is neither a list nor NULL
    if(!is.null(exclude)) {
      if(!is.list(exclude)) {stop("\nItems in \"exclude\" must be specified as a list\n")}
    }

  
    ##if exclude is list  
    if(is.list(exclude)) {

      ##determine number of elements in exclude
      nexcl <- length(exclude)

      ##check each formula for presence of exclude variable extracted with formula( )  
      not.include <- lapply(cand.set, FUN=formula)

      ##set up a new list with model formula
      forms <- list()
      for (i in 1:nmods) {
        form.tmp <- strsplit(as.character(not.include[i]), split="~")[[1]][-1]
        if(attr(regexpr("\\+", form.tmp), "match.length")==-1) {  #this line causes problems if intercept is removed from model
          forms[i] <- form.tmp
        } else {forms[i] <- strsplit(form.tmp, split=" \\+ ")}    #this line causes problems if intercept is removed from model
      }

      ##additional check to see whether some variable names include "+"
      check.forms <- unlist(lapply(forms, FUN=function(i) any(attr(regexpr("\\+", i), "match.length")>0)[[1]]))
      if (any(check.forms==TRUE)) stop("\nPlease avoid \"+\" in variable names\n")

      ##additional check to determine if intercept was removed from models
      check.forms <- unlist(lapply(forms, FUN=function(i) any(attr(regexpr("\\- 1", i), "match.length")>0)[[1]]))
      if (any(check.forms==TRUE)) stop("\nModels without intercept are not supported in this version, please use alternative parameterization\n")

 
      ##search within formula for variables to exclude
      mod.exclude <- matrix(NA, nrow=nmods, ncol=nexcl)

      ##iterate over each element in exclude list
      for (var in 1:nexcl) {

        ##iterate over each formula in mod_formula list
        for (i in 1:nmods) {
          idents <- NULL
          form.excl <- forms[[i]]
                    
          ##iterate over each element of forms[[i]]
          for (j in 1:length(form.excl)) {
            idents[j] <- identical(exclude[var][[1]], gsub("(^ +)|( +$)", "", form.excl[j]))
          }
          mod.exclude[i,var] <- ifelse(any(idents==1), 1, 0)
        }    
        
      }
  
      ##determine outcome across all variables to exclude
      to.exclude <- rowSums(mod.exclude)
  
  
      ##exclude models following models from model averaging  
      include[which(to.exclude>=1)] <- 0
      
      
    }



    ##add a check to determine if include always == 0
    if (sum(include)==0) {stop("\nParameter not found in any of the candidate models\n") }

    new.cand.set <- cand.set[which(include==1)] #select models including a given parameter
    new.mod.name <- modnames[which(include==1)]    #update model names
    ##

    new_table <- aictab(cand.set = new.cand.set, modnames = new.mod.name,
                        second.ord = second.ord, nobs = nobs, sort = FALSE)  #recompute AIC table and associated measures
    new_table$Beta_est <- unlist(lapply(new.cand.set, FUN = function(i) coef(i)[paste(parm)])) #extract beta estimate for parm
    new_table$SE <- unlist(lapply(new.cand.set, FUN = function(i) sqrt(diag(vcov(i)))[paste(parm)]))

   
    
    ##AICc
    ##compute model-averaged estimates, unconditional SE, and 95% CL
    if(second.ord == TRUE) {
      Modavg_beta <- sum(new_table$AICcWt*new_table$Beta_est)

      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Uncond_SE <- sum(new_table$AICcWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
      }
      
      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Uncond_SE <- sqrt(sum(new_table$AICcWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
      }
    }

    
    
    ##AIC
    if(second.ord == FALSE) {
      Modavg_beta <- sum(new_table$AICWt*new_table$Beta_est)
      
      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Uncond_SE <- sum(new_table$AICWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
      }
      
      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Uncond_SE <- sqrt(sum(new_table$AICWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
      }
    }

    
    zcrit <- qnorm(p=(1-conf.level)/2, lower.tail=FALSE)
    Lower_CL <- Modavg_beta-zcrit*Uncond_SE
    Upper_CL <- Modavg_beta+zcrit*Uncond_SE
    out.modavg <- list("Parameter"=paste(parm), "Mod.avg.table" = new_table, "Mod.avg.beta" = Modavg_beta,
                       "Uncond.SE" = Uncond_SE, "Conf.level" = conf.level, "Lower.CL"= Lower_CL, "Upper.CL" = Upper_CL)
  
    
    class(out.modavg) <- c("modavg", "list")
    return(out.modavg)

  }



##betareg
modavg.AICbetareg <-
  function(cand.set, parm, modnames = NULL, second.ord = TRUE,
           nobs = NULL,  uncond.se = "revised", conf.level = 0.95, 
           exclude = NULL, warn = TRUE, ...){

    ##check if named list if modnames are not supplied
    if(is.null(modnames)) {
      if(is.null(names(cand.set))) {
        modnames <- paste("Mod", 1:length(cand.set), sep = "")
        warning("\nModel names have been supplied automatically in the table\n")
      } else {
        modnames <- names(cand.set)
      }
    }


#####MODIFICATIONS BEGIN#######
    ##remove all leading and trailing white space and within parm
    parm <- gsub('[[:space:]]+', "", parm)
     
    ##reverse parm - problematic for parameters on "(phi)_temp:batch4" vs "batch4:(phi)_temp"
    reversed.parm <- reverse.parm(parm)
    exclude <- reverse.exclude(exclude = exclude)
#####MODIFICATIONS END######

    ##extract labels
    ##determine if parameter is on mean or phi
    if(regexpr(pattern = "\\(phi\\)_", parm) == "-1") {
      parm.phi <- NULL
      ##extract model formula for each model in cand.set
      mod_formula <- lapply(cand.set, FUN=function(i) rownames(summary(i)$coefficients$mean))
    } else {
      ##replace parm
      parm.phi <- gsub(pattern = "\\(phi\\)_", "", parm)
      if(regexpr(pattern = ":", parm) != "-1") {
        warning(cat("\nthis function does not yet support interaction terms on phi:\n",
                    "use 'modavgCustom' instead\n"))
      }
      mod_formula <- lapply(cand.set, FUN=function(i) rownames(summary(i)$coefficients$precision))
    }

    
    nmods <- length(cand.set)
  
    ##setup matrix to indicate presence of parm in the model
    include <- matrix(NA, nrow=nmods, ncol=1)
    ##add a check for multiple instances of same variable in given model (i.e., interactions)
    include.check <- matrix(NA, nrow=nmods, ncol=1)

    
    ##iterate over each formula in mod_formula list
    ##if parameters on mean
    if(is.null(parm.phi)) {
      for (i in 1:nmods) {
        idents <- NULL
        idents.check <- NULL
        form <- mod_formula[[i]]

######################################################################################################
######################################################################################################
###MODIFICATIONS BEGIN
        ##iterate over each element of formula[[i]] in list
        if(is.null(reversed.parm)) {
          for (j in 1:length(form)) {
            idents[j] <- identical(parm, form[j])
            idents.check[j] <- ifelse(attr(regexpr(parm, form[j], fixed=TRUE), "match.length")== "-1" , 0, 1)  
          }
        } else {
          for (j in 1:length(form)) {
            idents[j] <- identical(parm, form[j]) | identical(reversed.parm, form[j])
            idents.check[j] <- ifelse(attr(regexpr(parm, form[j], fixed=TRUE), "match.length")=="-1" & attr(regexpr(reversed.parm, form[j],
                                                                    fixed=TRUE), "match.length")=="-1" , 0, 1)  
          }
        }
###MODIFICATIONS END
######################################################################################################
######################################################################################################


        include[i] <- ifelse(any(idents==1), 1, 0)
        include.check[i] <- ifelse(sum(idents.check)>1, "duplicates", "OK")
      }
    }



    
    ##if parameters on phi
    if(!is.null(parm.phi)) {
      for (i in 1:nmods) {
        idents <- NULL
        idents.check <- NULL
        form <- mod_formula[[i]]

######################################################################################################
######################################################################################################
###MODIFICATIONS BEGIN
        ##iterate over each element of formula[[i]] in list
        #if(is.null(reversed.parm)) {
        ##do not consider reversed parm here because of "(phi)_" prefix in coefficients
          for (j in 1:length(form)) {
            idents[j] <- identical(parm.phi, form[j])
            idents.check[j] <- ifelse(attr(regexpr(parm.phi, form[j], fixed=TRUE), "match.length")== "-1" , 0, 1)  
          }
#        } else {
#          for (j in 1:length(form)) {
#            idents[j] <- identical(parm.phi, form[j]) | identical(reversed.parm, form[j])
#            idents.check[j] <- ifelse(attr(regexpr(parm, form[j], fixed=TRUE), "match.length")=="-1" & attr(regexpr(reversed.parm, form[j],
#                                                                    fixed=TRUE), "match.length")=="-1" , 0, 1)  
#          }
#        }
###MODIFICATIONS END
######################################################################################################
######################################################################################################

        
        include[i] <- ifelse(any(idents==1), 1, 0)
        include.check[i] <- ifelse(sum(idents.check)>1, "duplicates", "OK")
      }
    }


    

#####################################################
    ##exclude == NULL; warn=TRUE:  warn that duplicates occur and stop
    if(is.null(exclude) && identical(warn, TRUE)) {
      if(any(include.check == "duplicates")) {
        stop("\nSome models possibly include more than one instance of the parameter of interest.\n",
             "This may be due to the presence of interaction/polynomial terms, or variables\n",
             "with similar names:\n",
             "\tsee \"?modavg\" for details on variable specification and \"exclude\" argument\n")
      }

    }

    ##exclude == NULL; warn=FALSE:  compute model-averaged beta estimate from models including variable of interest,
    ##assuming that the variable is not involved in interaction or higher order polynomial (x^2, x^3, etc...),
    ##warn that models were not excluded
    if(is.null(exclude) && identical(warn, FALSE)) {
      if(any(include.check == "duplicates")) {
        warning("Multiple instances of parameter of interest in given model is presumably\n",
                "not due to interaction or polynomial terms - these models will not be\n",
                "excluded from the computation of model-averaged estimate\n")
      }
      
    }



    ##warn if exclude is neither a list nor NULL
    if(!is.null(exclude)) {
      if(!is.list(exclude)) {stop("\nItems in \"exclude\" must be specified as a list\n")}
    }

  
    ##if exclude is list  
    if(is.list(exclude)) {

      ##determine number of elements in exclude
      nexcl <- length(exclude)

      ##check each formula for presence of exclude variable extracted with formula( )  
      not.include <- lapply(cand.set, FUN=formula)

      ##set up a new list with model formula
      forms <- list()
      for (i in 1:nmods) {
        form.tmp <- strsplit(as.character(not.include[i]), split="~")[[1]][-1]
        if(attr(regexpr("\\+", form.tmp), "match.length")==-1) {  #this line causes problems if intercept is removed from model
          forms[i] <- form.tmp
        } else {forms[i] <- strsplit(form.tmp, split=" \\+ ")}    #this line causes problems if intercept is removed from model
      }

      ##additional check to see whether some variable names include "+"
      check.forms <- unlist(lapply(forms, FUN=function(i) any(attr(regexpr("\\+", i), "match.length")>0)[[1]]))
      if (any(check.forms==TRUE)) stop("\nPlease avoid \"+\" in variable names\n")

      ##additional check to determine if intercept was removed from models
      check.forms <- unlist(lapply(forms, FUN=function(i) any(attr(regexpr("\\- 1", i), "match.length")>0)[[1]]))
      if (any(check.forms==TRUE)) stop("\nModels without intercept are not supported in this version, please use alternative parameterization\n")

 
      ##search within formula for variables to exclude
      mod.exclude <- matrix(NA, nrow=nmods, ncol=nexcl)

      ##iterate over each element in exclude list
      for (var in 1:nexcl) {

        ##iterate over each formula in mod_formula list
        for (i in 1:nmods) {
          idents <- NULL
          form.excl <- forms[[i]]
                    
          ##iterate over each element of forms[[i]]
          for (j in 1:length(form.excl)) {
            idents[j] <- identical(exclude[var][[1]], gsub("(^ +)|( +$)", "", form.excl[j]))
          }
          mod.exclude[i,var] <- ifelse(any(idents==1), 1, 0)
        }    
        
      }
  
      ##determine outcome across all variables to exclude
      to.exclude <- rowSums(mod.exclude)
  
  
      ##exclude models following models from model averaging  
      include[which(to.exclude>=1)] <- 0
      
      
    }



    ##add a check to determine if include always == 0
    if (sum(include)==0) {stop("\nParameter not found in any of the candidate models\n") }

    new.cand.set <- cand.set[which(include==1)] #select models including a given parameter
    new.mod.name <- modnames[which(include==1)]    #update model names
    ##

    new_table <- aictab(cand.set = new.cand.set, modnames = new.mod.name,
                        second.ord = second.ord, nobs = nobs, sort = FALSE)  #recompute AIC table and associated measures
    new_table$Beta_est <- unlist(lapply(new.cand.set, FUN = function(i) coef(i)[paste(parm)])) #extract beta estimate for parm
    new_table$SE <- unlist(lapply(new.cand.set, FUN = function(i) sqrt(diag(vcov(i)))[paste(parm)]))

   
    
    ##AICc
    ##compute model-averaged estimates, unconditional SE, and 95% CL
    if(second.ord == TRUE) {
      Modavg_beta <- sum(new_table$AICcWt*new_table$Beta_est)

      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Uncond_SE <- sum(new_table$AICcWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
      }
      
      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Uncond_SE <- sqrt(sum(new_table$AICcWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
      }
    }

    
    
    ##AIC
    if(second.ord == FALSE) {
      Modavg_beta <- sum(new_table$AICWt*new_table$Beta_est)
      
      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Uncond_SE <- sum(new_table$AICWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
      }
      
      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Uncond_SE <- sqrt(sum(new_table$AICWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
      }
    }

    
    zcrit <- qnorm(p=(1-conf.level)/2, lower.tail=FALSE)
    Lower_CL <- Modavg_beta-zcrit*Uncond_SE
    Upper_CL <- Modavg_beta+zcrit*Uncond_SE
    out.modavg <- list("Parameter"=paste(parm), "Mod.avg.table" = new_table, "Mod.avg.beta" = Modavg_beta,
                       "Uncond.SE" = Uncond_SE, "Conf.level" = conf.level, "Lower.CL"= Lower_CL, "Upper.CL" = Upper_CL)
  
    
    class(out.modavg) <- c("modavg", "list")
    return(out.modavg)

  }



##clm
modavg.AICsclm.clm <-
  function(cand.set, parm, modnames = NULL, second.ord = TRUE,
           nobs = NULL, uncond.se = "revised", conf.level = 0.95,
           exclude = NULL, warn = TRUE, ...){

    ##check if named list if modnames are not supplied
    if(is.null(modnames)) {
      if(is.null(names(cand.set))) {
        modnames <- paste("Mod", 1:length(cand.set), sep = "")
        warning("\nModel names have been supplied automatically in the table\n")
      } else {
        modnames <- names(cand.set)
      }
    }
    

    ##check that link function is the same for all models
    check.link <- unlist(lapply(X = cand.set, FUN=function(i) i$link))
    unique.link <- unique(x = check.link)
    if(length(unique.link) > 1) stop("\nIt is not appropriate to compute a model-averaged beta estimate\n",
                                         "from models using different link functions\n")
  
    
#####MODIFICATIONS BEGIN#######
    ##remove all leading and trailing white space and within parm
    parm <- gsub('[[:space:]]+', "", parm)
  
    ##reverse parm
    reversed.parm <- reverse.parm(parm)
    exclude <- reverse.exclude(exclude = exclude)
#####MODIFICATIONS END######

  
    ##extract model formula for each model in cand.set
    mod_formula <- lapply(cand.set, FUN = function(i) rownames(summary(i)$coefficients))
       
    nmods <- length(cand.set)
  
    ##setup matrix to indicate presence of parms in the model
    include <- matrix(NA, nrow = nmods, ncol = 1)
    ##add a check for multiple instances of same variable in given model (i.e., interactions)
    include.check <- matrix(NA, nrow = nmods, ncol = 1)
    
    ##iterate over each formula in mod_formula list
    for (i in 1:nmods) {
      idents <- NULL
      idents.check <- NULL
      form <- mod_formula[[i]]

######################################################################################################
######################################################################################################
###MODIFICATIONS BEGIN
      ##iterate over each element of formula[[i]] in list
      if(is.null(reversed.parm)) {
        for (j in 1:length(form)) {
          idents[j] <- identical(parm, form[j])
          idents.check[j] <- ifelse(attr(regexpr(parm, form[j], fixed=TRUE), "match.length")== "-1" , 0, 1)  
        }
      } else {
        for (j in 1:length(form)) {
          idents[j] <- identical(parm, form[j]) | identical(reversed.parm, form[j])
          idents.check[j] <- ifelse(attr(regexpr(parm, form[j], fixed=TRUE), "match.length")=="-1" & attr(regexpr(reversed.parm, form[j],
                                                                  fixed=TRUE), "match.length")=="-1" , 0, 1)  
        }
      }
###MODIFICATIONS END
######################################################################################################
######################################################################################################

  
      include[i] <- ifelse(any(idents==1), 1, 0)
      include.check[i] <- ifelse(sum(idents.check)>1, "duplicates", "OK")
    }

  #####################################################
    ##exclude == NULL; warn=TRUE:  warn that duplicates occur and stop
    if(is.null(exclude) && identical(warn, TRUE)) {
      ##check for duplicates in same model
      if(any(include.check == "duplicates")) {
        stop("\nSome models include more than one instance of the parameter of interest. \n",
             "This may be due to the presence of interaction/polynomial terms, or variables\n",
             "with similar names:\n",
             "\tsee \"?modavg\" for details on variable specification and \"exclude\" argument\n")
      }

    }

    ##exclude == NULL; warn=FALSE:  compute model-averaged beta estimate from models including variable of interest,
    ##assuming that the variable is not involved in interaction or higher order polynomial (x^2, x^3, etc...),
    ##warn that models were not excluded
    if(is.null(exclude) && identical(warn, FALSE)) {
      if(any(include.check == "duplicates")) {
        warning("\nMultiple instances of parameter of interest in given model is presumably\n",
                "not due to interaction or polynomial terms - these models will not be\n",
                "excluded from the computation of model-averaged estimate\n")
      }
      
    }

    ##warn if exclude is neither a list nor NULL
    if(!is.null(exclude)) {
      if(!is.list(exclude)) {stop("\nItems in \"exclude\" must be specified as a list\n")}
    }


    ##if exclude is list  
    if(is.list(exclude)) {
      
      ##determine number of elements in exclude
      nexcl <- length(exclude)
      
      ##check each formula for presence of exclude variable extracted with formula( )  
      not.include <- lapply(cand.set, FUN = formula)

      ##set up a new list with model formula
      forms <- list()
      for (i in 1:nmods) {
        form.tmp <- as.character(not.include[[i]])[3]
        if(attr(regexpr("\\+", form.tmp), "match.length")==-1) {
          forms[i] <- form.tmp
        } else {forms[i] <- strsplit(form.tmp, split=" \\+ ")}
      }

      ##additional check to see whether some variable names include "+"
      check.forms <- unlist(lapply(forms, FUN=function(i) any(attr(regexpr("\\+", i), "match.length")>0)[[1]]))
      if (any(check.forms==TRUE)) stop("\nPlease avoid \"+\" in variable names\n")

      ##additional check to determine if intercept was removed from models
      check.forms <- unlist(lapply(forms, FUN=function(i) any(attr(regexpr("\\- 1", i), "match.length")>0)[[1]]))
      if (any(check.forms==TRUE)) stop("\nModels without intercept are not supported in this version, please use alternative parameterization\n")
         
      
      ##search within formula for variables to exclude
      mod.exclude <- matrix(NA, nrow=nmods, ncol=nexcl)
      
      ##iterate over each element in exclude list
      for (var in 1:nexcl) {

        ##iterate over each formula in mod_formula list
        for (i in 1:nmods) {
          idents <- NULL
          form.excl <- forms[[i]]
                    
          ##iterate over each element of forms[[i]]
          for (j in 1:length(form.excl)) {
            idents[j] <- identical(exclude[var][[1]], gsub("(^ +)|( +$)", "", form.excl[j]))
          }
          mod.exclude[i,var] <- ifelse(any(idents==1), 1, 0)
        }    
        
      }
  
      ##determine outcome across all variables to exclude
      to.exclude <- rowSums(mod.exclude)
         
      
      ##exclude models following models from model averaging  
      include[which(to.exclude>=1)] <- 0
      
      
    }

       
       
    ##add a check to determine if include always == 0
    if (sum(include)==0) {stop("\nParameter not found in any of the candidate models\n") }
    
    new.cand.set <- cand.set[which(include==1)] #select models including a given parameter
    new.mod.name <- modnames[which(include==1)]    #update model names
    
    new_table <- aictab(cand.set=new.cand.set, modnames=new.mod.name,
                        second.ord=second.ord, nobs=nobs, sort=FALSE)  #recompute AIC table and associated measures
    new_table$Beta_est <- unlist(lapply(new.cand.set, FUN=function(i) coef(i)[paste(parm)])) #extract beta estimate for parm
    new_table$SE <- unlist(lapply(new.cand.set, FUN=function(i) sqrt(diag(vcov(i)))[paste(parm)]))
  
    ##compute model-averaged estimates, unconditional SE, and 95% CL
    if(second.ord == TRUE) {
      Modavg_beta <- sum(new_table$AICcWt*new_table$Beta_est)
         
      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Uncond_SE <- sum(new_table$AICcWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
      }

      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Uncond_SE <- sqrt(sum(new_table$AICcWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
      } 


    } else {
      Modavg_beta <- sum(new_table$AICWt*new_table$Beta_est)
      
      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Uncond_SE <- sum(new_table$AICWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
      }

      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Uncond_SE <- sqrt(sum(new_table$AICWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
      }
    
      
    }
     
    zcrit <- qnorm(p=(1-conf.level)/2, lower.tail=FALSE)
    Lower_CL <- Modavg_beta - zcrit*Uncond_SE
    Upper_CL <- Modavg_beta + zcrit*Uncond_SE
    out.modavg <- list("Parameter"=paste(parm), "Mod.avg.table" = new_table, "Mod.avg.beta" = Modavg_beta,
                       "Uncond.SE" = Uncond_SE, "Conf.level" = conf.level, "Lower.CL" = Lower_CL,
                       "Upper.CL" = Upper_CL)
    
    class(out.modavg) <- c("modavg", "list")
    return(out.modavg)
  }



##clmm
modavg.AICclmm <-
  function(cand.set, parm, modnames = NULL, second.ord = TRUE,
           nobs = NULL, uncond.se = "revised", conf.level = 0.95,
           exclude = NULL, warn = TRUE, ...){

    ##check if named list if modnames are not supplied
    if(is.null(modnames)) {
      if(is.null(names(cand.set))) {
        modnames <- paste("Mod", 1:length(cand.set), sep = "")
        warning("\nModel names have been supplied automatically in the table\n")
      } else {
        modnames <- names(cand.set)
      }
    }
    
    
    ##check that link function is the same for all models
    check.link <- unlist(lapply(X = cand.set, FUN=function(i) i$link))
    unique.link <- unique(x = check.link)
    if(length(unique.link) > 1) stop("\nIt is not appropriate to compute a model-averaged beta estimate\n",
                                         "from models using different link functions\n")
  
    
#####MODIFICATIONS BEGIN#######
    ##remove all leading and trailing white space and within parm
    parm <- gsub('[[:space:]]+', "", parm)
  
    ##reverse parm
    reversed.parm <- reverse.parm(parm)
    exclude <- reverse.exclude(exclude = exclude)
#####MODIFICATIONS END######
    
  
    ##extract model formula for each model in cand.set
    mod_formula <- lapply(cand.set, FUN = function(i) rownames(summary(i)$coefficients))
       
    nmods <- length(cand.set)
  
    ##setup matrix to indicate presence of parms in the model
    include <- matrix(NA, nrow = nmods, ncol = 1)
    ##add a check for multiple instances of same variable in given model (i.e., interactions)
    include.check <- matrix(NA, nrow = nmods, ncol = 1)

    ##iterate over each formula in mod_formula list
    for (i in 1:nmods) {
      idents <- NULL
      idents.check <- NULL
      form <- mod_formula[[i]]

######################################################################################################
######################################################################################################
###MODIFICATIONS BEGIN
      ##iterate over each element of formula[[i]] in list
      if(is.null(reversed.parm)) {
        for (j in 1:length(form)) {
          idents[j] <- identical(parm, form[j])
          idents.check[j] <- ifelse(attr(regexpr(parm, form[j], fixed=TRUE), "match.length")== "-1" , 0, 1)  
        }
      } else {
        for (j in 1:length(form)) {
          idents[j] <- identical(parm, form[j]) | identical(reversed.parm, form[j])
          idents.check[j] <- ifelse(attr(regexpr(parm, form[j], fixed=TRUE), "match.length")=="-1" & attr(regexpr(reversed.parm, form[j],
                                                                  fixed=TRUE), "match.length")=="-1" , 0, 1)  
        }
      }
###MODIFICATIONS END
######################################################################################################
######################################################################################################

  
      include[i] <- ifelse(any(idents==1), 1, 0)
      include.check[i] <- ifelse(sum(idents.check)>1, "duplicates", "OK")
    }

  #####################################################
    ##exclude == NULL; warn=TRUE:  warn that duplicates occur and stop
    if(is.null(exclude) && identical(warn, TRUE)) {
      ##check for duplicates in same model
      if(any(include.check == "duplicates")) {
        stop("\nSome models include more than one instance of the parameter of interest. \n",
             "This may be due to the presence of interaction/polynomial terms, or variables\n",
             "with similar names:\n",
             "\tsee \"?modavg\" for details on variable specification and \"exclude\" argument\n")
      }

    }

    ##exclude == NULL; warn=FALSE:  compute model-averaged beta estimate from models including variable of interest,
    ##assuming that the variable is not involved in interaction or higher order polynomial (x^2, x^3, etc...),
    ##warn that models were not excluded
    if(is.null(exclude) && identical(warn, FALSE)) {
      if(any(include.check == "duplicates")) {
        warning("\nMultiple instances of parameter of interest in given model is presumably\n",
                "not due to interaction or polynomial terms - these models will not be\n",
                "excluded from the computation of model-averaged estimate\n")
      }
      
    }

    ##warn if exclude is neither a list nor NULL
    if(!is.null(exclude)) {
      if(!is.list(exclude)) {stop("\nItems in \"exclude\" must be specified as a list\n")}
    }


    ##if exclude is list  
    if(is.list(exclude)) {
      
      ##determine number of elements in exclude
      nexcl <- length(exclude)
      
      ##check each formula for presence of exclude variable extracted with formula( )  
      not.include <- lapply(cand.set, FUN = formula)

      ##set up a new list with model formula
      forms <- list()
      for (i in 1:nmods) {
        form.tmp <- as.character(not.include[[i]])[3]
        if(attr(regexpr("\\+", form.tmp), "match.length")==-1) {
          forms[i] <- form.tmp
        } else {forms[i] <- strsplit(form.tmp, split=" \\+ ")}
      }

      ##additional check to see whether some variable names include "+"
      check.forms <- unlist(lapply(forms, FUN=function(i) any(attr(regexpr("\\+", i), "match.length")>0)[[1]]))
      if (any(check.forms==TRUE)) stop("\nPlease avoid \"+\" in variable names\n")

      ##additional check to determine if intercept was removed from models
      check.forms <- unlist(lapply(forms, FUN=function(i) any(attr(regexpr("\\- 1", i), "match.length")>0)[[1]]))
      if (any(check.forms==TRUE)) stop("\nModels without intercept are not supported in this version, please use alternative parameterization\n")
         
      
      ##search within formula for variables to exclude
      mod.exclude <- matrix(NA, nrow=nmods, ncol=nexcl)
      
      ##iterate over each element in exclude list
      for (var in 1:nexcl) {

        ##iterate over each formula in mod_formula list
        for (i in 1:nmods) {
          idents <- NULL
          form.excl <- forms[[i]]
                    
          ##iterate over each element of forms[[i]]
          for (j in 1:length(form.excl)) {
            idents[j] <- identical(exclude[var][[1]], gsub("(^ +)|( +$)", "", form.excl[j]))
          }
          mod.exclude[i,var] <- ifelse(any(idents==1), 1, 0)
        }    
        
      }
  
      ##determine outcome across all variables to exclude
      to.exclude <- rowSums(mod.exclude)
         
      
      ##exclude models following models from model averaging  
      include[which(to.exclude>=1)] <- 0
      
      
    }

       
       
    ##add a check to determine if include always == 0
    if (sum(include)==0) {stop("\nParameter not found in any of the candidate models\n") }
    
    new.cand.set <- cand.set[which(include==1)] #select models including a given parameter
    new.mod.name <- modnames[which(include==1)]    #update model names
    
    new_table <- aictab(cand.set=new.cand.set, modnames=new.mod.name,
                        second.ord=second.ord, nobs=nobs, sort=FALSE)  #recompute AIC table and associated measures
    new_table$Beta_est <- unlist(lapply(new.cand.set, FUN=function(i) coef(i)[paste(parm)])) #extract beta estimate for parm
    new_table$SE <- unlist(lapply(new.cand.set, FUN=function(i) sqrt(diag(vcov(i)))[paste(parm)]))
  
    ##compute model-averaged estimates, unconditional SE, and 95% CL
    if(second.ord == TRUE) {
      Modavg_beta <- sum(new_table$AICcWt*new_table$Beta_est)
         
      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Uncond_SE <- sum(new_table$AICcWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
      }

      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Uncond_SE <- sqrt(sum(new_table$AICcWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
      } 


    } else {
      Modavg_beta <- sum(new_table$AICWt*new_table$Beta_est)
      
      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Uncond_SE <- sum(new_table$AICWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
      }

      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Uncond_SE <- sqrt(sum(new_table$AICWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
      }
    
      
    }
     
    zcrit <- qnorm(p=(1-conf.level)/2, lower.tail=FALSE)
    Lower_CL <- Modavg_beta - zcrit*Uncond_SE
    Upper_CL <- Modavg_beta + zcrit*Uncond_SE
    out.modavg <- list("Parameter"=paste(parm), "Mod.avg.table" = new_table, "Mod.avg.beta" = Modavg_beta,
                       "Uncond.SE" = Uncond_SE, "Conf.level" = conf.level, "Lower.CL" = Lower_CL,
                       "Upper.CL" = Upper_CL)
    
    class(out.modavg) <- c("modavg", "list")
    return(out.modavg)
  }



##coxme
modavg.AICcoxme <- function(cand.set, parm, modnames = NULL, second.ord = TRUE,
                               nobs = NULL, uncond.se = "revised", conf.level = 0.95,
                               exclude = NULL, warn = TRUE, ...){

  ##check if named list if modnames are not supplied
  if(is.null(modnames)) {
    if(is.null(names(cand.set))) {
      modnames <- paste("Mod", 1:length(cand.set), sep = "")
      warning("\nModel names have been supplied automatically in the table\n")
    } else {
      modnames <- names(cand.set)
    }
  }
    

    
#####MODIFICATIONS BEGIN#######
  ##remove all leading and trailing white space and within parm
  parm <- gsub('[[:space:]]+', "", parm)
    
  ##reverse parm
  reversed.parm <- reverse.parm(parm)
  exclude <- reverse.exclude(exclude = exclude)
#####MODIFICATIONS END######

      
  ##extract model formula for each model in cand.set
  mod_formula <- lapply(cand.set, FUN=function(i) names(fixef(i)))
    
  nmods <- length(cand.set)
    
  ##setup matrix to indicate presence of parms in the model
  include <- matrix(NA, nrow=nmods, ncol=1)
  ##add a check for multiple instances of same variable in given model (i.e., interactions)
  include.check <- matrix(NA, nrow=nmods, ncol=1)
  
  ##iterate over each formula in mod_formula list
  for (i in 1:nmods) {
    idents <- NULL
    idents.check <- NULL
    form <- mod_formula[[i]]
    
######################################################################################################
######################################################################################################
###MODIFICATIONS BEGIN
    ##iterate over each element of formula[[i]] in list
    if(is.null(reversed.parm)) {
      for (j in 1:length(form)) {
        idents[j] <- identical(parm, form[j])
        idents.check[j] <- ifelse(attr(regexpr(parm, form[j], fixed=TRUE), "match.length")== "-1" , 0, 1)  
      }
    } else {
      for (j in 1:length(form)) {
        idents[j] <- identical(parm, form[j]) | identical(reversed.parm, form[j])
        idents.check[j] <- ifelse(attr(regexpr(parm, form[j], fixed=TRUE), "match.length")=="-1" & attr(regexpr(reversed.parm, form[j],
                                                                fixed=TRUE), "match.length")=="-1" , 0, 1)  
      }
    }
###MODIFICATIONS END
######################################################################################################
######################################################################################################
    
    include[i] <- ifelse(any(idents==1), 1, 0)
    include.check[i] <- ifelse(sum(idents.check)>1, "duplicates", "OK")
  }
  
#####################################################
  ##exclude == NULL; warn=TRUE:  warn that duplicates occur and stop
  if(is.null(exclude) && identical(warn, TRUE)) {
    ##check for duplicates in same model
    if(any(include.check == "duplicates")) {
      stop("\nSome models include more than one instance of the parameter of interest. \n",
           "This may be due to the presence of interaction/polynomial terms, or variables\n",
           "with similar names:\n",
           "\tsee \"?modavg\" for details on variable specification and \"exclude\" argument\n")
    }
    
  }

  ##exclude == NULL; warn=FALSE:  compute model-averaged beta estimate from models including variable of interest,
  ##assuming that the variable is not involved in interaction or higher order polynomial (x^2, x^3, etc...),
  ##warn that models were not excluded
  if(is.null(exclude) && identical(warn, FALSE)) {
    if(any(include.check == "duplicates")) {
      warning("\nMultiple instances of parameter of interest in given model is presumably\n",
              "not due to interaction or polynomial terms - these models will not be\n",
              "excluded from the computation of model-averaged estimate\n")
    }
    
  }
  
  ##warn if exclude is neither a list nor NULL
  if(!is.null(exclude)) {
    if(!is.list(exclude)) {stop("\nItems in \"exclude\" must be specified as a list\n")}
  }
  
  
  ##if exclude is list  
  if(is.list(exclude)) {
    
    ##determine number of elements in exclude
    nexcl <- length(exclude)
    
    ##check each formula for presence of exclude variable extracted with formula( )  
    not.include <- lapply(cand.set, FUN = function(i) formula(i)$fixed)
    
    ##set up a new list with model formula
    forms <- list( )
    for (i in 1:nmods) {
      form.tmp <- strsplit(as.character(not.include[i]), split="~")[[1]][-1]
      if(attr(regexpr("\\+", form.tmp), "match.length")==-1) {
        forms[i] <- form.tmp
      } else {forms[i] <- strsplit(form.tmp, split=" \\+ ")}
    }
    
    ##additional check to see whether some variable names include "+"
    check.forms <- unlist(lapply(forms, FUN=function(i) any(attr(regexpr("\\+", i), "match.length")>0)[[1]]))
    if (any(check.forms==TRUE)) stop("\nPlease avoid \"+\" in variable names\n")
    
    ##search within formula for variables to exclude
    mod.exclude <- matrix(NA, nrow=nmods, ncol=nexcl)
    
    ##iterate over each element in exclude list
    for (var in 1:nexcl) {
      
      ##iterate over each formula in mod_formula list
      for (i in 1:nmods) {
        idents <- NULL
        form.excl <- forms[[i]]
        
        ##iterate over each element of forms[[i]]
        for (j in 1:length(form.excl)) {
          idents[j] <- identical(exclude[var][[1]], gsub("(^ +)|( +$)", "", form.excl[j]))
        }
        mod.exclude[i,var] <- ifelse(any(idents==1), 1, 0)
      }            
    }
    
    ##determine outcome across all variables to exclude
    to.exclude <- rowSums(mod.exclude)
        
      
    ##exclude models following models from model averaging  
    include[which(to.exclude>=1)] <- 0
    
    
  }

      
    
  ##add a check to determine if include always == 0
  if (sum(include)==0) {stop("\nParameter not found in any of the candidate models\n") }
  
  new.cand.set <- cand.set[which(include==1)] #select models including a given parameter
  new.mod.name <- modnames[which(include==1)]    #update model names
  
  new_table <- aictab(cand.set=new.cand.set, modnames=new.mod.name,
                      second.ord=second.ord, nobs=nobs, sort=FALSE)  #recompute AIC table and associated measures
  new_table$Beta_est <- unlist(lapply(new.cand.set, FUN=function(i) coef(i)[paste(parm)])) #extract beta estimate for parm
  new_table$SE <- unlist(lapply(new.cand.set, FUN=function(i) extractSE(i)[paste(parm)]))
  
  ##compute model-averaged estimates, unconditional SE, and 95% CL
  if(second.ord==TRUE) {
    Modavg_beta <- sum(new_table$AICcWt*new_table$Beta_est)
        
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Uncond_SE <- sum(new_table$AICcWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
    }

    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Uncond_SE <- sqrt(sum(new_table$AICcWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
    } 


  } else {
    Modavg_beta <- sum(new_table$AICWt*new_table$Beta_est)
    
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Uncond_SE <- sum(new_table$AICWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
    }

    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Uncond_SE <- sqrt(sum(new_table$AICWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
    }

                      
  }
  
  zcrit <- qnorm(p = (1 - conf.level)/2, lower.tail = FALSE)
  Lower_CL <- Modavg_beta - zcrit*Uncond_SE
  Upper_CL <- Modavg_beta + zcrit*Uncond_SE
  out.modavg <- list("Parameter" = paste(parm), "Mod.avg.table" = new_table, "Mod.avg.beta" = Modavg_beta,
                     "Uncond.SE" = Uncond_SE, "Conf.level" = conf.level, "Lower.CL" = Lower_CL,
                     "Upper.CL" = Upper_CL)
  
  class(out.modavg) <- c("modavg", "list")
  return(out.modavg)
}



##coxph and clogit
modavg.AICcoxph <- function(cand.set, parm, modnames = NULL, second.ord = TRUE,
                               nobs = NULL, uncond.se = "revised", conf.level = 0.95,
                               exclude = NULL, warn = TRUE, ...){

  ##check if named list if modnames are not supplied
  if(is.null(modnames)) {
    if(is.null(names(cand.set))) {
      modnames <- paste("Mod", 1:length(cand.set), sep = "")
      warning("\nModel names have been supplied automatically in the table\n")
    } else {
      modnames <- names(cand.set)
    }
  }
    

#####MODIFICATIONS BEGIN#######
    ##remove all leading and trailing white space and within parm
    parm <- gsub('[[:space:]]+', "", parm)
    
    ##reverse parm
    reversed.parm <- reverse.parm(parm)
    exclude <- reverse.exclude(exclude = exclude)
#####MODIFICATIONS END######

      
    ##extract model formula for each model in cand.set
    mod_formula <- lapply(cand.set, FUN=function(i) rownames(summary(i)$coefficients))
    
    nmods <- length(cand.set)
    
    ##setup matrix to indicate presence of parms in the model
    include <- matrix(NA, nrow=nmods, ncol=1)
    ##add a check for multiple instances of same variable in given model (i.e., interactions)
    include.check <- matrix(NA, nrow=nmods, ncol=1)

    ##iterate over each formula in mod_formula list
    for (i in 1:nmods) {
      idents <- NULL
      idents.check <- NULL
      form <- mod_formula[[i]]

######################################################################################################
######################################################################################################
###MODIFICATIONS BEGIN
      ##iterate over each element of formula[[i]] in list
      if(is.null(reversed.parm)) {
        for (j in 1:length(form)) {
          idents[j] <- identical(parm, form[j])
          idents.check[j] <- ifelse(attr(regexpr(parm, form[j], fixed=TRUE), "match.length")== "-1" , 0, 1)  
        }
      } else {
        for (j in 1:length(form)) {
          idents[j] <- identical(parm, form[j]) | identical(reversed.parm, form[j])
          idents.check[j] <- ifelse(attr(regexpr(parm, form[j], fixed=TRUE), "match.length")=="-1" & attr(regexpr(reversed.parm, form[j],
                                                                  fixed=TRUE), "match.length")=="-1" , 0, 1)  
        }
      }
###MODIFICATIONS END
######################################################################################################
######################################################################################################
  
      include[i] <- ifelse(any(idents==1), 1, 0)
      include.check[i] <- ifelse(sum(idents.check)>1, "duplicates", "OK")
    }

#####################################################
    ##exclude == NULL; warn=TRUE:  warn that duplicates occur and stop
    if(is.null(exclude) && identical(warn, TRUE)) {
      ##check for duplicates in same model
      if(any(include.check == "duplicates")) {
        stop("\nSome models include more than one instance of the parameter of interest. \n",
             "This may be due to the presence of interaction/polynomial terms, or variables\n",
             "with similar names:\n",
             "\tsee \"?modavg\" for details on variable specification and \"exclude\" argument\n")
      }
        
    }

    ##exclude == NULL; warn=FALSE:  compute model-averaged beta estimate from models including variable of interest,
    ##assuming that the variable is not involved in interaction or higher order polynomial (x^2, x^3, etc...),
    ##warn that models were not excluded
    if(is.null(exclude) && identical(warn, FALSE)) {
      if(any(include.check == "duplicates")) {
        warning("\nMultiple instances of parameter of interest in given model is presumably\n",
                "not due to interaction or polynomial terms - these models will not be\n",
                "excluded from the computation of model-averaged estimate\n")
      }
      
    }
    
    ##warn if exclude is neither a list nor NULL
    if(!is.null(exclude)) {
      if(!is.list(exclude)) {stop("\nItems in \"exclude\" must be specified as a list\n")}
    }


    ##if exclude is list  
    if(is.list(exclude)) {

      ##determine number of elements in exclude
      nexcl <- length(exclude)

      ##check each formula for presence of exclude variable extracted with formula( )  
      not.include <- lapply(cand.set, FUN=formula)

      ##set up a new list with model formula
      forms <- list( )
      for (i in 1:nmods) {
        form.tmp <- strsplit(as.character(not.include[i]), split="~")[[1]][-1]
        if(attr(regexpr("\\+", form.tmp), "match.length")==-1) {
          forms[i] <- form.tmp
        } else {forms[i] <- strsplit(form.tmp, split=" \\+ ")}
      }

      ##additional check to see whether some variable names include "+"
      check.forms <- unlist(lapply(forms, FUN=function(i) any(attr(regexpr("\\+", i), "match.length")>0)[[1]]))
      if (any(check.forms==TRUE)) stop("\nPlease avoid \"+\" in variable names\n")
      
      ##search within formula for variables to exclude
      mod.exclude <- matrix(NA, nrow=nmods, ncol=nexcl)

      ##iterate over each element in exclude list
      for (var in 1:nexcl) {

        ##iterate over each formula in mod_formula list
        for (i in 1:nmods) {
          idents <- NULL
          form.excl <- forms[[i]]
                    
          ##iterate over each element of forms[[i]]
          for (j in 1:length(form.excl)) {
            idents[j] <- identical(exclude[var][[1]], gsub("(^ +)|( +$)", "", form.excl[j]))
          }
          mod.exclude[i,var] <- ifelse(any(idents==1), 1, 0)
        }            
      }
  
      ##determine outcome across all variables to exclude
      to.exclude <- rowSums(mod.exclude)
        
      
      ##exclude models following models from model averaging  
      include[which(to.exclude>=1)] <- 0
      
        
    }

      
    
  ##add a check to determine if include always == 0
  if (sum(include)==0) {stop("\nParameter not found in any of the candidate models\n") }
    
  new.cand.set <- cand.set[which(include==1)] #select models including a given parameter
  new.mod.name <- modnames[which(include==1)]    #update model names

  new_table <- aictab(cand.set=new.cand.set, modnames=new.mod.name,
                      second.ord=second.ord, nobs=nobs, sort=FALSE)  #recompute AIC table and associated measures
  new_table$Beta_est <- unlist(lapply(new.cand.set, FUN=function(i) coef(i)[paste(parm)])) #extract beta estimate for parm
  new_table$SE <- unlist(lapply(new.cand.set, FUN=function(i) sqrt(diag(vcov(i)))[paste(parm)]))
  
  ##compute model-averaged estimates, unconditional SE, and 95% CL
  if(second.ord==TRUE) {
    Modavg_beta <- sum(new_table$AICcWt*new_table$Beta_est)
    
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Uncond_SE <- sum(new_table$AICcWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
    }
    
    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Uncond_SE <- sqrt(sum(new_table$AICcWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
    } 


  } else {
    Modavg_beta <- sum(new_table$AICWt*new_table$Beta_est)
    
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Uncond_SE <- sum(new_table$AICWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
    }

    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Uncond_SE <- sqrt(sum(new_table$AICWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
    }
    
                      
  }
  
  zcrit <- qnorm(p = (1 - conf.level)/2, lower.tail = FALSE)
  Lower_CL <- Modavg_beta - zcrit*Uncond_SE
  Upper_CL <- Modavg_beta + zcrit*Uncond_SE
  out.modavg <- list("Parameter" = paste(parm), "Mod.avg.table" = new_table, "Mod.avg.beta" = Modavg_beta,
                     "Uncond.SE" = Uncond_SE, "Conf.level" = conf.level, "Lower.CL" = Lower_CL,
                     "Upper.CL" = Upper_CL)
  
  class(out.modavg) <- c("modavg", "list")
  return(out.modavg)
}



##glm
modavg.AICglm.lm <-
  function(cand.set, parm, modnames = NULL, second.ord = TRUE,
           nobs = NULL,  uncond.se = "revised", conf.level = 0.95, 
           exclude = NULL, warn = TRUE, c.hat = 1, gamdisp = NULL, ...){

    ##check if named list if modnames are not supplied
    if(is.null(modnames)) {
      if(is.null(names(cand.set))) {
        modnames <- paste("Mod", 1:length(cand.set), sep = "")
        warning("\nModel names have been supplied automatically in the table\n")
      } else {
        modnames <- names(cand.set)
      }
    }

    ##check that link function is the same for all models
    check.link <- unlist(lapply(X = cand.set, FUN = function(i) i$family$link))
    unique.link <- unique(x = check.link)
    if(length(unique.link) > 1) stop("\nIt is not appropriate to compute a model-averaged beta estimate\n",
                                     "from models using different link functions\n")
  


    ##check family of glm to avoid problems when requesting predictions with argument 'dispersion'
    fam.type <- unlist(lapply(cand.set, FUN=function(i) family(i)$family))
    fam.unique <- unique(fam.type)
    if(identical(fam.unique, "gaussian")) {disp <- NULL} else{disp <- 1}
    ##poisson and binomial defaults to 1 (no separate parameter for variance)

    ##for negative binomial - reset to NULL
    if(any(regexpr("Negative Binomial", fam.type) != -1)) {
      disp <- NULL
      ##check for mixture of negative binomial and other
      ##number of models with negative binomial
      negbin.num <- sum(regexpr("Negative Binomial", fam.type) != -1)
      if(negbin.num < length(fam.type)) {
        stop("\nFunction does not support mixture of negative binomial with other distributions in model set\n")
      }
    }
    ##gamma is treated separately

#####MODIFICATIONS BEGIN#######
    ##remove all leading and trailing white space and within parm
    parm <- gsub('[[:space:]]+', "", parm)
     
    ##reverse parm
    reversed.parm <- reverse.parm(parm)
    exclude <- reverse.exclude(exclude = exclude)
#####MODIFICATIONS END######

    ##extract model formula for each model in cand.set
    mod_formula <- lapply(cand.set, FUN = function(i) rownames(summary(i)$coefficients)) #extract model formula for each model in cand.set
    
    nmods <- length(cand.set)
  
    ##setup matrix to indicate presence of parm in the model
    include <- matrix(NA, nrow=nmods, ncol=1)
    ##add a check for multiple instances of same variable in given model (i.e., interactions)
    include.check <- matrix(NA, nrow=nmods, ncol=1)

    ##iterate over each formula in mod_formula list
    for (i in 1:nmods) {
      idents <- NULL
      idents.check <- NULL
      form <- mod_formula[[i]]

######################################################################################################
######################################################################################################
###MODIFICATIONS BEGIN
      ##iterate over each element of formula[[i]] in list
      if(is.null(reversed.parm)) {
        for (j in 1:length(form)) {
          idents[j] <- identical(parm, form[j])
          idents.check[j] <- ifelse(attr(regexpr(parm, form[j], fixed=TRUE), "match.length")== "-1" , 0, 1)  
        }
      } else {
        for (j in 1:length(form)) {
          idents[j] <- identical(parm, form[j]) | identical(reversed.parm, form[j])
          idents.check[j] <- ifelse(attr(regexpr(parm, form[j], fixed=TRUE), "match.length")=="-1" & attr(regexpr(reversed.parm, form[j],
                                                                  fixed=TRUE), "match.length")=="-1" , 0, 1)  
        }
      }
###MODIFICATIONS END
######################################################################################################
######################################################################################################


      
      include[i] <- ifelse(any(idents==1), 1, 0)
      include.check[i] <- ifelse(sum(idents.check)>1, "duplicates", "OK")
    }

#####################################################
    ##exclude == NULL; warn=TRUE:  warn that duplicates occur and stop
    if(is.null(exclude) && identical(warn, TRUE)) {
      if(any(include.check == "duplicates")) {
        stop("\nSome models possibly include more than one instance of the parameter of interest.\n",
             "This may be due to the presence of interaction/polynomial terms, or variables\n",
             "with similar names:\n",
             "\tsee \"?modavg\" for details on variable specification and \"exclude\" argument\n")
      }

    }

    ##exclude == NULL; warn=FALSE:  compute model-averaged beta estimate from models including variable of interest,
    ##assuming that the variable is not involved in interaction or higher order polynomial (x^2, x^3, etc...),
    ##warn that models were not excluded
    if(is.null(exclude) && identical(warn, FALSE)) {
      if(any(include.check == "duplicates")) {
        warning("Multiple instances of parameter of interest in given model is presumably\n",
                "not due to interaction or polynomial terms - these models will not be\n",
                "excluded from the computation of model-averaged estimate\n")
      }
      
    }



    ##warn if exclude is neither a list nor NULL
    if(!is.null(exclude)) {
      if(!is.list(exclude)) {stop("\nItems in \"exclude\" must be specified as a list\n")}
    }

  
    ##if exclude is list  
    if(is.list(exclude)) {

      ##determine number of elements in exclude
      nexcl <- length(exclude)

      ##check each formula for presence of exclude variable extracted with formula( )  
      not.include <- lapply(cand.set, FUN=formula)

      ##set up a new list with model formula
      forms <- list()
      for (i in 1:nmods) {
        form.tmp <- strsplit(as.character(not.include[i]), split="~")[[1]][-1]
        if(attr(regexpr("\\+", form.tmp), "match.length")==-1) {  #this line causes problems if intercept is removed from model
          forms[i] <- form.tmp
        } else {forms[i] <- strsplit(form.tmp, split=" \\+ ")}    #this line causes problems if intercept is removed from model
      }

      ##additional check to see whether some variable names include "+"
      check.forms <- unlist(lapply(forms, FUN=function(i) any(attr(regexpr("\\+", i), "match.length")>0)[[1]]))
      if (any(check.forms==TRUE)) stop("\nPlease avoid \"+\" in variable names\n")

      ##additional check to determine if intercept was removed from models
      check.forms <- unlist(lapply(forms, FUN=function(i) any(attr(regexpr("\\- 1", i), "match.length")>0)[[1]]))
      if (any(check.forms==TRUE)) stop("\nModels without intercept are not supported in this version, please use alternative parameterization\n")

 
      ##search within formula for variables to exclude
      mod.exclude <- matrix(NA, nrow=nmods, ncol=nexcl)

      ##iterate over each element in exclude list
      for (var in 1:nexcl) {

        ##iterate over each formula in mod_formula list
        for (i in 1:nmods) {
          idents <- NULL
          form.excl <- forms[[i]]
                    
          ##iterate over each element of forms[[i]]
          for (j in 1:length(form.excl)) {
            idents[j] <- identical(exclude[var][[1]], gsub("(^ +)|( +$)", "", form.excl[j]))
          }
          mod.exclude[i,var] <- ifelse(any(idents==1), 1, 0)
        }    
        
      }
  
      ##determine outcome across all variables to exclude
      to.exclude <- rowSums(mod.exclude)
  
  
      ##exclude models following models from model averaging  
      include[which(to.exclude>=1)] <- 0
      
      
    }



    ##add a check to determine if include always == 0
    if (sum(include)==0) {stop("\nParameter not found in any of the candidate models\n") }

    new.cand.set <- cand.set[which(include==1)] #select models including a given parameter
    new.mod.name <- modnames[which(include==1)]    #update model names
    ##

    new_table <- aictab(cand.set = new.cand.set, modnames = new.mod.name,
                        second.ord = second.ord, nobs = nobs, sort = FALSE,
                        c.hat = c.hat)  #recompute AIC table and associated measures
    new_table$Beta_est <- unlist(lapply(new.cand.set, FUN = function(i) coef(i)[paste(parm)])) #extract beta estimate for parm
    new_table$SE <- unlist(lapply(new.cand.set, FUN = function(i) sqrt(diag(vcov(i, dispersion = disp)))[paste(parm)]))

   
    ##if c-hat is estimated adjust the SE's by multiplying with sqrt of c-hat
    if(c.hat > 1) {
      new_table$SE <- new_table$SE*sqrt(c.hat)
    } 

    gam1<-unlist(lapply(new.cand.set, FUN=function(i) family(i)$family[1]=="Gamma")) #check for gamma regression models
    ##correct SE's for estimates of gamma regressions
    if(any(gam1)==TRUE)  {
      ##check for specification of gamdisp argument
      if(is.null(gamdisp)) stop("\nYou must specify a gamma dispersion parameter with gamma generalized linear models\n")
      new_table$SE <- unlist(lapply(new.cand.set,
                                  FUN=function(i) sqrt(diag(vcov(i, dispersion=gamdisp)))[paste(parm)]))
    } 



    ##AICc
    ##compute model-averaged estimates, unconditional SE, and 95% CL
    if(c.hat == 1 && second.ord == TRUE) {
      Modavg_beta <- sum(new_table$AICcWt*new_table$Beta_est)

      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Uncond_SE <- sum(new_table$AICcWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
      }
      
      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Uncond_SE <- sqrt(sum(new_table$AICcWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
      }
    }

    
    
    ##QAICc
    if(c.hat > 1 && second.ord == TRUE) {
      Modavg_beta <- sum(new_table$QAICcWt*new_table$Beta_est)
      
      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {      
        Uncond_SE <- sum(new_table$QAICcWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
      }
      
      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Uncond_SE <- sqrt(sum(new_table$QAICcWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
      }
    }     



    ##AIC
    if(c.hat == 1 && second.ord == FALSE) {
      Modavg_beta <- sum(new_table$AICWt*new_table$Beta_est)
      
      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Uncond_SE <- sum(new_table$AICWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
      }
      
      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Uncond_SE <- sqrt(sum(new_table$AICWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
      }
    }



    ##QAIC
    if(c.hat > 1 && second.ord == FALSE) {
      Modavg_beta <- sum(new_table$QAICWt*new_table$Beta_est)

      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Uncond_SE <- sum(new_table$QAICWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
      }

      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Uncond_SE <- sqrt(sum(new_table$QAICWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
      }  
    }     

    zcrit <- qnorm(p=(1-conf.level)/2, lower.tail=FALSE)
    Lower_CL <- Modavg_beta-zcrit*Uncond_SE
    Upper_CL <- Modavg_beta+zcrit*Uncond_SE
    out.modavg <- list("Parameter"=paste(parm), "Mod.avg.table" = new_table, "Mod.avg.beta" = Modavg_beta,
                       "Uncond.SE" = Uncond_SE, "Conf.level" = conf.level, "Lower.CL"= Lower_CL, "Upper.CL" = Upper_CL)
  
    
    class(out.modavg) <- c("modavg", "list")
    return(out.modavg)

  }



##gls
modavg.AICgls <-
  function(cand.set, parm, modnames = NULL, second.ord = TRUE,
           nobs = NULL, uncond.se = "revised", conf.level = 0.95,
           exclude = NULL, warn = TRUE, ...){

    ##check if named list if modnames are not supplied
    if(is.null(modnames)) {
      if(is.null(names(cand.set))) {
        modnames <- paste("Mod", 1:length(cand.set), sep = "")
        warning("\nModel names have been supplied automatically in the table\n")
      } else {
        modnames <- names(cand.set)
      }
    }

    
#####MODIFICATIONS BEGIN#######
    ##remove all leading and trailing white space and within parm
    parm <- gsub('[[:space:]]+', "", parm)
    
    ##reverse parm
    reversed.parm <- reverse.parm(parm)
    exclude <- reverse.exclude(exclude = exclude)
#####MODIFICATIONS END######

    
    ##extract model formula for each model in cand.set
    mod_formula <- lapply(cand.set, FUN=function(i) labels(summary(i)$coefficients))
    
    nmods <- length(cand.set)
  
    ##setup matrix to indicate presence of parms in the model
    include <- matrix(NA, nrow=nmods, ncol=1)
    ##add a check for multiple instances of same variable in given model (i.e., interactions)
    include.check <- matrix(NA, nrow=nmods, ncol=1)
    
    ##iterate over each formula in mod_formula list
    for (i in 1:nmods) {
      idents <- NULL
      idents.check <- NULL
      form <- mod_formula[[i]]
      
######################################################################################################
######################################################################################################
###MODIFICATIONS BEGIN
      ##iterate over each element of formula[[i]] in list
      if(is.null(reversed.parm)) {
        for (j in 1:length(form)) {
          idents[j] <- identical(parm, form[j])
          idents.check[j] <- ifelse(attr(regexpr(parm, form[j], fixed=TRUE), "match.length")== "-1" , 0, 1)  
        }
      } else {
        for (j in 1:length(form)) {
          idents[j] <- identical(parm, form[j]) | identical(reversed.parm, form[j])
          idents.check[j] <- ifelse(attr(regexpr(parm, form[j], fixed=TRUE), "match.length")=="-1" & attr(regexpr(reversed.parm, form[j],
                                                                  fixed=TRUE), "match.length")=="-1" , 0, 1)  
        }
      }
###MODIFICATIONS END
######################################################################################################
######################################################################################################
  
      include[i] <- ifelse(any(idents==1), 1, 0)
      include.check[i] <- ifelse(sum(idents.check)>1, "duplicates", "OK")
    }

#####################################################
    ##exclude == NULL; warn=TRUE:  warn that duplicates occur and stop
    if(is.null(exclude) && identical(warn, TRUE)) {
      ##check for duplicates in same model
      if(any(include.check == "duplicates")) {
        stop("\nSome models include more than one instance of the parameter of interest. \n",
             "This may be due to the presence of interaction/polynomial terms, or variables\n",
             "with similar names:\n",
             "\tsee \"?modavg\" for details on variable specification and \"exclude\" argument\n")
      }
      
    }

    ##exclude == NULL; warn=FALSE:  compute model-averaged beta estimate from models including variable of interest,
    ##assuming that the variable is not involved in interaction or higher order polynomial (x^2, x^3, etc...),
    ##warn that models were not excluded
    if(is.null(exclude) && identical(warn, FALSE)) {
      if(any(include.check == "duplicates")) {
        warning("\nMultiple instances of parameter of interest in given model is presumably\n",
                "not due to interaction or polynomial terms - these models will not be\n",
                "excluded from the computation of model-averaged estimate\n")
      }
      
    }
    
    ##warn if exclude is neither a list nor NULL
    if(!is.null(exclude)) {
      if(!is.list(exclude)) {stop("\nItems in \"exclude\" must be specified as a list\n")}
    }

    
    ##if exclude is list  
    if(is.list(exclude)) {

      ##determine number of elements in exclude
      nexcl <- length(exclude)

      ##check each formula for presence of exclude variable extracted with formula( )  
      not.include <- lapply(cand.set, FUN=formula)
      
      ##set up a new list with model formula
      forms <- list()
      for (i in 1:nmods) {
        form.tmp <- strsplit(as.character(not.include[i]), split="~")[[1]][-1]
        if(attr(regexpr("\\+", form.tmp), "match.length")==-1) {
          forms[i] <- form.tmp
        } else {forms[i] <- strsplit(form.tmp, split=" \\+ ")}
      }
      
      ##additional check to see whether some variable names include "+"
      check.forms <- unlist(lapply(forms, FUN=function(i) any(attr(regexpr("\\+", i), "match.length")>0)[[1]]))
      if (any(check.forms==TRUE)) stop("\nPlease avoid \"+\" in variable names\n")
      
      ##additional check to determine if intercept was removed from models
      check.forms <- unlist(lapply(forms, FUN=function(i) any(attr(regexpr("\\- 1", i), "match.length")>0)[[1]]))
      if (any(check.forms==TRUE)) stop("\nModels without intercept are not supported in this version, please use alternative parameterization\n")
      
      
      ##search within formula for variables to exclude
      mod.exclude <- matrix(NA, nrow=nmods, ncol=nexcl)
      
      ##iterate over each element in exclude list
      for (var in 1:nexcl) {
        
        ##iterate over each formula in mod_formula list
        for (i in 1:nmods) {
          idents <- NULL
          form.excl <- forms[[i]]
                    
          ##iterate over each element of forms[[i]]
          for (j in 1:length(form.excl)) {
            idents[j] <- identical(exclude[var][[1]], gsub("(^ +)|( +$)", "", form.excl[j]))
          }
          mod.exclude[i,var] <- ifelse(any(idents==1), 1, 0)
        }            
      }
  
      ##determine outcome across all variables to exclude
      to.exclude <- rowSums(mod.exclude)
      
      
      ##exclude models following models from model averaging  
      include[which(to.exclude>=1)] <- 0
      
      
    }
    

 
    ##add a check to determine if include always == 0
    if (sum(include)==0) {stop("\nParameter not found in any of the candidate models\n") }
    
    new.cand.set <- cand.set[which(include==1)] #select models including a given parameter
    new.mod.name <- modnames[which(include==1)]    #update model names
    
    new_table <- aictab(cand.set=new.cand.set, modnames=new.mod.name, second.ord=second.ord,
                        nobs=nobs, sort=FALSE)  #recompute AIC table and associated measures
    new_table$Beta_est <- unlist(lapply(new.cand.set, FUN=function(i) coef(i)[paste(parm)])) #extract beta estimate for parm
    new_table$SE <- unlist(lapply(new.cand.set, FUN=function(i) sqrt(diag(vcov(i)))[paste(parm)]))
  
    ##compute model-averaged estimates, unconditional SE, and 95% CL
    if(second.ord==TRUE) {
      Modavg_beta <- sum(new_table$AICcWt*new_table$Beta_est)
      
      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Uncond_SE <- sum(new_table$AICcWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
      }

      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Uncond_SE <- sqrt(sum(new_table$AICcWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
      } 


    } else {
      Modavg_beta <- sum(new_table$AICWt*new_table$Beta_est)

      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Uncond_SE <- sum(new_table$AICWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
      }

      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Uncond_SE <- sqrt(sum(new_table$AICWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
      }

                      
    }
  
    zcrit <- qnorm(p=(1-conf.level)/2, lower.tail=FALSE)
    Lower_CL <- Modavg_beta-zcrit*Uncond_SE
    Upper_CL <- Modavg_beta+zcrit*Uncond_SE
    out.modavg <- list("Parameter"=paste(parm), "Mod.avg.table" = new_table, "Mod.avg.beta" = Modavg_beta,
                       "Uncond.SE" = Uncond_SE, "Conf.level" = conf.level, "Lower.CL" = Lower_CL,
                       "Upper.CL" = Upper_CL)
  
    class(out.modavg) <- c("modavg", "list")
    return(out.modavg)
  }



##hurdle
modavg.AIChurdle <-
function(cand.set, parm, modnames = NULL, second.ord = TRUE,
         nobs = NULL, uncond.se = "revised", conf.level = 0.95,
         exclude = NULL, warn = TRUE, ...){

  ##check if named list if modnames are not supplied
  if(is.null(modnames)) {
    if(is.null(names(cand.set))) {
      modnames <- paste("Mod", 1:length(cand.set), sep = "")
      warning("\nModel names have been supplied automatically in the table\n")
    } else {
      modnames <- names(cand.set)
    }
  }
    
  
  ##check that link function is the same for all models
  check.link <- unlist(lapply(X = cand.set, FUN=function(i) i$link))
  unique.link <- unique(x=check.link)
  if(length(unique.link) > 1) stop("\nIt is not appropriate to compute a model-averaged beta estimate\n",
"from models using different link functions\n")
  
#####MODIFICATIONS BEGIN#######
    ##remove all leading and trailing white space and within parm
    parm <- gsub('[[:space:]]+', "", parm)
     
    ##reverse parm
    reversed.parm <- reverse.parm(parm)
    exclude <- reverse.exclude(exclude = exclude)
#####MODIFICATIONS END######

    ##extract model formula for each model in cand.set
    mod_formula <- lapply(cand.set, FUN=function(i) labels(coefficients(i))) #extract model formula for each model in cand.set

    nmods <- length(cand.set)
  
    ##setup matrix to indicate presence of parm in the model
    include <- matrix(NA, nrow=nmods, ncol=1)
    ##add a check for multiple instances of same variable in given model (i.e., interactions)
    include.check <- matrix(NA, nrow=nmods, ncol=1)

    ##iterate over each formula in mod_formula list
    for (i in 1:nmods) {
      idents <- NULL
      idents.check <- NULL
      form <- mod_formula[[i]]

######################################################################################################
######################################################################################################
###MODIFICATIONS BEGIN
      ##iterate over each element of formula[[i]] in list
      if(is.null(reversed.parm)) {
        for (j in 1:length(form)) {
          idents[j] <- identical(parm, form[j])
          idents.check[j] <- ifelse(attr(regexpr(parm, form[j], fixed=TRUE), "match.length")== "-1" , 0, 1)  
        }
      } else {
        for (j in 1:length(form)) {
          idents[j] <- identical(parm, form[j]) | identical(reversed.parm, form[j])
          idents.check[j] <- ifelse(attr(regexpr(parm, form[j], fixed=TRUE), "match.length")=="-1" & attr(regexpr(reversed.parm, form[j],
                                                                  fixed=TRUE), "match.length")=="-1" , 0, 1)  
        }
      }
###MODIFICATIONS END
######################################################################################################
######################################################################################################


      
      include[i] <- ifelse(any(idents==1), 1, 0)
      include.check[i] <- ifelse(sum(idents.check)>1, "duplicates", "OK")
    }

#####################################################
    ##exclude == NULL; warn=TRUE:  warn that duplicates occur and stop
    if(is.null(exclude) && identical(warn, TRUE)) {
      if(any(include.check == "duplicates")) {
        stop("\nSome models possibly include more than one instance of the parameter of interest.\n",
             "This may be due to the presence of interaction/polynomial terms, or variables\n",
             "with similar names:\n",
             "\tsee \"?modavg\" for details on variable specification and \"exclude\" argument\n")
      }

    }

    ##exclude == NULL; warn=FALSE:  compute model-averaged beta estimate from models including variable of interest,
    ##assuming that the variable is not involved in interaction or higher order polynomial (x^2, x^3, etc...),
    ##warn that models were not excluded
    if(is.null(exclude) && identical(warn, FALSE)) {
      if(any(include.check == "duplicates")) {
        warning("Multiple instances of parameter of interest in given model is presumably\n",
                "not due to interaction or polynomial terms - these models will not be\n",
                "excluded from the computation of model-averaged estimate\n")
      }
      
    }



    ##warn if exclude is neither a list nor NULL
    if(!is.null(exclude)) {
      if(!is.list(exclude)) {stop("\nItems in \"exclude\" must be specified as a list\n")}
    }

  
    ##if exclude is list  
    if(is.list(exclude)) {

      ##determine number of elements in exclude
      nexcl <- length(exclude)

      ##check each formula for presence of exclude variable extracted with formula( )  
      not.include <- lapply(cand.set, FUN=formula)

      ##set up a new list with model formula
      forms <- list()
      for (i in 1:nmods) {
        form.tmp <- strsplit(as.character(not.include[i]), split="~")[[1]][-1]
        if(attr(regexpr("\\+", form.tmp), "match.length")==-1) {  #this line causes problems if intercept is removed from model
          forms[i] <- form.tmp
        } else {forms[i] <- strsplit(form.tmp, split=" \\+ ")}    #this line causes problems if intercept is removed from model
      }

      ##additional check to see whether some variable names include "+"
      check.forms <- unlist(lapply(forms, FUN=function(i) any(attr(regexpr("\\+", i), "match.length")>0)[[1]]))
      if (any(check.forms==TRUE)) stop("\nPlease avoid \"+\" in variable names\n")

      ##additional check to determine if intercept was removed from models
      check.forms <- unlist(lapply(forms, FUN=function(i) any(attr(regexpr("\\- 1", i), "match.length")>0)[[1]]))
      if (any(check.forms==TRUE)) stop("\nModels without intercept are not supported in this version, please use alternative parameterization\n")

 
      ##search within formula for variables to exclude
      mod.exclude <- matrix(NA, nrow=nmods, ncol=nexcl)

      ##iterate over each element in exclude list
      for (var in 1:nexcl) {

        ##iterate over each formula in mod_formula list
        for (i in 1:nmods) {
          idents <- NULL
          form.excl <- forms[[i]]
                    
          ##iterate over each element of forms[[i]]
          for (j in 1:length(form.excl)) {
            idents[j] <- identical(exclude[var][[1]], gsub("(^ +)|( +$)", "", form.excl[j]))
          }
          mod.exclude[i,var] <- ifelse(any(idents==1), 1, 0)
        }    
        
      }
  
      ##determine outcome across all variables to exclude
      to.exclude <- rowSums(mod.exclude)
  
  
      ##exclude models following models from model averaging  
      include[which(to.exclude>=1)] <- 0
      
      
    }



    ##add a check to determine if include always == 0
    if (sum(include)==0) {stop("\nParameter not found in any of the candidate models\n") }

    new.cand.set <- cand.set[which(include==1)] #select models including a given parameter
    new.mod.name <- modnames[which(include==1)]    #update model names
    ##

    new_table <- aictab(cand.set = new.cand.set, modnames = new.mod.name, 
                        second.ord = second.ord, nobs = nobs, sort = FALSE)  #recompute AIC table and associated measures
    new_table$Beta_est <- unlist(lapply(new.cand.set, FUN=function(i) coefficients(i)[paste(parm)])) #extract beta estimate for parm
    new_table$SE <- unlist(lapply(new.cand.set, FUN=function(i) sqrt(diag(vcov(i)))[paste(parm)]))

   
    ##AICc
    ##compute model-averaged estimates, unconditional SE, and 95% CL
    if(second.ord == TRUE) {
      Modavg_beta <- sum(new_table$AICcWt*new_table$Beta_est)

      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Uncond_SE <- sum(new_table$AICcWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
      }
      
      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Uncond_SE <- sqrt(sum(new_table$AICcWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
      }
    }

    
    
    ##AIC
    if(second.ord == FALSE) {
      Modavg_beta <- sum(new_table$AICWt*new_table$Beta_est)
      
      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Uncond_SE <- sum(new_table$AICWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
      }
      
      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Uncond_SE <- sqrt(sum(new_table$AICWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
      }
    }


    zcrit <- qnorm(p=(1-conf.level)/2, lower.tail=FALSE)
    Lower_CL <- Modavg_beta-zcrit*Uncond_SE
    Upper_CL <- Modavg_beta+zcrit*Uncond_SE
    out.modavg <- list("Parameter"=paste(parm), "Mod.avg.table" = new_table, "Mod.avg.beta" = Modavg_beta,
                       "Uncond.SE" = Uncond_SE, "Conf.level" = conf.level, "Lower.CL"= Lower_CL, "Upper.CL" = Upper_CL)
    
  class(out.modavg) <- c("modavg", "list")
  return(out.modavg)

}



##lm
modavg.AIClm <-
  function(cand.set, parm, modnames = NULL, second.ord = TRUE,
           nobs = NULL,  uncond.se = "revised", conf.level = 0.95, 
           exclude = NULL, warn = TRUE, ...){

    ##check if named list if modnames are not supplied
    if(is.null(modnames)) {
      if(is.null(names(cand.set))) {
        modnames <- paste("Mod", 1:length(cand.set), sep = "")
        warning("\nModel names have been supplied automatically in the table\n")
      } else {
        modnames <- names(cand.set)
      }
    }


#####MODIFICATIONS BEGIN#######
    ##remove all leading and trailing white space and within parm
    parm <- gsub('[[:space:]]+', "", parm)
     
    ##reverse parm
    reversed.parm <- reverse.parm(parm)
    exclude <- reverse.exclude(exclude = exclude)
#####MODIFICATIONS END######

    ##extract model formula for each model in cand.set
    mod_formula <- lapply(cand.set, FUN = function(i) rownames(summary(i)$coefficients)) #extract model formula for each model in cand.set
    
    nmods <- length(cand.set)
  
    ##setup matrix to indicate presence of parm in the model
    include <- matrix(NA, nrow=nmods, ncol=1)
    ##add a check for multiple instances of same variable in given model (i.e., interactions)
    include.check <- matrix(NA, nrow=nmods, ncol=1)

    ##iterate over each formula in mod_formula list
    for (i in 1:nmods) {
      idents <- NULL
      idents.check <- NULL
      form <- mod_formula[[i]]

######################################################################################################
######################################################################################################
###MODIFICATIONS BEGIN
      ##iterate over each element of formula[[i]] in list
      if(is.null(reversed.parm)) {
        for (j in 1:length(form)) {
          idents[j] <- identical(parm, form[j])
          idents.check[j] <- ifelse(attr(regexpr(parm, form[j], fixed=TRUE), "match.length")== "-1" , 0, 1)  
        }
      } else {
        for (j in 1:length(form)) {
          idents[j] <- identical(parm, form[j]) | identical(reversed.parm, form[j])
          idents.check[j] <- ifelse(attr(regexpr(parm, form[j], fixed=TRUE), "match.length")=="-1" & attr(regexpr(reversed.parm, form[j],
                                                                  fixed=TRUE), "match.length")=="-1" , 0, 1)  
        }
      }
###MODIFICATIONS END
######################################################################################################
######################################################################################################


      include[i] <- ifelse(any(idents==1), 1, 0)
      include.check[i] <- ifelse(sum(idents.check)>1, "duplicates", "OK")
    }

#####################################################
    ##exclude == NULL; warn=TRUE:  warn that duplicates occur and stop
    if(is.null(exclude) && identical(warn, TRUE)) {
      if(any(include.check == "duplicates")) {
        stop("\nSome models possibly include more than one instance of the parameter of interest.\n",
             "This may be due to the presence of interaction/polynomial terms, or variables\n",
             "with similar names:\n",
             "\tsee \"?modavg\" for details on variable specification and \"exclude\" argument\n")
      }

    }

    ##exclude == NULL; warn=FALSE:  compute model-averaged beta estimate from models including variable of interest,
    ##assuming that the variable is not involved in interaction or higher order polynomial (x^2, x^3, etc...),
    ##warn that models were not excluded
    if(is.null(exclude) && identical(warn, FALSE)) {
      if(any(include.check == "duplicates")) {
        warning("Multiple instances of parameter of interest in given model is presumably\n",
                "not due to interaction or polynomial terms - these models will not be\n",
                "excluded from the computation of model-averaged estimate\n")
      }
      
    }



    ##warn if exclude is neither a list nor NULL
    if(!is.null(exclude)) {
      if(!is.list(exclude)) {stop("\nItems in \"exclude\" must be specified as a list\n")}
    }

  
    ##if exclude is list  
    if(is.list(exclude)) {

      ##determine number of elements in exclude
      nexcl <- length(exclude)

      ##check each formula for presence of exclude variable extracted with formula( )  
      not.include <- lapply(cand.set, FUN=formula)

      ##set up a new list with model formula
      forms <- list()
      for (i in 1:nmods) {
        form.tmp <- strsplit(as.character(not.include[i]), split="~")[[1]][-1]
        if(attr(regexpr("\\+", form.tmp), "match.length")==-1) {  #this line causes problems if intercept is removed from model
          forms[i] <- form.tmp
        } else {forms[i] <- strsplit(form.tmp, split=" \\+ ")}    #this line causes problems if intercept is removed from model
      }

      ##additional check to see whether some variable names include "+"
      check.forms <- unlist(lapply(forms, FUN=function(i) any(attr(regexpr("\\+", i), "match.length")>0)[[1]]))
      if (any(check.forms==TRUE)) stop("\nPlease avoid \"+\" in variable names\n")

      ##additional check to determine if intercept was removed from models
      check.forms <- unlist(lapply(forms, FUN=function(i) any(attr(regexpr("\\- 1", i), "match.length")>0)[[1]]))
      if (any(check.forms==TRUE)) stop("\nModels without intercept are not supported in this version, please use alternative parameterization\n")

 
      ##search within formula for variables to exclude
      mod.exclude <- matrix(NA, nrow=nmods, ncol=nexcl)

      ##iterate over each element in exclude list
      for (var in 1:nexcl) {

        ##iterate over each formula in mod_formula list
        for (i in 1:nmods) {
          idents <- NULL
          form.excl <- forms[[i]]
                    
          ##iterate over each element of forms[[i]]
          for (j in 1:length(form.excl)) {
            idents[j] <- identical(exclude[var][[1]], gsub("(^ +)|( +$)", "", form.excl[j]))
          }
          mod.exclude[i,var] <- ifelse(any(idents==1), 1, 0)
        }    
        
      }
  
      ##determine outcome across all variables to exclude
      to.exclude <- rowSums(mod.exclude)
  
  
      ##exclude models following models from model averaging  
      include[which(to.exclude>=1)] <- 0
      
      
    }



    ##add a check to determine if include always == 0
    if (sum(include)==0) {stop("\nParameter not found in any of the candidate models\n") }

    new.cand.set <- cand.set[which(include==1)] #select models including a given parameter
    new.mod.name <- modnames[which(include==1)]    #update model names
    ##

    new_table <- aictab(cand.set = new.cand.set, modnames = new.mod.name,
                        second.ord = second.ord, nobs = nobs, sort = FALSE)  #recompute AIC table and associated measures
    new_table$Beta_est <- unlist(lapply(new.cand.set, FUN = function(i) coef(i)[paste(parm)])) #extract beta estimate for parm
    new_table$SE <- unlist(lapply(new.cand.set, FUN = function(i) sqrt(diag(vcov(i)))[paste(parm)]))

   
    
    ##AICc
    ##compute model-averaged estimates, unconditional SE, and 95% CL
    if(second.ord == TRUE) {
      Modavg_beta <- sum(new_table$AICcWt*new_table$Beta_est)

      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Uncond_SE <- sum(new_table$AICcWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
      }
      
      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Uncond_SE <- sqrt(sum(new_table$AICcWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
      }
    }

    
    
    ##AIC
    if(second.ord == FALSE) {
      Modavg_beta <- sum(new_table$AICWt*new_table$Beta_est)
      
      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Uncond_SE <- sum(new_table$AICWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
      }
      
      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Uncond_SE <- sqrt(sum(new_table$AICWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
      }
    }

    
    zcrit <- qnorm(p=(1-conf.level)/2, lower.tail=FALSE)
    Lower_CL <- Modavg_beta-zcrit*Uncond_SE
    Upper_CL <- Modavg_beta+zcrit*Uncond_SE
    out.modavg <- list("Parameter"=paste(parm), "Mod.avg.table" = new_table, "Mod.avg.beta" = Modavg_beta,
                       "Uncond.SE" = Uncond_SE, "Conf.level" = conf.level, "Lower.CL"= Lower_CL, "Upper.CL" = Upper_CL)
  
    
    class(out.modavg) <- c("modavg", "list")
    return(out.modavg)

  }


##lme
modavg.AIClme <-
  function(cand.set, parm, modnames = NULL, second.ord = TRUE,
           nobs = NULL, uncond.se = "revised", conf.level = 0.95,
           exclude = NULL, warn = TRUE, ...){

    ##check if named list if modnames are not supplied
    if(is.null(modnames)) {
      if(is.null(names(cand.set))) {
        modnames <- paste("Mod", 1:length(cand.set), sep = "")
        warning("\nModel names have been supplied automatically in the table\n")
      } else {
        modnames <- names(cand.set)
      }
    }

    
#####MODIFICATIONS BEGIN#######
    ##remove all leading and trailing white space and within parm
    parm <- gsub('[[:space:]]+', "", parm)
  
    ##reverse parm
    reversed.parm <- reverse.parm(parm)
    exclude <- reverse.exclude(exclude = exclude)
#####MODIFICATIONS END######

  
    ##extract model formula for each model in cand.set
    mod_formula <- lapply(cand.set, FUN=function(i) labels(summary(i)$coefficients$fixed))
    
    nmods <- length(cand.set)
  
    ##setup matrix to indicate presence of parms in the model
    include <- matrix(NA, nrow = nmods, ncol = 1)
    ##add a check for multiple instances of same variable in given model (i.e., interactions)
    include.check <- matrix(NA, nrow = nmods, ncol = 1)

    ##iterate over each formula in mod_formula list
    for (i in 1:nmods) {
      idents <- NULL
      idents.check <- NULL
      form <- mod_formula[[i]]

######################################################################################################
######################################################################################################
###MODIFICATIONS BEGIN
      ##iterate over each element of formula[[i]] in list
      if(is.null(reversed.parm)) {
        for (j in 1:length(form)) {
          idents[j] <- identical(parm, form[j])
          idents.check[j] <- ifelse(attr(regexpr(parm, form[j], fixed=TRUE), "match.length")== "-1" , 0, 1)  
        }
      } else {
        for (j in 1:length(form)) {
          idents[j] <- identical(parm, form[j]) | identical(reversed.parm, form[j])
          idents.check[j] <- ifelse(attr(regexpr(parm, form[j], fixed=TRUE), "match.length")=="-1" & attr(regexpr(reversed.parm, form[j],
                                                                  fixed=TRUE), "match.length")=="-1" , 0, 1)  
        }
      }
###MODIFICATIONS END
######################################################################################################
######################################################################################################

  
      include[i] <- ifelse(any(idents==1), 1, 0)
      include.check[i] <- ifelse(sum(idents.check)>1, "duplicates", "OK")
    }

#####################################################
    ##exclude == NULL; warn=TRUE:  warn that duplicates occur and stop
    if(is.null(exclude) && identical(warn, TRUE)) {
      ##check for duplicates in same model
      if(any(include.check == "duplicates")) {
        stop("\nSome models include more than one instance of the parameter of interest. \n",
             "This may be due to the presence of interaction/polynomial terms, or variables\n",
             "with similar names:\n",
             "\tsee \"?modavg\" for details on variable specification and \"exclude\" argument\n")
      }

    }

    ##exclude == NULL; warn=FALSE:  compute model-averaged beta estimate from models including variable of interest,
    ##assuming that the variable is not involved in interaction or higher order polynomial (x^2, x^3, etc...),
    ##warn that models were not excluded
    if(is.null(exclude) && identical(warn, FALSE)) {
      if(any(include.check == "duplicates")) {
        warning("\nMultiple instances of parameter of interest in given model is presumably\n",
                "not due to interaction or polynomial terms - these models will not be\n",
                "excluded from the computation of model-averaged estimate\n")
      }
    
    }

    ##warn if exclude is neither a list nor NULL
    if(!is.null(exclude)) {
      if(!is.list(exclude)) {stop("\nItems in \"exclude\" must be specified as a list\n")}
    }


    ##if exclude is list  
    if(is.list(exclude)) {
      
      ##determine number of elements in exclude
      nexcl <- length(exclude)

      ##check each formula for presence of exclude variable extracted with formula( )  
      not.include <- lapply(cand.set, FUN=formula)

      ##set up a new list with model formula
      forms <- list()
      for (i in 1:nmods) {
        form.tmp <- strsplit(as.character(not.include[i]), split="~")[[1]][-1]
        if(attr(regexpr("\\+", form.tmp), "match.length")==-1) {
          forms[i] <- form.tmp
        } else {forms[i] <- strsplit(form.tmp, split=" \\+ ")}
      }

      ##additional check to see whether some variable names include "+"
      check.forms <- unlist(lapply(forms, FUN=function(i) any(attr(regexpr("\\+", i), "match.length")>0)[[1]]))
      if (any(check.forms==TRUE)) stop("\nPlease avoid \"+\" in variable names\n")

      ##additional check to determine if intercept was removed from models
      check.forms <- unlist(lapply(forms, FUN=function(i) any(attr(regexpr("\\- 1", i), "match.length")>0)[[1]]))
      if (any(check.forms==TRUE)) stop("\nModels without intercept are not supported in this version, please use alternative parameterization\n")

      
      ##search within formula for variables to exclude
      mod.exclude <- matrix(NA, nrow=nmods, ncol=nexcl)

      ##iterate over each element in exclude list
      for (var in 1:nexcl) {

        ##iterate over each formula in mod_formula list
        for (i in 1:nmods) {
          idents <- NULL
          form.excl <- forms[[i]]
                    
          ##iterate over each element of forms[[i]]
          for (j in 1:length(form.excl)) {
            idents[j] <- identical(exclude[var][[1]], gsub("(^ +)|( +$)", "", form.excl[j]))
          }
          mod.exclude[i,var] <- ifelse(any(idents==1), 1, 0)
        }    
        
      }
  
      ##determine outcome across all variables to exclude
      to.exclude <- rowSums(mod.exclude)
    
      
      ##exclude models following models from model averaging  
      include[which(to.exclude>=1)] <- 0
      
      
    }


    
    ##add a check to determine if include always == 0
    if (sum(include)==0) {stop("\nParameter not found in any of the candidate models\n") }
    
    new.cand.set <- cand.set[which(include==1)] #select models including a given parameter
    new.mod.name <- modnames[which(include==1)]    #update model names

    new_table <- aictab(cand.set = new.cand.set, modnames = new.mod.name,
                        second.ord = second.ord, nobs = nobs, sort = FALSE)  #recompute AIC table and associated measures
    new_table$Beta_est <- unlist(lapply(new.cand.set, FUN = function(i) fixef(i)[paste(parm)])) #extract beta estimate for parm
    new_table$SE <- unlist(lapply(new.cand.set, FUN = function(i) sqrt(diag(vcov(i)))[paste(parm)]))
  
    ##compute model-averaged estimates, unconditional SE, and 95% CL
    if(second.ord == TRUE) {
      Modavg_beta <- sum(new_table$AICcWt*new_table$Beta_est)

      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Uncond_SE <- sum(new_table$AICcWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
      }

      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Uncond_SE <- sqrt(sum(new_table$AICcWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
      } 


    } else {
      Modavg_beta <- sum(new_table$AICWt*new_table$Beta_est)

      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Uncond_SE <- sum(new_table$AICWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
      }

      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Uncond_SE <- sqrt(sum(new_table$AICWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
      }
    
                      
    }
  
    zcrit <- qnorm(p=(1-conf.level)/2, lower.tail=FALSE)
    Lower_CL<-Modavg_beta-zcrit*Uncond_SE
    Upper_CL<-Modavg_beta+zcrit*Uncond_SE
    out.modavg <- list("Parameter"=paste(parm), "Mod.avg.table" = new_table, "Mod.avg.beta" = Modavg_beta,
                       "Uncond.SE" = Uncond_SE, "Conf.level" = conf.level, "Lower.CL" = Lower_CL,
                       "Upper.CL" = Upper_CL)
    
    class(out.modavg) <- c("modavg", "list")
    return(out.modavg)
  }



##lmekin
modavg.AIClmekin <-
function(cand.set, parm, modnames = NULL, second.ord = TRUE,
         nobs = NULL, uncond.se = "revised", conf.level = 0.95,
         exclude = NULL, warn = TRUE, ...){

  ##check if named list if modnames are not supplied
  if(is.null(modnames)) {
    if(is.null(names(cand.set))) {
      modnames <- paste("Mod", 1:length(cand.set), sep = "")
      warning("\nModel names have been supplied automatically in the table\n")
    } else {
      modnames <- names(cand.set)
    }
  }


#####MODIFICATIONS BEGIN#######
  ##remove all leading and trailing white space and within parm
  parm <- gsub('[[:space:]]+', "", parm)
  
  ##reverse parm
  reversed.parm <- reverse.parm(parm)
  exclude <- reverse.exclude(exclude = exclude)
#####MODIFICATIONS END######

  
  ##extract model formula for each model in cand.set
  mod_formula <- lapply(cand.set, FUN=function(i) labels(fixef(i)))

  nmods <- length(cand.set)
  
  ##setup matrix to indicate presence of parms in the model
  include <- matrix(NA, nrow = nmods, ncol = 1)
  ##add a check for multiple instances of same variable in given model (i.e., interactions)
  include.check <- matrix(NA, nrow = nmods, ncol = 1)

  ##iterate over each formula in mod_formula list
  for (i in 1:nmods) {
    idents <- NULL
    idents.check <- NULL
    form <- mod_formula[[i]]

######################################################################################################
######################################################################################################
###MODIFICATIONS BEGIN
    ##iterate over each element of formula[[i]] in list
    if(is.null(reversed.parm)) {
      for (j in 1:length(form)) {
        idents[j] <- identical(parm, form[j])
        idents.check[j] <- ifelse(attr(regexpr(parm, form[j], fixed=TRUE), "match.length")== "-1" , 0, 1)  
      }
    } else {
      for (j in 1:length(form)) {
        idents[j] <- identical(parm, form[j]) | identical(reversed.parm, form[j])
        idents.check[j] <- ifelse(attr(regexpr(parm, form[j], fixed=TRUE), "match.length")=="-1" & attr(regexpr(reversed.parm, form[j],
                                                                fixed=TRUE), "match.length")=="-1" , 0, 1)  
      }
    }
###MODIFICATIONS END
######################################################################################################
######################################################################################################

  
    include[i] <- ifelse(any(idents==1), 1, 0)
    include.check[i] <- ifelse(sum(idents.check)>1, "duplicates", "OK")
  }

  #####################################################
  ##exclude == NULL; warn=TRUE:  warn that duplicates occur and stop
  if(is.null(exclude) && identical(warn, TRUE)) {
    ##check for duplicates in same model
    if(any(include.check == "duplicates")) {
      stop("\nSome models include more than one instance of the parameter of interest. \n",
           "This may be due to the presence of interaction/polynomial terms, or variables\n",
           "with similar names:\n",
           "\tsee \"?modavg\" for details on variable specification and \"exclude\" argument\n")
    }

  }

  ##exclude == NULL; warn=FALSE:  compute model-averaged beta estimate from models including variable of interest,
  ##assuming that the variable is not involved in interaction or higher order polynomial (x^2, x^3, etc...),
  ##warn that models were not excluded
  if(is.null(exclude) && identical(warn, FALSE)) {
    if(any(include.check == "duplicates")) {
      warning("\nMultiple instances of parameter of interest in given model is presumably\n",
              "not due to interaction or polynomial terms - these models will not be\n",
              "excluded from the computation of model-averaged estimate\n")
    }
    
  }

  ##warn if exclude is neither a list nor NULL
  if(!is.null(exclude)) {
    if(!is.list(exclude)) {stop("\nItems in \"exclude\" must be specified as a list\n")}
  }


  ##if exclude is list  
  if(is.list(exclude)) {

    ##determine number of elements in exclude
    nexcl <- length(exclude)

    ##check each formula for presence of exclude variable extracted with formula( )  
    not.include <- lapply(cand.set, FUN=formula)

    ##set up a new list with model formula
    forms <- list()
    for (i in 1:nmods) {
      form.tmp <- strsplit(as.character(not.include[i]), split="~")[[1]][-1]
      if(attr(regexpr("\\+", form.tmp), "match.length")==-1) {
        forms[i] <- form.tmp
      } else {forms[i] <- strsplit(form.tmp, split=" \\+ ")}
    }

    ##additional check to see whether some variable names include "+"
    check.forms <- unlist(lapply(forms, FUN=function(i) any(attr(regexpr("\\+", i), "match.length")>0)[[1]]))
    if (any(check.forms==TRUE)) stop("\nPlease avoid \"+\" in variable names\n")

    ##additional check to determine if intercept was removed from models
    check.forms <- unlist(lapply(forms, FUN=function(i) any(attr(regexpr("\\- 1", i), "match.length")>0)[[1]]))
    if (any(check.forms==TRUE)) stop("\nModels without intercept are not supported in this version, please use alternative parameterization\n")

 
    ##search within formula for variables to exclude
    mod.exclude <- matrix(NA, nrow=nmods, ncol=nexcl)

    ##iterate over each element in exclude list
    for (var in 1:nexcl) {

      ##iterate over each formula in mod_formula list
      for (i in 1:nmods) {
        idents <- NULL
        form.excl <- forms[[i]]
                    
        ##iterate over each element of forms[[i]]
        for (j in 1:length(form.excl)) {
          idents[j] <- identical(exclude[var][[1]], gsub("(^ +)|( +$)", "", form.excl[j]))
        }
        mod.exclude[i,var] <- ifelse(any(idents==1), 1, 0)
      }    
        
    }
  
    ##determine outcome across all variables to exclude
    to.exclude <- rowSums(mod.exclude)
    
  
    ##exclude models following models from model averaging  
    include[which(to.exclude>=1)] <- 0
      
      
  }

  ##add a check to determine if include always == 0
  if (sum(include)==0) {stop("\nParameter not found in any of the candidate models\n") }
  
  new.cand.set <- cand.set[which(include==1)] #select models including a given parameter
  new.mod.name <- modnames[which(include==1)]    #update model names

  ##compute table
  new_table <- aictab(cand.set = new.cand.set, modnames = new.mod.name,
                      second.ord = second.ord, nobs = nobs, sort = FALSE)  #recompute AIC table and associated measures
  new_table$Beta_est <- unlist(lapply(new.cand.set, FUN = function(i) fixef(i)[paste(parm)])) #extract beta estimate for parm
  new_table$SE <- unlist(lapply(new.cand.set, FUN = function(i) extractSE(i)[paste(parm)]))               

  ##compute model-averaged estimates, unconditional SE, and 95% CL
  if(second.ord == TRUE) {
    Modavg_beta <- sum(new_table$AICcWt*new_table$Beta_est)

    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Uncond_SE <- sum(new_table$AICcWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
    }

    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Uncond_SE <- sqrt(sum(new_table$AICcWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
    } 


  } else {
    Modavg_beta <- sum(new_table$AICWt*new_table$Beta_est)

    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Uncond_SE <- sum(new_table$AICWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
    }

    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Uncond_SE <- sqrt(sum(new_table$AICWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
    }
    
                      
  }
  
  zcrit <- qnorm(p = (1 - conf.level)/2, lower.tail = FALSE)
  Lower_CL <- Modavg_beta - zcrit*Uncond_SE
  Upper_CL <- Modavg_beta + zcrit*Uncond_SE
  out.modavg <- list("Parameter" = paste(parm), "Mod.avg.table" = new_table, "Mod.avg.beta" = Modavg_beta,
                     "Uncond.SE" = Uncond_SE, "Conf.level" = conf.level, "Lower.CL" = Lower_CL,
                     "Upper.CL" = Upper_CL)
  
  class(out.modavg) <- c("modavg", "list")
  return(out.modavg)
}



##maxlike
modavg.AICmaxlikeFit.list <-
  function(cand.set, parm, modnames = NULL, second.ord = TRUE,
           nobs = NULL, uncond.se = "revised", conf.level = 0.95,
           exclude = NULL, warn = TRUE, c.hat = 1, ...){

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

    ##check that link function is the same for all models
    check.link <- unlist(lapply(X = cand.set, FUN=function(i) i$link))
    unique.link <- unique(x = check.link)
    if(length(unique.link) > 1) stop("\nIt is not appropriate to compute a model-averaged beta estimate\n",
                                         "from models using different link functions\n")
  
    
#####MODIFICATIONS BEGIN#######
    ##remove all leading and trailing white space and within parm
    parm <- gsub('[[:space:]]+', "", parm)
  
    ##reverse parm
    reversed.parm <- reverse.parm(parm)
    exclude <- reverse.exclude(exclude = exclude)
#####MODIFICATIONS END######

  
    ##extract model formula for each model in cand.set
    mod_formula <- lapply(cand.set, FUN = function(i) names(coef(i)))
       
    nmods <- length(cand.set)
  
    ##setup matrix to indicate presence of parms in the model
    include <- matrix(NA, nrow = nmods, ncol = 1)
    ##add a check for multiple instances of same variable in given model (i.e., interactions)
    include.check <- matrix(NA, nrow = nmods, ncol = 1)

    ##iterate over each formula in mod_formula list
    for (i in 1:nmods) {
      idents <- NULL
      idents.check <- NULL
      form <- mod_formula[[i]]

######################################################################################################
######################################################################################################
###MODIFICATIONS BEGIN
      ##iterate over each element of formula[[i]] in list
      if(is.null(reversed.parm)) {
        for (j in 1:length(form)) {
          idents[j] <- identical(parm, form[j])
          idents.check[j] <- ifelse(attr(regexpr(parm, form[j], fixed=TRUE), "match.length")== "-1" , 0, 1)  
        }
      } else {
        for (j in 1:length(form)) {
          idents[j] <- identical(parm, form[j]) | identical(reversed.parm, form[j])
          idents.check[j] <- ifelse(attr(regexpr(parm, form[j], fixed=TRUE), "match.length")=="-1" & attr(regexpr(reversed.parm, form[j],
                                                                  fixed=TRUE), "match.length")=="-1" , 0, 1)  
        }
      }
###MODIFICATIONS END
######################################################################################################
######################################################################################################

  
      include[i] <- ifelse(any(idents==1), 1, 0)
      include.check[i] <- ifelse(sum(idents.check)>1, "duplicates", "OK")
    }

  #####################################################
    ##exclude == NULL; warn=TRUE:  warn that duplicates occur and stop
    if(is.null(exclude) && identical(warn, TRUE)) {
      ##check for duplicates in same model
      if(any(include.check == "duplicates")) {
        stop("\nSome models include more than one instance of the parameter of interest. \n",
             "This may be due to the presence of interaction/polynomial terms, or variables\n",
             "with similar names:\n",
             "\tsee \"?modavg\" for details on variable specification and \"exclude\" argument\n")
      }

    }

    ##exclude == NULL; warn=FALSE:  compute model-averaged beta estimate from models including variable of interest,
    ##assuming that the variable is not involved in interaction or higher order polynomial (x^2, x^3, etc...),
    ##warn that models were not excluded
    if(is.null(exclude) && identical(warn, FALSE)) {
      if(any(include.check == "duplicates")) {
        warning("\nMultiple instances of parameter of interest in given model is presumably\n",
                "not due to interaction or polynomial terms - these models will not be\n",
                "excluded from the computation of model-averaged estimate\n")
      }
      
    }

    ##warn if exclude is neither a list nor NULL
    if(!is.null(exclude)) {
      if(!is.list(exclude)) {stop("\nItems in \"exclude\" must be specified as a list\n")}
    }


    ##if exclude is list  
    if(is.list(exclude)) {
      
      ##determine number of elements in exclude
      nexcl <- length(exclude)
      
      ##check each formula for presence of exclude variable extracted with formula( )  
      not.include <- lapply(cand.set, FUN = formula)

      ##set up a new list with model formula
      forms <- list()
      for (i in 1:nmods) {
        form.tmp <- as.character(not.include[[i]])[3]
        if(attr(regexpr("\\+", form.tmp), "match.length")==-1) {
          forms[i] <- form.tmp
        } else {forms[i] <- strsplit(form.tmp, split=" \\+ ")}
      }

      ##additional check to see whether some variable names include "+"
      check.forms <- unlist(lapply(forms, FUN=function(i) any(attr(regexpr("\\+", i), "match.length")>0)[[1]]))
      if (any(check.forms==TRUE)) stop("\nPlease avoid \"+\" in variable names\n")

      ##additional check to determine if intercept was removed from models
      check.forms <- unlist(lapply(forms, FUN=function(i) any(attr(regexpr("\\- 1", i), "match.length")>0)[[1]]))
      if (any(check.forms==TRUE)) stop("\nModels without intercept are not supported in this version, please use alternative parameterization\n")
         
      
      ##search within formula for variables to exclude
      mod.exclude <- matrix(NA, nrow=nmods, ncol=nexcl)
      
      ##iterate over each element in exclude list
      for (var in 1:nexcl) {

        ##iterate over each formula in mod_formula list
        for (i in 1:nmods) {
          idents <- NULL
          form.excl <- forms[[i]]
                    
          ##iterate over each element of forms[[i]]
          for (j in 1:length(form.excl)) {
            idents[j] <- identical(exclude[var][[1]], gsub("(^ +)|( +$)", "", form.excl[j]))
          }
          mod.exclude[i,var] <- ifelse(any(idents==1), 1, 0)
        }    
        
      }
  
      ##determine outcome across all variables to exclude
      to.exclude <- rowSums(mod.exclude)
         
      
      ##exclude models following models from model averaging  
      include[which(to.exclude>=1)] <- 0
      
      
    }

       
       
    ##add a check to determine if include always == 0
    if (sum(include)==0) {stop("\nParameter not found in any of the candidate models\n") }
    
    new.cand.set <- cand.set[which(include==1)] #select models including a given parameter
    new.mod.name <- modnames[which(include==1)]    #update model names
    
    new_table <- aictab(cand.set=new.cand.set, modnames=new.mod.name, second.ord=second.ord,
                        nobs=nobs, sort=FALSE, c.hat = 1)  #recompute AIC table and associated measures
    new_table$Beta_est <- unlist(lapply(new.cand.set, FUN=function(i) coef(i)[paste(parm)])) #extract beta estimate for parm
    new_table$SE <- unlist(lapply(new.cand.set, FUN=function(i) sqrt(diag(vcov(i)))[paste(parm)]))
  
    ##compute model-averaged estimates, unconditional SE, and 95% CL
    if(second.ord == TRUE) {
      Modavg_beta <- sum(new_table$AICcWt*new_table$Beta_est)
         
      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Uncond_SE <- sum(new_table$AICcWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
      }

      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Uncond_SE <- sqrt(sum(new_table$AICcWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
      } 


    } else {
      Modavg_beta <- sum(new_table$AICWt*new_table$Beta_est)
      
      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Uncond_SE <- sum(new_table$AICWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
      }

      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Uncond_SE <- sqrt(sum(new_table$AICWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
      }
    
      
    }
     
    zcrit <- qnorm(p=(1-conf.level)/2, lower.tail=FALSE)
    Lower_CL <- Modavg_beta - zcrit*Uncond_SE
    Upper_CL <- Modavg_beta + zcrit*Uncond_SE
    out.modavg <- list("Parameter"=paste(parm), "Mod.avg.table" = new_table, "Mod.avg.beta" = Modavg_beta,
                       "Uncond.SE" = Uncond_SE, "Conf.level" = conf.level, "Lower.CL" = Lower_CL,
                       "Upper.CL" = Upper_CL)
    
    class(out.modavg) <- c("modavg", "list")
    return(out.modavg)
  }



##mer
modavg.AICmer <-
function(cand.set, parm, modnames = NULL, second.ord = TRUE,
         nobs = NULL, uncond.se = "revised", conf.level = 0.95,
         exclude = NULL, warn = TRUE, ...){

  ##check if named list if modnames are not supplied
  if(is.null(modnames)) {
    if(is.null(names(cand.set))) {
      modnames <- paste("Mod", 1:length(cand.set), sep = "")
      warning("\nModel names have been supplied automatically in the table\n")
    } else {
      modnames <- names(cand.set)
    }
  }
    
  
#####MODIFICATIONS BEGIN#######
  ##remove all leading and trailing white space and within parm
  parm <- gsub('[[:space:]]+', "", parm)
  
  ##reverse parm
  reversed.parm <- reverse.parm(parm)
  exclude <- reverse.exclude(exclude = exclude)
#####MODIFICATIONS END######

##check for nlmer models
if(any(unlist(lapply(X = cand.set, FUN = function(i) as.character(i@call)[1])) == "nlmer")) warning("\nNon-linear models are part of the candidate model set: model-averaged estimates may not be meaningful\n")

###################
  ##determine families of model
  fam.list <- unlist(lapply(X = cand.set, FUN = function(i) fam.link.mer(i)$family))
  check.fam <- unique(fam.list)
  if(length(check.fam) > 1) stop("\nIt is not appropriate to compute a model-averaged beta estimate\n",
                                       "from models using different families of distributions\n")
  ##determine link functions
  link.list <- unlist(lapply(X = cand.set, FUN = function(i) fam.link.mer(i)$link))
  check.link <- unique(link.list)
  if(length(check.link) > 1) stop("\nIt is not appropriate to compute a model-averaged beta estimate\n",
                                       "from models using different link functions\n")
###################       

  
  ##extract model formula for each model in cand.set
  mod_formula <- lapply(cand.set, FUN=function(i) labels(fixef(i)))

  nmods <- length(cand.set)
  
  ##setup matrix to indicate presence of parms in the model
  include <- matrix(NA, nrow=nmods, ncol=1)
  ##add a check for multiple instances of same variable in given model (i.e., interactions)
  include.check <- matrix(NA, nrow=nmods, ncol=1)

  ##iterate over each formula in mod_formula list
  for (i in 1:nmods) {
    idents <- NULL
    idents.check <- NULL
    form <- mod_formula[[i]]

######################################################################################################
######################################################################################################
###MODIFICATIONS BEGIN
    ##iterate over each element of formula[[i]] in list
      if(is.null(reversed.parm)) {
        for (j in 1:length(form)) {
          idents[j] <- identical(parm, form[j])
          idents.check[j] <- ifelse(attr(regexpr(parm, form[j], fixed=TRUE), "match.length")== "-1" , 0, 1)  
        }
      } else {
        for (j in 1:length(form)) {
          idents[j] <- identical(parm, form[j]) | identical(reversed.parm, form[j])
          idents.check[j] <- ifelse(attr(regexpr(parm, form[j], fixed=TRUE), "match.length")=="-1" & attr(regexpr(reversed.parm, form[j],
                                                                  fixed=TRUE), "match.length")=="-1" , 0, 1)  
        }
      }
###MODIFICATIONS END
######################################################################################################
######################################################################################################

    include[i] <- ifelse(any(idents==1), 1, 0)
    include.check[i] <- ifelse(sum(idents.check)>1, "duplicates", "OK")
  }

#####################################################
###exclude == NULL; warn=TRUE:  warn that duplicates occur and stop
  if(is.null(exclude) && identical(warn, TRUE)) {
    ##check for duplicates in same model
    if(any(include.check == "duplicates")) {
      stop("\nSome models include more than one instance of the parameter of interest. \n",
           "This may be due to the presence of interaction/polynomial terms, or variables\n",
           "with similar names:\n",
           "\tsee \"?modavg\" for details on variable specification and \"exclude\" argument\n")
    }
    
  }
  
  ##exclude == NULL; warn=FALSE:  compute model-averaged beta estimate from models including variable of interest,
  ##assuming that the variable is not involved in interaction or higher order polynomial (x^2, x^3, etc...),
  ##warn that models were not excluded
  if(is.null(exclude) && identical(warn, FALSE)) {
    if(any(include.check == "duplicates")) {
      warning("\nMultiple instances of parameter of interest in given model is presumably\n",
              "not due to interaction or polynomial terms - these models will not be\n",
              "excluded from the computation of model-averaged estimate\n")
    }
    
  }
  
  ##warn if exclude is neither a list nor NULL
  if(!is.null(exclude)) {
    if(!is.list(exclude)) {stop("\nItems in \"exclude\" must be specified as a list\n")}
  }
  

  ##if exclude is list  
  if(is.list(exclude)) {
    
    ##determine number of elements in exclude
    nexcl <- length(exclude)

    ##check each formula for presence of exclude variable extracted with formula( )  
    not.include <- lapply(cand.set, FUN=formula)  #random effect portion is returned within parentheses
    ##because matching uses identical( ) to check fixed effects against formula( ),
    ##should not be problematic for variables included in random effects
    
    
    ##set up a new list with model formula
    forms <- list()
    for (i in 1:nmods) {
      form.tmp <- strsplit(as.character(not.include[i]), split="~")[[1]][-1]
      if(attr(regexpr("\\+", form.tmp), "match.length")==-1) {
        forms[i] <- form.tmp
      } else {forms[i] <- strsplit(form.tmp, split=" \\+ ")}
    }

    ######################################
    ##remove leading and trailing spaces as well as spaces within string
    ##forms <- lapply(forms.space, FUN = function(b) gsub('[[:space:]]+', "", b)) 
    ######################################
    
    ##additional check to see whether some variable names include "+"
    check.forms <- unlist(lapply(forms, FUN=function(i) any(attr(regexpr("\\+", i), "match.length")>0)[[1]]))
    if (any(check.forms==TRUE)) stop("\nPlease avoid \"+\" in variable names\n")
    
    ##additional check to determine if intercept was removed from models
    check.forms <- unlist(lapply(forms, FUN=function(i) any(attr(regexpr("\\- 1", i), "match.length")>0)[[1]]))
    if (any(check.forms==TRUE)) stop("\nModels without intercept are not supported in this version, please use alternative parameterization\n")

 
    ##search within formula for variables to exclude
    mod.exclude <- matrix(NA, nrow=nmods, ncol=nexcl)

    ##iterate over each element in exclude list
    for (var in 1:nexcl) {
      
      ##iterate over each formula in mod_formula list
      for (i in 1:nmods) {
        idents <- NULL
        form.excl <- forms[[i]]
                    
        ##iterate over each element of forms[[i]]
        for (j in 1:length(form.excl)) {
          idents[j] <- identical(exclude[var][[1]], gsub("(^ +)|( +$)", "", form.excl[j]))
        }
        mod.exclude[i,var] <- ifelse(any(idents==1), 1, 0)
      }
      
    }
  
    ##determine outcome across all variables to exclude
    to.exclude <- rowSums(mod.exclude)
    
    
    ##exclude models following models from model averaging  
    include[which(to.exclude>=1)] <- 0
    
    
  }


  
  ##add a check to determine if include always == 0
  if (sum(include)==0) {stop("\nParameter not found in any of the candidate models\n") }
  
  new.cand.set <- cand.set[which(include==1)] #select models including a given parameter
  new.mod.name <- modnames[which(include==1)]    #update model names
  
  new_table <- aictab(cand.set=new.cand.set, modnames=new.mod.name,
                      second.ord=second.ord, nobs=nobs, sort=FALSE)  #recompute AIC table and associated measures
  new_table$Beta_est <- unlist(lapply(new.cand.set, FUN=function(i) fixef(i)[paste(parm)])) #extract beta estimate for parm
  new_table$SE <- unlist(lapply(new.cand.set, FUN=function(i) extractSE(i)[paste(parm)]))
  
  ##compute model-averaged estimates, unconditional SE, and 95% CL
  if(second.ord==TRUE) {
    Modavg_beta <- sum(new_table$AICcWt*new_table$Beta_est)
    
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Uncond_SE <- sum(new_table$AICcWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
    }

    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Uncond_SE <- sqrt(sum(new_table$AICcWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
    } 


  } else {
    Modavg_beta<-sum(new_table$AICWt*new_table$Beta_est)
    
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Uncond_SE<-sum(new_table$AICWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
    }
    
    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Uncond_SE<-sqrt(sum(new_table$AICWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
    }
    
    
  }
  
  zcrit <- qnorm(p=(1-conf.level)/2, lower.tail=FALSE)
  Lower_CL<-Modavg_beta-zcrit*Uncond_SE
  Upper_CL<-Modavg_beta+zcrit*Uncond_SE
  out.modavg <- list("Parameter"=paste(parm), "Mod.avg.table" = new_table, "Mod.avg.beta" = Modavg_beta,
                     "Uncond.SE" = Uncond_SE, "Conf.level" = conf.level, "Lower.CL" = Lower_CL,
                     "Upper.CL" = Upper_CL)
  
  class(out.modavg) <- c("modavg", "list")
  return(out.modavg)
}



##lmerMod
modavg.AIClmerMod <-
  function(cand.set, parm, modnames = NULL, second.ord = TRUE,
           nobs = NULL, uncond.se = "revised", conf.level = 0.95,
           exclude = NULL, warn = TRUE, ...){

    ##check if named list if modnames are not supplied
    if(is.null(modnames)) {
      if(is.null(names(cand.set))) {
        modnames <- paste("Mod", 1:length(cand.set), sep = "")
        warning("\nModel names have been supplied automatically in the table\n")
      } else {
        modnames <- names(cand.set)
      }
    }
  
    
#####MODIFICATIONS BEGIN#######
      ##remove all leading and trailing white space and within parm
      parm <- gsub('[[:space:]]+', "", parm)
      
      ##reverse parm
      reversed.parm <- reverse.parm(parm)
      exclude <- reverse.exclude(exclude = exclude)
#####MODIFICATIONS END######

      ##check for nlmer models
                                        #if(any(unlist(lapply(X = cand.set, FUN = function(i) as.character(i@call)[1])) == "nlmer")) warning("\nNon-linear models are part of the candidate model set: model-averaged estimates may not be meaningful\n")
      
###################
      
###################       

  
      ##extract model formula for each model in cand.set
      mod_formula <- lapply(cand.set, FUN=function(i) labels(fixef(i)))
      
      nmods <- length(cand.set)
  
      ##setup matrix to indicate presence of parms in the model
      include <- matrix(NA, nrow=nmods, ncol=1)
      ##add a check for multiple instances of same variable in given model (i.e., interactions)
      include.check <- matrix(NA, nrow=nmods, ncol=1)

      ##iterate over each formula in mod_formula list
      for (i in 1:nmods) {
        idents <- NULL
        idents.check <- NULL
        form <- mod_formula[[i]]

######################################################################################################
######################################################################################################
###MODIFICATIONS BEGIN
        ##iterate over each element of formula[[i]] in list
        if(is.null(reversed.parm)) {
          for (j in 1:length(form)) {
            idents[j] <- identical(parm, form[j])
            idents.check[j] <- ifelse(attr(regexpr(parm, form[j], fixed=TRUE), "match.length")== "-1" , 0, 1)  
          }
        } else {
          for (j in 1:length(form)) {
            idents[j] <- identical(parm, form[j]) | identical(reversed.parm, form[j])
            idents.check[j] <- ifelse(attr(regexpr(parm, form[j], fixed=TRUE), "match.length")=="-1" & attr(regexpr(reversed.parm, form[j],
                                                                    fixed=TRUE), "match.length")=="-1" , 0, 1)  
          }
        }
###MODIFICATIONS END
######################################################################################################
######################################################################################################
        
        include[i] <- ifelse(any(idents==1), 1, 0)
        include.check[i] <- ifelse(sum(idents.check)>1, "duplicates", "OK")
      }

#####################################################
###exclude == NULL; warn=TRUE:  warn that duplicates occur and stop
      if(is.null(exclude) && identical(warn, TRUE)) {
        ##check for duplicates in same model
        if(any(include.check == "duplicates")) {
          stop("\nSome models include more than one instance of the parameter of interest. \n",
               "This may be due to the presence of interaction/polynomial terms, or variables\n",
               "with similar names:\n",
               "\tsee \"?modavg\" for details on variable specification and \"exclude\" argument\n")
        }
        
      }
  
      ##exclude == NULL; warn=FALSE:  compute model-averaged beta estimate from models including variable of interest,
      ##assuming that the variable is not involved in interaction or higher order polynomial (x^2, x^3, etc...),
      ##warn that models were not excluded
      if(is.null(exclude) && identical(warn, FALSE)) {
        if(any(include.check == "duplicates")) {
          warning("\nMultiple instances of parameter of interest in given model is presumably\n",
                  "not due to interaction or polynomial terms - these models will not be\n",
                  "excluded from the computation of model-averaged estimate\n")
        }
    
      }
  
      ##warn if exclude is neither a list nor NULL
      if(!is.null(exclude)) {
        if(!is.list(exclude)) {stop("\nItems in \"exclude\" must be specified as a list\n")}
      }
  

      ##if exclude is list  
      if(is.list(exclude)) {
        
        ##determine number of elements in exclude
        nexcl <- length(exclude)

        ##check each formula for presence of exclude variable extracted with formula( )  
        not.include <- lapply(cand.set, FUN = formula)  #random effect portion is returned within parentheses
        ##because matching uses identical( ) to check fixed effects against formula( ),
        ##should not be problematic for variables included in random effects
        
        ##set up a new list with model formula
        forms <- list( )
        for (i in 1:nmods) {
          form.tmp <- strsplit(as.character(not.include[i]), split="~")[[1]][-1]
          if(attr(regexpr("\\+", form.tmp), "match.length")==-1) {
            forms[i] <- form.tmp
          } else {forms[i] <- strsplit(form.tmp, split=" \\+ ")}
        }

######################################
        ##remove leading and trailing spaces as well as spaces within string
        ##forms <- lapply(forms.space, FUN = function(b) gsub('[[:space:]]+', "", b)) 
######################################
    
        ##additional check to see whether some variable names include "+"
        check.forms <- unlist(lapply(forms, FUN=function(i) any(attr(regexpr("\\+", i), "match.length")>0)[[1]]))
        if (any(check.forms==TRUE)) stop("\nPlease avoid \"+\" in variable names\n")
        
        ##additional check to determine if intercept was removed from models
        check.forms <- unlist(lapply(forms, FUN=function(i) any(attr(regexpr("\\- 1", i), "match.length")>0)[[1]]))
        if (any(check.forms==TRUE)) stop("\nModels without intercept are not supported in this version, please use alternative parameterization\n")

        
        ##search within formula for variables to exclude
        mod.exclude <- matrix(NA, nrow=nmods, ncol=nexcl)

        ##iterate over each element in exclude list
        for (var in 1:nexcl) {
      
          ##iterate over each formula in mod_formula list
          for (i in 1:nmods) {
            idents <- NULL
            form.excl <- forms[[i]]
            
            ##iterate over each element of forms[[i]]
            for (j in 1:length(form.excl)) {
              idents[j] <- identical(exclude[var][[1]], gsub("(^ +)|( +$)", "", form.excl[j]))
            }
            mod.exclude[i,var] <- ifelse(any(idents==1), 1, 0)
          }
          
        }
  
        ##determine outcome across all variables to exclude
        to.exclude <- rowSums(mod.exclude)
    
        ##exclude models following models from model averaging  
        include[which(to.exclude>=1)] <- 0
    
      }
      
      ##add a check to determine if include always == 0
      if (sum(include)==0) {stop("\nParameter not found in any of the candidate models\n") }
  
      new.cand.set <- cand.set[which(include==1)] #select models including a given parameter
      new.mod.name <- modnames[which(include==1)]    #update model names
      
      new_table <- aictab(cand.set=new.cand.set, modnames=new.mod.name,
                          second.ord=second.ord, nobs=nobs, sort=FALSE)  #recompute AIC table and associated measures
      new_table$Beta_est <- unlist(lapply(new.cand.set, FUN=function(i) fixef(i)[paste(parm)])) #extract beta estimate for parm
      new_table$SE <- unlist(lapply(new.cand.set, FUN=function(i) extractSE(i)[paste(parm)]))
  
      ##compute model-averaged estimates, unconditional SE, and 95% CL
      if(second.ord==TRUE) {
        Modavg_beta <- sum(new_table$AICcWt*new_table$Beta_est)
        
        ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
        if(identical(uncond.se, "old")) {
          Uncond_SE <- sum(new_table$AICcWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
        }

        ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
        if(identical(uncond.se, "revised")) {
          Uncond_SE <- sqrt(sum(new_table$AICcWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
        } 


      } else {
        Modavg_beta<-sum(new_table$AICWt*new_table$Beta_est)
    
        ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
        if(identical(uncond.se, "old")) {
          Uncond_SE<-sum(new_table$AICWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
        }
    
        ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
        if(identical(uncond.se, "revised")) {
          Uncond_SE<-sqrt(sum(new_table$AICWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
        }
    
    
      }
  
      zcrit <- qnorm(p=(1-conf.level)/2, lower.tail=FALSE)
      Lower_CL <- Modavg_beta-zcrit*Uncond_SE
      Upper_CL <- Modavg_beta+zcrit*Uncond_SE
      out.modavg <- list("Parameter"=paste(parm), "Mod.avg.table" = new_table, "Mod.avg.beta" = Modavg_beta,
                         "Uncond.SE" = Uncond_SE, "Conf.level" = conf.level, "Lower.CL" = Lower_CL,
                         "Upper.CL" = Upper_CL)
  
      class(out.modavg) <- c("modavg", "list")
      return(out.modavg)
      
    }




##glmerMod
modavg.AICglmerMod <-
  function(cand.set, parm, modnames = NULL, second.ord = TRUE,
           nobs = NULL, uncond.se = "revised", conf.level = 0.95,
           exclude = NULL, warn = TRUE, ...){

    ##check if named list if modnames are not supplied
    if(is.null(modnames)) {
      if(is.null(names(cand.set))) {
        modnames <- paste("Mod", 1:length(cand.set), sep = "")
        warning("\nModel names have been supplied automatically in the table\n")
      } else {
        modnames <- names(cand.set)
      }
    }
  
    
#####MODIFICATIONS BEGIN#######
      ##remove all leading and trailing white space and within parm
      parm <- gsub('[[:space:]]+', "", parm)
      
      ##reverse parm
      reversed.parm <- reverse.parm(parm)
      exclude <- reverse.exclude(exclude = exclude)
#####MODIFICATIONS END######

      ##check for nlmer models
                                        #if(any(unlist(lapply(X = cand.set, FUN = function(i) as.character(i@call)[1])) == "nlmer")) warning("\nNon-linear models are part of the candidate model set: model-averaged estimates may not be meaningful\n")
      
###################
      ##determine families of model
      fam.list <- unlist(lapply(X = cand.set, FUN = function(i) fam.link.mer(i)$family))
      check.fam <- unique(fam.list)
      if(length(check.fam) > 1) stop("\nIt is not appropriate to compute a model-averaged beta estimate\n",
                                     "from models using different families of distributions\n")
      ##determine link functions
      link.list <- unlist(lapply(X = cand.set, FUN = function(i) fam.link.mer(i)$link))
      check.link <- unique(link.list)
      if(length(check.link) > 1) stop("\nIt is not appropriate to compute a model-averaged beta estimate\n",
                                      "from models using different link functions\n")
###################       

  
      ##extract model formula for each model in cand.set
      mod_formula <- lapply(cand.set, FUN=function(i) labels(fixef(i)))
      
      nmods <- length(cand.set)
  
      ##setup matrix to indicate presence of parms in the model
      include <- matrix(NA, nrow=nmods, ncol=1)
      ##add a check for multiple instances of same variable in given model (i.e., interactions)
      include.check <- matrix(NA, nrow=nmods, ncol=1)

      ##iterate over each formula in mod_formula list
      for (i in 1:nmods) {
        idents <- NULL
        idents.check <- NULL
        form <- mod_formula[[i]]

######################################################################################################
######################################################################################################
###MODIFICATIONS BEGIN
        ##iterate over each element of formula[[i]] in list
        if(is.null(reversed.parm)) {
          for (j in 1:length(form)) {
            idents[j] <- identical(parm, form[j])
            idents.check[j] <- ifelse(attr(regexpr(parm, form[j], fixed=TRUE), "match.length")== "-1" , 0, 1)  
          }
        } else {
          for (j in 1:length(form)) {
            idents[j] <- identical(parm, form[j]) | identical(reversed.parm, form[j])
            idents.check[j] <- ifelse(attr(regexpr(parm, form[j], fixed=TRUE), "match.length")=="-1" & attr(regexpr(reversed.parm, form[j],
                                                                    fixed=TRUE), "match.length")=="-1" , 0, 1)  
          }
        }
###MODIFICATIONS END
######################################################################################################
######################################################################################################
        
        include[i] <- ifelse(any(idents==1), 1, 0)
        include.check[i] <- ifelse(sum(idents.check)>1, "duplicates", "OK")
      }

#####################################################
###exclude == NULL; warn=TRUE:  warn that duplicates occur and stop
      if(is.null(exclude) && identical(warn, TRUE)) {
        ##check for duplicates in same model
        if(any(include.check == "duplicates")) {
          stop("\nSome models include more than one instance of the parameter of interest. \n",
               "This may be due to the presence of interaction/polynomial terms, or variables\n",
               "with similar names:\n",
               "\tsee \"?modavg\" for details on variable specification and \"exclude\" argument\n")
        }
        
      }
  
      ##exclude == NULL; warn=FALSE:  compute model-averaged beta estimate from models including variable of interest,
      ##assuming that the variable is not involved in interaction or higher order polynomial (x^2, x^3, etc...),
      ##warn that models were not excluded
      if(is.null(exclude) && identical(warn, FALSE)) {
        if(any(include.check == "duplicates")) {
          warning("\nMultiple instances of parameter of interest in given model is presumably\n",
                  "not due to interaction or polynomial terms - these models will not be\n",
                  "excluded from the computation of model-averaged estimate\n")
        }
    
      }
  
      ##warn if exclude is neither a list nor NULL
      if(!is.null(exclude)) {
        if(!is.list(exclude)) {stop("\nItems in \"exclude\" must be specified as a list\n")}
      }
  

      ##if exclude is list  
      if(is.list(exclude)) {
        
        ##determine number of elements in exclude
        nexcl <- length(exclude)

        ##check each formula for presence of exclude variable extracted with formula( )  
        not.include <- lapply(cand.set, FUN = formula)  #random effect portion is returned within parentheses
        ##because matching uses identical( ) to check fixed effects against formula( ),
        ##should not be problematic for variables included in random effects
        
        ##set up a new list with model formula
        forms <- list( )
        for (i in 1:nmods) {
          form.tmp <- strsplit(as.character(not.include[i]), split="~")[[1]][-1]
          if(attr(regexpr("\\+", form.tmp), "match.length")==-1) {
            forms[i] <- form.tmp
          } else {forms[i] <- strsplit(form.tmp, split=" \\+ ")}
        }

######################################
        ##remove leading and trailing spaces as well as spaces within string
        ##forms <- lapply(forms.space, FUN = function(b) gsub('[[:space:]]+', "", b)) 
######################################
    
        ##additional check to see whether some variable names include "+"
        check.forms <- unlist(lapply(forms, FUN=function(i) any(attr(regexpr("\\+", i), "match.length")>0)[[1]]))
        if (any(check.forms==TRUE)) stop("\nPlease avoid \"+\" in variable names\n")
        
        ##additional check to determine if intercept was removed from models
        check.forms <- unlist(lapply(forms, FUN=function(i) any(attr(regexpr("\\- 1", i), "match.length")>0)[[1]]))
        if (any(check.forms==TRUE)) stop("\nModels without intercept are not supported in this version, please use alternative parameterization\n")

        
        ##search within formula for variables to exclude
        mod.exclude <- matrix(NA, nrow=nmods, ncol=nexcl)

        ##iterate over each element in exclude list
        for (var in 1:nexcl) {
      
          ##iterate over each formula in mod_formula list
          for (i in 1:nmods) {
            idents <- NULL
            form.excl <- forms[[i]]
            
            ##iterate over each element of forms[[i]]
            for (j in 1:length(form.excl)) {
              idents[j] <- identical(exclude[var][[1]], gsub("(^ +)|( +$)", "", form.excl[j]))
            }
            mod.exclude[i,var] <- ifelse(any(idents==1), 1, 0)
          }
          
        }
  
        ##determine outcome across all variables to exclude
        to.exclude <- rowSums(mod.exclude)
    
        ##exclude models following models from model averaging  
        include[which(to.exclude>=1)] <- 0
    
      }
      
      ##add a check to determine if include always == 0
      if (sum(include)==0) {stop("\nParameter not found in any of the candidate models\n") }
  
      new.cand.set <- cand.set[which(include==1)] #select models including a given parameter
      new.mod.name <- modnames[which(include==1)]    #update model names
      
      new_table <- aictab(cand.set=new.cand.set, modnames=new.mod.name,
                          second.ord=second.ord, nobs=nobs, sort=FALSE)  #recompute AIC table and associated measures
      new_table$Beta_est <- unlist(lapply(new.cand.set, FUN=function(i) fixef(i)[paste(parm)])) #extract beta estimate for parm
      new_table$SE <- unlist(lapply(new.cand.set, FUN=function(i) extractSE(i)[paste(parm)]))
  
      ##compute model-averaged estimates, unconditional SE, and 95% CL
      if(second.ord==TRUE) {
        Modavg_beta <- sum(new_table$AICcWt*new_table$Beta_est)
        
        ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
        if(identical(uncond.se, "old")) {
          Uncond_SE <- sum(new_table$AICcWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
        }

        ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
        if(identical(uncond.se, "revised")) {
          Uncond_SE <- sqrt(sum(new_table$AICcWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
        } 


      } else {
        Modavg_beta<-sum(new_table$AICWt*new_table$Beta_est)
    
        ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
        if(identical(uncond.se, "old")) {
          Uncond_SE<-sum(new_table$AICWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
        }
    
        ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
        if(identical(uncond.se, "revised")) {
          Uncond_SE<-sqrt(sum(new_table$AICWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
        }
    
    
      }
  
      zcrit <- qnorm(p=(1-conf.level)/2, lower.tail=FALSE)
      Lower_CL <- Modavg_beta-zcrit*Uncond_SE
      Upper_CL <- Modavg_beta+zcrit*Uncond_SE
      out.modavg <- list("Parameter"=paste(parm), "Mod.avg.table" = new_table, "Mod.avg.beta" = Modavg_beta,
                         "Uncond.SE" = Uncond_SE, "Conf.level" = conf.level, "Lower.CL" = Lower_CL,
                         "Upper.CL" = Upper_CL)
  
      class(out.modavg) <- c("modavg", "list")
      return(out.modavg)
      
    }




##multinom
modavg.AICmultinom.nnet <-
  function(cand.set, parm, modnames = NULL, second.ord = TRUE,
           nobs = NULL, uncond.se = "revised", conf.level = 0.95,
           exclude = NULL, warn = TRUE, c.hat = 1, ...){

    ##check if named list if modnames are not supplied
    if(is.null(modnames)) {
      if(is.null(names(cand.set))) {
        modnames <- paste("Mod", 1:length(cand.set), sep = "")
        warning("\nModel names have been supplied automatically in the table\n")
      } else {
        modnames <- names(cand.set)
      }
}
    
  
#####MODIFICATIONS BEGIN#######
  ##remove all leading and trailing white space and within parm
  parm <- gsub('[[:space:]]+', "", parm)
  
  ##reverse parm
  reversed.parm <- reverse.parm(parm)
  exclude <- reverse.exclude(exclude = exclude)
#####MODIFICATIONS END######

  
  ##extract model formula for each model in cand.set    
  mod_formula<-lapply(cand.set, FUN=function(i) colnames(summary(i)$coefficients)) 

nmods <- length(cand.set)
  
  ##setup matrix to indicate presence of parm in the model
  include <- matrix(NA, nrow=nmods, ncol=1)
  ##add a check for multiple instances of same variable in given model (i.e., interactions)
  include.check <- matrix(NA, nrow=nmods, ncol=1)

  ##iterate over each formula in mod_formula list
  for (i in 1:nmods) {
    idents <- NULL
    idents.check <- NULL
    form <- mod_formula[[i]]

######################################################################################################
######################################################################################################
###MODIFICATIONS BEGIN
    ##iterate over each element of formula[[i]] in list
    if(is.null(reversed.parm)) {
      for (j in 1:length(form)) {
        idents[j] <- identical(parm, form[j])
        idents.check[j] <- ifelse(attr(regexpr(parm, form[j], fixed=TRUE), "match.length")== "-1" , 0, 1)  
      }
    } else {
      for (j in 1:length(form)) {
        idents[j] <- identical(parm, form[j]) | identical(reversed.parm, form[j])
        idents.check[j] <- ifelse(attr(regexpr(parm, form[j], fixed=TRUE), "match.length")=="-1" & attr(regexpr(reversed.parm, form[j],
                                                                fixed=TRUE), "match.length")=="-1" , 0, 1)  
      }
    }
###MODIFICATIONS END
######################################################################################################
######################################################################################################

    
    include[i] <- ifelse(any(idents==1), 1, 0)
    include.check[i] <- ifelse(sum(idents.check)>1, "duplicates", "OK")
  }

#####################################################
##exclude == NULL; warn=TRUE:  warn that duplicates occur and stop
  if(is.null(exclude) && identical(warn, TRUE)) {
    if(any(include.check == "duplicates")) {
      stop("\nSome models possibly include more than one instance of the parameter of interest.\n",
           "This may be due to the presence of interaction/polynomial terms, or variables\n",
           "with similar names:\n",
           "\tsee \"?modavg\" for details on variable specification and \"exclude\" argument\n")
    }

  }
  


  ##exclude == NULL; warn=FALSE:  compute model-averaged beta estimate from models including variable of interest,
  ##assuming that the variable is not involved in interaction or higher order polynomial (x^2, x^3, etc...),
  ##warn that models were not excluded
  if(is.null(exclude) && identical(warn, FALSE)) {
    if(any(include.check == "duplicates")) {
      warning("\nMultiple instances of parameter of interest in given model is presumably\n",
              "not due to interaction or polynomial terms - these models will not be\n",
              "excluded from the computation of model-averaged estimate\n")
    }
    
  }

  

  ##warn if exclude is neither a list nor NULL
  if(!is.null(exclude)) {
    if(!is.list(exclude)) {stop("\nItems in \"exclude\" must be specified as a list\n")}
  }


  
  ##if exclude is list  
  if(is.list(exclude)) {
    
    ##determine number of elements in exclude
    nexcl <- length(exclude)

    ##check each formula for presence of exclude variable extracted with formula( )  - in multinom( ) must be extracted from call   
    not.include <- lapply(cand.set, FUN=function(i) formula(i$call))

    ##set up a new list with model formula
    forms <- list()
    for (i in 1:nmods) {
      form.tmp <- strsplit(as.character(not.include[i]), split="~")[[1]][-1]
      if(attr(regexpr("\\+", form.tmp), "match.length")==-1) {
        forms[i] <- form.tmp
      } else {forms[i] <- strsplit(form.tmp, split=" \\+ ")}
    }

    ##additional check to see whether some variable names include "+"
    check.forms <- unlist(lapply(forms, FUN=function(i) any(attr(regexpr("\\+", i), "match.length")>0)[[1]]))
    if (any(check.forms==TRUE)) stop("\nPlease avoid \"+\" in variable names\n")

    ##additional check to determine if intercept was removed from models
    check.forms <- unlist(lapply(forms, FUN=function(i) any(attr(regexpr("\\- 1", i), "match.length")>0)[[1]]))
    if (any(check.forms==TRUE)) stop("\nModels without intercept are not supported in this version, please use alternative parameterization\n")

 
    ##search within formula for variables to exclude
    mod.exclude <- matrix(NA, nrow=nmods, ncol=nexcl)
    
    ##iterate over each element in exclude list
    for (var in 1:nexcl) {

      ##iterate over each formula in mod_formula list
      for (i in 1:nmods) {
        idents <- NULL
        form.excl <- forms[[i]]
                    
        ##iterate over each element of forms[[i]]
        for (j in 1:length(form.excl)) {
          idents[j] <- identical(exclude[var][[1]], gsub("(^ +)|( +$)", "", form.excl[j]))
        }
        mod.exclude[i,var] <- ifelse(any(idents==1), 1, 0)
      }

    }
  
    ##determine outcome across all variables to exclude
    to.exclude <- rowSums(mod.exclude)
  
  
    ##exclude models following models from model averaging  
    include[which(to.exclude>=1)] <- 0
      
    
  }
  

  
  ##add a check to determine if include always == 0
  if (sum(include)==0) {stop("\nParameter not found in any of the candidate models\n") }

  new.cand.set<-cand.set[which(include==1)] #select models including a given parameter
  new.mod.name<-modnames[which(include==1)]    #update model names
  ##


  ##determine number of levels - 1
  mod.levels <- lapply(cand.set, FUN=function(i) rownames(summary(i)$coefficients)) #extract level of response variable 
  check.levels <- unlist(unique(mod.levels))


  ##recompute AIC table and associated measures
  new_table<-aictab(cand.set=new.cand.set, modnames=new.mod.name,
                    second.ord=second.ord, nobs=nobs, sort=FALSE, c.hat=c.hat) 
  
  ##create object to store model-averaged estimate and SE's of k - 1 level of response
  out.est <- matrix(data=NA, nrow=length(check.levels), ncol=4)
  colnames(out.est) <- c("Mod.avg.est", "Uncond.SE", "Lower.CL", "Upper.CL")
  rownames(out.est) <- check.levels

  ##iterate over levels of response variable
  for (g in 1:length(check.levels)) {
    ##extract beta estimate for parm
    new_table$Beta_est <- unlist(lapply(new.cand.set,
                                        FUN=function(i) coef(i)[check.levels[g], paste(parm)]))
    ##extract SE of estimate for parm
    new_table$SE <- unlist(lapply(new.cand.set,
                                  FUN=function(i) sqrt(diag(vcov(i)))[paste(check.levels[g], ":",
                                    parm, sep="")]))

    ##if c-hat is estimated adjust the SE's by multiplying with sqrt of c-hat
    if(c.hat > 1) {new_table$SE<-new_table$SE*sqrt(c.hat)} 

    ##compute model-averaged estimates, unconditional SE, and 95% CL
    ##AICc
    if(c.hat == 1 && second.ord == TRUE) {
      Modavg_beta <- sum(new_table$AICcWt*new_table$Beta_est)
      
      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Uncond_SE <- sum(new_table$AICcWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
      }
      
      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Uncond_SE <- sqrt(sum(new_table$AICcWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
      }
    }


    ##QAICc
    ##if c-hat is estimated compute values accordingly and adjust table names
    if(c.hat > 1 && second.ord == TRUE) {
      Modavg_beta <- sum(new_table$QAICcWt*new_table$Beta_est)

      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Uncond_SE <- sum(new_table$QAICcWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
      }

      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Uncond_SE <- sqrt(sum(new_table$QAICcWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
      } 
    }

    
    ##AIC
    if(c.hat == 1 && second.ord == FALSE) {
      Modavg_beta <- sum(new_table$AICWt*new_table$Beta_est)

      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Uncond_SE <- sum(new_table$AICWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
      }

      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Uncond_SE <- sqrt(sum(new_table$AICWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
      }
    }
    

    ##QAIC
    ##if c-hat is estimated compute values accordingly and adjust table names  
    if(c.hat > 1 && second.ord == FALSE) {
      Modavg_beta <- sum(new_table$QAICWt*new_table$Beta_est)

      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Uncond_SE <- sum(new_table$QAICWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
      }

      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Uncond_SE <- sqrt(sum(new_table$QAICWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
      } 
    }
    

    out.est[g, 1] <- Modavg_beta
    out.est[g, 2] <- Uncond_SE
  }
     
  zcrit <- qnorm(p=(1-conf.level)/2, lower.tail=FALSE)
  out.est[,3] <- out.est[,1] - zcrit*out.est[,2]
  out.est[,4] <- out.est[,1] + zcrit*out.est[,2]
  out.modavg <- list("Parameter"=paste(parm), "Mod.avg.table" = new_table, "Mod.avg.beta" = out.est[,1],
                     "Uncond.SE" = out.est[,2], "Conf.level" = conf.level, "Lower.CL"= out.est[,3],
                     "Upper.CL" = out.est[,4])
  
  class(out.modavg) <- c("modavg", "list")
  return(out.modavg)
  
}



##polr
modavg.AICpolr <-
function(cand.set, parm, modnames = NULL, second.ord = TRUE,
         nobs = NULL, uncond.se = "revised", conf.level = 0.95, 
         exclude = NULL, warn = TRUE, ...){

  ##check if named list if modnames are not supplied
  if(is.null(modnames)) {
    if(is.null(names(cand.set))) {
      modnames <- paste("Mod", 1:length(cand.set), sep = "")
      warning("\nModel names have been supplied automatically in the table\n")
    } else {
      modnames <- names(cand.set)
    }
  }
    

#####MODIFICATIONS BEGIN#######
  ##remove all leading and trailing white space and within parm
  parm <- gsub('[[:space:]]+', "", parm)
  
  ##reverse parm
  reversed.parm <- reverse.parm(parm)
  exclude <- reverse.exclude(exclude = exclude)
#####MODIFICATIONS END######

  
  ##extract model formula for each model in cand.set    
  mod_formula<-lapply(cand.set, FUN=function(i) rownames(summary(i)$coefficients)) 

  nmods <- length(cand.set)

  ##setup matrix to indicate presence of parm in the model
  include <- matrix(NA, nrow=nmods, ncol=1)
  ##add a check for multiple instances of same variable in given model (i.e., interactions)
  include.check <- matrix(NA, nrow=nmods, ncol=1)

  ##iterate over each formula in mod_formula list
  for (i in 1:nmods) {
    idents <- NULL
    idents.check <- NULL
    form <- mod_formula[[i]]

######################################################################################################
######################################################################################################
###MODIFICATIONS BEGIN
    ##iterate over each element of formula[[i]] in list
      if(is.null(reversed.parm)) {
        for (j in 1:length(form)) {
          idents[j] <- identical(parm, form[j])
          idents.check[j] <- ifelse(attr(regexpr(parm, form[j], fixed=TRUE), "match.length")== "-1" , 0, 1)  
        }
      } else {
        for (j in 1:length(form)) {
          idents[j] <- identical(parm, form[j]) | identical(reversed.parm, form[j])
          idents.check[j] <- ifelse(attr(regexpr(parm, form[j], fixed=TRUE), "match.length")=="-1" & attr(regexpr(reversed.parm, form[j],
                                                                  fixed=TRUE), "match.length")=="-1" , 0, 1)  
        }
      }
###MODIFICATIONS END
######################################################################################################
######################################################################################################

    
    include[i] <- ifelse(any(idents==1), 1, 0)
    include.check[i] <- ifelse(sum(idents.check)>1, "duplicates", "OK")
  }

  #####################################################
  ##exclude == NULL; warn=TRUE:  warn that duplicates occur and stop
  if(is.null(exclude) && identical(warn, TRUE)) {
    if(any(include.check == "duplicates")) {
      stop("\nSome models possibly include more than one instance of the parameter of interest.\n",
           "This may be due to the presence of interaction/polynomial terms, or variables\n",
           "with similar names:\n",
           "\tsee \"?modavg\" for details on variable specification and \"exclude\" argument\n")
    }

  }

  ##exclude == NULL; warn=FALSE:  compute model-averaged beta estimate from models including variable of interest,
  ##assuming that the variable is not involved in interaction or higher order polynomial (x^2, x^3, etc...),
  ##warn that models were not excluded
  if(is.null(exclude) && identical(warn, FALSE)) {
    if(any(include.check == "duplicates")) {
      warning("Multiple instances of parameter of interest in given model is presumably\n",
              "not due to interaction or polynomial terms - these models will not be\n",
              "excluded from the computation of model-averaged estimate\n")
    }
    
  }



  ##warn if exclude is neither a list nor NULL
    if(!is.null(exclude)) {
      if(!is.list(exclude)) {stop("\nItems in \"exclude\" must be specified as a list\n")}
    }

  
  ##if exclude is list  
  if(is.list(exclude)) {

    ##determine number of elements in exclude
    nexcl <- length(exclude)

    ##check each formula for presence of exclude variable extracted with formula( ) - in polr( ) must be extracted from call  
    not.include <- lapply(cand.set, FUN=function(i) formula(i$call))

    
    ##set up a new list with model formula
    forms <- list()
    for (i in 1:nmods) {
      form.tmp <- strsplit(as.character(not.include[i]), split="~")[[1]][-1]
      if(attr(regexpr("\\+", form.tmp), "match.length")==-1) {
        forms[i] <- form.tmp
      } else {forms[i] <- strsplit(form.tmp, split=" \\+ ")}
    }

    ##additional check to see whether some variable names include "+"
    check.forms <- unlist(lapply(forms, FUN=function(i) any(attr(regexpr("\\+", i), "match.length")>0)[[1]]))
    if (any(check.forms==TRUE)) stop("\nPlease avoid \"+\" in variable names\n")

    ##additional check to determine if intercept was removed from models
    check.forms <- unlist(lapply(forms, FUN=function(i) any(attr(regexpr("\\- 1", i), "match.length")>0)[[1]]))
    if (any(check.forms==TRUE)) stop("\nModels without intercept are not supported in this version, please use alternative parameterization\n")
 

    ##search within formula for variables to exclude
    mod.exclude <- matrix(NA, nrow=nmods, ncol=nexcl)

    ##iterate over each element in exclude list
    for (var in 1:nexcl) {

      ##iterate over each formula in mod_formula list
      for (i in 1:nmods) {
        idents <- NULL
        form.excl <- forms[[i]]
                    
        ##iterate over each element of forms[[i]]
        for (j in 1:length(form.excl)) {
          idents[j] <- identical(exclude[var][[1]], gsub("(^ +)|( +$)", "", form.excl[j]))
        }
        mod.exclude[i,var] <- ifelse(any(idents==1), 1, 0)
      }
      
    }
  
    ##determine outcome across all variables to exclude
    to.exclude <- rowSums(mod.exclude)
  
  
    ##exclude models following models from model averaging  
    include[which(to.exclude>=1)] <- 0
      
    
  }


    
  ##add a check to determine if include always == 0
  if (sum(include)==0) {stop("\nParameter not found in any of the candidate models\n") }

  new.cand.set<-cand.set[which(include==1)] #select models including a given parameter
  new.mod.name<-modnames[which(include==1)]    #update model names


  new_table<-aictab(cand.set=new.cand.set, modnames=new.mod.name, 
                    second.ord=second.ord, nobs=nobs, sort=FALSE)  #recompute AIC table and associated measures
  
  ##add logical test to distinguish between intercepts and other coefs
  if(attr(regexpr(pattern = "\\|", text = parm), "match.length")==-1) {
  new_table$Beta_est<-unlist(lapply(new.cand.set, FUN=function(i) coef(i)[paste(parm)]))
} else {new_table$Beta_est <- unlist(lapply(new.cand.set, FUN=function(i) (i)$zeta[paste(parm)])) }
        
  ##extract beta estimate for parm
  new_table$SE<-unlist(lapply(new.cand.set, FUN=function(i) sqrt(diag(vcov(i)))[paste(parm)]))

  ##AICc
  ##compute model-averaged estimates, unconditional SE, and 95% CL based on AICc
  if(second.ord == TRUE) {
    Modavg_beta<-sum(new_table$AICcWt*new_table$Beta_est)
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Uncond_SE<-sum(new_table$AICcWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
    }

    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Uncond_SE<-sqrt(sum(new_table$AICcWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
    }   
  }

  
  ##AICc  
  ##compute model-averaged estimates, unconditional SE, and 95% CL based on AIC
  if(second.ord == FALSE) {
    Modavg_beta<-sum(new_table$AICWt*new_table$Beta_est)

    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Uncond_SE<-sum(new_table$AICWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
    }

    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Uncond_SE<-sqrt(sum(new_table$AICWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
    }
  }

  
  zcrit <- qnorm(p=(1-conf.level)/2, lower.tail=FALSE)
  Lower_CL<-Modavg_beta-zcrit*Uncond_SE
  Upper_CL<-Modavg_beta+zcrit*Uncond_SE
  out.modavg <- list("Parameter"=paste(parm), "Mod.avg.table" = new_table, "Mod.avg.beta" = Modavg_beta,
                     "Uncond.SE" = Uncond_SE, "Conf.level" = conf.level, "Lower.CL"= Lower_CL,
                     "Upper.CL" = Upper_CL)
  class(out.modavg) <- c("modavg", "list")
  return(out.modavg)
}



##rlm
modavg.AICrlm.lm <-
  function(cand.set, parm, modnames = NULL, second.ord = TRUE,
           nobs = NULL, uncond.se = "revised", conf.level = 0.95,
           exclude = NULL, warn = TRUE, ...){

    ##check if named list if modnames are not supplied
    if(is.null(modnames)) {
      if(is.null(names(cand.set))) {
        modnames <- paste("Mod", 1:length(cand.set), sep = "")
        warning("\nModel names have been supplied automatically in the table\n")
      } else {
        modnames <- names(cand.set)
      }
    }

    
#####MODIFICATIONS BEGIN#######
    ##remove all leading and trailing white space and within parm
    parm <- gsub('[[:space:]]+', "", parm)
  
    ##reverse parm
    reversed.parm <- reverse.parm(parm)
    exclude <- reverse.exclude(exclude = exclude)
#####MODIFICATIONS END######

  
    ##extract model formula for each model in cand.set
    mod_formula <- lapply(cand.set, FUN = function(i) rownames(summary(i)$coefficients))
       
    nmods <- length(cand.set)
  
    ##setup matrix to indicate presence of parms in the model
    include <- matrix(NA, nrow = nmods, ncol = 1)
    ##add a check for multiple instances of same variable in given model (i.e., interactions)
    include.check <- matrix(NA, nrow = nmods, ncol = 1)

    ##iterate over each formula in mod_formula list
    for (i in 1:nmods) {
      idents <- NULL
      idents.check <- NULL
      form <- mod_formula[[i]]

######################################################################################################
######################################################################################################
###MODIFICATIONS BEGIN
      ##iterate over each element of formula[[i]] in list
      if(is.null(reversed.parm)) {
        for (j in 1:length(form)) {
          idents[j] <- identical(parm, form[j])
          idents.check[j] <- ifelse(attr(regexpr(parm, form[j], fixed=TRUE), "match.length")== "-1" , 0, 1)  
        }
      } else {
        for (j in 1:length(form)) {
          idents[j] <- identical(parm, form[j]) | identical(reversed.parm, form[j])
          idents.check[j] <- ifelse(attr(regexpr(parm, form[j], fixed=TRUE), "match.length")=="-1" & attr(regexpr(reversed.parm, form[j],
                                                                  fixed=TRUE), "match.length")=="-1" , 0, 1)  
        }
      }
###MODIFICATIONS END
######################################################################################################
######################################################################################################

  
      include[i] <- ifelse(any(idents==1), 1, 0)
      include.check[i] <- ifelse(sum(idents.check)>1, "duplicates", "OK")
    }

  #####################################################
    ##exclude == NULL; warn=TRUE:  warn that duplicates occur and stop
    if(is.null(exclude) && identical(warn, TRUE)) {
      ##check for duplicates in same model
      if(any(include.check == "duplicates")) {
        stop("\nSome models include more than one instance of the parameter of interest. \n",
             "This may be due to the presence of interaction/polynomial terms, or variables\n",
             "with similar names:\n",
             "\tsee \"?modavg\" for details on variable specification and \"exclude\" argument\n")
      }

    }

    ##exclude == NULL; warn=FALSE:  compute model-averaged beta estimate from models including variable of interest,
    ##assuming that the variable is not involved in interaction or higher order polynomial (x^2, x^3, etc...),
    ##warn that models were not excluded
    if(is.null(exclude) && identical(warn, FALSE)) {
      if(any(include.check == "duplicates")) {
        warning("\nMultiple instances of parameter of interest in given model is presumably\n",
                "not due to interaction or polynomial terms - these models will not be\n",
                "excluded from the computation of model-averaged estimate\n")
      }
      
    }

    ##warn if exclude is neither a list nor NULL
    if(!is.null(exclude)) {
      if(!is.list(exclude)) {stop("\nItems in \"exclude\" must be specified as a list\n")}
    }


    ##if exclude is list  
    if(is.list(exclude)) {
      
      ##determine number of elements in exclude
      nexcl <- length(exclude)
      
      ##check each formula for presence of exclude variable extracted with formula( )  
      not.include <- lapply(cand.set, FUN = formula)

      ##set up a new list with model formula
      forms <- list()
      for (i in 1:nmods) {
        form.tmp <- strsplit(as.character(not.include[i]), split="~")[[1]][-1]
        if(attr(regexpr("\\+", form.tmp), "match.length")==-1) {
          forms[i] <- form.tmp
        } else {forms[i] <- strsplit(form.tmp, split=" \\+ ")}
      }

      ##additional check to see whether some variable names include "+"
      check.forms <- unlist(lapply(forms, FUN=function(i) any(attr(regexpr("\\+", i), "match.length")>0)[[1]]))
      if (any(check.forms==TRUE)) stop("\nPlease avoid \"+\" in variable names\n")

      ##additional check to determine if intercept was removed from models
      check.forms <- unlist(lapply(forms, FUN=function(i) any(attr(regexpr("\\- 1", i), "match.length")>0)[[1]]))
      if (any(check.forms==TRUE)) stop("\nModels without intercept are not supported in this version, please use alternative parameterization\n")
         
      
      ##search within formula for variables to exclude
      mod.exclude <- matrix(NA, nrow=nmods, ncol=nexcl)
      
      ##iterate over each element in exclude list
      for (var in 1:nexcl) {

        ##iterate over each formula in mod_formula list
        for (i in 1:nmods) {
          idents <- NULL
          form.excl <- forms[[i]]
                    
          ##iterate over each element of forms[[i]]
          for (j in 1:length(form.excl)) {
            idents[j] <- identical(exclude[var][[1]], gsub("(^ +)|( +$)", "", form.excl[j]))
          }
          mod.exclude[i,var] <- ifelse(any(idents==1), 1, 0)
        }    
        
      }
  
      ##determine outcome across all variables to exclude
      to.exclude <- rowSums(mod.exclude)
         
      
      ##exclude models following models from model averaging  
      include[which(to.exclude>=1)] <- 0
      
      
    }

       
       
    ##add a check to determine if include always == 0
    if (sum(include)==0) {stop("\nParameter not found in any of the candidate models\n") }
    
    new.cand.set <- cand.set[which(include==1)] #select models including a given parameter
    new.mod.name <- modnames[which(include==1)]    #update model names
    
    new_table <- aictab(cand.set=new.cand.set, modnames=new.mod.name, 
                        second.ord=second.ord, nobs=nobs, sort=FALSE)  #recompute AIC table and associated measures
    new_table$Beta_est <- unlist(lapply(new.cand.set, FUN=function(i) coef(i)[paste(parm)])) #extract beta estimate for parm
    new_table$SE <- unlist(lapply(new.cand.set, FUN=function(i) sqrt(diag(vcov(i)))[paste(parm)]))
  
    ##compute model-averaged estimates, unconditional SE, and 95% CL
    if(second.ord == TRUE) {
      Modavg_beta <- sum(new_table$AICcWt*new_table$Beta_est)
         
      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Uncond_SE <- sum(new_table$AICcWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
      }

      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Uncond_SE <- sqrt(sum(new_table$AICcWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
      } 


    } else {
      Modavg_beta <- sum(new_table$AICWt*new_table$Beta_est)
      
      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Uncond_SE <- sum(new_table$AICWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
      }

      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Uncond_SE <- sqrt(sum(new_table$AICWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
      }
    
      
    }
     
    zcrit <- qnorm(p=(1-conf.level)/2, lower.tail=FALSE)
    Lower_CL <- Modavg_beta - zcrit*Uncond_SE
    Upper_CL <- Modavg_beta + zcrit*Uncond_SE
    out.modavg <- list("Parameter"=paste(parm), "Mod.avg.table" = new_table, "Mod.avg.beta" = Modavg_beta,
                       "Uncond.SE" = Uncond_SE, "Conf.level" = conf.level, "Lower.CL" = Lower_CL,
                       "Upper.CL" = Upper_CL)
    
    class(out.modavg) <- c("modavg", "list")
    return(out.modavg)
  }



##survreg
modavg.AICsurvreg <-
  function(cand.set, parm, modnames = NULL, second.ord = TRUE,
           nobs = NULL,  uncond.se = "revised", conf.level = 0.95, 
           exclude = NULL, warn = TRUE, ...){

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
    check.dist <- sapply(X = cand.set, FUN = function(i) i$dist)
    unique.dist <- unique(x = check.dist)
    if(length(unique.dist) > 1) stop("\nFunction does not support model-averaging estimates from different distributions\n")
  

    

#####MODIFICATIONS BEGIN#######
    ##remove all leading and trailing white space and within parm
    parm <- gsub('[[:space:]]+', "", parm)
     
    ##reverse parm
    reversed.parm <- reverse.parm(parm)
    exclude <- reverse.exclude(exclude = exclude)
#####MODIFICATIONS END######

    ##extract model formula for each model in cand.set
    mod_formula <- lapply(cand.set, FUN = function(i) names(summary(i)$coefficients)) #extract model formula for each model in cand.set
    
    nmods <- length(cand.set)
  
    ##setup matrix to indicate presence of parm in the model
    include <- matrix(NA, nrow=nmods, ncol=1)
    ##add a check for multiple instances of same variable in given model (i.e., interactions)
    include.check <- matrix(NA, nrow=nmods, ncol=1)

    ##iterate over each formula in mod_formula list
    for (i in 1:nmods) {
      idents <- NULL
      idents.check <- NULL
      form <- mod_formula[[i]]

######################################################################################################
######################################################################################################
###MODIFICATIONS BEGIN
      ##iterate over each element of formula[[i]] in list
      if(is.null(reversed.parm)) {
        for (j in 1:length(form)) {
          idents[j] <- identical(parm, form[j])
          idents.check[j] <- ifelse(attr(regexpr(parm, form[j], fixed=TRUE), "match.length")== "-1" , 0, 1)  
        }
      } else {
        for (j in 1:length(form)) {
          idents[j] <- identical(parm, form[j]) | identical(reversed.parm, form[j])
          idents.check[j] <- ifelse(attr(regexpr(parm, form[j], fixed=TRUE), "match.length")=="-1" & attr(regexpr(reversed.parm, form[j],
                                                                  fixed=TRUE), "match.length")=="-1" , 0, 1)  
        }
      }
###MODIFICATIONS END
######################################################################################################
######################################################################################################


      include[i] <- ifelse(any(idents==1), 1, 0)
      include.check[i] <- ifelse(sum(idents.check)>1, "duplicates", "OK")
    }

#####################################################
    ##exclude == NULL; warn=TRUE:  warn that duplicates occur and stop
    if(is.null(exclude) && identical(warn, TRUE)) {
      if(any(include.check == "duplicates")) {
        stop("\nSome models possibly include more than one instance of the parameter of interest.\n",
             "This may be due to the presence of interaction/polynomial terms, or variables\n",
             "with similar names:\n",
             "\tsee \"?modavg\" for details on variable specification and \"exclude\" argument\n")
      }

    }

    ##exclude == NULL; warn=FALSE:  compute model-averaged beta estimate from models including variable of interest,
    ##assuming that the variable is not involved in interaction or higher order polynomial (x^2, x^3, etc...),
    ##warn that models were not excluded
    if(is.null(exclude) && identical(warn, FALSE)) {
      if(any(include.check == "duplicates")) {
        warning("Multiple instances of parameter of interest in given model is presumably\n",
                "not due to interaction or polynomial terms - these models will not be\n",
                "excluded from the computation of model-averaged estimate\n")
      }
      
    }



    ##warn if exclude is neither a list nor NULL
    if(!is.null(exclude)) {
      if(!is.list(exclude)) {stop("\nItems in \"exclude\" must be specified as a list\n")}
    }

  
    ##if exclude is list  
    if(is.list(exclude)) {

      ##determine number of elements in exclude
      nexcl <- length(exclude)

      ##check each formula for presence of exclude variable extracted with formula( )  
      not.include <- lapply(cand.set, FUN=formula)

      ##set up a new list with model formula
      forms <- list()
      for (i in 1:nmods) {
        form.tmp <- strsplit(as.character(not.include[i]), split="~")[[1]][-1]
        if(attr(regexpr("\\+", form.tmp), "match.length")==-1) {  #this line causes problems if intercept is removed from model
          forms[i] <- form.tmp
        } else {forms[i] <- strsplit(form.tmp, split=" \\+ ")}    #this line causes problems if intercept is removed from model
      }

      ##additional check to see whether some variable names include "+"
      check.forms <- unlist(lapply(forms, FUN=function(i) any(attr(regexpr("\\+", i), "match.length")>0)[[1]]))
      if (any(check.forms==TRUE)) stop("\nPlease avoid \"+\" in variable names\n")

      ##additional check to determine if intercept was removed from models
      check.forms <- unlist(lapply(forms, FUN=function(i) any(attr(regexpr("\\- 1", i), "match.length")>0)[[1]]))
      if (any(check.forms==TRUE)) stop("\nModels without intercept are not supported in this version, please use alternative parameterization\n")

 
      ##search within formula for variables to exclude
      mod.exclude <- matrix(NA, nrow=nmods, ncol=nexcl)

      ##iterate over each element in exclude list
      for (var in 1:nexcl) {

        ##iterate over each formula in mod_formula list
        for (i in 1:nmods) {
          idents <- NULL
          form.excl <- forms[[i]]
                    
          ##iterate over each element of forms[[i]]
          for (j in 1:length(form.excl)) {
            idents[j] <- identical(exclude[var][[1]], gsub("(^ +)|( +$)", "", form.excl[j]))
          }
          mod.exclude[i,var] <- ifelse(any(idents==1), 1, 0)
        }    
        
      }
  
      ##determine outcome across all variables to exclude
      to.exclude <- rowSums(mod.exclude)
  
  
      ##exclude models following models from model averaging  
      include[which(to.exclude>=1)] <- 0
      
      
    }



    ##add a check to determine if include always == 0
    if (sum(include)==0) {stop("\nParameter not found in any of the candidate models\n") }

    new.cand.set <- cand.set[which(include==1)] #select models including a given parameter
    new.mod.name <- modnames[which(include==1)]    #update model names
    ##

    new_table <- aictab(cand.set = new.cand.set, modnames = new.mod.name,
                        second.ord = second.ord, nobs = nobs, sort = FALSE)  #recompute AIC table and associated measures
    new_table$Beta_est <- unlist(lapply(new.cand.set, FUN = function(i) coef(i)[paste(parm)])) #extract beta estimate for parm
    new_table$SE <- unlist(lapply(new.cand.set, FUN = function(i) sqrt(diag(vcov(i)))[paste(parm)]))

   
    
    ##AICc
    ##compute model-averaged estimates, unconditional SE, and 95% CL
    if(second.ord == TRUE) {
      Modavg_beta <- sum(new_table$AICcWt*new_table$Beta_est)

      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Uncond_SE <- sum(new_table$AICcWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
      }
      
      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Uncond_SE <- sqrt(sum(new_table$AICcWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
      }
    }

    
    
    ##AIC
    if(second.ord == FALSE) {
      Modavg_beta <- sum(new_table$AICWt*new_table$Beta_est)
      
      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Uncond_SE <- sum(new_table$AICWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
      }
      
      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Uncond_SE <- sqrt(sum(new_table$AICWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
      }
    }

    
    zcrit <- qnorm(p=(1-conf.level)/2, lower.tail=FALSE)
    Lower_CL <- Modavg_beta-zcrit*Uncond_SE
    Upper_CL <- Modavg_beta+zcrit*Uncond_SE
    out.modavg <- list("Parameter"=paste(parm), "Mod.avg.table" = new_table, "Mod.avg.beta" = Modavg_beta,
                       "Uncond.SE" = Uncond_SE, "Conf.level" = conf.level, "Lower.CL"= Lower_CL, "Upper.CL" = Upper_CL)
  
    
    class(out.modavg) <- c("modavg", "list")
    return(out.modavg)

  }



##vglm
modavg.AICvglm <-
function(cand.set, parm, modnames = NULL, second.ord = TRUE,
         nobs = NULL, uncond.se = "revised", conf.level = 0.95,
         exclude = NULL, warn = TRUE, c.hat = 1, ...){

  ##check if named list if modnames are not supplied
  if(is.null(modnames)) {
    if(is.null(names(cand.set))) {
      modnames <- paste("Mod", 1:length(cand.set), sep = "")
      warning("\nModel names have been supplied automatically in the table\n")
    } else {
      modnames <- names(cand.set)
    }
  }
    
  ##check that link function is the same for all models
  check.link <- unlist(lapply(X = cand.set, FUN=function(i) i@family@blurb[3])) 
    unique.link <- unique(x=check.link)
    if(length(unique.link) > 1) stop("\nIt is not appropriate to compute a model-averaged beta estimate\n",
                                         "from models using different link functions\n")
    
    ##check family of vglm to avoid problems when requesting predictions with argument 'dispersion'
    fam.type <- unlist(lapply(cand.set, FUN=function(i) i@family@vfamily))
    fam.unique <- unique(fam.type)
    if(identical(fam.unique, "gaussianff")) {disp <- NULL} else{disp <- 1}
    if(identical(fam.unique, "gammaff")) stop("\nGamma distribution is not supported yet\n")
  ##poisson and binomial defaults to 1 (no separate parameter for variance)
  ##for negative binomial - reset to NULL
  if(identical(fam.unique, "negbinomial")) {disp <- NULL}


#####MODIFICATIONS BEGIN#######
    ##remove all leading and trailing white space and within parm
    parm <- gsub('[[:space:]]+', "", parm)
     
    ##reverse parm
    reversed.parm <- reverse.parm(parm)
    exclude <- reverse.exclude(exclude = exclude)
#####MODIFICATIONS END######

    ##extract model formula for each model in cand.set
    mod_formula <- lapply(cand.set, FUN=function(i) labels(coefficients(i))) #extract model formula for each model in cand.set

    nmods <- length(cand.set)

    ##check whether parm is involved in interaction or if label changes for some models - e.g., ZIP models
    ##if : not already included
    if(regexpr(":", parm, fixed = TRUE) == -1){

      ##if : not included
      parm.inter <- c(paste(parm, ":", sep = ""), paste(":", parm, sep = ""))
      inter.check <- ifelse(attr(regexpr(parm.inter[1], mod_formula, fixed = TRUE), "match.length") == "-1" & attr(regexpr(parm.inter[2],
                                                                       mod_formula, fixed = TRUE), "match.length") == "-1", 0, 1)
      if(sum(inter.check) > 0) warning("\nLabel of parameter of interest seems to change across models:\n",
                                       "check model syntax for possible problems\n")
    } else {
      ##if : already included
      ##remove : from parm
      simple.parm <- unlist(strsplit(parm, split = ":"))[1]

      ##search for simple.parm and parm in model formulae
      no.colon <- sum(ifelse(attr(regexpr(simple.parm, mod_formula, fixed = TRUE), "match.length") != "-1", 1, 0))
      with.colon <- sum(ifelse(attr(regexpr(parm, mod_formula, fixed = TRUE), "match.length") != "-1", 0, 1))

      ##check if both are > 0
      if(no.colon > 0 && with.colon > 0) warning("\nLabel of parameter of interest seems to change across models:\n",
                                                 "check model syntax for possible problems\n")
    }

    
    ##setup matrix to indicate presence of parm in the model
    include <- matrix(NA, nrow=nmods, ncol=1)
    ##add a check for multiple instances of same variable in given model (i.e., interactions)
    include.check <- matrix(NA, nrow=nmods, ncol=1)

    ##iterate over each formula in mod_formula list
    for (i in 1:nmods) {
      idents <- NULL
      idents.check <- NULL
      form <- mod_formula[[i]]

######################################################################################################
######################################################################################################
###MODIFICATIONS BEGIN
      ##iterate over each element of formula[[i]] in list
      if(is.null(reversed.parm)) {
        for (j in 1:length(form)) {
          idents[j] <- identical(parm, form[j])
          idents.check[j] <- ifelse(attr(regexpr(parm, form[j], fixed=TRUE), "match.length")== "-1" , 0, 1)  
        }
      } else {
        for (j in 1:length(form)) {
          idents[j] <- identical(parm, form[j]) | identical(reversed.parm, form[j])
          idents.check[j] <- ifelse(attr(regexpr(parm, form[j], fixed=TRUE), "match.length")=="-1" & attr(regexpr(reversed.parm, form[j],
                                                                  fixed=TRUE), "match.length")=="-1" , 0, 1)  
        }
      }
###MODIFICATIONS END
######################################################################################################
######################################################################################################


      
      include[i] <- ifelse(any(idents==1), 1, 0)
      include.check[i] <- ifelse(sum(idents.check)>1, "duplicates", "OK")
    }

#####################################################
    ##exclude == NULL; warn=TRUE:  warn that duplicates occur and stop
    if(is.null(exclude) && identical(warn, TRUE)) {
      if(any(include.check == "duplicates")) {
        stop("\nSome models possibly include more than one instance of the parameter of interest.\n",
             "This may be due to the presence of interaction/polynomial terms, or variables\n",
             "with similar names:\n",
             "\tsee \"?modavg\" for details on variable specification and \"exclude\" argument\n")
      }

    }

    ##exclude == NULL; warn=FALSE:  compute model-averaged beta estimate from models including variable of interest,
    ##assuming that the variable is not involved in interaction or higher order polynomial (x^2, x^3, etc...),
    ##warn that models were not excluded
    if(is.null(exclude) && identical(warn, FALSE)) {
      if(any(include.check == "duplicates")) {
        warning("Multiple instances of parameter of interest in given model is presumably\n",
                "not due to interaction or polynomial terms - these models will not be\n",
                "excluded from the computation of model-averaged estimate\n")
      }
      
    }



    ##warn if exclude is neither a list nor NULL
    if(!is.null(exclude)) {
      if(!is.list(exclude)) {stop("\nItems in \"exclude\" must be specified as a list\n")}
    }

  
    ##if exclude is list  
    if(is.list(exclude)) {

      ##determine number of elements in exclude
      nexcl <- length(exclude)

      ##check each formula for presence of exclude variable extracted with formula( )  
      not.include <- lapply(cand.set, FUN = formula)

      ##set up a new list with model formula
      forms <- list()
      for (i in 1:nmods) {
        form.tmp <- strsplit(as.character(not.include[i]), split="~")[[1]][-1]
        if(attr(regexpr("\\+", form.tmp), "match.length")==-1) {  #this line causes problems if intercept is removed from model
          forms[i] <- form.tmp
        } else {forms[i] <- strsplit(form.tmp, split=" \\+ ")}    #this line causes problems if intercept is removed from model
      }

      ##additional check to see whether some variable names include "+"
      check.forms <- unlist(lapply(forms, FUN=function(i) any(attr(regexpr("\\+", i), "match.length")>0)[[1]]))
      if (any(check.forms==TRUE)) stop("\nPlease avoid \"+\" in variable names\n")

      ##additional check to determine if intercept was removed from models
      check.forms <- unlist(lapply(forms, FUN=function(i) any(attr(regexpr("\\- 1", i), "match.length")>0)[[1]]))
      if (any(check.forms==TRUE)) stop("\nModels without intercept are not supported in this version, please use alternative parameterization\n")

 
      ##search within formula for variables to exclude
      mod.exclude <- matrix(NA, nrow=nmods, ncol=nexcl)

      ##iterate over each element in exclude list
      for (var in 1:nexcl) {

        ##iterate over each formula in mod_formula list
        for (i in 1:nmods) {
          idents <- NULL
          form.excl <- forms[[i]]
                    
          ##iterate over each element of forms[[i]]
          for (j in 1:length(form.excl)) {
            idents[j] <- identical(exclude[var][[1]], gsub("(^ +)|( +$)", "", form.excl[j]))
          }
          mod.exclude[i,var] <- ifelse(any(idents==1), 1, 0)
        }    
        
      }
  
      ##determine outcome across all variables to exclude
      to.exclude <- rowSums(mod.exclude)
  
  
      ##exclude models following models from model averaging  
      include[which(to.exclude>=1)] <- 0
      
      
    }



    ##add a check to determine if include always == 0
    if (sum(include)==0) {stop("\nParameter not found in any of the candidate models\n") }

    new.cand.set <- cand.set[which(include==1)] #select models including a given parameter
    new.mod.name <- modnames[which(include==1)]    #update model names
    ##

    new_table <- aictab(cand.set = new.cand.set, modnames = new.mod.name,
                        second.ord = second.ord, nobs = nobs, sort = FALSE, c.hat = c.hat)  #recompute AIC table and associated measures
    new_table$Beta_est <- unlist(lapply(new.cand.set, FUN=function(i) coefficients(i)[paste(parm)])) #extract beta estimate for parm
    new_table$SE <- unlist(lapply(new.cand.set, FUN=function(i) sqrt(diag(vcov(i, dispersion = disp)))[paste(parm)]))

   
    ##if c-hat is estimated adjust the SE's by multiplying with sqrt of c-hat
    if(c.hat > 1) {
      new_table$SE <- new_table$SE*sqrt(c.hat)
    } 

    

    ##AICc
    ##compute model-averaged estimates, unconditional SE, and 95% CL
    if(c.hat == 1 && second.ord == TRUE) {
      Modavg_beta <- sum(new_table$AICcWt*new_table$Beta_est)

      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Uncond_SE <- sum(new_table$AICcWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
      }
      
      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Uncond_SE <- sqrt(sum(new_table$AICcWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
      }
    }

    
    
    ##QAICc
    if(c.hat > 1 && second.ord == TRUE) {
      Modavg_beta <- sum(new_table$QAICcWt*new_table$Beta_est)
      
      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {      
        Uncond_SE <- sum(new_table$QAICcWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
      }
      
      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Uncond_SE <- sqrt(sum(new_table$QAICcWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
      }
    }     



    ##AIC
    if(c.hat == 1 && second.ord == FALSE) {
      Modavg_beta <- sum(new_table$AICWt*new_table$Beta_est)
      
      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Uncond_SE <- sum(new_table$AICWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
      }
      
      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Uncond_SE <- sqrt(sum(new_table$AICWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
      }
    }



    ##QAIC
    if(c.hat > 1 && second.ord == FALSE) {
      Modavg_beta <- sum(new_table$QAICWt*new_table$Beta_est)

      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Uncond_SE <- sum(new_table$QAICWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
      }

      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Uncond_SE <- sqrt(sum(new_table$QAICWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
      }  
    }     

    zcrit <- qnorm(p=(1-conf.level)/2, lower.tail=FALSE)
    Lower_CL <- Modavg_beta-zcrit*Uncond_SE
    Upper_CL <- Modavg_beta+zcrit*Uncond_SE
    out.modavg <- list("Parameter"=paste(parm), "Mod.avg.table" = new_table, "Mod.avg.beta" = Modavg_beta,
                       "Uncond.SE" = Uncond_SE, "Conf.level" = conf.level, "Lower.CL"= Lower_CL, "Upper.CL" = Upper_CL)

  class(out.modavg) <- c("modavg", "list")
  return(out.modavg)

}




##zeroinfl
modavg.AICzeroinfl <-
function(cand.set, parm, modnames = NULL, second.ord = TRUE,
         nobs = NULL, uncond.se = "revised", conf.level = 0.95,
         exclude = NULL, warn = TRUE, ...){

  ##check if named list if modnames are not supplied
  if(is.null(modnames)) {
    if(is.null(names(cand.set))) {
      modnames <- paste("Mod", 1:length(cand.set), sep = "")
      warning("\nModel names have been supplied automatically in the table\n")
    } else {
      modnames <- names(cand.set)
    }
  }
    
  
  ##check that link function is the same for all models
  check.link <- unlist(lapply(X = cand.set, FUN=function(i) i$link))
  unique.link <- unique(x=check.link)
  if(length(unique.link) > 1) stop("\nIt is not appropriate to compute a model-averaged beta estimate\n",
"from models using different link functions\n")
  
#####MODIFICATIONS BEGIN#######
    ##remove all leading and trailing white space and within parm
    parm <- gsub('[[:space:]]+', "", parm)
     
    ##reverse parm
    reversed.parm <- reverse.parm(parm)
    exclude <- reverse.exclude(exclude = exclude)
#####MODIFICATIONS END######

    ##extract model formula for each model in cand.set
    mod_formula <- lapply(cand.set, FUN=function(i) labels(coefficients(i))) #extract model formula for each model in cand.set

    nmods <- length(cand.set)
  
    ##setup matrix to indicate presence of parm in the model
    include <- matrix(NA, nrow=nmods, ncol=1)
    ##add a check for multiple instances of same variable in given model (i.e., interactions)
    include.check <- matrix(NA, nrow=nmods, ncol=1)

    ##iterate over each formula in mod_formula list
    for (i in 1:nmods) {
      idents <- NULL
      idents.check <- NULL
      form <- mod_formula[[i]]

######################################################################################################
######################################################################################################
###MODIFICATIONS BEGIN
      ##iterate over each element of formula[[i]] in list
      if(is.null(reversed.parm)) {
        for (j in 1:length(form)) {
          idents[j] <- identical(parm, form[j])
          idents.check[j] <- ifelse(attr(regexpr(parm, form[j], fixed=TRUE), "match.length")== "-1" , 0, 1)  
        }
      } else {
        for (j in 1:length(form)) {
          idents[j] <- identical(parm, form[j]) | identical(reversed.parm, form[j])
          idents.check[j] <- ifelse(attr(regexpr(parm, form[j], fixed=TRUE), "match.length")=="-1" & attr(regexpr(reversed.parm, form[j],
                                                                  fixed=TRUE), "match.length")=="-1" , 0, 1)  
        }
      }
###MODIFICATIONS END
######################################################################################################
######################################################################################################


      
      include[i] <- ifelse(any(idents==1), 1, 0)
      include.check[i] <- ifelse(sum(idents.check)>1, "duplicates", "OK")
    }

#####################################################
    ##exclude == NULL; warn=TRUE:  warn that duplicates occur and stop
    if(is.null(exclude) && identical(warn, TRUE)) {
      if(any(include.check == "duplicates")) {
        stop("\nSome models possibly include more than one instance of the parameter of interest.\n",
             "This may be due to the presence of interaction/polynomial terms, or variables\n",
             "with similar names:\n",
             "\tsee \"?modavg\" for details on variable specification and \"exclude\" argument\n")
      }

    }

    ##exclude == NULL; warn=FALSE:  compute model-averaged beta estimate from models including variable of interest,
    ##assuming that the variable is not involved in interaction or higher order polynomial (x^2, x^3, etc...),
    ##warn that models were not excluded
    if(is.null(exclude) && identical(warn, FALSE)) {
      if(any(include.check == "duplicates")) {
        warning("Multiple instances of parameter of interest in given model is presumably\n",
                "not due to interaction or polynomial terms - these models will not be\n",
                "excluded from the computation of model-averaged estimate\n")
      }
      
    }



    ##warn if exclude is neither a list nor NULL
    if(!is.null(exclude)) {
      if(!is.list(exclude)) {stop("\nItems in \"exclude\" must be specified as a list\n")}
    }

  
    ##if exclude is list  
    if(is.list(exclude)) {

      ##determine number of elements in exclude
      nexcl <- length(exclude)

      ##check each formula for presence of exclude variable extracted with formula( )  
      not.include <- lapply(cand.set, FUN=formula)

      ##set up a new list with model formula
      forms <- list()
      for (i in 1:nmods) {
        form.tmp <- strsplit(as.character(not.include[i]), split="~")[[1]][-1]
        if(attr(regexpr("\\+", form.tmp), "match.length")==-1) {  #this line causes problems if intercept is removed from model
          forms[i] <- form.tmp
        } else {forms[i] <- strsplit(form.tmp, split=" \\+ ")}    #this line causes problems if intercept is removed from model
      }

      ##additional check to see whether some variable names include "+"
      check.forms <- unlist(lapply(forms, FUN=function(i) any(attr(regexpr("\\+", i), "match.length")>0)[[1]]))
      if (any(check.forms==TRUE)) stop("\nPlease avoid \"+\" in variable names\n")

      ##additional check to determine if intercept was removed from models
      check.forms <- unlist(lapply(forms, FUN=function(i) any(attr(regexpr("\\- 1", i), "match.length")>0)[[1]]))
      if (any(check.forms==TRUE)) stop("\nModels without intercept are not supported in this version, please use alternative parameterization\n")

 
      ##search within formula for variables to exclude
      mod.exclude <- matrix(NA, nrow=nmods, ncol=nexcl)

      ##iterate over each element in exclude list
      for (var in 1:nexcl) {

        ##iterate over each formula in mod_formula list
        for (i in 1:nmods) {
          idents <- NULL
          form.excl <- forms[[i]]
                    
          ##iterate over each element of forms[[i]]
          for (j in 1:length(form.excl)) {
            idents[j] <- identical(exclude[var][[1]], gsub("(^ +)|( +$)", "", form.excl[j]))
          }
          mod.exclude[i,var] <- ifelse(any(idents==1), 1, 0)
        }    
        
      }
  
      ##determine outcome across all variables to exclude
      to.exclude <- rowSums(mod.exclude)
  
  
      ##exclude models following models from model averaging  
      include[which(to.exclude>=1)] <- 0
      
      
    }



    ##add a check to determine if include always == 0
    if (sum(include)==0) {stop("\nParameter not found in any of the candidate models\n") }

    new.cand.set <- cand.set[which(include==1)] #select models including a given parameter
    new.mod.name <- modnames[which(include==1)]    #update model names
    ##

    new_table <- aictab(cand.set = new.cand.set, modnames = new.mod.name, 
                        second.ord = second.ord, nobs = nobs, sort = FALSE)  #recompute AIC table and associated measures
    new_table$Beta_est <- unlist(lapply(new.cand.set, FUN=function(i) coefficients(i)[paste(parm)])) #extract beta estimate for parm
    new_table$SE <- unlist(lapply(new.cand.set, FUN=function(i) sqrt(diag(vcov(i)))[paste(parm)]))

   
    ##AICc
    ##compute model-averaged estimates, unconditional SE, and 95% CL
    if(second.ord == TRUE) {
      Modavg_beta <- sum(new_table$AICcWt*new_table$Beta_est)

      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Uncond_SE <- sum(new_table$AICcWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
      }
      
      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Uncond_SE <- sqrt(sum(new_table$AICcWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
      }
    }

    
    
    ##AIC
    if(second.ord == FALSE) {
      Modavg_beta <- sum(new_table$AICWt*new_table$Beta_est)
      
      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Uncond_SE <- sum(new_table$AICWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
      }
      
      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Uncond_SE <- sqrt(sum(new_table$AICWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
      }
    }


    zcrit <- qnorm(p=(1-conf.level)/2, lower.tail=FALSE)
    Lower_CL <- Modavg_beta-zcrit*Uncond_SE
    Upper_CL <- Modavg_beta+zcrit*Uncond_SE
    out.modavg <- list("Parameter"=paste(parm), "Mod.avg.table" = new_table, "Mod.avg.beta" = Modavg_beta,
                       "Uncond.SE" = Uncond_SE, "Conf.level" = conf.level, "Lower.CL"= Lower_CL, "Upper.CL" = Upper_CL)
    
  class(out.modavg) <- c("modavg", "list")
  return(out.modavg)

}




####added functionality for reversing parameters
##added additional argument parm.type = "psi", "gamma", "epsilon", "lambda", "omega", "detect"
##model type:  parameters labeled in unmarked - parameters labeled in AICcmodavg.unmarked
##single season: state, det - USE psi, detect
##multiseason model:  psi, col, ext, det - USE psi, gamma, epsilon, detect
##RN heterogeneity model: state, det - USE lambda, detect
##N-mixture: state, det - USE lambda, detect
##Open N-mixture: lambda, gamma, omega, det - USE lambda, gamma, omega, detect
##distsamp: state, det - USE lambda, detect
##gdistsamp: state, det, phi - USE lambda, detect, phi
##false-positive occupancy: state, det, fp - USE psi, detect, fp
##gpcount: lambda, phi, det - USE lambda, phi, detect
##gmultmix: lambda, phi, det - USE lambda, phi, detect
##multinomPois: state, det - USE lambda, detect

##occu
modavg.AICunmarkedFitOccu <-
function(cand.set, parm, modnames = NULL, second.ord = TRUE, nobs = NULL, 
         uncond.se = "revised", conf.level = 0.95, exclude = NULL, warn = TRUE,
         c.hat = 1, parm.type = NULL, ...){

  ##note that parameter is referenced differently from unmarked object - see labels( )

  ##check if named list if modnames are not supplied
  if(is.null(modnames)) {
    if(is.null(names(cand.set))) {
      modnames <- paste("Mod", 1:length(cand.set), sep = "")
      warning("\nModel names have been supplied automatically in the table\n")
    } else {
      modnames <- names(cand.set)
    }
  }
    
  
  ##check for parm.type and stop if NULL
  if(is.null(parm.type)) {stop("\n'parm.type' must be specified for this model type, see ?modavg for details\n")}
  
#####MODIFICATIONS BEGIN#######
  ##remove all leading and trailing white space and within parm
  parm <- gsub('[[:space:]]+', "", parm)
  parm.strip <- parm #to use later
  
  ##if (Intercept) is chosen assign (Int) - for compatibility
  if(identical(parm, "(Intercept)")) parm <- "Int"

  ##reverse parm
  reversed.parm <- reverse.parm(parm)
  reversed.parm.strip <- reversed.parm #to use later
  exclude <- reverse.exclude(exclude = exclude)
#####MODIFICATIONS END######
  
    ##single-season occupancy model
    ##psi
    if(identical(parm.type, "psi")) {
      ##extract model formula for each model in cand.set
      mod_formula <- lapply(cand.set, FUN = function(i) labels(coef(i@estimates@estimates$state)))
      parm <- paste("psi", "(", parm, ")", sep="")
      if(!is.null(reversed.parm)) {reversed.parm <- paste("psi", "(", reversed.parm, ")", sep="")}
      not.include <- lapply(cand.set, FUN = function(i) i@formula[[3]])      
    }
  
    ##detect
    if(identical(parm.type, "detect")) {
      mod_formula <- lapply(cand.set, FUN = function(i) labels(coef(i@estimates@estimates$det)))
      parm <- paste("p", "(", parm, ")", sep="")
      if(!is.null(reversed.parm)) {reversed.parm <- paste("p", "(", reversed.parm, ")", sep="")}
      not.include <- lapply(cand.set, FUN = function(i) i@formula[[2]])
    }
  
  nmods <- length(cand.set)
  
  ##setup matrix to indicate presence of parms in the model
  include <- matrix(NA, nrow=nmods, ncol=1)
  ##add a check for multiple instances of same variable in given model (i.e., interactions)
  include.check <- matrix(NA, nrow=nmods, ncol=1)

##################################
##################################
###ADDED A NEW OBJECT TO STRIP AWAY lam( ) from parm on line 35
###to enable search with regexpr( ) 
  
  
  ##iterate over each formula in mod_formula list
  for (i in 1:nmods) {
    idents <- NULL
    idents.check <- NULL
    form <- mod_formula[[i]]


######################################################################################################
######################################################################################################
###MODIFICATIONS BEGIN
      ##iterate over each element of formula[[i]] in list
      if(is.null(reversed.parm)) {
        for (j in 1:length(form)) {
          idents[j] <- identical(parm, form[j])
          ##added parm.strip here for regexpr( )
          idents.check[j] <- ifelse(attr(regexpr(parm.strip, form[j], fixed=TRUE), "match.length")== "-1" , 0, 1)  
        }
      } else {
        for (j in 1:length(form)) {
          idents[j] <- identical(parm, form[j]) | identical(reversed.parm, form[j])
          idents.check[j] <- ifelse(attr(regexpr(parm.strip, form[j], fixed=TRUE), "match.length")=="-1" & attr(regexpr(reversed.parm.strip, form[j],
                                                                  fixed=TRUE), "match.length")=="-1" , 0, 1)  
        }
      }
###MODIFICATIONS END
######################################################################################################
######################################################################################################

    
    include[i] <- ifelse(any(idents==1), 1, 0)
    include.check[i] <- ifelse(sum(idents.check)>1, "duplicates", "OK")
  }

  #####################################################
  #exclude == NULL; warn=TRUE:  warn that duplicates occur and stop
  if(is.null(exclude) && identical(warn, TRUE)) {
    #check for duplicates in same model
    if(any(include.check == "duplicates")) {
      stop("\nSome models include more than one instance of the parameter of interest. \n",
           "This may be due to the presence of interaction/polynomial terms, or variables\n",
             "with similar names:\n",
             "\tsee \"?modavg\" for details on variable specification and \"exclude\" argument\n")
      }

    }

    #exclude == NULL; warn=FALSE:  compute model-averaged beta estimate from models including variable of interest,
    #assuming that the variable is not involved in interaction or higher order polynomial (x^2, x^3, etc...),
    #warn that models were not excluded
    if(is.null(exclude) && identical(warn, FALSE)) {
      if(any(include.check == "duplicates")) {
        warning("\nMultiple instances of parameter of interest in given model is presumably\n",
                "not due to interaction or polynomial terms - these models will not be\n",
                "excluded from the computation of model-averaged estimate\n")
      }
      
    }

    #warn if exclude is neither a list nor NULL
    if(!is.null(exclude)) {
      if(!is.list(exclude)) {stop("\nItems in \"exclude\" must be specified as a list")}
    }


    #if exclude is list  
    if(is.list(exclude)) {

    #determine number of elements in exclude
      nexcl <- length(exclude)

      #check each formula for presence of exclude variable in not.include list
      #not.include <- lapply(cand.set, FUN = formula)

      #set up a new list with model formula
      forms <- list()
      for (i in 1:nmods) {
        form.tmp <- as.character(not.include[i]) #changed from other versions as formula returned is of different structure for unmarked objects
        if(attr(regexpr("\\+", form.tmp), "match.length")==-1) {
          forms[i] <- form.tmp
        } else {forms[i] <- strsplit(form.tmp, split=" \\+ ")}
      }

      #additional check to see whether some variable names include "+"
      check.forms <- unlist(lapply(forms, FUN=function(i) any(attr(regexpr("\\+", i), "match.length")>0)[[1]]))
      if (any(check.forms==TRUE)) stop("\nPlease avoid \"+\" in variable names\n")

      ##additional check to determine if intercept was removed from models
      check.forms <- unlist(lapply(forms, FUN=function(i) any(attr(regexpr("\\- 1", i), "match.length")>0)[[1]]))
      if (any(check.forms==TRUE)) stop("\nModels without intercept are not supported in this version, please use alternative parameterization\n")

 
      #search within formula for variables to exclude
      mod.exclude <- matrix(NA, nrow=nmods, ncol=nexcl)

      #iterate over each element in exclude list
      for (var in 1:nexcl) {

      #iterate over each formula in mod_formula list
        for (i in 1:nmods) {
          idents <- NULL
          form.excl <- forms[[i]]

          #iterate over each element of forms[[i]]
          for (j in 1:length(form.excl)) {
            idents[j] <- identical(exclude[var][[1]], form.excl[j])
          }
          mod.exclude[i,var] <- ifelse(any(idents == 1), 1, 0)
        }    
        
      }
  
      #determine outcome across all variables to exclude
      to.exclude <- rowSums(mod.exclude)
  
  
      #exclude models following models from model averaging  
      include[which(to.exclude >= 1)] <- 0
      
      
    }


 
  ##add a check to determine if include always == 0
  if (sum(include) == 0) {stop("\nParameter not found in any of the candidate models\n") }
  
  new.cand.set <- cand.set[which(include == 1)] #select models including a given parameter
  new.mod.name <- modnames[which(include == 1)]    #update model names

  new_table <- aictab(cand.set = new.cand.set, modnames = new.mod.name,
                      second.ord = second.ord, nobs = nobs, sort = FALSE, c.hat = c.hat)  #recompute AIC table and associated measures
  new_table$Beta_est <- unlist(lapply(new.cand.set, FUN = function(i) coef(i)[paste(parm)])) #extract beta estimate for parm
  ##if reversed.parm is not null and varies across models, potentially check for it here
  new_table$SE <- unlist(lapply(new.cand.set, FUN = function(i) sqrt(diag(vcov(i)))[paste(parm)]))
  ##if reversed.parm is not null and varies across models, potentially check for it here

  ##if c-hat is estimated adjust the SE's by multiplying with sqrt of c-hat
  if(c.hat > 1) {
    new_table$SE <- new_table$SE*sqrt(c.hat)
  } 

  ##AICc
  ##compute model-averaged estimates, unconditional SE, and 95% CL
  if(c.hat == 1 && second.ord == TRUE) {
    Modavg_beta <- sum(new_table$AICcWt*new_table$Beta_est)

    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Uncond_SE <- sum(new_table$AICcWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
    }
      
    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Uncond_SE <- sqrt(sum(new_table$AICcWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
    }
  }

  ##QAICc
  if(c.hat > 1 && second.ord == TRUE) {
    Modavg_beta <- sum(new_table$QAICcWt*new_table$Beta_est)
      
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {      
      Uncond_SE <- sum(new_table$QAICcWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
    }
    
    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Uncond_SE <- sqrt(sum(new_table$QAICcWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
    }
  }     


  ##AIC
  if(c.hat == 1 && second.ord == FALSE) {
    Modavg_beta <- sum(new_table$AICWt*new_table$Beta_est)
      
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Uncond_SE <- sum(new_table$AICWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
    }
      
    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Uncond_SE <- sqrt(sum(new_table$AICWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
    }
  }


  ##QAIC
  if(c.hat > 1 && second.ord == FALSE) {
    Modavg_beta <- sum(new_table$QAICWt*new_table$Beta_est)

    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Uncond_SE <- sum(new_table$QAICWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
    }

    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Uncond_SE <- sqrt(sum(new_table$QAICWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
    }  
  }     

  
  zcrit <- qnorm(p = (1 - conf.level)/2, lower.tail = FALSE)
  Lower_CL <- Modavg_beta - zcrit*Uncond_SE
  Upper_CL <- Modavg_beta + zcrit*Uncond_SE
  out.modavg <- list("Parameter"=paste(parm), "Mod.avg.table" = new_table, "Mod.avg.beta" = Modavg_beta,
                     "Uncond.SE" = Uncond_SE, "Conf.level" = conf.level, "Lower.CL" = Lower_CL,
                     "Upper.CL" = Upper_CL)
  
  class(out.modavg) <- c("modavg", "list")
  return(out.modavg)
}



##colext
modavg.AICunmarkedFitColExt <-
  function(cand.set, parm, modnames = NULL, second.ord = TRUE, nobs = NULL, 
           uncond.se = "revised", conf.level = 0.95, exclude = NULL, warn = TRUE,
           c.hat = 1, parm.type = NULL, ...){

  ##note that parameter is referenced differently from unmarked object - see labels( )

  ##check if named list if modnames are not supplied
  if(is.null(modnames)) {
    if(is.null(names(cand.set))) {
      modnames <- paste("Mod", 1:length(cand.set), sep = "")
      warning("\nModel names have been supplied automatically in the table\n")
    } else {
      modnames <- names(cand.set)
    }
  }
    
  
  ##check for parm.type and stop if NULL
  if(is.null(parm.type)) {stop("\n'parm.type' must be specified for this model type, see ?modavg for details\n")}

  
#####MODIFICATIONS BEGIN#######
  ##remove all leading and trailing white space and within parm
  parm <- gsub('[[:space:]]+', "", parm)
  parm.strip <- parm #to use later
  
  ##if (Intercept) is chosen assign (Int) - for compatibility
  if(identical(parm, "(Intercept)")) parm <- "Int"

  ##reverse parm
  reversed.parm <- reverse.parm(parm)
  reversed.parm.strip <- reversed.parm #to use later
  exclude <- reverse.exclude(exclude = exclude)
#####MODIFICATIONS END######
  
  

    ##multiseason occupancy model
    ##psi - initial occupancy
    if(identical(parm.type, "psi")) {
      ##extract model formula for each model in cand.set
      mod_formula <- lapply(cand.set, FUN = function(i) labels(coef(i@estimates@estimates$psi)))
      ##create label for parm
      parm <- paste("psi", "(", parm, ")", sep="")
      if(!is.null(reversed.parm)) {reversed.parm <- paste("psi", "(", reversed.parm, ")", sep="")}
      not.include <- lapply(cand.set, FUN = function(i) i@psiformula)
    }
    ##gamma - extinction
    if(identical(parm.type, "gamma")) {
      ##extract model formula for each model in cand.set
      mod_formula <- lapply(cand.set, FUN = function(i) labels(coef(i@estimates@estimates$col)))
      ##create label for parm
      parm <- paste("col", "(", parm, ")", sep="")
      if(!is.null(reversed.parm)) {reversed.parm <- paste("col", "(", reversed.parm, ")", sep="")}
      not.include <- lapply(cand.set, FUN = function(i) i@gamformula)
    }
    ##epsilon - extinction
    if(identical(parm.type, "epsilon")) {
      ##extract model formula for each model in cand.set
      mod_formula <- lapply(cand.set, FUN = function(i) labels(coef(i@estimates@estimates$ext)))
      ##create label for parm
      parm <- paste("ext", "(", parm, ")", sep="")
      if(!is.null(reversed.parm)) {reversed.parm <- paste("ext", "(", reversed.parm, ")", sep="")}
      ##for epsilon
      not.include <- lapply(cand.set, FUN = function(i) i@epsformula) 
    }
    ##detect
    if(identical(parm.type, "detect")) {
      mod_formula <- lapply(cand.set, FUN = function(i) labels(coef(i@estimates@estimates$det)))
      parm <- paste("p", "(", parm, ")", sep="")
      if(!is.null(reversed.parm)) {reversed.parm <- paste("p", "(", reversed.parm, ")", sep="")}
      not.include <- lapply(cand.set, FUN = function(i) i@detformula)
    }

  nmods <- length(cand.set)
  
  ##setup matrix to indicate presence of parms in the model
  include <- matrix(NA, nrow=nmods, ncol=1)
  ##add a check for multiple instances of same variable in given model (i.e., interactions)
  include.check <- matrix(NA, nrow=nmods, ncol=1)

##################################
##################################
###ADDED A NEW OBJECT TO STRIP AWAY lam( ) from parm on line 35
###to enable search with regexpr( ) 
  
  
  ##iterate over each formula in mod_formula list
  for (i in 1:nmods) {
    idents <- NULL
    idents.check <- NULL
    form <- mod_formula[[i]]


######################################################################################################
######################################################################################################
###MODIFICATIONS BEGIN
      ##iterate over each element of formula[[i]] in list
      if(is.null(reversed.parm)) {
        for (j in 1:length(form)) {
          idents[j] <- identical(parm, form[j])
          ##added parm.strip here for regexpr( )
          idents.check[j] <- ifelse(attr(regexpr(parm.strip, form[j], fixed=TRUE), "match.length")== "-1" , 0, 1)  
        }
      } else {
        for (j in 1:length(form)) {
          idents[j] <- identical(parm, form[j]) | identical(reversed.parm, form[j])
          idents.check[j] <- ifelse(attr(regexpr(parm.strip, form[j], fixed=TRUE), "match.length")=="-1" & attr(regexpr(reversed.parm.strip, form[j],
                                                                  fixed=TRUE), "match.length")=="-1" , 0, 1)  
        }
      }
###MODIFICATIONS END
######################################################################################################
######################################################################################################

    
    include[i] <- ifelse(any(idents==1), 1, 0)
    include.check[i] <- ifelse(sum(idents.check)>1, "duplicates", "OK")
  }

  #####################################################
  #exclude == NULL; warn=TRUE:  warn that duplicates occur and stop
  if(is.null(exclude) && identical(warn, TRUE)) {
    #check for duplicates in same model
    if(any(include.check == "duplicates")) {
      stop("\nSome models include more than one instance of the parameter of interest. \n",
           "This may be due to the presence of interaction/polynomial terms, or variables\n",
             "with similar names:\n",
             "\tsee \"?modavg\" for details on variable specification and \"exclude\" argument\n")
      }

    }

    #exclude == NULL; warn=FALSE:  compute model-averaged beta estimate from models including variable of interest,
    #assuming that the variable is not involved in interaction or higher order polynomial (x^2, x^3, etc...),
    #warn that models were not excluded
    if(is.null(exclude) && identical(warn, FALSE)) {
      if(any(include.check == "duplicates")) {
        warning("\nMultiple instances of parameter of interest in given model is presumably\n",
                "not due to interaction or polynomial terms - these models will not be\n",
                "excluded from the computation of model-averaged estimate\n")
      }
      
    }

    #warn if exclude is neither a list nor NULL
    if(!is.null(exclude)) {
      if(!is.list(exclude)) {stop("\nItems in \"exclude\" must be specified as a list")}
    }


    #if exclude is list  
    if(is.list(exclude)) {

    #determine number of elements in exclude
      nexcl <- length(exclude)

      #check each formula for presence of exclude variable in not.include list
      #not.include <- lapply(cand.set, FUN = formula)

      #set up a new list with model formula
      forms <- list()
      for (i in 1:nmods) {
        form.tmp <- as.character(not.include[i]) #changed from other versions as formula returned is of different structure for unmarked objects
        if(attr(regexpr("\\+", form.tmp), "match.length")==-1) {
          forms[i] <- form.tmp
        } else {forms[i] <- strsplit(form.tmp, split=" \\+ ")}
      }

      #additional check to see whether some variable names include "+"
      check.forms <- unlist(lapply(forms, FUN=function(i) any(attr(regexpr("\\+", i), "match.length")>0)[[1]]))
      if (any(check.forms==TRUE)) stop("\nPlease avoid \"+\" in variable names\n")

      ##additional check to determine if intercept was removed from models
      check.forms <- unlist(lapply(forms, FUN=function(i) any(attr(regexpr("\\- 1", i), "match.length")>0)[[1]]))
      if (any(check.forms==TRUE)) stop("\nModels without intercept are not supported in this version, please use alternative parameterization\n")

 
      #search within formula for variables to exclude
      mod.exclude <- matrix(NA, nrow=nmods, ncol=nexcl)

      #iterate over each element in exclude list
      for (var in 1:nexcl) {

      #iterate over each formula in mod_formula list
        for (i in 1:nmods) {
          idents <- NULL
          form.excl <- forms[[i]]

          #iterate over each element of forms[[i]]
          for (j in 1:length(form.excl)) {
            idents[j] <- identical(exclude[var][[1]], form.excl[j])
          }
          mod.exclude[i,var] <- ifelse(any(idents == 1), 1, 0)
        }    
        
      }
  
      #determine outcome across all variables to exclude
      to.exclude <- rowSums(mod.exclude)
  
  
      #exclude models following models from model averaging  
      include[which(to.exclude >= 1)] <- 0
      
      
    }


 
  ##add a check to determine if include always == 0
  if (sum(include) == 0) {stop("\nParameter not found in any of the candidate models\n") }
  
  new.cand.set <- cand.set[which(include == 1)] #select models including a given parameter
  new.mod.name <- modnames[which(include == 1)]    #update model names

  new_table <- aictab(cand.set = new.cand.set, modnames = new.mod.name,
                      second.ord = second.ord, nobs = nobs, sort = FALSE, c.hat = c.hat) #recompute AIC table and associated measures
  new_table$Beta_est <- unlist(lapply(new.cand.set, FUN = function(i) coef(i)[paste(parm)])) #extract beta estimate for parm
  ##if reversed.parm is not null and varies across models, potentially check for it here
  new_table$SE <- unlist(lapply(new.cand.set, FUN = function(i) sqrt(diag(vcov(i)))[paste(parm)]))
  ##if reversed.parm is not null and varies across models, potentially check for it here

  ##if c-hat is estimated adjust the SE's by multiplying with sqrt of c-hat
  if(c.hat > 1) {
    new_table$SE <- new_table$SE*sqrt(c.hat)
  } 

  ##AICc
  ##compute model-averaged estimates, unconditional SE, and 95% CL
  if(c.hat == 1 && second.ord == TRUE) {
    Modavg_beta <- sum(new_table$AICcWt*new_table$Beta_est)

    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Uncond_SE <- sum(new_table$AICcWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
    }
      
    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Uncond_SE <- sqrt(sum(new_table$AICcWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
    }
  }

  ##QAICc
  if(c.hat > 1 && second.ord == TRUE) {
    Modavg_beta <- sum(new_table$QAICcWt*new_table$Beta_est)
      
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {      
      Uncond_SE <- sum(new_table$QAICcWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
    }
    
    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Uncond_SE <- sqrt(sum(new_table$QAICcWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
    }
  }     


  ##AIC
  if(c.hat == 1 && second.ord == FALSE) {
    Modavg_beta <- sum(new_table$AICWt*new_table$Beta_est)
      
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Uncond_SE <- sum(new_table$AICWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
    }
      
    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Uncond_SE <- sqrt(sum(new_table$AICWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
    }
  }


  ##QAIC
  if(c.hat > 1 && second.ord == FALSE) {
    Modavg_beta <- sum(new_table$QAICWt*new_table$Beta_est)

    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Uncond_SE <- sum(new_table$QAICWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
    }

    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Uncond_SE <- sqrt(sum(new_table$QAICWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
    }  
  }     

  
  zcrit <- qnorm(p = (1 - conf.level)/2, lower.tail = FALSE)
  Lower_CL <- Modavg_beta - zcrit*Uncond_SE
  Upper_CL <- Modavg_beta + zcrit*Uncond_SE
  out.modavg <- list("Parameter"=paste(parm), "Mod.avg.table" = new_table, "Mod.avg.beta" = Modavg_beta,
                     "Uncond.SE" = Uncond_SE, "Conf.level" = conf.level, "Lower.CL" = Lower_CL,
                     "Upper.CL" = Upper_CL)
  
  class(out.modavg) <- c("modavg", "list")
  return(out.modavg)
}



##occuRN
modavg.AICunmarkedFitOccuRN <-
function(cand.set, parm, modnames = NULL, second.ord = TRUE, nobs = NULL, 
         uncond.se = "revised", conf.level = 0.95, exclude = NULL, warn = TRUE,
         c.hat = 1, parm.type = NULL, ...){

  ##note that parameter is referenced differently from unmarked object - see labels( )

  ##check if named list if modnames are not supplied
  if(is.null(modnames)) {
    if(is.null(names(cand.set))) {
      modnames <- paste("Mod", 1:length(cand.set), sep = "")
      warning("\nModel names have been supplied automatically in the table\n")
    } else {
      modnames <- names(cand.set)
    }
  }
    
  
  ##check for parm.type and stop if NULL
  if(is.null(parm.type)) {stop("\n'parm.type' must be specified for this model type, see ?modavg for details\n")}

  
#####MODIFICATIONS BEGIN#######
  ##remove all leading and trailing white space and within parm
  parm <- gsub('[[:space:]]+', "", parm)
  parm.strip <- parm #to use later
  
  ##if (Intercept) is chosen assign (Int) - for compatibility
  if(identical(parm, "(Intercept)")) parm <- "Int"

  ##reverse parm
  reversed.parm <- reverse.parm(parm)
  reversed.parm.strip <- reversed.parm #to use later
  exclude <- reverse.exclude(exclude = exclude)
#####MODIFICATIONS END######
  
  ##Royle-Nichols heterogeneity model
  ##lambda - abundance
  if(identical(parm.type, "lambda")) {
    ##extract model formula for each model in cand.set
    mod_formula <- lapply(cand.set, FUN = function(i) labels(coef(i@estimates@estimates$state)))
    parm <- paste("lam", "(", parm, ")", sep="")
    if(!is.null(reversed.parm)) {reversed.parm <- paste("lam", "(", reversed.parm, ")", sep="")}
    not.include <- lapply(cand.set, FUN = function(i) i@formula[[3]])
  }
  ##detect
  if(identical(parm.type, "detect")) {
    mod_formula<-lapply(cand.set, FUN = function(i) labels(coef(i@estimates@estimates$det)))
    parm <- paste("p", "(", parm, ")", sep="")
    if(!is.null(reversed.parm)) {reversed.parm <- paste("p", "(", reversed.parm, ")", sep="")}
    not.include <- lapply(cand.set, FUN = function(i) i@formula[[2]])
  }

  
  nmods <- length(cand.set)
  
  ##setup matrix to indicate presence of parms in the model
  include <- matrix(NA, nrow=nmods, ncol=1)
  ##add a check for multiple instances of same variable in given model (i.e., interactions)
  include.check <- matrix(NA, nrow=nmods, ncol=1)

##################################
##################################
###ADDED A NEW OBJECT TO STRIP AWAY lam( ) from parm on line 35
###to enable search with regexpr( ) 
  
  
  ##iterate over each formula in mod_formula list
  for (i in 1:nmods) {
    idents <- NULL
    idents.check <- NULL
    form <- mod_formula[[i]]


######################################################################################################
######################################################################################################
###MODIFICATIONS BEGIN
      ##iterate over each element of formula[[i]] in list
      if(is.null(reversed.parm)) {
        for (j in 1:length(form)) {
          idents[j] <- identical(parm, form[j])
          ##added parm.strip here for regexpr( )
          idents.check[j] <- ifelse(attr(regexpr(parm.strip, form[j], fixed=TRUE), "match.length")== "-1" , 0, 1)  
        }
      } else {
        for (j in 1:length(form)) {
          idents[j] <- identical(parm, form[j]) | identical(reversed.parm, form[j])
          idents.check[j] <- ifelse(attr(regexpr(parm.strip, form[j], fixed=TRUE), "match.length")=="-1" & attr(regexpr(reversed.parm.strip, form[j],
                                                                  fixed=TRUE), "match.length")=="-1" , 0, 1)  
        }
      }
###MODIFICATIONS END
######################################################################################################
######################################################################################################

    
    include[i] <- ifelse(any(idents==1), 1, 0)
    include.check[i] <- ifelse(sum(idents.check)>1, "duplicates", "OK")
  }

  #####################################################
  #exclude == NULL; warn=TRUE:  warn that duplicates occur and stop
  if(is.null(exclude) && identical(warn, TRUE)) {
    #check for duplicates in same model
    if(any(include.check == "duplicates")) {
      stop("\nSome models include more than one instance of the parameter of interest. \n",
           "This may be due to the presence of interaction/polynomial terms, or variables\n",
             "with similar names:\n",
             "\tsee \"?modavg\" for details on variable specification and \"exclude\" argument\n")
      }

    }

    #exclude == NULL; warn=FALSE:  compute model-averaged beta estimate from models including variable of interest,
    #assuming that the variable is not involved in interaction or higher order polynomial (x^2, x^3, etc...),
    #warn that models were not excluded
    if(is.null(exclude) && identical(warn, FALSE)) {
      if(any(include.check == "duplicates")) {
        warning("\nMultiple instances of parameter of interest in given model is presumably\n",
                "not due to interaction or polynomial terms - these models will not be\n",
                "excluded from the computation of model-averaged estimate\n")
      }
      
    }

    #warn if exclude is neither a list nor NULL
    if(!is.null(exclude)) {
      if(!is.list(exclude)) {stop("\nItems in \"exclude\" must be specified as a list")}
    }


    #if exclude is list  
    if(is.list(exclude)) {

    #determine number of elements in exclude
      nexcl <- length(exclude)

      #check each formula for presence of exclude variable in not.include list
      #not.include <- lapply(cand.set, FUN = formula)

      #set up a new list with model formula
      forms <- list()
      for (i in 1:nmods) {
        form.tmp <- as.character(not.include[i]) #changed from other versions as formula returned is of different structure for unmarked objects
        if(attr(regexpr("\\+", form.tmp), "match.length")==-1) {
          forms[i] <- form.tmp
        } else {forms[i] <- strsplit(form.tmp, split=" \\+ ")}
      }

      #additional check to see whether some variable names include "+"
      check.forms <- unlist(lapply(forms, FUN=function(i) any(attr(regexpr("\\+", i), "match.length")>0)[[1]]))
      if (any(check.forms==TRUE)) stop("\nPlease avoid \"+\" in variable names\n")

      ##additional check to determine if intercept was removed from models
      check.forms <- unlist(lapply(forms, FUN=function(i) any(attr(regexpr("\\- 1", i), "match.length")>0)[[1]]))
      if (any(check.forms==TRUE)) stop("\nModels without intercept are not supported in this version, please use alternative parameterization\n")

 
      #search within formula for variables to exclude
      mod.exclude <- matrix(NA, nrow=nmods, ncol=nexcl)

      #iterate over each element in exclude list
      for (var in 1:nexcl) {

      #iterate over each formula in mod_formula list
        for (i in 1:nmods) {
          idents <- NULL
          form.excl <- forms[[i]]

          #iterate over each element of forms[[i]]
          for (j in 1:length(form.excl)) {
            idents[j] <- identical(exclude[var][[1]], form.excl[j])
          }
          mod.exclude[i,var] <- ifelse(any(idents == 1), 1, 0)
        }    
        
      }
  
      #determine outcome across all variables to exclude
      to.exclude <- rowSums(mod.exclude)
  
  
      #exclude models following models from model averaging  
      include[which(to.exclude >= 1)] <- 0
      
      
    }


 
  ##add a check to determine if include always == 0
  if (sum(include) == 0) {stop("\nParameter not found in any of the candidate models\n") }
  
  new.cand.set <- cand.set[which(include == 1)] #select models including a given parameter
  new.mod.name <- modnames[which(include == 1)]    #update model names

  new_table <- aictab(cand.set = new.cand.set, modnames = new.mod.name,
                      second.ord = second.ord, nobs = nobs, sort = FALSE, c.hat = c.hat) #recompute AIC table and associated measures
  new_table$Beta_est <- unlist(lapply(new.cand.set, FUN = function(i) coef(i)[paste(parm)])) #extract beta estimate for parm
  ##if reversed.parm is not null and varies across models, potentially check for it here
  new_table$SE <- unlist(lapply(new.cand.set, FUN = function(i) sqrt(diag(vcov(i)))[paste(parm)]))
  ##if reversed.parm is not null and varies across models, potentially check for it here

  ##if c-hat is estimated adjust the SE's by multiplying with sqrt of c-hat
  if(c.hat > 1) {
    new_table$SE <- new_table$SE*sqrt(c.hat)
  } 

  ##AICc
  ##compute model-averaged estimates, unconditional SE, and 95% CL
  if(c.hat == 1 && second.ord == TRUE) {
    Modavg_beta <- sum(new_table$AICcWt*new_table$Beta_est)

    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Uncond_SE <- sum(new_table$AICcWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
    }
      
    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Uncond_SE <- sqrt(sum(new_table$AICcWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
    }
  }

  ##QAICc
  if(c.hat > 1 && second.ord == TRUE) {
    Modavg_beta <- sum(new_table$QAICcWt*new_table$Beta_est)
      
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {      
      Uncond_SE <- sum(new_table$QAICcWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
    }
    
    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Uncond_SE <- sqrt(sum(new_table$QAICcWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
    }
  }     


  ##AIC
  if(c.hat == 1 && second.ord == FALSE) {
    Modavg_beta <- sum(new_table$AICWt*new_table$Beta_est)
      
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Uncond_SE <- sum(new_table$AICWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
    }
      
    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Uncond_SE <- sqrt(sum(new_table$AICWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
    }
  }


  ##QAIC
  if(c.hat > 1 && second.ord == FALSE) {
    Modavg_beta <- sum(new_table$QAICWt*new_table$Beta_est)

    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Uncond_SE <- sum(new_table$QAICWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
    }

    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Uncond_SE <- sqrt(sum(new_table$QAICWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
    }  
  }     

  
  zcrit <- qnorm(p = (1 - conf.level)/2, lower.tail = FALSE)
  Lower_CL <- Modavg_beta - zcrit*Uncond_SE
  Upper_CL <- Modavg_beta + zcrit*Uncond_SE
  out.modavg <- list("Parameter"=paste(parm), "Mod.avg.table" = new_table, "Mod.avg.beta" = Modavg_beta,
                     "Uncond.SE" = Uncond_SE, "Conf.level" = conf.level, "Lower.CL" = Lower_CL,
                     "Upper.CL" = Upper_CL)
  
  class(out.modavg) <- c("modavg", "list")
  return(out.modavg)
}



##pcount
modavg.AICunmarkedFitPCount <-
function(cand.set, parm, modnames = NULL, second.ord = TRUE, nobs = NULL, 
         uncond.se = "revised", conf.level = 0.95, exclude = NULL, warn = TRUE,
         c.hat = 1, parm.type = NULL, ...){

  ##note that parameter is referenced differently from unmarked object - see labels( )

  ##check if named list if modnames are not supplied
  if(is.null(modnames)) {
    if(is.null(names(cand.set))) {
      modnames <- paste("Mod", 1:length(cand.set), sep = "")
      warning("\nModel names have been supplied automatically in the table\n")
    } else {
      modnames <- names(cand.set)
    }
  }
    
  
  ##check for parm.type and stop if NULL
  if(is.null(parm.type)) {stop("\n'parm.type' must be specified for this model type, see ?modavg for details\n")}

  
#####MODIFICATIONS BEGIN#######
  ##remove all leading and trailing white space and within parm
  parm <- gsub('[[:space:]]+', "", parm)
  parm.strip <- parm #to use later
  
  ##if (Intercept) is chosen assign (Int) - for compatibility
  if(identical(parm, "(Intercept)")) parm <- "Int"

  ##reverse parm
  reversed.parm <- reverse.parm(parm)
  reversed.parm.strip <- reversed.parm #to use later
  exclude <- reverse.exclude(exclude = exclude)
#####MODIFICATIONS END######

  ##single season N-mixture model
  ##lambda - abundance
  if(identical(parm.type, "lambda")) {
    ##extract model formula for each model in cand.set
    mod_formula <- lapply(cand.set, FUN = function(i) labels(coef(i@estimates@estimates$state)))
    ##create label for parm
    parm <- paste("lam", "(", parm, ")", sep="")
    if(!is.null(reversed.parm)) {reversed.parm <- paste("lam", "(", reversed.parm, ")", sep="")}
    not.include <- lapply(cand.set, FUN = function(i) i@formula[[3]])
  }
  ##detect
  if(identical(parm.type, "detect")) {
    mod_formula <- lapply(cand.set, FUN = function(i) labels(coef(i@estimates@estimates$det)))
    parm <- paste("p", "(", parm, ")", sep="")
    if(!is.null(reversed.parm)) {reversed.parm <- paste("p", "(", reversed.parm, ")", sep="")}
    not.include <- lapply(cand.set, FUN = function(i) i@formula[[2]])
  }
  
    
  nmods <- length(cand.set)
  
  ##setup matrix to indicate presence of parms in the model
  include <- matrix(NA, nrow=nmods, ncol=1)
  ##add a check for multiple instances of same variable in given model (i.e., interactions)
  include.check <- matrix(NA, nrow=nmods, ncol=1)

##################################
##################################
###ADDED A NEW OBJECT TO STRIP AWAY lam( ) from parm on line 35
###to enable search with regexpr( ) 
  
  
  ##iterate over each formula in mod_formula list
  for (i in 1:nmods) {
    idents <- NULL
    idents.check <- NULL
    form <- mod_formula[[i]]


######################################################################################################
######################################################################################################
###MODIFICATIONS BEGIN
      ##iterate over each element of formula[[i]] in list
      if(is.null(reversed.parm)) {
        for (j in 1:length(form)) {
          idents[j] <- identical(parm, form[j])
          ##added parm.strip here for regexpr( )
          idents.check[j] <- ifelse(attr(regexpr(parm.strip, form[j], fixed=TRUE), "match.length")== "-1" , 0, 1)  
        }
      } else {
        for (j in 1:length(form)) {
          idents[j] <- identical(parm, form[j]) | identical(reversed.parm, form[j])
          idents.check[j] <- ifelse(attr(regexpr(parm.strip, form[j], fixed=TRUE), "match.length")=="-1" & attr(regexpr(reversed.parm.strip, form[j],
                                                                  fixed=TRUE), "match.length")=="-1" , 0, 1)  
        }
      }
###MODIFICATIONS END
######################################################################################################
######################################################################################################

    
    include[i] <- ifelse(any(idents==1), 1, 0)
    include.check[i] <- ifelse(sum(idents.check)>1, "duplicates", "OK")
  }

  #####################################################
  #exclude == NULL; warn=TRUE:  warn that duplicates occur and stop
  if(is.null(exclude) && identical(warn, TRUE)) {
    #check for duplicates in same model
    if(any(include.check == "duplicates")) {
      stop("\nSome models include more than one instance of the parameter of interest. \n",
           "This may be due to the presence of interaction/polynomial terms, or variables\n",
             "with similar names:\n",
             "\tsee \"?modavg\" for details on variable specification and \"exclude\" argument\n")
      }

    }

    #exclude == NULL; warn=FALSE:  compute model-averaged beta estimate from models including variable of interest,
    #assuming that the variable is not involved in interaction or higher order polynomial (x^2, x^3, etc...),
    #warn that models were not excluded
    if(is.null(exclude) && identical(warn, FALSE)) {
      if(any(include.check == "duplicates")) {
        warning("\nMultiple instances of parameter of interest in given model is presumably\n",
                "not due to interaction or polynomial terms - these models will not be\n",
                "excluded from the computation of model-averaged estimate\n")
      }
      
    }

    #warn if exclude is neither a list nor NULL
    if(!is.null(exclude)) {
      if(!is.list(exclude)) {stop("\nItems in \"exclude\" must be specified as a list")}
    }


    #if exclude is list  
    if(is.list(exclude)) {

    #determine number of elements in exclude
      nexcl <- length(exclude)

      #check each formula for presence of exclude variable in not.include list
      #not.include <- lapply(cand.set, FUN = formula)

      #set up a new list with model formula
      forms <- list()
      for (i in 1:nmods) {
        form.tmp <- as.character(not.include[i]) #changed from other versions as formula returned is of different structure for unmarked objects
        if(attr(regexpr("\\+", form.tmp), "match.length")==-1) {
          forms[i] <- form.tmp
        } else {forms[i] <- strsplit(form.tmp, split=" \\+ ")}
      }

      #additional check to see whether some variable names include "+"
      check.forms <- unlist(lapply(forms, FUN=function(i) any(attr(regexpr("\\+", i), "match.length")>0)[[1]]))
      if (any(check.forms==TRUE)) stop("\nPlease avoid \"+\" in variable names\n")

      ##additional check to determine if intercept was removed from models
      check.forms <- unlist(lapply(forms, FUN=function(i) any(attr(regexpr("\\- 1", i), "match.length")>0)[[1]]))
      if (any(check.forms==TRUE)) stop("\nModels without intercept are not supported in this version, please use alternative parameterization\n")

 
      #search within formula for variables to exclude
      mod.exclude <- matrix(NA, nrow=nmods, ncol=nexcl)

      #iterate over each element in exclude list
      for (var in 1:nexcl) {

      #iterate over each formula in mod_formula list
        for (i in 1:nmods) {
          idents <- NULL
          form.excl <- forms[[i]]

          #iterate over each element of forms[[i]]
          for (j in 1:length(form.excl)) {
            idents[j] <- identical(exclude[var][[1]], form.excl[j])
          }
          mod.exclude[i,var] <- ifelse(any(idents == 1), 1, 0)
        }    
        
      }
  
      #determine outcome across all variables to exclude
      to.exclude <- rowSums(mod.exclude)
  
  
      #exclude models following models from model averaging  
      include[which(to.exclude >= 1)] <- 0
      
      
    }


 
  ##add a check to determine if include always == 0
  if (sum(include) == 0) {stop("\nParameter not found in any of the candidate models\n") }
  
  new.cand.set <- cand.set[which(include == 1)] #select models including a given parameter
  new.mod.name <- modnames[which(include == 1)]    #update model names

  new_table <- aictab(cand.set = new.cand.set, modnames = new.mod.name,
                      second.ord = second.ord, nobs = nobs, sort = FALSE, c.hat = c.hat) #recompute AIC table and associated measures
  new_table$Beta_est <- unlist(lapply(new.cand.set, FUN = function(i) coef(i)[paste(parm)])) #extract beta estimate for parm
  ##if reversed.parm is not null and varies across models, potentially check for it here
  new_table$SE <- unlist(lapply(new.cand.set, FUN = function(i) sqrt(diag(vcov(i)))[paste(parm)]))
  ##if reversed.parm is not null and varies across models, potentially check for it here

  ##if c-hat is estimated adjust the SE's by multiplying with sqrt of c-hat
  if(c.hat > 1) {
    new_table$SE <- new_table$SE*sqrt(c.hat)
  } 

  ##AICc
  ##compute model-averaged estimates, unconditional SE, and 95% CL
  if(c.hat == 1 && second.ord == TRUE) {
    Modavg_beta <- sum(new_table$AICcWt*new_table$Beta_est)

    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Uncond_SE <- sum(new_table$AICcWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
    }
      
    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Uncond_SE <- sqrt(sum(new_table$AICcWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
    }
  }

  ##QAICc
  if(c.hat > 1 && second.ord == TRUE) {
    Modavg_beta <- sum(new_table$QAICcWt*new_table$Beta_est)
      
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {      
      Uncond_SE <- sum(new_table$QAICcWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
    }
    
    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Uncond_SE <- sqrt(sum(new_table$QAICcWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
    }
  }     


  ##AIC
  if(c.hat == 1 && second.ord == FALSE) {
    Modavg_beta <- sum(new_table$AICWt*new_table$Beta_est)
      
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Uncond_SE <- sum(new_table$AICWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
    }
      
    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Uncond_SE <- sqrt(sum(new_table$AICWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
    }
  }


  ##QAIC
  if(c.hat > 1 && second.ord == FALSE) {
    Modavg_beta <- sum(new_table$QAICWt*new_table$Beta_est)

    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Uncond_SE <- sum(new_table$QAICWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
    }

    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Uncond_SE <- sqrt(sum(new_table$QAICWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
    }  
  }     

  
  zcrit <- qnorm(p = (1 - conf.level)/2, lower.tail = FALSE)
  Lower_CL <- Modavg_beta - zcrit*Uncond_SE
  Upper_CL <- Modavg_beta + zcrit*Uncond_SE
  out.modavg <- list("Parameter"=paste(parm), "Mod.avg.table" = new_table, "Mod.avg.beta" = Modavg_beta,
                     "Uncond.SE" = Uncond_SE, "Conf.level" = conf.level, "Lower.CL" = Lower_CL,
                     "Upper.CL" = Upper_CL)
  
  class(out.modavg) <- c("modavg", "list")
  return(out.modavg)
}



##pcountOpen
modavg.AICunmarkedFitPCO <-
  function(cand.set, parm, modnames = NULL, second.ord = TRUE, nobs = NULL, 
           uncond.se = "revised", conf.level = 0.95, exclude = NULL, warn = TRUE,
           c.hat = 1, parm.type = NULL, ...){

  ##note that parameter is referenced differently from unmarked object - see labels( )

  ##check if named list if modnames are not supplied
  if(is.null(modnames)) {
    if(is.null(names(cand.set))) {
      modnames <- paste("Mod", 1:length(cand.set), sep = "")
      warning("\nModel names have been supplied automatically in the table\n")
    } else {
      modnames <- names(cand.set)
    }
  }
    
  
  ##check for parm.type and stop if NULL
  if(is.null(parm.type)) {stop("\n'parm.type' must be specified for this model type, see ?modavg for details\n")}

  
#####MODIFICATIONS BEGIN#######
  ##remove all leading and trailing white space and within parm
  parm <- gsub('[[:space:]]+', "", parm)
  parm.strip <- parm #to use later
  
  ##if (Intercept) is chosen assign (Int) - for compatibility
  if(identical(parm, "(Intercept)")) parm <- "Int"

  ##reverse parm
  reversed.parm <- reverse.parm(parm)
  reversed.parm.strip <- reversed.parm #to use later
  exclude <- reverse.exclude(exclude = exclude)
#####MODIFICATIONS END######
  
  
    
  ##open version of N-mixture model
  ##lambda - abundance
  if(identical(parm.type, "lambda")) {
    ##extract model formula for each model in cand.set
    mod_formula <- lapply(cand.set, FUN = function(i) labels(coef(i@estimates@estimates$lambda)))
    parm <- paste("lam", "(", parm, ")", sep="")
    if(!is.null(reversed.parm)) {reversed.parm <- paste("lam", "(", reversed.parm, ")", sep="")}
    not.include <- lapply(cand.set, FUN = function(i) i@formlist$lambdaformula)
  }
  ##gamma - recruitment
  if(identical(parm.type, "gamma")) {
    ##extract model formula for each model in cand.set
    mod_formula <- lapply(cand.set, FUN = function(i) labels(coef(i@estimates@estimates$gamma)))
    
    ##determine if same H0 on gamma (gamConst, gamAR, gamTrend)
    strip.gam <- sapply(mod_formula, FUN = function(i) unlist(strsplit(i, "\\("))[[1]])
    unique.gam <- unique(strip.gam)
    if(length(unique.gam) > 1) stop("\nDifferent formulations of gamma parameter occur among models:\n
beta estimates cannot be model-averaged\n")
    
    ##create label for parm
    parm <- paste(unique.gam, "(", parm, ")", sep="")
    if(!is.null(reversed.parm)) {reversed.parm <- paste(unique.gam, "(", reversed.parm, ")", sep="")}
    not.include <- lapply(cand.set, FUN = function(i) i@formlist$gammaformula)
  }
  ##omega - apparent survival
  if(identical(parm.type, "omega")) {
    ##extract model formula for each model in cand.set
    mod_formula <- lapply(cand.set, FUN = function(i) labels(coef(i@estimates@estimates$omega)))   
    ##create label for parm
    parm <- paste("omega", "(", parm, ")", sep="")
    if(!is.null(reversed.parm)) {reversed.parm <- paste("omega", "(", reversed.parm, ")", sep="")}
    not.include <- lapply(cand.set, FUN = function(i) i@formlist$omegaformula)
  }
  ##detect
  if(identical(parm.type, "detect")) {
    mod_formula<-lapply(cand.set, FUN = function(i) labels(coef(i@estimates@estimates$det)))
    parm <- paste("p", "(", parm, ")", sep="")
    if(!is.null(reversed.parm)) {reversed.parm <- paste("p", "(", reversed.parm, ")", sep="")}
    not.include <- lapply(cand.set, FUN = function(i) i@formlist$pformula)
  }
  
  nmods <- length(cand.set)
  
  ##setup matrix to indicate presence of parms in the model
  include <- matrix(NA, nrow=nmods, ncol=1)
  ##add a check for multiple instances of same variable in given model (i.e., interactions)
  include.check <- matrix(NA, nrow=nmods, ncol=1)

##################################
##################################
###ADDED A NEW OBJECT TO STRIP AWAY lam( ) from parm on line 35
###to enable search with regexpr( ) 
  
  
  ##iterate over each formula in mod_formula list
  for (i in 1:nmods) {
    idents <- NULL
    idents.check <- NULL
    form <- mod_formula[[i]]


######################################################################################################
######################################################################################################
###MODIFICATIONS BEGIN
      ##iterate over each element of formula[[i]] in list
      if(is.null(reversed.parm)) {
        for (j in 1:length(form)) {
          idents[j] <- identical(parm, form[j])
          ##added parm.strip here for regexpr( )
          idents.check[j] <- ifelse(attr(regexpr(parm.strip, form[j], fixed=TRUE), "match.length")== "-1" , 0, 1)  
        }
      } else {
        for (j in 1:length(form)) {
          idents[j] <- identical(parm, form[j]) | identical(reversed.parm, form[j])
          idents.check[j] <- ifelse(attr(regexpr(parm.strip, form[j], fixed=TRUE), "match.length")=="-1" & attr(regexpr(reversed.parm.strip, form[j],
                                                                  fixed=TRUE), "match.length")=="-1" , 0, 1)  
        }
      }
###MODIFICATIONS END
######################################################################################################
######################################################################################################

    
    include[i] <- ifelse(any(idents==1), 1, 0)
    include.check[i] <- ifelse(sum(idents.check)>1, "duplicates", "OK")
  }

  #####################################################
  #exclude == NULL; warn=TRUE:  warn that duplicates occur and stop
  if(is.null(exclude) && identical(warn, TRUE)) {
    #check for duplicates in same model
    if(any(include.check == "duplicates")) {
      stop("\nSome models include more than one instance of the parameter of interest. \n",
           "This may be due to the presence of interaction/polynomial terms, or variables\n",
             "with similar names:\n",
             "\tsee \"?modavg\" for details on variable specification and \"exclude\" argument\n")
      }

    }

    #exclude == NULL; warn=FALSE:  compute model-averaged beta estimate from models including variable of interest,
    #assuming that the variable is not involved in interaction or higher order polynomial (x^2, x^3, etc...),
    #warn that models were not excluded
    if(is.null(exclude) && identical(warn, FALSE)) {
      if(any(include.check == "duplicates")) {
        warning("\nMultiple instances of parameter of interest in given model is presumably\n",
                "not due to interaction or polynomial terms - these models will not be\n",
                "excluded from the computation of model-averaged estimate\n")
      }
      
    }

    #warn if exclude is neither a list nor NULL
    if(!is.null(exclude)) {
      if(!is.list(exclude)) {stop("\nItems in \"exclude\" must be specified as a list")}
    }


    #if exclude is list  
    if(is.list(exclude)) {

    #determine number of elements in exclude
      nexcl <- length(exclude)

      #check each formula for presence of exclude variable in not.include list
      #not.include <- lapply(cand.set, FUN = formula)

      #set up a new list with model formula
      forms <- list()
      for (i in 1:nmods) {
        form.tmp <- as.character(not.include[i]) #changed from other versions as formula returned is of different structure for unmarked objects
        if(attr(regexpr("\\+", form.tmp), "match.length")==-1) {
          forms[i] <- form.tmp
        } else {forms[i] <- strsplit(form.tmp, split=" \\+ ")}
      }

      #additional check to see whether some variable names include "+"
      check.forms <- unlist(lapply(forms, FUN=function(i) any(attr(regexpr("\\+", i), "match.length")>0)[[1]]))
      if (any(check.forms==TRUE)) stop("\nPlease avoid \"+\" in variable names\n")

      ##additional check to determine if intercept was removed from models
      check.forms <- unlist(lapply(forms, FUN=function(i) any(attr(regexpr("\\- 1", i), "match.length")>0)[[1]]))
      if (any(check.forms==TRUE)) stop("\nModels without intercept are not supported in this version, please use alternative parameterization\n")

 
      #search within formula for variables to exclude
      mod.exclude <- matrix(NA, nrow=nmods, ncol=nexcl)

      #iterate over each element in exclude list
      for (var in 1:nexcl) {

      #iterate over each formula in mod_formula list
        for (i in 1:nmods) {
          idents <- NULL
          form.excl <- forms[[i]]

          #iterate over each element of forms[[i]]
          for (j in 1:length(form.excl)) {
            idents[j] <- identical(exclude[var][[1]], form.excl[j])
          }
          mod.exclude[i,var] <- ifelse(any(idents == 1), 1, 0)
        }    
        
      }
  
      #determine outcome across all variables to exclude
      to.exclude <- rowSums(mod.exclude)
  
  
      #exclude models following models from model averaging  
      include[which(to.exclude >= 1)] <- 0
      
      
    }


 
  ##add a check to determine if include always == 0
  if (sum(include) == 0) {stop("\nParameter not found in any of the candidate models\n") }
  
  new.cand.set <- cand.set[which(include == 1)] #select models including a given parameter
  new.mod.name <- modnames[which(include == 1)]    #update model names

  new_table <- aictab(cand.set = new.cand.set, modnames = new.mod.name,
                      second.ord = second.ord, nobs = nobs, sort = FALSE, c.hat = c.hat) #recompute AIC table and associated measures
  new_table$Beta_est <- unlist(lapply(new.cand.set, FUN = function(i) coef(i)[paste(parm)])) #extract beta estimate for parm
  ##if reversed.parm is not null and varies across models, potentially check for it here
  new_table$SE <- unlist(lapply(new.cand.set, FUN = function(i) sqrt(diag(vcov(i)))[paste(parm)]))
  ##if reversed.parm is not null and varies across models, potentially check for it here

  ##if c-hat is estimated adjust the SE's by multiplying with sqrt of c-hat
  if(c.hat > 1) {
    new_table$SE <- new_table$SE*sqrt(c.hat)
  } 

  ##AICc
  ##compute model-averaged estimates, unconditional SE, and 95% CL
  if(c.hat == 1 && second.ord == TRUE) {
    Modavg_beta <- sum(new_table$AICcWt*new_table$Beta_est)

    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Uncond_SE <- sum(new_table$AICcWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
    }
      
    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Uncond_SE <- sqrt(sum(new_table$AICcWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
    }
  }

  ##QAICc
  if(c.hat > 1 && second.ord == TRUE) {
    Modavg_beta <- sum(new_table$QAICcWt*new_table$Beta_est)
      
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {      
      Uncond_SE <- sum(new_table$QAICcWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
    }
    
    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Uncond_SE <- sqrt(sum(new_table$QAICcWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
    }
  }     


  ##AIC
  if(c.hat == 1 && second.ord == FALSE) {
    Modavg_beta <- sum(new_table$AICWt*new_table$Beta_est)
      
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Uncond_SE <- sum(new_table$AICWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
    }
      
    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Uncond_SE <- sqrt(sum(new_table$AICWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
    }
  }


  ##QAIC
  if(c.hat > 1 && second.ord == FALSE) {
    Modavg_beta <- sum(new_table$QAICWt*new_table$Beta_est)

    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Uncond_SE <- sum(new_table$QAICWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
    }

    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Uncond_SE <- sqrt(sum(new_table$QAICWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
    }  
  }     

  
  zcrit <- qnorm(p = (1 - conf.level)/2, lower.tail = FALSE)
  Lower_CL <- Modavg_beta - zcrit*Uncond_SE
  Upper_CL <- Modavg_beta + zcrit*Uncond_SE
  out.modavg <- list("Parameter"=paste(parm), "Mod.avg.table" = new_table, "Mod.avg.beta" = Modavg_beta,
                     "Uncond.SE" = Uncond_SE, "Conf.level" = conf.level, "Lower.CL" = Lower_CL,
                     "Upper.CL" = Upper_CL)
  
  class(out.modavg) <- c("modavg", "list")
  return(out.modavg)
}



##distsamp
modavg.AICunmarkedFitDS <-
  function(cand.set, parm, modnames = NULL, second.ord = TRUE, nobs = NULL, 
           uncond.se = "revised", conf.level = 0.95, exclude = NULL, warn = TRUE,
           c.hat = 1, parm.type = NULL, ...){

    ##note that parameter is referenced differently from unmarked object - see labels( )

    ##check if named list if modnames are not supplied
    if(is.null(modnames)) {
      if(is.null(names(cand.set))) {
        modnames <- paste("Mod", 1:length(cand.set), sep = "")
        warning("\nModel names have been supplied automatically in the table\n")
      } else {
        modnames <- names(cand.set)
      }
    }
    
  
  ##check for parm.type and stop if NULL
  if(is.null(parm.type)) {stop("\n'parm.type' must be specified for this model type, see ?modavg for details\n")}

  
#####MODIFICATIONS BEGIN#######
  ##remove all leading and trailing white space and within parm
  parm <- gsub('[[:space:]]+', "", parm)
  parm.strip <- parm #to use later
  
  ##if (Intercept) is chosen assign (Int) - for compatibility
  if(identical(parm, "(Intercept)")) parm <- "Int"

  ##reverse parm
  reversed.parm <- reverse.parm(parm)
  reversed.parm.strip <- reversed.parm #to use later
  exclude <- reverse.exclude(exclude = exclude)
#####MODIFICATIONS END######
  
  

    ##Distance sampling model
    ##lambda - abundance
    if(identical(parm.type, "lambda")) {
      ##extract model formula for each model in cand.set
      mod_formula <- lapply(cand.set, FUN = function(i) labels(coef(i@estimates@estimates$state)))
      parm <- paste("lam", "(", parm, ")", sep="")
      if(!is.null(reversed.parm)) {reversed.parm <- paste("lam", "(", reversed.parm, ")", sep="")}
      not.include <- lapply(cand.set, FUN = function(i) i@formula[[3]])
    }
    ##detect
    if(identical(parm.type, "detect")) {
      if(identical(parm.type, "detect")) {
        stop("\nModel-averaging estimates of detection covariates not yet supported for unmarkedFitDS class\n")
      }
    }
    
    nmods <- length(cand.set)
  
    ##setup matrix to indicate presence of parms in the model
    include <- matrix(NA, nrow=nmods, ncol=1)
    ##add a check for multiple instances of same variable in given model (i.e., interactions)
    include.check <- matrix(NA, nrow=nmods, ncol=1)

##################################
##################################
###ADDED A NEW OBJECT TO STRIP AWAY lam( ) from parm on line 35
###to enable search with regexpr( ) 
  
  
    ##iterate over each formula in mod_formula list
    for (i in 1:nmods) {
      idents <- NULL
      idents.check <- NULL
      form <- mod_formula[[i]]
      

######################################################################################################
######################################################################################################
###MODIFICATIONS BEGIN
      ##iterate over each element of formula[[i]] in list
      if(is.null(reversed.parm)) {
        for (j in 1:length(form)) {
          idents[j] <- identical(parm, form[j])
          ##added parm.strip here for regexpr( )
          idents.check[j] <- ifelse(attr(regexpr(parm.strip, form[j], fixed=TRUE), "match.length")== "-1" , 0, 1)  
        }
      } else {
        for (j in 1:length(form)) {
          idents[j] <- identical(parm, form[j]) | identical(reversed.parm, form[j])
          idents.check[j] <- ifelse(attr(regexpr(parm.strip, form[j], fixed=TRUE), "match.length")=="-1" & attr(regexpr(reversed.parm.strip, form[j],
                                                                  fixed=TRUE), "match.length")=="-1" , 0, 1)  
        }
      }
###MODIFICATIONS END
######################################################################################################
######################################################################################################

    
    include[i] <- ifelse(any(idents==1), 1, 0)
    include.check[i] <- ifelse(sum(idents.check)>1, "duplicates", "OK")
  }

  #####################################################
  #exclude == NULL; warn=TRUE:  warn that duplicates occur and stop
  if(is.null(exclude) && identical(warn, TRUE)) {
    #check for duplicates in same model
    if(any(include.check == "duplicates")) {
      stop("\nSome models include more than one instance of the parameter of interest. \n",
           "This may be due to the presence of interaction/polynomial terms, or variables\n",
             "with similar names:\n",
             "\tsee \"?modavg\" for details on variable specification and \"exclude\" argument\n")
      }

    }

    #exclude == NULL; warn=FALSE:  compute model-averaged beta estimate from models including variable of interest,
    #assuming that the variable is not involved in interaction or higher order polynomial (x^2, x^3, etc...),
    #warn that models were not excluded
    if(is.null(exclude) && identical(warn, FALSE)) {
      if(any(include.check == "duplicates")) {
        warning("\nMultiple instances of parameter of interest in given model is presumably\n",
                "not due to interaction or polynomial terms - these models will not be\n",
                "excluded from the computation of model-averaged estimate\n")
      }
      
    }

    #warn if exclude is neither a list nor NULL
    if(!is.null(exclude)) {
      if(!is.list(exclude)) {stop("\nItems in \"exclude\" must be specified as a list")}
    }


    #if exclude is list  
    if(is.list(exclude)) {

    #determine number of elements in exclude
      nexcl <- length(exclude)

      #check each formula for presence of exclude variable in not.include list
      #not.include <- lapply(cand.set, FUN = formula)

      #set up a new list with model formula
      forms <- list()
      for (i in 1:nmods) {
        form.tmp <- as.character(not.include[i]) #changed from other versions as formula returned is of different structure for unmarked objects
        if(attr(regexpr("\\+", form.tmp), "match.length")==-1) {
          forms[i] <- form.tmp
        } else {forms[i] <- strsplit(form.tmp, split=" \\+ ")}
      }

      #additional check to see whether some variable names include "+"
      check.forms <- unlist(lapply(forms, FUN=function(i) any(attr(regexpr("\\+", i), "match.length")>0)[[1]]))
      if (any(check.forms==TRUE)) stop("\nPlease avoid \"+\" in variable names\n")

      ##additional check to determine if intercept was removed from models
      check.forms <- unlist(lapply(forms, FUN=function(i) any(attr(regexpr("\\- 1", i), "match.length")>0)[[1]]))
      if (any(check.forms==TRUE)) stop("\nModels without intercept are not supported in this version, please use alternative parameterization\n")

 
      #search within formula for variables to exclude
      mod.exclude <- matrix(NA, nrow=nmods, ncol=nexcl)

      #iterate over each element in exclude list
      for (var in 1:nexcl) {

      #iterate over each formula in mod_formula list
        for (i in 1:nmods) {
          idents <- NULL
          form.excl <- forms[[i]]

          #iterate over each element of forms[[i]]
          for (j in 1:length(form.excl)) {
            idents[j] <- identical(exclude[var][[1]], form.excl[j])
          }
          mod.exclude[i,var] <- ifelse(any(idents == 1), 1, 0)
        }    
        
      }
  
      #determine outcome across all variables to exclude
      to.exclude <- rowSums(mod.exclude)
  
  
      #exclude models following models from model averaging  
      include[which(to.exclude >= 1)] <- 0
      
      
    }


 
  ##add a check to determine if include always == 0
  if (sum(include) == 0) {stop("\nParameter not found in any of the candidate models\n") }
  
  new.cand.set <- cand.set[which(include == 1)] #select models including a given parameter
  new.mod.name <- modnames[which(include == 1)]    #update model names

  new_table <- aictab(cand.set = new.cand.set, modnames = new.mod.name,
                      second.ord = second.ord, nobs = nobs, sort = FALSE, c.hat = c.hat) #recompute AIC table and associated measures
  new_table$Beta_est <- unlist(lapply(new.cand.set, FUN = function(i) coef(i)[paste(parm)])) #extract beta estimate for parm
  ##if reversed.parm is not null and varies across models, potentially check for it here
  new_table$SE <- unlist(lapply(new.cand.set, FUN = function(i) sqrt(diag(vcov(i)))[paste(parm)]))
  ##if reversed.parm is not null and varies across models, potentially check for it here

  ##if c-hat is estimated adjust the SE's by multiplying with sqrt of c-hat
  if(c.hat > 1) {
    new_table$SE <- new_table$SE*sqrt(c.hat)
  } 

  ##AICc
  ##compute model-averaged estimates, unconditional SE, and 95% CL
  if(c.hat == 1 && second.ord == TRUE) {
    Modavg_beta <- sum(new_table$AICcWt*new_table$Beta_est)

    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Uncond_SE <- sum(new_table$AICcWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
    }
      
    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Uncond_SE <- sqrt(sum(new_table$AICcWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
    }
  }

  ##QAICc
  if(c.hat > 1 && second.ord == TRUE) {
    Modavg_beta <- sum(new_table$QAICcWt*new_table$Beta_est)
      
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {      
      Uncond_SE <- sum(new_table$QAICcWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
    }
    
    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Uncond_SE <- sqrt(sum(new_table$QAICcWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
    }
  }     


  ##AIC
  if(c.hat == 1 && second.ord == FALSE) {
    Modavg_beta <- sum(new_table$AICWt*new_table$Beta_est)
      
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Uncond_SE <- sum(new_table$AICWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
    }
      
    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Uncond_SE <- sqrt(sum(new_table$AICWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
    }
  }


  ##QAIC
  if(c.hat > 1 && second.ord == FALSE) {
    Modavg_beta <- sum(new_table$QAICWt*new_table$Beta_est)

    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Uncond_SE <- sum(new_table$QAICWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
    }

    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Uncond_SE <- sqrt(sum(new_table$QAICWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
    }  
  }     

  
  zcrit <- qnorm(p = (1 - conf.level)/2, lower.tail = FALSE)
  Lower_CL <- Modavg_beta - zcrit*Uncond_SE
  Upper_CL <- Modavg_beta + zcrit*Uncond_SE
  out.modavg <- list("Parameter"=paste(parm), "Mod.avg.table" = new_table, "Mod.avg.beta" = Modavg_beta,
                     "Uncond.SE" = Uncond_SE, "Conf.level" = conf.level, "Lower.CL" = Lower_CL,
                     "Upper.CL" = Upper_CL)
  
  class(out.modavg) <- c("modavg", "list")
  return(out.modavg)
}



##gdistsamp
modavg.AICunmarkedFitGDS <-
  function(cand.set, parm, modnames = NULL, second.ord = TRUE, nobs = NULL, 
           uncond.se = "revised", conf.level = 0.95, exclude = NULL, warn = TRUE,
           c.hat = 1, parm.type = NULL, ...){

  ##note that parameter is referenced differently from unmarked object - see labels( )

  ##check if named list if modnames are not supplied
  if(is.null(modnames)) {
    if(is.null(names(cand.set))) {
      modnames <- paste("Mod", 1:length(cand.set), sep = "")
      warning("\nModel names have been supplied automatically in the table\n")
    } else {
      modnames <- names(cand.set)
    }
  }
    
  
  ##check for parm.type and stop if NULL
  if(is.null(parm.type)) {stop("\n'parm.type' must be specified for this model type, see ?modavg for details\n")}

  
#####MODIFICATIONS BEGIN#######
  ##remove all leading and trailing white space and within parm
  parm <- gsub('[[:space:]]+', "", parm)
  parm.strip <- parm #to use later
  
  ##if (Intercept) is chosen assign (Int) - for compatibility
  if(identical(parm, "(Intercept)")) parm <- "Int"

  ##reverse parm
  reversed.parm <- reverse.parm(parm)
  reversed.parm.strip <- reversed.parm #to use later
  exclude <- reverse.exclude(exclude = exclude)
#####MODIFICATIONS END######
  
  ##Distance sampling model with availability
  ##lambda - abundance
  if(identical(parm.type, "lambda")) {
    ##extract model formula for each model in cand.set
    mod_formula <- lapply(cand.set, FUN = function(i) labels(coef(i@estimates@estimates$lambda)))
    parm <- paste("lambda", "(", parm, ")", sep="")
    if(!is.null(reversed.parm)) {reversed.parm <- paste("lam", "(", reversed.parm, ")", sep="")}
    not.include <- lapply(cand.set, FUN = function(i) i@formlist$lambdaformula)
  }
  ##detect
  if(identical(parm.type, "detect")) {
    if(identical(parm.type, "detect")) {
      stop("\nModel-averaging estimates of detection covariates not yet supported for unmarkedFitGDS class\n")
    }
    }
  ##availability
  if(identical(parm.type, "phi")) {
    ##extract model formula for each model in cand.set
    mod_formula <- lapply(cand.set, FUN = function(i) labels(coef(i@estimates@estimates$phi)))
    parm <- paste("phi", "(", parm, ")", sep="")
    if(!is.null(reversed.parm)) {reversed.parm <- paste("phi", "(", reversed.parm, ")", sep="")}
    not.include <- lapply(cand.set, FUN = function(i) i@formlist$phiformula)
  }
  
    
  nmods <- length(cand.set)
  
  ##setup matrix to indicate presence of parms in the model
  include <- matrix(NA, nrow=nmods, ncol=1)
  ##add a check for multiple instances of same variable in given model (i.e., interactions)
  include.check <- matrix(NA, nrow=nmods, ncol=1)

##################################
##################################
###ADDED A NEW OBJECT TO STRIP AWAY lam( ) from parm on line 35
###to enable search with regexpr( ) 
  
  
  ##iterate over each formula in mod_formula list
  for (i in 1:nmods) {
    idents <- NULL
    idents.check <- NULL
    form <- mod_formula[[i]]


######################################################################################################
######################################################################################################
###MODIFICATIONS BEGIN
      ##iterate over each element of formula[[i]] in list
      if(is.null(reversed.parm)) {
        for (j in 1:length(form)) {
          idents[j] <- identical(parm, form[j])
          ##added parm.strip here for regexpr( )
          idents.check[j] <- ifelse(attr(regexpr(parm.strip, form[j], fixed=TRUE), "match.length")== "-1" , 0, 1)  
        }
      } else {
        for (j in 1:length(form)) {
          idents[j] <- identical(parm, form[j]) | identical(reversed.parm, form[j])
          idents.check[j] <- ifelse(attr(regexpr(parm.strip, form[j], fixed=TRUE), "match.length")=="-1" & attr(regexpr(reversed.parm.strip, form[j],
                                                                  fixed=TRUE), "match.length")=="-1" , 0, 1)  
        }
      }
###MODIFICATIONS END
######################################################################################################
######################################################################################################

    
    include[i] <- ifelse(any(idents==1), 1, 0)
    include.check[i] <- ifelse(sum(idents.check)>1, "duplicates", "OK")
  }

  #####################################################
  #exclude == NULL; warn=TRUE:  warn that duplicates occur and stop
  if(is.null(exclude) && identical(warn, TRUE)) {
    #check for duplicates in same model
    if(any(include.check == "duplicates")) {
      stop("\nSome models include more than one instance of the parameter of interest. \n",
           "This may be due to the presence of interaction/polynomial terms, or variables\n",
             "with similar names:\n",
             "\tsee \"?modavg\" for details on variable specification and \"exclude\" argument\n")
      }

    }

    #exclude == NULL; warn=FALSE:  compute model-averaged beta estimate from models including variable of interest,
    #assuming that the variable is not involved in interaction or higher order polynomial (x^2, x^3, etc...),
    #warn that models were not excluded
    if(is.null(exclude) && identical(warn, FALSE)) {
      if(any(include.check == "duplicates")) {
        warning("\nMultiple instances of parameter of interest in given model is presumably\n",
                "not due to interaction or polynomial terms - these models will not be\n",
                "excluded from the computation of model-averaged estimate\n")
      }
      
    }

    #warn if exclude is neither a list nor NULL
    if(!is.null(exclude)) {
      if(!is.list(exclude)) {stop("\nItems in \"exclude\" must be specified as a list")}
    }


    #if exclude is list  
    if(is.list(exclude)) {

    #determine number of elements in exclude
      nexcl <- length(exclude)

      #check each formula for presence of exclude variable in not.include list
      #not.include <- lapply(cand.set, FUN = formula)

      #set up a new list with model formula
      forms <- list()
      for (i in 1:nmods) {
        form.tmp <- as.character(not.include[i]) #changed from other versions as formula returned is of different structure for unmarked objects
        if(attr(regexpr("\\+", form.tmp), "match.length")==-1) {
          forms[i] <- form.tmp
        } else {forms[i] <- strsplit(form.tmp, split=" \\+ ")}
      }

      #additional check to see whether some variable names include "+"
      check.forms <- unlist(lapply(forms, FUN=function(i) any(attr(regexpr("\\+", i), "match.length")>0)[[1]]))
      if (any(check.forms==TRUE)) stop("\nPlease avoid \"+\" in variable names\n")

      ##additional check to determine if intercept was removed from models
      check.forms <- unlist(lapply(forms, FUN=function(i) any(attr(regexpr("\\- 1", i), "match.length")>0)[[1]]))
      if (any(check.forms==TRUE)) stop("\nModels without intercept are not supported in this version, please use alternative parameterization\n")

 
      #search within formula for variables to exclude
      mod.exclude <- matrix(NA, nrow=nmods, ncol=nexcl)

      #iterate over each element in exclude list
      for (var in 1:nexcl) {

      #iterate over each formula in mod_formula list
        for (i in 1:nmods) {
          idents <- NULL
          form.excl <- forms[[i]]

          #iterate over each element of forms[[i]]
          for (j in 1:length(form.excl)) {
            idents[j] <- identical(exclude[var][[1]], form.excl[j])
          }
          mod.exclude[i,var] <- ifelse(any(idents == 1), 1, 0)
        }    
        
      }
  
      #determine outcome across all variables to exclude
      to.exclude <- rowSums(mod.exclude)
  
  
      #exclude models following models from model averaging  
      include[which(to.exclude >= 1)] <- 0
      
      
    }


 
  ##add a check to determine if include always == 0
  if (sum(include) == 0) {stop("\nParameter not found in any of the candidate models\n") }
  
  new.cand.set <- cand.set[which(include == 1)] #select models including a given parameter
  new.mod.name <- modnames[which(include == 1)]    #update model names

  new_table <- aictab(cand.set = new.cand.set, modnames = new.mod.name,
                      second.ord = second.ord, nobs = nobs, sort = FALSE, c.hat = c.hat) #recompute AIC table and associated measures
  new_table$Beta_est <- unlist(lapply(new.cand.set, FUN = function(i) coef(i)[paste(parm)])) #extract beta estimate for parm
  ##if reversed.parm is not null and varies across models, potentially check for it here
  new_table$SE <- unlist(lapply(new.cand.set, FUN = function(i) sqrt(diag(vcov(i)))[paste(parm)]))
  ##if reversed.parm is not null and varies across models, potentially check for it here

  ##if c-hat is estimated adjust the SE's by multiplying with sqrt of c-hat
  if(c.hat > 1) {
    new_table$SE <- new_table$SE*sqrt(c.hat)
  } 

  ##AICc
  ##compute model-averaged estimates, unconditional SE, and 95% CL
  if(c.hat == 1 && second.ord == TRUE) {
    Modavg_beta <- sum(new_table$AICcWt*new_table$Beta_est)

    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Uncond_SE <- sum(new_table$AICcWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
    }
      
    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Uncond_SE <- sqrt(sum(new_table$AICcWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
    }
  }

  ##QAICc
  if(c.hat > 1 && second.ord == TRUE) {
    Modavg_beta <- sum(new_table$QAICcWt*new_table$Beta_est)
      
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {      
      Uncond_SE <- sum(new_table$QAICcWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
    }
    
    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Uncond_SE <- sqrt(sum(new_table$QAICcWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
    }
  }     


  ##AIC
  if(c.hat == 1 && second.ord == FALSE) {
    Modavg_beta <- sum(new_table$AICWt*new_table$Beta_est)
      
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Uncond_SE <- sum(new_table$AICWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
    }
      
    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Uncond_SE <- sqrt(sum(new_table$AICWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
    }
  }


  ##QAIC
  if(c.hat > 1 && second.ord == FALSE) {
    Modavg_beta <- sum(new_table$QAICWt*new_table$Beta_est)

    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Uncond_SE <- sum(new_table$QAICWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
    }

    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Uncond_SE <- sqrt(sum(new_table$QAICWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
    }  
  }     

  
  zcrit <- qnorm(p = (1 - conf.level)/2, lower.tail = FALSE)
  Lower_CL <- Modavg_beta - zcrit*Uncond_SE
  Upper_CL <- Modavg_beta + zcrit*Uncond_SE
  out.modavg <- list("Parameter"=paste(parm), "Mod.avg.table" = new_table, "Mod.avg.beta" = Modavg_beta,
                     "Uncond.SE" = Uncond_SE, "Conf.level" = conf.level, "Lower.CL" = Lower_CL,
                     "Upper.CL" = Upper_CL)
  
  class(out.modavg) <- c("modavg", "list")
  return(out.modavg)
}



##occuFP
modavg.AICunmarkedFitOccuFP <-
  function(cand.set, parm, modnames = NULL, second.ord = TRUE, nobs = NULL, 
           uncond.se = "revised", conf.level = 0.95, exclude = NULL, warn = TRUE,
           c.hat = 1, parm.type = NULL, ...){

    ##note that parameter is referenced differently from unmarked object - see labels( )

    ##check if named list if modnames are not supplied
    if(is.null(modnames)) {
      if(is.null(names(cand.set))) {
        modnames <- paste("Mod", 1:length(cand.set), sep = "")
        warning("\nModel names have been supplied automatically in the table\n")
      } else {
        modnames <- names(cand.set)
      }
    }
    
  
    ##check for parm.type and stop if NULL
    if(is.null(parm.type)) {stop("\n'parm.type' must be specified for this model type, see ?modavg for details\n")}
    
  
#####MODIFICATIONS BEGIN#######
    ##remove all leading and trailing white space and within parm
    parm <- gsub('[[:space:]]+', "", parm)
    parm.strip <- parm #to use later
  
    ##if (Intercept) is chosen assign (Int) - for compatibility
    if(identical(parm, "(Intercept)")) parm <- "Int"

    ##reverse parm
    reversed.parm <- reverse.parm(parm)
    reversed.parm.strip <- reversed.parm #to use later
    exclude <- reverse.exclude(exclude = exclude)
#####MODIFICATIONS END######
  
  
    ##single-season false-positive occupancy model
    ##psi
    if(identical(parm.type, "psi")) {
      ##extract model formula for each model in cand.set
      mod_formula <- lapply(cand.set, FUN = function(i) labels(coef(i@estimates@estimates$state)))
      parm <- paste("psi", "(", parm, ")", sep="")
      if(!is.null(reversed.parm)) {reversed.parm <- paste("psi", "(", reversed.parm, ")", sep="")}
      not.include <- lapply(cand.set, FUN = function(i) i@stateformula)      
    }
    ##detect
    if(identical(parm.type, "detect")) {
      mod_formula <- lapply(cand.set, FUN = function(i) labels(coef(i@estimates@estimates$det)))
      parm <- paste("p", "(", parm, ")", sep="")
      if(!is.null(reversed.parm)) {reversed.parm <- paste("p", "(", reversed.parm, ")", sep="")}
      not.include <- lapply(cand.set, FUN = function(i) i@detformula)
    }
    ##false positives - fp
    if(identical(parm.type, "fp")) {
      mod_formula <- lapply(cand.set, FUN = function(i) labels(coef(i@estimates@estimates$fp)))
      parm <- paste("fp", "(", parm, ")", sep="")
      if(!is.null(reversed.parm)) {reversed.parm <- paste("fp", "(", reversed.parm, ")", sep="")}
      not.include <- lapply(cand.set, FUN = function(i) i@FPformula)
    }
  
  
  nmods <- length(cand.set)
  
  ##setup matrix to indicate presence of parms in the model
  include <- matrix(NA, nrow=nmods, ncol=1)
  ##add a check for multiple instances of same variable in given model (i.e., interactions)
  include.check <- matrix(NA, nrow=nmods, ncol=1)

##################################
##################################
###ADDED A NEW OBJECT TO STRIP AWAY lam( ) from parm on line 35
###to enable search with regexpr( ) 
  
  
  ##iterate over each formula in mod_formula list
  for (i in 1:nmods) {
    idents <- NULL
    idents.check <- NULL
    form <- mod_formula[[i]]


######################################################################################################
######################################################################################################
###MODIFICATIONS BEGIN
      ##iterate over each element of formula[[i]] in list
      if(is.null(reversed.parm)) {
        for (j in 1:length(form)) {
          idents[j] <- identical(parm, form[j])
          ##added parm.strip here for regexpr( )
          idents.check[j] <- ifelse(attr(regexpr(parm.strip, form[j], fixed=TRUE), "match.length")== "-1" , 0, 1)  
        }
      } else {
        for (j in 1:length(form)) {
          idents[j] <- identical(parm, form[j]) | identical(reversed.parm, form[j])
          idents.check[j] <- ifelse(attr(regexpr(parm.strip, form[j], fixed=TRUE), "match.length")=="-1" & attr(regexpr(reversed.parm.strip, form[j],
                                                                  fixed=TRUE), "match.length")=="-1" , 0, 1)  
        }
      }
###MODIFICATIONS END
######################################################################################################
######################################################################################################

    
    include[i] <- ifelse(any(idents==1), 1, 0)
    include.check[i] <- ifelse(sum(idents.check)>1, "duplicates", "OK")
  }

  #####################################################
  #exclude == NULL; warn=TRUE:  warn that duplicates occur and stop
  if(is.null(exclude) && identical(warn, TRUE)) {
    #check for duplicates in same model
    if(any(include.check == "duplicates")) {
      stop("\nSome models include more than one instance of the parameter of interest. \n",
           "This may be due to the presence of interaction/polynomial terms, or variables\n",
             "with similar names:\n",
             "\tsee \"?modavg\" for details on variable specification and \"exclude\" argument\n")
      }

    }

    #exclude == NULL; warn=FALSE:  compute model-averaged beta estimate from models including variable of interest,
    #assuming that the variable is not involved in interaction or higher order polynomial (x^2, x^3, etc...),
    #warn that models were not excluded
    if(is.null(exclude) && identical(warn, FALSE)) {
      if(any(include.check == "duplicates")) {
        warning("\nMultiple instances of parameter of interest in given model is presumably\n",
                "not due to interaction or polynomial terms - these models will not be\n",
                "excluded from the computation of model-averaged estimate\n")
      }
      
    }

    #warn if exclude is neither a list nor NULL
    if(!is.null(exclude)) {
      if(!is.list(exclude)) {stop("\nItems in \"exclude\" must be specified as a list")}
    }


    #if exclude is list  
    if(is.list(exclude)) {

    #determine number of elements in exclude
      nexcl <- length(exclude)

      #check each formula for presence of exclude variable in not.include list
      #not.include <- lapply(cand.set, FUN = formula)

      #set up a new list with model formula
      forms <- list()
      for (i in 1:nmods) {
        form.tmp <- as.character(not.include[i]) #changed from other versions as formula returned is of different structure for unmarked objects
        if(attr(regexpr("\\+", form.tmp), "match.length")==-1) {
          forms[i] <- form.tmp
        } else {forms[i] <- strsplit(form.tmp, split=" \\+ ")}
      }

      #additional check to see whether some variable names include "+"
      check.forms <- unlist(lapply(forms, FUN=function(i) any(attr(regexpr("\\+", i), "match.length")>0)[[1]]))
      if (any(check.forms==TRUE)) stop("\nPlease avoid \"+\" in variable names\n")

      ##additional check to determine if intercept was removed from models
      check.forms <- unlist(lapply(forms, FUN=function(i) any(attr(regexpr("\\- 1", i), "match.length")>0)[[1]]))
      if (any(check.forms==TRUE)) stop("\nModels without intercept are not supported in this version, please use alternative parameterization\n")

 
      #search within formula for variables to exclude
      mod.exclude <- matrix(NA, nrow=nmods, ncol=nexcl)

      #iterate over each element in exclude list
      for (var in 1:nexcl) {

      #iterate over each formula in mod_formula list
        for (i in 1:nmods) {
          idents <- NULL
          form.excl <- forms[[i]]

          #iterate over each element of forms[[i]]
          for (j in 1:length(form.excl)) {
            idents[j] <- identical(exclude[var][[1]], form.excl[j])
          }
          mod.exclude[i,var] <- ifelse(any(idents == 1), 1, 0)
        }    
        
      }
  
      #determine outcome across all variables to exclude
      to.exclude <- rowSums(mod.exclude)
  
  
      #exclude models following models from model averaging  
      include[which(to.exclude >= 1)] <- 0
      
      
    }


 
  ##add a check to determine if include always == 0
  if (sum(include) == 0) {stop("\nParameter not found in any of the candidate models\n") }
  
  new.cand.set <- cand.set[which(include == 1)] #select models including a given parameter
  new.mod.name <- modnames[which(include == 1)]    #update model names

  new_table <- aictab(cand.set = new.cand.set, modnames = new.mod.name,
                      second.ord = second.ord, nobs = nobs, sort = FALSE, c.hat = c.hat) #recompute AIC table and associated measures
  new_table$Beta_est <- unlist(lapply(new.cand.set, FUN = function(i) coef(i)[paste(parm)])) #extract beta estimate for parm
  ##if reversed.parm is not null and varies across models, potentially check for it here
  new_table$SE <- unlist(lapply(new.cand.set, FUN = function(i) sqrt(diag(vcov(i)))[paste(parm)]))
  ##if reversed.parm is not null and varies across models, potentially check for it here

  ##if c-hat is estimated adjust the SE's by multiplying with sqrt of c-hat
  if(c.hat > 1) {
    new_table$SE <- new_table$SE*sqrt(c.hat)
  } 

  ##AICc
  ##compute model-averaged estimates, unconditional SE, and 95% CL
  if(c.hat == 1 && second.ord == TRUE) {
    Modavg_beta <- sum(new_table$AICcWt*new_table$Beta_est)

    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Uncond_SE <- sum(new_table$AICcWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
    }
      
    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Uncond_SE <- sqrt(sum(new_table$AICcWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
    }
  }

  ##QAICc
  if(c.hat > 1 && second.ord == TRUE) {
    Modavg_beta <- sum(new_table$QAICcWt*new_table$Beta_est)
      
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {      
      Uncond_SE <- sum(new_table$QAICcWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
    }
    
    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Uncond_SE <- sqrt(sum(new_table$QAICcWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
    }
  }     


  ##AIC
  if(c.hat == 1 && second.ord == FALSE) {
    Modavg_beta <- sum(new_table$AICWt*new_table$Beta_est)
      
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Uncond_SE <- sum(new_table$AICWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
    }
      
    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Uncond_SE <- sqrt(sum(new_table$AICWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
    }
  }


  ##QAIC
  if(c.hat > 1 && second.ord == FALSE) {
    Modavg_beta <- sum(new_table$QAICWt*new_table$Beta_est)

    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Uncond_SE <- sum(new_table$QAICWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
    }

    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Uncond_SE <- sqrt(sum(new_table$QAICWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
    }  
  }     

  
  zcrit <- qnorm(p = (1 - conf.level)/2, lower.tail = FALSE)
  Lower_CL <- Modavg_beta - zcrit*Uncond_SE
  Upper_CL <- Modavg_beta + zcrit*Uncond_SE
  out.modavg <- list("Parameter"=paste(parm), "Mod.avg.table" = new_table, "Mod.avg.beta" = Modavg_beta,
                     "Uncond.SE" = Uncond_SE, "Conf.level" = conf.level, "Lower.CL" = Lower_CL,
                     "Upper.CL" = Upper_CL)
  
  class(out.modavg) <- c("modavg", "list")
  return(out.modavg)
}



##multinomPois
modavg.AICunmarkedFitMPois <-
  function(cand.set, parm, modnames = NULL, second.ord = TRUE, nobs = NULL, 
           uncond.se = "revised", conf.level = 0.95, exclude = NULL, warn = TRUE,
           c.hat = 1, parm.type = NULL, ...){

    ##note that parameter is referenced differently from unmarked object - see labels( )

    ##check if named list if modnames are not supplied
    if(is.null(modnames)) {
    if(is.null(names(cand.set))) {
      modnames <- paste("Mod", 1:length(cand.set), sep = "")
      warning("\nModel names have been supplied automatically in the table\n")
    } else {
      modnames <- names(cand.set)
    }
  }
    
  
    ##check for parm.type and stop if NULL
    if(is.null(parm.type)) {stop("\n'parm.type' must be specified for this model type, see ?modavg for details\n")}

  
#####MODIFICATIONS BEGIN#######
    ##remove all leading and trailing white space and within parm
    parm <- gsub('[[:space:]]+', "", parm)
    parm.strip <- parm #to use later
    
    ##if (Intercept) is chosen assign (Int) - for compatibility
    if(identical(parm, "(Intercept)")) parm <- "Int"
    
    ##reverse parm
    reversed.parm <- reverse.parm(parm)
    reversed.parm.strip <- reversed.parm #to use later
    exclude <- reverse.exclude(exclude = exclude)
#####MODIFICATIONS END######
  
  
    ##multinomPois model
    ##lambda - abundance
    if(identical(parm.type, "lambda")) {
      ##extract model formula for each model in cand.set
      mod_formula <- lapply(cand.set, FUN = function(i) labels(coef(i@estimates@estimates$state)))
      ##create label for parm
      parm <- paste("lambda", "(", parm, ")", sep="")
      if(!is.null(reversed.parm)) {reversed.parm <- paste("lambda", "(", reversed.parm, ")", sep="")}
      not.include <- lapply(cand.set, FUN = function(i) i@formula[[3]])
    }
    ##detect
    if(identical(parm.type, "detect")) {
      mod_formula <- lapply(cand.set, FUN = function(i) labels(coef(i@estimates@estimates$det)))
      parm <- paste("p", "(", parm, ")", sep="")
      if(!is.null(reversed.parm)) {reversed.parm <- paste("p", "(", reversed.parm, ")", sep="")}
      not.include <- lapply(cand.set, FUN = function(i) i@formula[[2]])
    }
  
  
  nmods <- length(cand.set)
  
  ##setup matrix to indicate presence of parms in the model
  include <- matrix(NA, nrow=nmods, ncol=1)
  ##add a check for multiple instances of same variable in given model (i.e., interactions)
  include.check <- matrix(NA, nrow=nmods, ncol=1)

##################################
##################################
###ADDED A NEW OBJECT TO STRIP AWAY lam( ) from parm on line 35
###to enable search with regexpr( ) 
  
  
  ##iterate over each formula in mod_formula list
  for (i in 1:nmods) {
    idents <- NULL
    idents.check <- NULL
    form <- mod_formula[[i]]


######################################################################################################
######################################################################################################
###MODIFICATIONS BEGIN
      ##iterate over each element of formula[[i]] in list
      if(is.null(reversed.parm)) {
        for (j in 1:length(form)) {
          idents[j] <- identical(parm, form[j])
          ##added parm.strip here for regexpr( )
          idents.check[j] <- ifelse(attr(regexpr(parm.strip, form[j], fixed=TRUE), "match.length")== "-1" , 0, 1)  
        }
      } else {
        for (j in 1:length(form)) {
          idents[j] <- identical(parm, form[j]) | identical(reversed.parm, form[j])
          idents.check[j] <- ifelse(attr(regexpr(parm.strip, form[j], fixed=TRUE), "match.length")=="-1" & attr(regexpr(reversed.parm.strip, form[j],
                                                                  fixed=TRUE), "match.length")=="-1" , 0, 1)  
        }
      }
###MODIFICATIONS END
######################################################################################################
######################################################################################################

    
    include[i] <- ifelse(any(idents==1), 1, 0)
    include.check[i] <- ifelse(sum(idents.check)>1, "duplicates", "OK")
  }

  #####################################################
  #exclude == NULL; warn=TRUE:  warn that duplicates occur and stop
  if(is.null(exclude) && identical(warn, TRUE)) {
    #check for duplicates in same model
    if(any(include.check == "duplicates")) {
      stop("\nSome models include more than one instance of the parameter of interest. \n",
           "This may be due to the presence of interaction/polynomial terms, or variables\n",
             "with similar names:\n",
             "\tsee \"?modavg\" for details on variable specification and \"exclude\" argument\n")
      }

    }

    #exclude == NULL; warn=FALSE:  compute model-averaged beta estimate from models including variable of interest,
    #assuming that the variable is not involved in interaction or higher order polynomial (x^2, x^3, etc...),
    #warn that models were not excluded
    if(is.null(exclude) && identical(warn, FALSE)) {
      if(any(include.check == "duplicates")) {
        warning("\nMultiple instances of parameter of interest in given model is presumably\n",
                "not due to interaction or polynomial terms - these models will not be\n",
                "excluded from the computation of model-averaged estimate\n")
      }
      
    }

    #warn if exclude is neither a list nor NULL
    if(!is.null(exclude)) {
      if(!is.list(exclude)) {stop("\nItems in \"exclude\" must be specified as a list")}
    }


    #if exclude is list  
    if(is.list(exclude)) {

    #determine number of elements in exclude
      nexcl <- length(exclude)

      #check each formula for presence of exclude variable in not.include list
      #not.include <- lapply(cand.set, FUN = formula)

      #set up a new list with model formula
      forms <- list()
      for (i in 1:nmods) {
        form.tmp <- as.character(not.include[i]) #changed from other versions as formula returned is of different structure for unmarked objects
        if(attr(regexpr("\\+", form.tmp), "match.length")==-1) {
          forms[i] <- form.tmp
        } else {forms[i] <- strsplit(form.tmp, split=" \\+ ")}
      }

      #additional check to see whether some variable names include "+"
      check.forms <- unlist(lapply(forms, FUN=function(i) any(attr(regexpr("\\+", i), "match.length")>0)[[1]]))
      if (any(check.forms==TRUE)) stop("\nPlease avoid \"+\" in variable names\n")

      ##additional check to determine if intercept was removed from models
      check.forms <- unlist(lapply(forms, FUN=function(i) any(attr(regexpr("\\- 1", i), "match.length")>0)[[1]]))
      if (any(check.forms==TRUE)) stop("\nModels without intercept are not supported in this version, please use alternative parameterization\n")

 
      #search within formula for variables to exclude
      mod.exclude <- matrix(NA, nrow=nmods, ncol=nexcl)

      #iterate over each element in exclude list
      for (var in 1:nexcl) {

      #iterate over each formula in mod_formula list
        for (i in 1:nmods) {
          idents <- NULL
          form.excl <- forms[[i]]

          #iterate over each element of forms[[i]]
          for (j in 1:length(form.excl)) {
            idents[j] <- identical(exclude[var][[1]], form.excl[j])
          }
          mod.exclude[i,var] <- ifelse(any(idents == 1), 1, 0)
        }    
        
      }
  
      #determine outcome across all variables to exclude
      to.exclude <- rowSums(mod.exclude)
  
  
      #exclude models following models from model averaging  
      include[which(to.exclude >= 1)] <- 0
      
      
    }


 
  ##add a check to determine if include always == 0
  if (sum(include) == 0) {stop("\nParameter not found in any of the candidate models\n") }
  
  new.cand.set <- cand.set[which(include == 1)] #select models including a given parameter
  new.mod.name <- modnames[which(include == 1)]    #update model names

  new_table <- aictab(cand.set = new.cand.set, modnames = new.mod.name,
                      second.ord = second.ord, nobs = nobs, sort = FALSE, c.hat = c.hat) #recompute AIC table and associated measures
  new_table$Beta_est <- unlist(lapply(new.cand.set, FUN = function(i) coef(i)[paste(parm)])) #extract beta estimate for parm
  ##if reversed.parm is not null and varies across models, potentially check for it here
  new_table$SE <- unlist(lapply(new.cand.set, FUN = function(i) sqrt(diag(vcov(i)))[paste(parm)]))
  ##if reversed.parm is not null and varies across models, potentially check for it here

  ##if c-hat is estimated adjust the SE's by multiplying with sqrt of c-hat
  if(c.hat > 1) {
    new_table$SE <- new_table$SE*sqrt(c.hat)
  } 

  ##AICc
  ##compute model-averaged estimates, unconditional SE, and 95% CL
  if(c.hat == 1 && second.ord == TRUE) {
    Modavg_beta <- sum(new_table$AICcWt*new_table$Beta_est)

    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Uncond_SE <- sum(new_table$AICcWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
    }
      
    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Uncond_SE <- sqrt(sum(new_table$AICcWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
    }
  }

  ##QAICc
  if(c.hat > 1 && second.ord == TRUE) {
    Modavg_beta <- sum(new_table$QAICcWt*new_table$Beta_est)
      
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {      
      Uncond_SE <- sum(new_table$QAICcWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
    }
    
    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Uncond_SE <- sqrt(sum(new_table$QAICcWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
    }
  }     


  ##AIC
  if(c.hat == 1 && second.ord == FALSE) {
    Modavg_beta <- sum(new_table$AICWt*new_table$Beta_est)
      
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Uncond_SE <- sum(new_table$AICWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
    }
      
    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Uncond_SE <- sqrt(sum(new_table$AICWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
    }
  }


  ##QAIC
  if(c.hat > 1 && second.ord == FALSE) {
    Modavg_beta <- sum(new_table$QAICWt*new_table$Beta_est)

    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Uncond_SE <- sum(new_table$QAICWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
    }

    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Uncond_SE <- sqrt(sum(new_table$QAICWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
    }  
  }     

  
  zcrit <- qnorm(p = (1 - conf.level)/2, lower.tail = FALSE)
  Lower_CL <- Modavg_beta - zcrit*Uncond_SE
  Upper_CL <- Modavg_beta + zcrit*Uncond_SE
  out.modavg <- list("Parameter"=paste(parm), "Mod.avg.table" = new_table, "Mod.avg.beta" = Modavg_beta,
                     "Uncond.SE" = Uncond_SE, "Conf.level" = conf.level, "Lower.CL" = Lower_CL,
                     "Upper.CL" = Upper_CL)
  
  class(out.modavg) <- c("modavg", "list")
  return(out.modavg)
}



##gmultmix
modavg.AICunmarkedFitGMM <-
  function(cand.set, parm, modnames = NULL, second.ord = TRUE, nobs = NULL, 
           uncond.se = "revised", conf.level = 0.95, exclude = NULL, warn = TRUE,
           c.hat = 1, parm.type = NULL, ...){

    ##note that parameter is referenced differently from unmarked object - see labels( )
    
    ##check if named list if modnames are not supplied
    if(is.null(modnames)) {
      if(is.null(names(cand.set))) {
        modnames <- paste("Mod", 1:length(cand.set), sep = "")
        warning("\nModel names have been supplied automatically in the table\n")
      } else {
        modnames <- names(cand.set)
      }
    }
    
  
    ##check for parm.type and stop if NULL
    if(is.null(parm.type)) {stop("\n'parm.type' must be specified for this model type, see ?modavg for details\n")}
    
  
#####MODIFICATIONS BEGIN#######
    ##remove all leading and trailing white space and within parm
    parm <- gsub('[[:space:]]+', "", parm)
    parm.strip <- parm #to use later
  
    ##if (Intercept) is chosen assign (Int) - for compatibility
    if(identical(parm, "(Intercept)")) parm <- "Int"
    
    ##reverse parm
    reversed.parm <- reverse.parm(parm)
    reversed.parm.strip <- reversed.parm #to use later
    exclude <- reverse.exclude(exclude = exclude)
#####MODIFICATIONS END######
  
    ##gmultmix model
    ##lambda - abundance
    if(identical(parm.type, "lambda")) {
      ##extract model formula for each model in cand.set
      mod_formula <- lapply(cand.set, FUN = function(i) labels(coef(i@estimates@estimates$lambda)))
      parm <- paste("lambda", "(", parm, ")", sep="")
      if(!is.null(reversed.parm)) {reversed.parm <- paste("lambda", "(", reversed.parm, ")", sep="")}
      not.include <- lapply(cand.set, FUN = function(i) i@formlist$lambdaformula)
    }
    ##detect
    if(identical(parm.type, "detect")) {
      mod_formula <- lapply(cand.set, FUN = function(i) labels(coef(i@estimates@estimates$det)))
      parm <- paste("p", "(", parm, ")", sep="")
      if(!is.null(reversed.parm)) {reversed.parm <- paste("p", "(", reversed.parm, ")", sep="")}
      not.include <- lapply(cand.set, FUN = function(i) i@formlist$pformula)
    }
    ##availability
    if(identical(parm.type, "phi")) {
      mod_formula <- lapply(cand.set, FUN = function(i) labels(coef(i@estimates@estimates$phi)))
      parm <- paste("phi", "(", parm, ")", sep="")
      if(!is.null(reversed.parm)) {reversed.parm <- paste("phi", "(", reversed.parm, ")", sep="")}
      not.include <- lapply(cand.set, FUN = function(i) i@formlist$phiformula)
    }
  
  
  nmods <- length(cand.set)
  
  ##setup matrix to indicate presence of parms in the model
  include <- matrix(NA, nrow=nmods, ncol=1)
  ##add a check for multiple instances of same variable in given model (i.e., interactions)
  include.check <- matrix(NA, nrow=nmods, ncol=1)

##################################
##################################
###ADDED A NEW OBJECT TO STRIP AWAY lam( ) from parm on line 35
###to enable search with regexpr( ) 
  
  
  ##iterate over each formula in mod_formula list
  for (i in 1:nmods) {
    idents <- NULL
    idents.check <- NULL
    form <- mod_formula[[i]]


######################################################################################################
######################################################################################################
###MODIFICATIONS BEGIN
      ##iterate over each element of formula[[i]] in list
      if(is.null(reversed.parm)) {
        for (j in 1:length(form)) {
          idents[j] <- identical(parm, form[j])
          ##added parm.strip here for regexpr( )
          idents.check[j] <- ifelse(attr(regexpr(parm.strip, form[j], fixed=TRUE), "match.length")== "-1" , 0, 1)  
        }
      } else {
        for (j in 1:length(form)) {
          idents[j] <- identical(parm, form[j]) | identical(reversed.parm, form[j])
          idents.check[j] <- ifelse(attr(regexpr(parm.strip, form[j], fixed=TRUE), "match.length")=="-1" & attr(regexpr(reversed.parm.strip, form[j],
                                                                  fixed=TRUE), "match.length")=="-1" , 0, 1)  
        }
      }
###MODIFICATIONS END
######################################################################################################
######################################################################################################

    
    include[i] <- ifelse(any(idents==1), 1, 0)
    include.check[i] <- ifelse(sum(idents.check)>1, "duplicates", "OK")
  }

  #####################################################
  #exclude == NULL; warn=TRUE:  warn that duplicates occur and stop
  if(is.null(exclude) && identical(warn, TRUE)) {
    #check for duplicates in same model
    if(any(include.check == "duplicates")) {
      stop("\nSome models include more than one instance of the parameter of interest. \n",
           "This may be due to the presence of interaction/polynomial terms, or variables\n",
             "with similar names:\n",
             "\tsee \"?modavg\" for details on variable specification and \"exclude\" argument\n")
      }

    }

    #exclude == NULL; warn=FALSE:  compute model-averaged beta estimate from models including variable of interest,
    #assuming that the variable is not involved in interaction or higher order polynomial (x^2, x^3, etc...),
    #warn that models were not excluded
    if(is.null(exclude) && identical(warn, FALSE)) {
      if(any(include.check == "duplicates")) {
        warning("\nMultiple instances of parameter of interest in given model is presumably\n",
                "not due to interaction or polynomial terms - these models will not be\n",
                "excluded from the computation of model-averaged estimate\n")
      }
      
    }

    #warn if exclude is neither a list nor NULL
    if(!is.null(exclude)) {
      if(!is.list(exclude)) {stop("\nItems in \"exclude\" must be specified as a list")}
    }


    #if exclude is list  
    if(is.list(exclude)) {

    #determine number of elements in exclude
      nexcl <- length(exclude)

      #check each formula for presence of exclude variable in not.include list
      #not.include <- lapply(cand.set, FUN = formula)

      #set up a new list with model formula
      forms <- list()
      for (i in 1:nmods) {
        form.tmp <- as.character(not.include[i]) #changed from other versions as formula returned is of different structure for unmarked objects
        if(attr(regexpr("\\+", form.tmp), "match.length")==-1) {
          forms[i] <- form.tmp
        } else {forms[i] <- strsplit(form.tmp, split=" \\+ ")}
      }

      #additional check to see whether some variable names include "+"
      check.forms <- unlist(lapply(forms, FUN=function(i) any(attr(regexpr("\\+", i), "match.length")>0)[[1]]))
      if (any(check.forms==TRUE)) stop("\nPlease avoid \"+\" in variable names\n")

      ##additional check to determine if intercept was removed from models
      check.forms <- unlist(lapply(forms, FUN=function(i) any(attr(regexpr("\\- 1", i), "match.length")>0)[[1]]))
      if (any(check.forms==TRUE)) stop("\nModels without intercept are not supported in this version, please use alternative parameterization\n")

 
      #search within formula for variables to exclude
      mod.exclude <- matrix(NA, nrow=nmods, ncol=nexcl)

      #iterate over each element in exclude list
      for (var in 1:nexcl) {

      #iterate over each formula in mod_formula list
        for (i in 1:nmods) {
          idents <- NULL
          form.excl <- forms[[i]]

          #iterate over each element of forms[[i]]
          for (j in 1:length(form.excl)) {
            idents[j] <- identical(exclude[var][[1]], form.excl[j])
          }
          mod.exclude[i,var] <- ifelse(any(idents == 1), 1, 0)
        }    
        
      }
  
      #determine outcome across all variables to exclude
      to.exclude <- rowSums(mod.exclude)
  
  
      #exclude models following models from model averaging  
      include[which(to.exclude >= 1)] <- 0
      
      
    }


 
  ##add a check to determine if include always == 0
  if (sum(include) == 0) {stop("\nParameter not found in any of the candidate models\n") }
  
  new.cand.set <- cand.set[which(include == 1)] #select models including a given parameter
  new.mod.name <- modnames[which(include == 1)]    #update model names

  new_table <- aictab(cand.set = new.cand.set, modnames = new.mod.name,
                      second.ord = second.ord, nobs = nobs, sort = FALSE, c.hat = c.hat) #recompute AIC table and associated measures
  new_table$Beta_est <- unlist(lapply(new.cand.set, FUN = function(i) coef(i)[paste(parm)])) #extract beta estimate for parm
  ##if reversed.parm is not null and varies across models, potentially check for it here
  new_table$SE <- unlist(lapply(new.cand.set, FUN = function(i) sqrt(diag(vcov(i)))[paste(parm)]))
  ##if reversed.parm is not null and varies across models, potentially check for it here

  ##if c-hat is estimated adjust the SE's by multiplying with sqrt of c-hat
  if(c.hat > 1) {
    new_table$SE <- new_table$SE*sqrt(c.hat)
  } 

  ##AICc
  ##compute model-averaged estimates, unconditional SE, and 95% CL
  if(c.hat == 1 && second.ord == TRUE) {
    Modavg_beta <- sum(new_table$AICcWt*new_table$Beta_est)

    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Uncond_SE <- sum(new_table$AICcWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
    }
      
    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Uncond_SE <- sqrt(sum(new_table$AICcWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
    }
  }

  ##QAICc
  if(c.hat > 1 && second.ord == TRUE) {
    Modavg_beta <- sum(new_table$QAICcWt*new_table$Beta_est)
      
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {      
      Uncond_SE <- sum(new_table$QAICcWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
    }
    
    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Uncond_SE <- sqrt(sum(new_table$QAICcWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
    }
  }     


  ##AIC
  if(c.hat == 1 && second.ord == FALSE) {
    Modavg_beta <- sum(new_table$AICWt*new_table$Beta_est)
      
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Uncond_SE <- sum(new_table$AICWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
    }
      
    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Uncond_SE <- sqrt(sum(new_table$AICWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
    }
  }


  ##QAIC
  if(c.hat > 1 && second.ord == FALSE) {
    Modavg_beta <- sum(new_table$QAICWt*new_table$Beta_est)

    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Uncond_SE <- sum(new_table$QAICWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
    }

    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Uncond_SE <- sqrt(sum(new_table$QAICWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
    }  
  }     

  
  zcrit <- qnorm(p = (1 - conf.level)/2, lower.tail = FALSE)
  Lower_CL <- Modavg_beta - zcrit*Uncond_SE
  Upper_CL <- Modavg_beta + zcrit*Uncond_SE
  out.modavg <- list("Parameter"=paste(parm), "Mod.avg.table" = new_table, "Mod.avg.beta" = Modavg_beta,
                     "Uncond.SE" = Uncond_SE, "Conf.level" = conf.level, "Lower.CL" = Lower_CL,
                     "Upper.CL" = Upper_CL)
  
  class(out.modavg) <- c("modavg", "list")
  return(out.modavg)
}



##gpcount
modavg.AICunmarkedFitGPC <-
  function(cand.set, parm, modnames = NULL, second.ord = TRUE, nobs = NULL, 
           uncond.se = "revised", conf.level = 0.95, exclude = NULL, warn = TRUE,
           c.hat = 1, parm.type = NULL, ...){

    ##note that parameter is referenced differently from unmarked object - see labels( )
    
    ##check if named list if modnames are not supplied
    if(is.null(modnames)) {
      if(is.null(names(cand.set))) {
        modnames <- paste("Mod", 1:length(cand.set), sep = "")
        warning("\nModel names have been supplied automatically in the table\n")
      } else {
        modnames <- names(cand.set)
      }
    }
    
  
    ##check for parm.type and stop if NULL
    if(is.null(parm.type)) {stop("\n'parm.type' must be specified for this model type, see ?modavg for details\n")}

  
#####MODIFICATIONS BEGIN#######
    ##remove all leading and trailing white space and within parm
    parm <- gsub('[[:space:]]+', "", parm)
    parm.strip <- parm #to use later
  
    ##if (Intercept) is chosen assign (Int) - for compatibility
    if(identical(parm, "(Intercept)")) parm <- "Int"

    ##reverse parm
    reversed.parm <- reverse.parm(parm)
    reversed.parm.strip <- reversed.parm #to use later
    exclude <- reverse.exclude(exclude = exclude)
#####MODIFICATIONS END######
  
  
  ##gpcount
    ##lambda - abundance
    if(identical(parm.type, "lambda")) {
      ##extract model formula for each model in cand.set
      mod_formula <- lapply(cand.set, FUN = function(i) labels(coef(i@estimates@estimates$lambda)))
      parm <- paste("lambda", "(", parm, ")", sep="")
      if(!is.null(reversed.parm)) {reversed.parm <- paste("lambda", "(", reversed.parm, ")", sep="")}
      not.include <- lapply(cand.set, FUN = function(i) i@formlist$lambdaformula)
    }
    ##detect
    if(identical(parm.type, "detect")) {
      mod_formula <- lapply(cand.set, FUN = function(i) labels(coef(i@estimates@estimates$det)))
      parm <- paste("p", "(", parm, ")", sep="")
      if(!is.null(reversed.parm)) {reversed.parm <- paste("p", "(", reversed.parm, ")", sep="")}
      not.include <- lapply(cand.set, FUN = function(i) i@formlist$pformula)
    }
    ##availability
    if(identical(parm.type, "phi")) {
      mod_formula <- lapply(cand.set, FUN = function(i) labels(coef(i@estimates@estimates$phi)))
      parm <- paste("phi", "(", parm, ")", sep="")
      if(!is.null(reversed.parm)) {reversed.parm <- paste("phi", "(", reversed.parm, ")", sep="")}
      not.include <- lapply(cand.set, FUN = function(i) i@formlist$phiformula)
    }
    
  
  
  nmods <- length(cand.set)
  
  ##setup matrix to indicate presence of parms in the model
  include <- matrix(NA, nrow=nmods, ncol=1)
  ##add a check for multiple instances of same variable in given model (i.e., interactions)
  include.check <- matrix(NA, nrow=nmods, ncol=1)

##################################
##################################
###ADDED A NEW OBJECT TO STRIP AWAY lam( ) from parm on line 35
###to enable search with regexpr( ) 
  
  
  ##iterate over each formula in mod_formula list
  for (i in 1:nmods) {
    idents <- NULL
    idents.check <- NULL
    form <- mod_formula[[i]]


######################################################################################################
######################################################################################################
###MODIFICATIONS BEGIN
      ##iterate over each element of formula[[i]] in list
      if(is.null(reversed.parm)) {
        for (j in 1:length(form)) {
          idents[j] <- identical(parm, form[j])
          ##added parm.strip here for regexpr( )
          idents.check[j] <- ifelse(attr(regexpr(parm.strip, form[j], fixed=TRUE), "match.length")== "-1" , 0, 1)  
        }
      } else {
        for (j in 1:length(form)) {
          idents[j] <- identical(parm, form[j]) | identical(reversed.parm, form[j])
          idents.check[j] <- ifelse(attr(regexpr(parm.strip, form[j], fixed=TRUE), "match.length")=="-1" & attr(regexpr(reversed.parm.strip, form[j],
                                                                  fixed=TRUE), "match.length")=="-1" , 0, 1)  
        }
      }
###MODIFICATIONS END
######################################################################################################
######################################################################################################

    
    include[i] <- ifelse(any(idents==1), 1, 0)
    include.check[i] <- ifelse(sum(idents.check)>1, "duplicates", "OK")
  }

  #####################################################
  #exclude == NULL; warn=TRUE:  warn that duplicates occur and stop
  if(is.null(exclude) && identical(warn, TRUE)) {
    #check for duplicates in same model
    if(any(include.check == "duplicates")) {
      stop("\nSome models include more than one instance of the parameter of interest. \n",
           "This may be due to the presence of interaction/polynomial terms, or variables\n",
             "with similar names:\n",
             "\tsee \"?modavg\" for details on variable specification and \"exclude\" argument\n")
      }

    }

    #exclude == NULL; warn=FALSE:  compute model-averaged beta estimate from models including variable of interest,
    #assuming that the variable is not involved in interaction or higher order polynomial (x^2, x^3, etc...),
    #warn that models were not excluded
    if(is.null(exclude) && identical(warn, FALSE)) {
      if(any(include.check == "duplicates")) {
        warning("\nMultiple instances of parameter of interest in given model is presumably\n",
                "not due to interaction or polynomial terms - these models will not be\n",
                "excluded from the computation of model-averaged estimate\n")
      }
      
    }

    #warn if exclude is neither a list nor NULL
    if(!is.null(exclude)) {
      if(!is.list(exclude)) {stop("\nItems in \"exclude\" must be specified as a list")}
    }


    #if exclude is list  
    if(is.list(exclude)) {

    #determine number of elements in exclude
      nexcl <- length(exclude)

      #check each formula for presence of exclude variable in not.include list
      #not.include <- lapply(cand.set, FUN = formula)

      #set up a new list with model formula
      forms <- list()
      for (i in 1:nmods) {
        form.tmp <- as.character(not.include[i]) #changed from other versions as formula returned is of different structure for unmarked objects
        if(attr(regexpr("\\+", form.tmp), "match.length")==-1) {
          forms[i] <- form.tmp
        } else {forms[i] <- strsplit(form.tmp, split=" \\+ ")}
      }

      #additional check to see whether some variable names include "+"
      check.forms <- unlist(lapply(forms, FUN=function(i) any(attr(regexpr("\\+", i), "match.length")>0)[[1]]))
      if (any(check.forms==TRUE)) stop("\nPlease avoid \"+\" in variable names\n")

      ##additional check to determine if intercept was removed from models
      check.forms <- unlist(lapply(forms, FUN=function(i) any(attr(regexpr("\\- 1", i), "match.length")>0)[[1]]))
      if (any(check.forms==TRUE)) stop("\nModels without intercept are not supported in this version, please use alternative parameterization\n")

 
      #search within formula for variables to exclude
      mod.exclude <- matrix(NA, nrow=nmods, ncol=nexcl)

      #iterate over each element in exclude list
      for (var in 1:nexcl) {

      #iterate over each formula in mod_formula list
        for (i in 1:nmods) {
          idents <- NULL
          form.excl <- forms[[i]]

          #iterate over each element of forms[[i]]
          for (j in 1:length(form.excl)) {
            idents[j] <- identical(exclude[var][[1]], form.excl[j])
          }
          mod.exclude[i,var] <- ifelse(any(idents == 1), 1, 0)
        }    
        
      }
  
      #determine outcome across all variables to exclude
      to.exclude <- rowSums(mod.exclude)
  
  
      #exclude models following models from model averaging  
      include[which(to.exclude >= 1)] <- 0
      
      
    }


 
  ##add a check to determine if include always == 0
  if (sum(include) == 0) {stop("\nParameter not found in any of the candidate models\n") }
  
  new.cand.set <- cand.set[which(include == 1)] #select models including a given parameter
  new.mod.name <- modnames[which(include == 1)]    #update model names

  new_table <- aictab(cand.set = new.cand.set, modnames = new.mod.name,
                      second.ord = second.ord, nobs = nobs, sort = FALSE, c.hat = c.hat) #recompute AIC table and associated measures
  new_table$Beta_est <- unlist(lapply(new.cand.set, FUN = function(i) coef(i)[paste(parm)])) #extract beta estimate for parm
  ##if reversed.parm is not null and varies across models, potentially check for it here
  new_table$SE <- unlist(lapply(new.cand.set, FUN = function(i) sqrt(diag(vcov(i)))[paste(parm)]))
  ##if reversed.parm is not null and varies across models, potentially check for it here

  ##if c-hat is estimated adjust the SE's by multiplying with sqrt of c-hat
  if(c.hat > 1) {
    new_table$SE <- new_table$SE*sqrt(c.hat)
  } 

  ##AICc
  ##compute model-averaged estimates, unconditional SE, and 95% CL
  if(c.hat == 1 && second.ord == TRUE) {
    Modavg_beta <- sum(new_table$AICcWt*new_table$Beta_est)

    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Uncond_SE <- sum(new_table$AICcWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
    }
      
    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Uncond_SE <- sqrt(sum(new_table$AICcWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
    }
  }

  ##QAICc
  if(c.hat > 1 && second.ord == TRUE) {
    Modavg_beta <- sum(new_table$QAICcWt*new_table$Beta_est)
      
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {      
      Uncond_SE <- sum(new_table$QAICcWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
    }
    
    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Uncond_SE <- sqrt(sum(new_table$QAICcWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
    }
  }     


  ##AIC
  if(c.hat == 1 && second.ord == FALSE) {
    Modavg_beta <- sum(new_table$AICWt*new_table$Beta_est)
      
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Uncond_SE <- sum(new_table$AICWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
    }
      
    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Uncond_SE <- sqrt(sum(new_table$AICWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
    }
  }


  ##QAIC
  if(c.hat > 1 && second.ord == FALSE) {
    Modavg_beta <- sum(new_table$QAICWt*new_table$Beta_est)

    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Uncond_SE <- sum(new_table$QAICWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
    }

    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Uncond_SE <- sqrt(sum(new_table$QAICWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
    }  
  }     

  
  zcrit <- qnorm(p = (1 - conf.level)/2, lower.tail = FALSE)
  Lower_CL <- Modavg_beta - zcrit*Uncond_SE
  Upper_CL <- Modavg_beta + zcrit*Uncond_SE
  out.modavg <- list("Parameter"=paste(parm), "Mod.avg.table" = new_table, "Mod.avg.beta" = Modavg_beta,
                     "Uncond.SE" = Uncond_SE, "Conf.level" = conf.level, "Lower.CL" = Lower_CL,
                     "Upper.CL" = Upper_CL)
  
  class(out.modavg) <- c("modavg", "list")
  return(out.modavg)
}





##print method
print.modavg <-
  function(x, digits = 2, ...) {
    ic <- colnames(x$Mod.avg.table)[3]
    cat("\nMultimodel inference on \"", x$Parameter, "\" based on ", ic, "\n", sep = "")
    cat("\n", ic, " table used to obtain model-averaged estimate:\n", sep = "")
    oldtab <- x$Mod.avg.table
    if (any(names(oldtab)=="c_hat")) {cat("\t(c-hat estimate = ", oldtab$c_hat[1], ")\n", sep = "")}
    cat("\n")
    if (any(names(oldtab)=="c_hat")) {
      nice.tab <- cbind(oldtab[,2], oldtab[,3], oldtab[,4], oldtab[,6],
                        oldtab[,9], oldtab[,10])
    } else {nice.tab <- cbind(oldtab[,2], oldtab[,3], oldtab[,4], oldtab[,6],
                              oldtab[,8], oldtab[,9])
          }

##modify printing style if multinomial model is used  
    if(length(x$Mod.avg.beta)==1) {
      colnames(nice.tab) <- c(colnames(oldtab)[c(2,3,4,6)], "Estimate", "SE")
      rownames(nice.tab) <- oldtab[,1]
      print(round(nice.tab, digits=digits))
      cat("\nModel-averaged estimate:", eval(round(x$Mod.avg.beta, digits=digits)), "\n")
      cat("Unconditional SE:", eval(round(x$Uncond.SE, digits=digits)), "\n")
      cat("",x$Conf.level*100, "% Unconditional confidence interval: ", round(x$Lower.CL, digits=digits),
          ", ", round(x$Upper.CL, digits=digits), "\n\n", sep = "")
    } else {
      col.ns <- ncol(nice.tab)
      nice.tab <- nice.tab[,-c(col.ns-1,col.ns)]
      colnames(nice.tab) <- c(colnames(oldtab)[c(2,3,4,6)])
      rownames(nice.tab) <- oldtab[,1]
      print(round(nice.tab, digits=digits))
      cat("\n\nModel-averaged estimates for different levels of response variable:", "\n\n")
      resp.labels <- labels(x$Mod.avg.beta)
      mult.out <- matrix(NA, nrow=length(resp.labels), ncol=4)
      colnames(mult.out) <- c("Model-averaged estimate", "Uncond. SE", paste(x$Conf.level*100,"% lower CL", sep = ""),
                              paste(x$Conf.level*100, "% upper CL", sep = ""))
      rownames(mult.out) <- resp.labels
      mult.out[,1] <- round(x$Mod.avg.beta, digits=digits)
      mult.out[,2] <- round(x$Uncond.SE, digits=digits)
      mult.out[,3] <- round(x$Lower.CL, digits=digits)
      mult.out[,4] <- round(x$Upper.CL, digits=digits)
      print(mult.out)
      cat("\n")
    }
  }
