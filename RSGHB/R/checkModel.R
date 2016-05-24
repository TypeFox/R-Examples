checkModel <- function(nodiagnostics = FALSE, verbose = TRUE, env = parent.frame())
{
     
     passChecks <- TRUE

     # model checks to make sure things are coded properly
     
     # the user needs to specify some variables. lets make sure they exist.
     
     if (is.null(env$gDIST) & env$gNIV > 0)
     {
          passChecks <- FALSE
          cat("\n********FATAL ERROR: Variable - gDIST - is undefined.\n")          
     }
     if (is.null(env$gNCREP))
     {
          passChecks <- FALSE
          cat("\n********FATAL ERROR: Variable - gNCREP - is undefined.\n")
     }
     if (is.null(env$gNEREP))
     {
          passChecks <- FALSE
          cat("\n********FATAL ERROR: Variable - gNEREP - is undefined.\n")
     }
     if (is.null(env$gNSKIP))
     {
          passChecks <- FALSE
          cat("\n********FATAL ERROR: Variable - gNSKIP - is undefined.\n")
     }
     if (is.null(env$gINFOSKIP))
     {
          passChecks <- FALSE
          cat("\n********FATAL ERROR: Variable - gINFOSKIP - is undefined.\n")
     }
     if (is.null(env$likelihood))
     {
          passChecks <- FALSE
          cat("\n********FATAL ERROR: The likelihood function is undefined.\n")
     }
     
     if (env$gNIV + env$gFIV == 0)
     {
          passChecks <- FALSE
          cat("\n********FATAL ERROR: Please specify at least one coefficient to be estimated in either in gVarNamesNormal or gVarNamesFixed.\n")
     }     
     # check to make sure that the number of distributions specified equals the 
     # the number of random coefficients in the model
     if (length(env$gDIST) != env$gNIV)
     {
          passChecks <- FALSE
          cat("\n********FATAL ERROR: The number of distributions specified in gDist doesn't equal the number of random coefficients in the model.\n")
     }
     # check to the see if the distributions specified exists in the set of allowable distributions
     if (env$gNIV > 0)
     {
          for(d in env$gDIST)
          {
               if(d < 1 | d > length(env$distNames))
               {
                    passChecks <- FALSE
                    cat("\n********FATAL ERROR: The specified distributions ", d, " in gDist do not exist\n")          
               }     
          }
     }
     # check to see if we have enough starting values for both the random and fixed coefficients
     if (env$gNIV != length(env$svN))
     {
          passChecks <- FALSE    
          cat("\n********FATAL ERROR: There are too many/not enough starting values for the random coefficients. Check your sVN vector.\n")          
     }
     if (env$gFIV != length(env$FC))
     {
          passChecks <- FALSE    
          cat("\n********FATAL ERROR: There are too many/not enough starting values for the fixed coefficients. Check your FC vector.\n")          
     }
     # the software assumes that there exists a respondent identifier is called ID
     if (is.null(env$choicedata$ID))
     {
          passChecks <- FALSE
          cat("\n********FATAL ERROR: Expecting to find a respondent identifier column called - ID - in your dataset. None found.\n")          
     }         

     # the software needs the data sorted by ID
     if (sum(sort(env$choicedata$ID) == env$choicedata$ID) != length(env$choicedata$ID))
     {
          passChecks <- FALSE    
          cat("\n********FATAL ERROR: The choice data is not sorted by ID.\n")               
     }
     
     if ((!is.null(env$fixedA)) & length(env$fixedA) != length(env$gVarNamesNormal))
     {
          passChecks <- FALSE    
          cat("\n********FATAL ERROR: The fixedA vector is not of the same length as the gVarNamesNormal vector.\n") 
     }
     
     if ((!is.null(env$fixedD)) & length(env$fixedD) != length(env$gVarNamesNormal))
     {
          passChecks <- FALSE    
          cat("\n********FATAL ERROR: The fixedD vector is not of the same length as the gVarNamesNormal vector.\n") 
     }
     
     if (passChecks) prepareModel(env)
     
     if (passChecks & verbose)
     {
          
          cat(rep("\n", 128))
          cat("Diagnostic checks passed. Please review before proceeding\n")
          diagnostics <- data.frame(` ` = c("Number of Individuals:",
                                            "Number of Observations:",
                                            "Custom Prior Matrix Used:",
                                            "Prior variance:",
                                            "Target Acceptance (Fixed):",
                                            "Target Acceptance (Normal):",
                                            "Degrees of Freedom:",
                                            "Avg. Number of Observations per Individual:",
                                            "Initial Log-Likelihood:"),
                                    ` ` = as.character(rep(NA, 9)),
                                    check.names = FALSE, stringsAsFactors = FALSE)
          
          diagnostics[1,2] <- env$gNP
          diagnostics[2,2] <- env$gNOBS
          if (env$useCustomPVMatrix) {
               diagnostics[3,2] <- TRUE
          } else {
               diagnostics[4,2] <- signif(env$priorVariance, env$gSIGDIG)
          }
          if (env$gFIV > 0) diagnostics[5,2]  <- signif(env$targetAcceptanceFixed, env$gSIGDIG)
          if (env$gNIV > 0) diagnostics[6,2] <- signif(env$targetAcceptanceNormal, env$gSIGDIG)
          diagnostics[7,2] <- signif(env$degreesOfFreedom, env$gSIGDIG)
          diagnostics[8,2] <- signif(env$gNOBS / env$gNP, env$gSIGDIG)
          diagnostics[9,2] <- signif(sum(log(env$likelihood(env$FC, env$B, env))), env$gSIGDIG)
          cat("-----------------------------------------------------------\n")
          print(diagnostics[complete.cases(diagnostics), , drop = FALSE], row.names = FALSE)
          cat("\n-----------------------------------------------------------\n\n")
          
   	     if (env$gFIV > 0)
	     {
               print(data.frame(`Fixed Parameters` = env$gVarNamesFixed, Start = env$FC, check.names = FALSE), row.names = FALSE)
	  	     cat("\n-----------------------------------------------------------\n\n")
	     }
          
	     if (env$gNIV > 0)
	     {
	          print(data.frame(`Random Parameters` = env$gVarNamesNormal, Start = env$svN, `Dist.` = env$distNames[env$gDIST], check.names = FALSE), row.names = FALSE)
               cat("\n-----------------------------------------------------------\n\n")
          }
          
          if (!is.null(env$constraintsNorm))
          {
               cat("Constraints applied to random parameters:\n")
               
               diagconstraints <- data.frame(` ` = unlist(lapply(env$constraintsNorm, FUN = `[`, 1)),
                                             ` ` = unlist(lapply(env$constraintsNorm, FUN = `[`, 2)),
                                             ` ` = unlist(lapply(env$constraintsNorm, FUN = `[`, 3)),
                                             check.names = FALSE)
               
               diagconstraints[, 1] <- env$gVarNamesNormal[diagconstraints[, 1]]
               diagconstraints[, 2] <- env$constraintLabels[diagconstraints[, 2]]
               diagconstraints[diagconstraints[, 3] != 0, 3] <- env$gVarNamesNormal[diagconstraints[diagconstraints[, 3] != 0, 3]]
               print(diagconstraints, row.names = FALSE, right = FALSE)
               cat("\n-----------------------------------------------------------")
          }
          
          if (!is.null(env$Choice))
          {
               cat("\n", "Choice Matrix", "\n")
               choiceMatrix <- cbind(table(env$Choice), round(prop.table(table(env$Choice)), 2))
          
               dimnames(choiceMatrix)[[2]] <- c("Count", "%")
               print(choiceMatrix)
          }
          
          cat("\n\n\n")
          if (!nodiagnostics & verbose)
          {
               rl <- readline("Estimate Model? (Y/N): ")
               if(rl != "Y" & rl != "y") passChecks <- FALSE
          }
     }
     return(passChecks)
}
