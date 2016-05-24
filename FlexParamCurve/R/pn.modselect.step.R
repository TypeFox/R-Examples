pn.modselect.step <-
structure(function # Backwards Stepwise Selection of Positive-Negative Richards \eqn{nlslist} Models
                               (x,
                                ### a numeric vector of the primary predictor
                                y,
                                ### a numeric vector of the response variable
                                grp,
                                ### a factor of same length as x and y that distinguishes groups within
                                ### the dataset
          			pn.options,
          	                ### required character value for name of list object populated with starting 
          	                ### parameter estimates, fitting 
          	                ### options and bounds or destination for modpar to write a new list (see Details)
                                forcemod = 0,
                                ### optional numeric value to constrain model selection (see Details)
                                existing = FALSE,
                                ### optional logical value specifying whether some of the relevant models
                                ### have already been fitted
                                penaliz = "1/sqrt(n)",
                                ### optional character value to determine how models are ranked (see Details)
                                taper.ends = 0.45,
                                Envir = .GlobalEnv,
                                ...
                                ) {
##decription<< This function performs backawards stepwise model selection for \code{\link{nlsList}}
## models fitted using
## \code{\link{SSposnegRichards}}. 
## details<< First, whether parameter M should be fixed
## (see \code{\link{SSposnegRichards}}) is determined by fitting models 12 and 20 and comparing
## their perfomance using \code{\link{extraF}}.
## If model 12 provides superior performance (variable values of  M) then 16 models that estimate M
## are run
## (models 1 through 16), otherwise the models with fixed M are fitted (models 21 through 36).
## Model selection then proceeds by fitting the most general model (8-parameter, model 1 for variable M;
## 7-parameter, model 21 for fixed M). At each subsequent step a more reduced model is evaluated
	## by creating \code{\link{nlsList}} models through removal of a single parameter from the decreasing
	## section of the curve (i.e. RAsym, Rk, Ri or RM). This is repeated until all possible models with
	## one less parameter have been fitted and then these models are then ranked by modified pooled residual
	## standard error (see below) to determine which reduced parameter model provides the best fit.
	## The best reduced parameter model is then compared with the more general model retained from the
	## the previous step using the function \code{\link{extraF}} to determine whether the more general
	## model provides significant improvement over the best reduced model. The most appropriate model
	## is then retained to be used as the general model in the next step. This process continues
	## for up to six steps (all steps will be attempted even if the general model provides better
	## performance to allow for much more reduced models to also be evaluated). The most reduced model
	## possible to evaluate in this function contains only parameters for the positive section of the curve
	## (4-parameters for variable M, 3-parameters for fixed M).
	##
	## Fitting these \code{\link{nlsList}} models can be time-consuming (2-4 hours using the dataset
	## \code{\link{posneg.data}} that encompasses 100 individuals) and if several of the relevant
	## models are already fitted the option existing=TRUE can be used to avoid refitting models that
	## already exist (note that a model object in which no grouping levels were successfully
	## parameterized will be refitted, as will objects that are not of class nlsList).
	##
	## Specifying forcemod=3 will force model selection to only consider fixed M models and setting
	## forcemod=4 will force model selection to consider models with varying values of M only.
	## If fitting both models
	## 12 and 20 fails, fixed M models will be used by default.
	##
	## Models are ranked by modified pooled residual square error. By default residual standard error
	## is divided by the square root of sample size. This exponentially penalizes models for which very few
	## grouping levels (individuals) are successfully parameterized (the few individuals that are
	## parameterized in these models are fit unsuprisingly well) using a function based on the relationship
	## between standard error and sample size. However, different users may have different preferences
	## and these can be specified in the argument penaliz (which residual
	## standard error is multiplied by). This argument must be a character value
	## that contains the character n (sample size) and must be a valid right hand side (RHS) of a formula:
	## e.g. 1*(n), (n)^2. It cannot contain more than one n but could be a custom function, e.g. FUN(n).
	    pcklib<- FPCEnv
	    options( warn = -1)
	    pnoptm=NULL
	    pnoptnm <- as.character(pn.options)
	    checkpen <- try(unlist(strsplit(penaliz, "(n)")), silent = TRUE)
	    if (length(checkpen) != 2 | class(checkpen)[1] == "try-error") {
		stop("penaliz parameter is ill defined: see ?pn.mod.compare")
	    } else {
		checkpen <- try(eval(parse(text = sprintf("%s", paste(checkpen[1],
		    "1", checkpen[2], sep = "")))))
		if (class(checkpen)[1] == "try-error")
		    stop("penaliz parameter is ill defined: see ?pn.mod.compare")
	    }
	    datamerg <- data.frame(x, y, grp)
	    userdata <- groupedData(y ~ x | grp, outer = ~grp, data = datamerg)
	    assign("userdata", userdata, envir = Envir)
	    if(as.character(list(Envir)) != "<environment>") stop ("No such environment")
            if(exists("Envir", mode = "environment") == FALSE) stop ("No such environment")
   	    FPCEnv$env <- Envir
	    testbounds <- 1
	    testpar <- 1
	    is.na(testbounds) <- TRUE
	    is.na(testpar) <- TRUE
	    testbounds <- try(get(pnoptnm, envir = Envir)[16:32],
		silent = T)
	    testpar <- try(get(pnoptnm, envir = Envir)[1:15],
		silent = T)
	    if (class(testbounds)[1] == "try-error" | class(testpar)[1] ==
		"try-error" | is.na(testbounds[1]) == TRUE | is.na(testpar[1]) ==
		TRUE)
       	    try({
       	        FPCEnv$mod.sel <- TRUE
         	modpar(datamerg[,1], datamerg[,2], pn.options = pnoptnm, taper.ends = taper.ends,
        	verbose=FALSE, Envir = Envir, ...)
        	options(warn=-1)
        	try(rm("mod.sel", envir = FPCEnv), silent =T)
        	options(warn=0)
        	}, silent = FALSE)    
	    extraF <- try(get("extraF", pos = 1), silent = TRUE)
	    if (class(extraF)[1] == "try-error") {
		stop("cannot find function: extraF - please reload FlexParamCurve")
	    }
	    mostreducedmod<-1
	    print("checking fit of positive section of the curve for variable M*************************************")
	    richardsR12.lis <- try(FPCEnv$richardsR12.lis,
		silent = TRUE)
	    if (class(richardsR12.lis)[1] == "try-error" | existing ==
		FALSE)
		richardsR12.lis <- eval(parse(text=sprintf("%s",paste("try(nlsList(y ~ SSposnegRichards(x,
		            Asym = Asym, K = K, Infl = Infl, M = M, modno = 12, pn.options = ",pnoptnm, "), data = userdata),
            silent = TRUE)",sep=""))))
	    print("checking fit of positive section of the curve for fixed M*************************************")
	    pnmodelparams <- get(pnoptnm, envir = Envir)[1:15]
	    change.pnparameters <- try(get("change.pnparameters", pos = 1),
		silent = TRUE)
   	    chk <- try(unlist(summary(richardsR12.lis))["RSE"], silent = TRUE)
	    richardsR20.lis <- try(FPCEnv$richardsR20.lis,
		silent = TRUE)
	    if (class(richardsR20.lis)[1] == "try-error" | existing ==
		FALSE)
		richardsR20.lis <- eval(parse(text=sprintf("%s",paste("try(nlsList(y ~ SSposnegRichards(x,
		            Asym = Asym, K = K, Infl = Infl, modno = 20, pn.options = ",pnoptnm,"), data = userdata),
            silent = TRUE)",sep=""))))
	    chk1 <- try(unlist(summary(richardsR20.lis))["RSE"], silent = TRUE)     
	    if ((class(richardsR20.lis)[1]) == "try-error" | class(richardsR20.lis)[[1]] != "nlsList" 
	    	| class(chk1)[1] == "try-error") {
	        print("3 parameter positive richards model failed/not fitted*************************************")
                if(forcemod != 3) forcemod = 4
		richardsR20.lis <- 1
		} else {
		FPCEnv$richardsR20.lis <- richardsR20.lis
	    					}
	    if ((class(richardsR12.lis)[1]) == "try-error" | class(richardsR12.lis)[[1]] != "nlsList" 
	    	| class(chk)[1] == "try-error")
		{
                print("4 parameter positive richards model failed/not fitted*************************************")
                if(forcemod != 4)  forcemod = 3
		    richardsR12.lis <- 1
		} else {
		FPCEnv$richardsR12.lis <- richardsR12.lis
					}
	    currentmodel <- 1
	    testmod <- try(extraF(richardsR20.lis, richardsR12.lis, warn = F))
	    if (forcemod == 0) {
				if (class(testmod) == "try-error") {
		    modelsig = 0.1
		} else {
		    modelsig = testmod[4]
		    if ((testmod[4]) > 0.05 & sqrt(testmod[5]/(testmod[3]-testmod[2])) > sqrt(testmod[6]/testmod[3])) {
			currentmodel <- richardsR20.lis
			mostreducednm <- substr("richardsR20.lis", 10,
			  11)
		    } else {
			currentmodel <- richardsR12.lis
			mostreducednm <- substr("richardsR12.lis", 10,
			  11)
		    }
		}
	    }
            mostreducedmod <- currentmodel
            if (class(testmod) != "try-error") {
            mostreducednm <- substr("richardsR20.lis", 10,
                  11)
            mostreducedmod <- richardsR20.lis
            } else {
            mostreducednm <- "NONE"
            }
	    if (forcemod == 3)
		{
		    modelsig = 0.1
		}
	    if (forcemod == 4)
		{
		    modelsig = 0.04
		}
	    if (modelsig < 0.05) {
		print("Variable M models most appropriate*************************************")
	    } else {
		print("Fixed M models most appropriate*************************************")
	    }
	    options(warn = 0)
	    options(warn = -1)
	    rankmod <- function(model1 = 1, model2 = 1, model3 = 1, model4 = 1) {
		nm <- rep(0, 4)
		nm[1] <- (as.character(substitute(model1)))
		nm[2] <- (as.character(substitute(model2)))
		nm[3] <- (as.character(substitute(model3)))
		nm[4] <- (as.character(substitute(model4)))
		if (class(model1)[[1]] == "nlsList" & class(model1)[1] !=
		    "try-error") {
		    if (is.null(nrow(coef(model1))) == TRUE) {
			model1 <- 1
		    }
		}
		if (class(model2)[[1]] == "nlsList" & class(model2)[1] !=
		    "try-error") {
		    if (is.null(nrow(coef(model2))) == TRUE) {
			model2 <- 1
		    }
		}
		if (class(model3)[[1]] == "nlsList" & class(model3)[1] !=
		    "try-error") {
		    if (is.null(nrow(coef(model3))) == TRUE) {
			model3 <- 1
		    }
		}
		if (class(model4)[[1]] == "nlsList" & class(model4)[1] !=
		    "try-error") {
		    if (is.null(nrow(coef(model4))) == TRUE) {
			model4 <- 1
		    }
		}
		modrank <- data.frame(modno = c(1, 2, 3, 4), rank = rep(-999,
		    4))
		nomods <- 4
		RSEstr <- "RSE"
		dfstr <- "df"
		usefun <- unlist(strsplit(penaliz, "(n)"))
		if (class(model1)[[1]] == "nlsList" & class(model1)[1] !=
		    "try-error") {
        		evfun <- parse(text = sprintf("%s", paste("summary(model1)[['",
           		 RSEstr, "']]*(", usefun[1], "(1+sum( summary(model1)[['",
           		 dfstr, "']],na.rm=TRUE)))", usefun[2], sep = "")))
		    modrank[1, 2] <- eval(evfun)
		} else {
		    nomods = nomods - 1
		}
		if (class(model2)[[1]] == "nlsList" & class(model2)[1] !=
		    "try-error") {
      			  evfun <- parse(text = sprintf("%s", paste("summary(model2)[['",
         			RSEstr, "']]*(", usefun[1], "(1+sum( summary(model2)[['",
          		 	 dfstr, "']],na.rm=TRUE)))", usefun[2], sep = "")))
		    modrank[2, 2] <- eval(evfun)
		} else {
		    nomods = nomods - 1
		}
		if (class(model3)[[1]] == "nlsList" & class(model3)[1] !=
		    "try-error") {
		          evfun <- parse(text = sprintf("%s", paste("summary(model3)[['",
		               RSEstr, "']]*(", usefun[1], "(1+sum( summary(model3)[['",
           			dfstr, "']],na.rm=TRUE)))", usefun[2], sep = "")))
		} else {
		    nomods = nomods - 1
		}
		if (class(model4)[[1]] == "nlsList" & class(model4)[1] !=
		    "try-error") {
        		evfun <- parse(text = sprintf("%s", paste("summary(model4)[['",
         		   RSEstr, "']]*(", usefun[1], "(1+sum( summary(model4)[['",
          		  dfstr, "']],na.rm=TRUE)))", usefun[2], sep = "")))
            modrank[4, 2] <- eval(evfun)
		} else {
		    nomods = nomods - 1
		}
		if (nomods == 0) {
		    for (j in 1:4) {
			if (modrank[j, 2] == -999 & nm[j] != 1)
			  print(paste("Model", nm[j], "failed to converge for all individuals",
			    sep = " "))
		    }
		    return(print("*************************************no models to evaluate*************************************"))
		} else {
		    modnmsav = ""
		    for (j in 1:4) {
			if (modrank[j, 2] == -999 & nm[j] != 1)
			  print(paste("Model", nm[j], "failed to converge for all individuals",
			    sep = " "))
			if (j == 1)
			  if (modrank[j, 2] != -999)
			    modnmsav <- nm[j]
			if (j > 1)
			  if (modrank[j, 2] != -999)
			    modnmsav <- paste(modnmsav, nm[j], sep = " vs. ")
		    }
		    print(paste("**************Ranking this step's models (all have same # parameters): ",
			modnmsav, sep = ""))
		    modrank <- modrank[modrank[, 2] > -999, ]
		    modrank <- modrank[order(modrank[, 2], modrank[,
			1]), ]
		    if (modrank[1, 1] == 1) {
			model <- model1
			submod <- as.character(substitute(model1))
		    }
		    if (modrank[1, 1] == 2) {
			model <- model2
			submod <- as.character(substitute(model2))
		    }
		    if (modrank[1, 1] == 3) {
			model <- model3
			submod <- as.character(substitute(model3))
		    }
		    if (modrank[1, 1] == 4) {
			model <- model4
			submod <- as.character(substitute(model4))
		    }
		    submod <- substr(submod, 10, 11)
		    FPCEnv$tempparam.select <- submod
		    return(model)
		}
	    }
	    rncheck <- function(modname, existing = FALSE) {
		modname <- as.character(substitute(modname))
		FPCEnv$tempmodnm <- modname
		modobj <- try(get(as.character(substitute(modname)),envir = pcklib)
		, silent = TRUE)
		if (class(modobj)[1] == "try-error" | existing == FALSE |
		    class(modobj)[1] == "NULL") {
		    outp <- TRUE
		} else {
		    outp <- FALSE
		}
		return(outp)
	    }
	    rncheckfirst <- function(modname, existing = FALSE) {
		modname <- as.character(substitute(modname))
		modobj <- try(get(as.character(substitute(modname)),envir = pcklib)
		, silent = TRUE)
		if (class(modobj)[1] == "try-error" | existing == FALSE) {
		} else {
		    return(modobj)
		}
	    }
	    rnassign <- function() {
		modname <- parse(text = sprintf("%s", FPCEnv$tempmodnm))
		if (class(eval(modname)[1]) != "try-error") {
		    chk <- try(unlist(summary(eval(modname)))["RSE"],
			silent = TRUE)
		    if (class(chk)[1] == "try-error") {
			return(1)
		    } else {
			modnm <- sprintf("%s", FPCEnv$tempmodnm)
			FPCEnv$tempparam.select <- substr(modnm, 10,
			  11)
			assign(modnm,eval(modname), pcklib)
			return(eval(modname))
		    }
		} else {
		    return(1)
		}
	    }
	    tstmod <- function(modelsub, modelcurrent) {
		extraF <- try(get("extraF", pos = 1), silent = TRUE)
		if (is.na(extraF(modelsub, modelcurrent, warn = F)[4]) == FALSE) {
		    if ((extraF(modelsub, modelcurrent, warn = F)[4]) > 0.05 &
			 sign(extraF(modelsub, modelcurrent, warn = F)[1]) != -1 &
			FPCEnv$legitmodel[1] == "legitmodelreset") {
			currentmodel <- modelsub
			return(currentmodel)
		    } else {
		     if(sign(extraF(modelsub, modelcurrent, warn = F)[1]) != -1 &
		       FPCEnv$legitmodel[1] ==  "legitmodelreset") {
		    	  currentmodel <- modelcurrent
			  return(currentmodel)
		    }else{
			if (FPCEnv$legitmodel[1] == "legitmodelreset" ) {
			  currentmodel <- modelsub
			  return(currentmodel)
			} else {
			  currentmodel <- FPCEnv$legitmodel
			  return(currentmodel)
			}
		     }
		    }
		} else {
		    if (FPCEnv$legitmodel[1] == "legitmodelreset") {
			currentmodel <- modelcurrent
			return(currentmodel)
		    } else {
			currentmodel <- FPCEnv$legitmodel
			return(currentmodel)
		    }
		}
	    }
	    cnt <- 1
	    step1submod <- FALSE
	    step5submod <- FALSE
	    FPCEnv$tempparam.select <- "NONE"
	    while (cnt < 6) {
		if (modelsig < 0.05) {
		    if (cnt == 1) {
			print("Step 1 of a maximum of 6*********************************************************************")
			print("--ASSESSING MODEL: richardsR1.lis  --")
			richardsR1.lis <- rncheckfirst(richardsR1.lis,
			  existing = existing)
			if (rncheck(richardsR1.lis, existing = existing) ==
			  TRUE)
			  richardsR1.lis <- eval(parse(text=sprintf("%s",paste("try({
			    nlsList(y ~ SSposnegRichards(x, Asym = Asym,
			      K = K, Infl = Infl, M = M, RAsym = RAsym,
			      Rk = Rk, Ri = Ri, RM = RM, modno = 1, pn.options = ",pnoptnm,"),
			      data = userdata, ...)
			  }, silent = TRUE)",sep=""))))
			currentmodel <- rnassign()
			if (class(currentmodel)[1] != "numeric")
			  step1submod <- TRUE
		    }
		    if (cnt == 2) {
			print("Step 2 of a maximum of 6*********************************************************************")
			print("--ASSESSING MODEL: richardsR2.lis  --")
			richardsR2.lis <- rncheckfirst(richardsR2.lis,
			  existing = existing)
			if (rncheck(richardsR2.lis, existing = existing) ==
			  TRUE)
			  richardsR2.lis <- eval(parse(text=sprintf("%s",paste("try({
			    nlsList(y ~ SSposnegRichards(x, Asym = Asym,
			      K = K, Infl = Infl, M = M, RAsym = RAsym,
			      Rk = Rk, Ri = Ri, RM = 1, modno = 2, pn.options = ",pnoptnm,"), data = userdata, ...)
			  }, silent = TRUE)",sep=""))))
			dump <- rnassign()
			print("--ASSESSING MODEL: richardsR7.lis  --")
			richardsR7.lis <- rncheckfirst(richardsR7.lis,
			  existing = existing)
			if (rncheck(richardsR7.lis, existing = existing) ==
			  TRUE)
			  richardsR7.lis <- eval(parse(text=sprintf("%s",paste("try({
			    nlsList(y ~ SSposnegRichards(x, Asym = Asym,
			      K = K, Infl = Infl, M = M, RAsym = 1, Rk = Rk,
			      Ri = Ri, RM = RM, modno = 7, pn.options = ",pnoptnm,"), data = userdata, ...)
			  }, silent = TRUE)",sep=""))))
			dump <- rnassign()
			print("--ASSESSING MODEL: richardsR6.lis  --")
			richardsR6.lis <- rncheckfirst(richardsR6.lis,
			  existing = existing)
			if (rncheck(richardsR6.lis, existing = existing) ==
			  TRUE)
			  richardsR6.lis <- eval(parse(text=sprintf("%s",paste("try({
			    nlsList(y ~ SSposnegRichards(x, Asym = Asym,
			      K = K, Infl = Infl, M = M, RAsym = RAsym,
			      Rk = 1, Ri = Ri, RM = RM, modno = 6, pn.options = ",pnoptnm,"), data = userdata, ...)
			  }, silent = TRUE)",sep=""))))
			dump <- rnassign()
			print("--ASSESSING MODEL: richardsR8.lis  --")
			richardsR8.lis <- rncheckfirst(richardsR8.lis,
			  existing = existing)
			if (rncheck(richardsR8.lis, existing = existing) ==
			  TRUE)
			  richardsR8.lis <- eval(parse(text=sprintf("%s",paste("try({
			    nlsList(y ~ SSposnegRichards(x, Asym = Asym,
			      K = K, Infl = Infl, M = M, RAsym = RAsym,
			      Rk = Rk, Ri = 1, RM = RM, modno = 8, pn.options = ",pnoptnm,"), data = userdata, ...)
			  }, silent = TRUE)",sep=""))))
			dump <- rnassign()
			submodel <- rankmod(richardsR2.lis, richardsR7.lis,
			  richardsR6.lis, richardsR8.lis)
			step2stat <- extraF(submodel, currentmodel, warn = F)
			currentmodel <- tstmod(submodel, currentmodel)
		    }
		    if (cnt == 3) {
			print("Step 3 of a maximum of 6*********************************************************************")
			currentmodID3 <- FPCEnv$tempparam.select
			if (currentmodID3 == "NONE")
			  currentmodID3 = "2."
			if (currentmodID3 == "2.") {
			  print("--ASSESSING MODEL: richardsR14.lis  --")
			  richardsR14.lis <- rncheckfirst(richardsR14.lis,
			    existing = existing)
			  if (rncheck(richardsR14.lis, existing = existing) ==
			    TRUE)
			    richardsR14.lis <- eval(parse(text=sprintf("%s",paste("try({
			      nlsList(y ~ SSposnegRichards(x, Asym = Asym,
				K = K, Infl = Infl, M = M, RAsym = 1,
				Rk = Rk, Ri = Ri, RM = 1, modno = 14, pn.options = ",pnoptnm,"),
				data = userdata, ...)
			    }, silent = TRUE)",sep=""))))
			  dump <- rnassign()
			  print("--ASSESSING MODEL: richardsR13.lis  --")
			  richardsR13.lis <- rncheckfirst(richardsR13.lis,
			    existing = existing)
			  if (rncheck(richardsR13.lis, existing = existing) ==
			    TRUE)
			    richardsR13.lis <- eval(parse(text=sprintf("%s",paste("try({
			      nlsList(y ~ SSposnegRichards(x, Asym = Asym,
				K = K, Infl = Infl, M = M, RAsym = RAsym,
				Rk = 1, Ri = Ri, RM = 1, modno = 13, pn.options = ",pnoptnm,"),
				data = userdata, ...)
			    }, silent = TRUE)",sep=""))))
			  dump <- rnassign()
			  print("--ASSESSING MODEL: richardsR15.lis  --")
			  richardsR15.lis <- rncheckfirst(richardsR15.lis,
			    existing = existing)
			  if (rncheck(richardsR15.lis, existing = existing) ==
			    TRUE)
			    richardsR15.lis <- eval(parse(text=sprintf("%s",paste("try({
			      nlsList(y ~ SSposnegRichards(x, Asym = Asym,
				K = K, Infl = Infl, M = M, RAsym = RAsym,
				Rk = Rk, Ri = 1, RM = 1, modno = 15, pn.options = ",pnoptnm,"),
				data = userdata, ...)
			    }, silent = TRUE)",sep=""))))
			  dump <- rnassign()
			  submodel <- rankmod(richardsR14.lis, richardsR13.lis,
			    richardsR15.lis)
			  step3stat <- extraF(submodel, currentmodel, warn = F)
			  currentmodel <- tstmod(submodel, currentmodel)
			}
			if (currentmodID3 == "NONE")
			  currentmodID3 = "7."
			if (currentmodID3 == "7.") {
			  print("--ASSESSING MODEL: richardsR14.lis  --")
			  richardsR14.lis <- rncheckfirst(richardsR14.lis,
			    existing = existing)
			  if (rncheck(richardsR14.lis, existing = existing) ==
			    TRUE)
			    richardsR14.lis <- eval(parse(text=sprintf("%s",paste("try({
			      nlsList(y ~ SSposnegRichards(x, Asym = Asym,
				K = K, Infl = Infl, M = M, RAsym = 1,
				Rk = Rk, Ri = Ri, RM = 1, modno = 14, pn.options = ",pnoptnm,"),
				data = userdata, ...)
			    }, silent = TRUE)",sep=""))))
			  dump <- rnassign()
			  print("--ASSESSING MODEL: richardsR3.lis  --")
			  richardsR3.lis <- rncheckfirst(richardsR3.lis,
			    existing = existing)
			  if (rncheck(richardsR3.lis, existing = existing) ==
			    TRUE)
			    richardsR3.lis <- eval(parse(text=sprintf("%s",paste("try({
			      nlsList(y ~ SSposnegRichards(x, Asym = Asym,
				K = K, Infl = Infl, M = M, RAsym = 1,
				Rk = 1, Ri = Ri, RM = RM, modno = 3, pn.options = ",pnoptnm,"),
				data = userdata, ...)
			    }, silent = TRUE)",sep=""))))
			  dump <- rnassign()
			  print("--ASSESSING MODEL: richardsR9.lis  --")
			  richardsR9.lis <- rncheckfirst(richardsR9.lis,
			    existing = existing)
			  if (rncheck(richardsR9.lis, existing = existing) ==
			    TRUE)
			    richardsR9.lis <- eval(parse(text=sprintf("%s",paste("try({
			      nlsList(y ~ SSposnegRichards(x, Asym = Asym,
				K = K, Infl = Infl, M = M, RAsym = 1,
				Rk = Rk, Ri = 1, RM = RM, modno = 9, pn.options = ",pnoptnm,"),
				data = userdata, ...)
			    }, silent = TRUE)",sep=""))))
			  dump <- rnassign()
			  submodel <- rankmod(richardsR14.lis, richardsR3.lis,
			    richardsR9.lis)
			  step3stat <- extraF(submodel, currentmodel, warn = F)
			  currentmodel <- tstmod(submodel, currentmodel)
			}
			if (currentmodID3 == "NONE")
			  currentmodID3 = "6."
			if (currentmodID3 == "6.") {
			  print("--ASSESSING MODEL: richardsR13.lis  --")
			  richardsR13.lis <- rncheckfirst(richardsR13.lis,
			    existing = existing)
			  if (rncheck(richardsR13.lis, existing = existing) ==
			    TRUE)
			    richardsR13.lis <- eval(parse(text=sprintf("%s",paste("try({
			      nlsList(y ~ SSposnegRichards(x, Asym = Asym,
				K = K, Infl = Infl, M = M, RAsym = RAsym,
				Rk = 1, Ri = Ri, RM = 1, modno = 13, pn.options = ",pnoptnm,"),
				data = userdata, ...)
			    }, silent = TRUE)",sep=""))))
			  dump <- rnassign()
			  print("--ASSESSING MODEL: richardsR3.lis  --")
			  richardsR3.lis <- rncheckfirst(richardsR3.lis,
			    existing = existing)
			  if (rncheck(richardsR3.lis, existing = existing) ==
			    TRUE)
			    richardsR3.lis <- eval(parse(text=sprintf("%s",paste("try({
			      nlsList(y ~ SSposnegRichards(x, Asym = Asym,
				K = K, Infl = Infl, M = M, RAsym = 1,
				Rk = 1, Ri = Ri, RM = RM, modno = 3, pn.options = ",pnoptnm,"),
				data = userdata, ...)
			    }, silent = TRUE)",sep=""))))
			  dump <- rnassign()
			  print("--ASSESSING MODEL: richardsR4.lis  --")
			  richardsR4.lis <- rncheckfirst(richardsR4.lis,
			    existing = existing)
			  if (rncheck(richardsR4.lis, existing = existing) ==
			    TRUE)
			    richardsR4.lis <- eval(parse(text=sprintf("%s",paste("try({
			      nlsList(y ~ SSposnegRichards(x, Asym = Asym,
				K = K, Infl = Infl, M = M, RAsym = RAsym,
				Rk = 1, Ri = 1, RM = RM, modno = 4, pn.options = ",pnoptnm,"),
				data = userdata, ...)
			    }, silent = TRUE)",sep=""))))
			  dump <- rnassign()
			  submodel <- rankmod(richardsR13.lis, richardsR3.lis,
			    richardsR4.lis)
			  step3stat <- extraF(submodel, currentmodel, warn = F)
			  currentmodel <- tstmod(submodel, currentmodel)
			}
			if (currentmodID3 == "NONE")
			  currentmodID3 = "8."
			if (currentmodID3 == "8.") {
			  print("--ASSESSING MODEL: richardsR15.lis  --")
			  richardsR15.lis <- rncheckfirst(richardsR15.lis,
			    existing = existing)
			  if (rncheck(richardsR15.lis, existing = existing) ==
			    TRUE)
			    richardsR15.lis <- eval(parse(text=sprintf("%s",paste("try({
			      nlsList(y ~ SSposnegRichards(x, Asym = Asym,
				K = K, Infl = Infl, M = M, RAsym = RAsym,
				Rk = Rk, Ri = 1, RM = 1, modno = 15, pn.options = ",pnoptnm,"),
				data = userdata, ...)
			    }, silent = TRUE)",sep=""))))
			  dump <- rnassign()
			  print("--ASSESSING MODEL: richardsR9.lis  --")
			  richardsR9.lis <- rncheckfirst(richardsR9.lis,
			    existing = existing)
			  if (rncheck(richardsR9.lis, existing = existing) ==
			    TRUE)
			    richardsR9.lis <- eval(parse(text=sprintf("%s",paste("try({
			      nlsList(y ~ SSposnegRichards(x, Asym = Asym,
				K = K, Infl = Infl, M = M, RAsym = 1,
				Rk = Rk, Ri = 1, RM = RM, modno = 9, pn.options = ",pnoptnm,"),
				data = userdata, ...)
			    }, silent = TRUE)",sep=""))))
			  dump <- rnassign()
			  print("--ASSESSING MODEL: richardsR4.lis  --")
			  richardsR4.lis <- rncheckfirst(richardsR4.lis,
			    existing = existing)
			  if (rncheck(richardsR4.lis, existing = existing) ==
			    TRUE)
			    richardsR4.lis <- eval(parse(text=sprintf("%s",paste("try({
			      nlsList(y ~ SSposnegRichards(x, Asym = Asym,
				K = K, Infl = Infl, M = M, RAsym = RAsym,
				Rk = 1, Ri = 1, RM = RM, modno = 4, pn.options = ",pnoptnm,"),
				data = userdata, ...)
			    }, silent = TRUE)",sep=""))))
			  dump <- rnassign()
			  submodel <- rankmod(richardsR15.lis, richardsR9.lis,
			    richardsR4.lis)
			  step3stat <- extraF(submodel, currentmodel, warn = F)
			  currentmodel <- tstmod(submodel, currentmodel)
			}
		    }
		    if (cnt == 4) {
			print("Step 4 of a maximum of 6*********************************************************************")
			currentmodID2 <- FPCEnv$tempparam.select
			if (currentmodID3 == "NONE" & currentmodID2 ==
			  "NONE") {
			  currentmodID3 = "2."
			  currentmodID2 = "14"
			}
			if (currentmodID3 == "2." & currentmodID2 ==
			  "14") {
			  print("--ASSESSING MODEL: richardsR10.lis  --")
			  richardsR10.lis <- rncheckfirst(richardsR10.lis,
			    existing = existing)
			  if (rncheck(richardsR10.lis, existing = existing) ==
			    TRUE)
			    richardsR10.lis <- eval(parse(text=sprintf("%s",paste("try({
			      nlsList(y ~ SSposnegRichards(x, Asym = Asym,
				K = K, Infl = Infl, M = M, RAsym = 1,
				Rk = 1, Ri = Ri, RM = 1, modno = 10, pn.options = ",pnoptnm,"),
				data = userdata, ...)
			    }, silent = TRUE)",sep=""))))
			  dump <- rnassign()
			  print("--ASSESSING MODEL: richardsR16.lis  --")
			  richardsR16.lis <- rncheckfirst(richardsR16.lis,
			    existing = existing)
			  if (rncheck(richardsR16.lis, existing = existing) ==
			    TRUE)
			    richardsR16.lis <- eval(parse(text=sprintf("%s",paste("try({
			      nlsList(y ~ SSposnegRichards(x, Asym = Asym,
				K = K, Infl = Infl, M = M, RAsym = 1,
				Rk = Rk, Ri = 1, RM = 1, modno = 16, pn.options = ",pnoptnm,"),
				data = userdata, ...)
			    }, silent = TRUE)",sep=""))))
			  dump <- rnassign()
			  submodel <- rankmod(richardsR10.lis, richardsR16.lis)
			  step4stat <- extraF(submodel, currentmodel, warn = F)
			  currentmodel <- tstmod(submodel, currentmodel)
			}
			if (currentmodID3 == "NONE" & currentmodID2 ==
			  "NONE") {
			  currentmodID3 = "2."
			  currentmodID2 = "13"
			}
			if (currentmodID3 == "2." & currentmodID2 ==
			  "13") {
			  print("--ASSESSING MODEL: richardsR10.lis  --")
			  richardsR10.lis <- rncheckfirst(richardsR10.lis,
			    existing = existing)
			  if (rncheck(richardsR10.lis, existing = existing) ==
			    TRUE)
			    richardsR10.lis <- eval(parse(text=sprintf("%s",paste("try({
			      nlsList(y ~ SSposnegRichards(x, Asym = Asym,
				K = K, Infl = Infl, M = M, RAsym = 1,
				Rk = 1, Ri = Ri, RM = 1, modno = 10, pn.options = ",pnoptnm,"),
				data = userdata, ...)
			    }, silent = TRUE)",sep=""))))
			  dump <- rnassign()
			  print("--ASSESSING MODEL: richardsR11.lis  --")
			  richardsR11.lis <- rncheckfirst(richardsR11.lis,
			    existing = existing)
			  if (rncheck(richardsR11.lis, existing = existing) ==
			    TRUE)
			    richardsR11.lis <- eval(parse(text=sprintf("%s",paste("try({
			      nlsList(y ~ SSposnegRichards(x, Asym = Asym,
				K = K, Infl = Infl, M = M, RAsym = RAsym,
				Rk = 1, Ri = 1, RM = 1, modno = 11, pn.options = ",pnoptnm,"),
				data = userdata, ...)
			    }, silent = TRUE)",sep=""))))
			  dump <- rnassign()
			  submodel <- rankmod(richardsR10.lis, richardsR11.lis)
			  step4stat <- extraF(submodel, currentmodel, warn = F)
			  currentmodel <- tstmod(submodel, currentmodel)
			}
			if (currentmodID3 == "NONE" & currentmodID2 ==
			  "NONE") {
			  currentmodID3 = "2."
			  currentmodID2 = "15"
			}
			if (currentmodID3 == "2." & currentmodID2 ==
			  "15") {
			  print("--ASSESSING MODEL: richardsR16.lis  --")
			  richardsR16.lis <- rncheckfirst(richardsR16.lis,
			    existing = existing)
			  if (rncheck(richardsR16.lis, existing = existing) ==
			    TRUE)
			    richardsR16.lis <- eval(parse(text=sprintf("%s",paste("try({
			      nlsList(y ~ SSposnegRichards(x, Asym = Asym,
				K = K, Infl = Infl, M = M, RAsym = 1,
				Rk = Rk, Ri = 1, RM = 1, modno = 16, pn.options = ",pnoptnm,"),
				data = userdata, ...)
			    }, silent = TRUE)",sep=""))))
			  dump <- rnassign()
			  print("--ASSESSING MODEL: richardsR11.lis  --")
			  richardsR11.lis <- rncheckfirst(richardsR11.lis,
			    existing = existing)
			  if (rncheck(richardsR11.lis, existing = existing) ==
			    TRUE)
			    richardsR11.lis <- eval(parse(text=sprintf("%s",paste("try({
			      nlsList(y ~ SSposnegRichards(x, Asym = Asym,
				K = K, Infl = Infl, M = M, RAsym = RAsym,
				Rk = 1, Ri = 1, RM = 1, modno = 11, pn.options = ",pnoptnm,"),
				data = userdata, ...)
			    }, silent = TRUE)",sep=""))))
			  dump <- rnassign()
			  submodel <- rankmod(richardsR16.lis, richardsR11.lis)
			  step4stat <- extraF(submodel, currentmodel, warn = F)
			  currentmodel <- tstmod(submodel, currentmodel)
			}
			if (currentmodID3 == "NONE" & currentmodID2 ==
			  "NONE") {
			  currentmodID3 = "7."
			  currentmodID2 = "14"
			}
			if (currentmodID3 == "7." & currentmodID2 ==
			  "14") {
			  print("--ASSESSING MODEL: richardsR10.lis  --")
			  richardsR10.lis <- rncheckfirst(richardsR10.lis,
			    existing = existing)
			  if (rncheck(richardsR10.lis, existing = existing) ==
			    TRUE)
			    richardsR10.lis <- eval(parse(text=sprintf("%s",paste("try({
			      nlsList(y ~ SSposnegRichards(x, Asym = Asym,
				K = K, Infl = Infl, M = M, RAsym = 1,
				Rk = 1, Ri = Ri, RM = 1, modno = 10, pn.options = ",pnoptnm,"),
				data = userdata, ...)
			    }, silent = TRUE)",sep=""))))
			  dump <- rnassign()
			  print("--ASSESSING MODEL: richardsR16.lis  --")
			  richardsR16.lis <- rncheckfirst(richardsR16.lis,
			    existing = existing)
			  if (rncheck(richardsR16.lis, existing = existing) ==
			    TRUE)
			    richardsR16.lis <- eval(parse(text=sprintf("%s",paste("try({
			      nlsList(y ~ SSposnegRichards(x, Asym = Asym,
				K = K, Infl = Infl, M = M, RAsym = 1,
				Rk = Rk, Ri = 1, RM = 1, modno = 16, pn.options = ",pnoptnm,"),
				data = userdata, ...)
			    }, silent = TRUE)",sep=""))))
			  dump <- rnassign()
			  submodel <- rankmod(richardsR10.lis, richardsR16.lis)
			  step4stat <- extraF(submodel, currentmodel, warn = F)
			  currentmodel <- tstmod(submodel, currentmodel)
			}
			if (currentmodID3 == "NONE" & currentmodID2 ==
			  "NONE") {
			  currentmodID3 = "7."
			  currentmodID2 = "3."
			}
			if (currentmodID3 == "7." & currentmodID2 ==
			  "3.") {
			  print("--ASSESSING MODEL: richardsR10.lis  --")
			  richardsR10.lis <- rncheckfirst(richardsR10.lis,
			    existing = existing)
			  if (rncheck(richardsR10.lis, existing = existing) ==
			    TRUE)
			    richardsR10.lis <- eval(parse(text=sprintf("%s",paste("try({
			      nlsList(y ~ SSposnegRichards(x, Asym = Asym,
				K = K, Infl = Infl, M = M, RAsym = 1,
				Rk = 1, Ri = Ri, RM = 1, modno = 10, pn.options = ",pnoptnm,"),
				data = userdata, ...)
			    }, silent = TRUE)",sep=""))))
			  dump <- rnassign()
			  print("--ASSESSING MODEL: richardsR5.lis  --")
			  richardsR5.lis <- rncheckfirst(richardsR5.lis,
			    existing = existing)
			  if (rncheck(richardsR5.lis, existing = existing) ==
			    TRUE)
			    richardsR5.lis <- eval(parse(text=sprintf("%s",paste("try({
			      nlsList(y ~ SSposnegRichards(x, Asym = Asym,
				K = K, Infl = Infl, M = M, RAsym = 1,
				Rk = 1, Ri = 1, RM = RM, modno = 5, pn.options = ",pnoptnm,"),
				data = userdata, ...)
			    }, silent = TRUE)",sep=""))))
			  dump <- rnassign()
			  submodel <- rankmod(richardsR10.lis, richardsR5.lis)
			  step4stat <- extraF(submodel, currentmodel, warn = F)
			  currentmodel <- tstmod(submodel, currentmodel)
			}
			if (currentmodID3 == "NONE" & currentmodID2 ==
			  "NONE") {
			  currentmodID3 = "7."
			  currentmodID2 = "9."
			}
			if (currentmodID3 == "7." & currentmodID2 ==
			  "9.") {
			  print("--ASSESSING MODEL: richardsR16.lis  --")
			  richardsR16.lis <- rncheckfirst(richardsR16.lis,
			    existing = existing)
			  if (rncheck(richardsR16.lis, existing = existing) ==
			    TRUE)
			    richardsR16.lis <- eval(parse(text=sprintf("%s",paste("try({
			      nlsList(y ~ SSposnegRichards(x, Asym = Asym,
				K = K, Infl = Infl, M = M, RAsym = 1,
				Rk = Rk, Ri = 1, RM = 1, modno = 16, pn.options = ",pnoptnm,"),
				data = userdata, ...)
			    }, silent = TRUE)",sep=""))))
			  dump <- rnassign()
			  print("--ASSESSING MODEL: richardsR5.lis  --")
			  richardsR5.lis <- rncheckfirst(richardsR5.lis,
			    existing = existing)
			  if (rncheck(richardsR5.lis, existing = existing) ==
			    TRUE)
			    richardsR5.lis <- eval(parse(text=sprintf("%s",paste("try({
			      nlsList(y ~ SSposnegRichards(x, Asym = Asym,
				K = K, Infl = Infl, M = M, RAsym = 1,
				Rk = 1, Ri = 1, RM = RM, modno = 5, pn.options = ",pnoptnm,"),
				data = userdata, ...)
			    }, silent = TRUE)",sep=""))))
			  dump <- rnassign()
			  submodel <- rankmod(richardsR16.lis, richardsR5.lis)
			  step4stat <- extraF(submodel, currentmodel, warn = F)
			  currentmodel <- tstmod(submodel, currentmodel)
			}
			if (currentmodID3 == "NONE" & currentmodID2 ==
			  "NONE") {
			  currentmodID3 = "6."
			  currentmodID2 = "13"
			}
			if (currentmodID3 == "6." & currentmodID2 ==
			  "13") {
			  print("--ASSESSING MODEL: richardsR10.lis  --")
			  richardsR10.lis <- rncheckfirst(richardsR10.lis,
			    existing = existing)
			  if (rncheck(richardsR10.lis, existing = existing) ==
			    TRUE)
			    richardsR10.lis <- eval(parse(text=sprintf("%s",paste("try({
			      nlsList(y ~ SSposnegRichards(x, Asym = Asym,
				K = K, Infl = Infl, M = M, RAsym = 1,
				Rk = 1, Ri = Ri, RM = 1, modno = 10, pn.options = ",pnoptnm,"),
				data = userdata, ...)
			    }, silent = TRUE)",sep=""))))
			  dump <- rnassign()
			  print("--ASSESSING MODEL: richardsR11.lis  --")
			  richardsR11.lis <- rncheckfirst(richardsR11.lis,
			    existing = existing)
			  if (rncheck(richardsR11.lis, existing = existing) ==
			    TRUE)
			    richardsR11.lis <- eval(parse(text=sprintf("%s",paste("try({
			      nlsList(y ~ SSposnegRichards(x, Asym = Asym,
				K = K, Infl = Infl, M = M, RAsym = RAsym,
				Rk = 1, Ri = 1, RM = 1, modno = 11, pn.options = ",pnoptnm,"),
				data = userdata, ...)
			    }, silent = TRUE)",sep=""))))
			  dump <- rnassign()
			  submodel <- rankmod(richardsR10.lis, richardsR11.lis)
			  step4stat <- extraF(submodel, currentmodel, warn = F)
			  currentmodel <- tstmod(submodel, currentmodel)
			}
			if (currentmodID3 == "NONE" & currentmodID2 ==
			  "NONE") {
			  currentmodID3 = "6."
			  currentmodID2 = "3."
			}
			if (currentmodID3 == "6." & currentmodID2 ==
			  "3.") {
			  print("--ASSESSING MODEL: richardsR10.lis  --")
			  richardsR10.lis <- rncheckfirst(richardsR10.lis,
			    existing = existing)
			  if (rncheck(richardsR10.lis, existing = existing) ==
			    TRUE)
			    richardsR10.lis <- eval(parse(text=sprintf("%s",paste("try({
			      nlsList(y ~ SSposnegRichards(x, Asym = Asym,
				K = K, Infl = Infl, M = M, RAsym = 1,
				Rk = 1, Ri = Ri, RM = 1, modno = 10, pn.options = ",pnoptnm,"),
				data = userdata, ...)
			    }, silent = TRUE)",sep=""))))
			  dump <- rnassign()
			  print("--ASSESSING MODEL: richardsR5.lis  --")
			  richardsR5.lis <- rncheckfirst(richardsR5.lis,
			    existing = existing)
			  if (rncheck(richardsR5.lis, existing = existing) ==
			    TRUE)
			    richardsR5.lis <- eval(parse(text=sprintf("%s",paste("try({
			      nlsList(y ~ SSposnegRichards(x, Asym = Asym,
				K = K, Infl = Infl, M = M, RAsym = 1,
				Rk = 1, Ri = 1, RM = RM, modno = 5, pn.options = ",pnoptnm,"),
				data = userdata, ...)
			    }, silent = TRUE)",sep=""))))
			  dump <- rnassign()
			  submodel <- rankmod(richardsR10.lis, richardsR5.lis)
			  step4stat <- extraF(submodel, currentmodel, warn = F)
			  currentmodel <- tstmod(submodel, currentmodel)
			}
			if (currentmodID3 == "NONE" & currentmodID2 ==
			  "NONE") {
			  currentmodID3 = "6."
			  currentmodID2 = "4."
			}
			if (currentmodID3 == "6." & currentmodID2 ==
			  "4.") {
			  print("--ASSESSING MODEL: richardsR11.lis  --")
			  richardsR11.lis <- rncheckfirst(richardsR11.lis,
			    existing = existing)
			  if (rncheck(richardsR11.lis, existing = existing) ==
			    TRUE)
			    richardsR11.lis <- eval(parse(text=sprintf("%s",paste("try({
			      nlsList(y ~ SSposnegRichards(x, Asym = Asym,
				K = K, Infl = Infl, M = M, RAsym = RAsym,
				Rk = 1, Ri = 1, RM = 1, modno = 11, pn.options = ",pnoptnm,"),
				data = userdata, ...)
			    }, silent = TRUE)",sep=""))))
			  dump <- rnassign()
			  print("--ASSESSING MODEL: richardsR5.lis  --")
			  richardsR5.lis <- rncheckfirst(richardsR5.lis,
			    existing = existing)
			  if (rncheck(richardsR5.lis, existing = existing) ==
			    TRUE)
			    richardsR5.lis <- eval(parse(text=sprintf("%s",paste("try({
			      nlsList(y ~ SSposnegRichards(x, Asym = Asym,
				K = K, Infl = Infl, M = M, RAsym = 1,
				Rk = 1, Ri = 1, RM = RM, modno = 5, pn.options = ",pnoptnm,"),
				data = userdata, ...)
			    }, silent = TRUE)",sep=""))))
			  dump <- rnassign()
			  submodel <- rankmod(richardsR11.lis, richardsR5.lis)
			  step4stat <- extraF(submodel, currentmodel, warn = F)
			  currentmodel <- tstmod(submodel, currentmodel)
			}
			if (currentmodID3 == "NONE" & currentmodID2 ==
			  "NONE") {
			  currentmodID3 = "8."
			  currentmodID2 = "15"
			}
			if (currentmodID3 == "8." & currentmodID2 ==
			  "15") {
			  print("--ASSESSING MODEL: richardsR16.lis  --")
			  richardsR16.lis <- rncheckfirst(richardsR16.lis,
			    existing = existing)
			  if (rncheck(richardsR16.lis, existing = existing) ==
			    TRUE)
			    richardsR16.lis <- eval(parse(text=sprintf("%s",paste("try({
			      nlsList(y ~ SSposnegRichards(x, Asym = Asym,
				K = K, Infl = Infl, M = M, RAsym = 1,
				Rk = Rk, Ri = 1, RM = 1, modno = 16, pn.options = ",pnoptnm,"),
				data = userdata, ...)
			    }, silent = TRUE)",sep=""))))
			  dump <- rnassign()
			  print("--ASSESSING MODEL: richardsR11.lis  --")
			  richardsR11.lis <- rncheckfirst(richardsR11.lis,
			    existing = existing)
			  if (rncheck(richardsR11.lis, existing = existing) ==
			    TRUE)
			    richardsR11.lis <- eval(parse(text=sprintf("%s",paste("try({
			      nlsList(y ~ SSposnegRichards(x, Asym = Asym,
				K = K, Infl = Infl, M = M, RAsym = RAsym,
				Rk = 1, Ri = 1, RM = 1, modno = 11, pn.options = ",pnoptnm,"),
				data = userdata, ...)
			    }, silent = TRUE)",sep=""))))
			  dump <- rnassign()
			  submodel <- rankmod(richardsR16.lis, richardsR11.lis)
			  step4stat <- extraF(submodel, currentmodel, warn = F)
			  currentmodel <- tstmod(submodel, currentmodel)
			}
			if (currentmodID3 == "NONE" & currentmodID2 ==
			  "NONE") {
			  currentmodID3 = "8."
			  currentmodID2 = "9."
			}
			if (currentmodID3 == "8." & currentmodID2 ==
			  "9.") {
			  print("--ASSESSING MODEL: richardsR16.lis  --")
			  richardsR16.lis <- rncheckfirst(richardsR16.lis,
			    existing = existing)
			  if (rncheck(richardsR16.lis, existing = existing) ==
			    TRUE)
			    richardsR16.lis <- eval(parse(text=sprintf("%s",paste("try({
			      nlsList(y ~ SSposnegRichards(x, Asym = Asym,
				K = K, Infl = Infl, M = M, RAsym = 1,
				Rk = Rk, Ri = 1, RM = 1, modno = 16, pn.options = ",pnoptnm,"),
				data = userdata, ...)
			    }, silent = TRUE)",sep=""))))
			  dump <- rnassign()
			  print("--ASSESSING MODEL: richardsR5.lis  --")
			  richardsR5.lis <- rncheckfirst(richardsR5.lis,
			    existing = existing)
			  if (rncheck(richardsR5.lis, existing = existing) ==
			    TRUE)
			    richardsR5.lis <- eval(parse(text=sprintf("%s",paste("try({
			      nlsList(y ~ SSposnegRichards(x, Asym = Asym,
				K = K, Infl = Infl, M = M, RAsym = 1,
				Rk = 1, Ri = 1, RM = RM, modno = 5, pn.options = ",pnoptnm,"),
				data = userdata, ...)
			    }, silent = TRUE)",sep=""))))
			  dump <- rnassign()
			  submodel <- rankmod(richardsR16.lis, richardsR5.lis)
			  step4stat <- extraF(submodel, currentmodel, warn = F)
			  currentmodel <- tstmod(submodel, currentmodel)
			}
			if (currentmodID3 == "NONE" & currentmodID2 ==
			  "NONE") {
			  currentmodID3 = "8."
			  currentmodID2 = "4."
			}
			if (currentmodID3 == "8." & currentmodID2 ==
			  "4.") {
			  print("--ASSESSING MODEL: richardsR11.lis  --")
			  richardsR11.lis <- rncheckfirst(richardsR11.lis,
			    existing = existing)
			  if (rncheck(richardsR11.lis, existing = existing) ==
			    TRUE)
			    richardsR11.lis <- eval(parse(text=sprintf("%s",paste("try({
			      nlsList(y ~ SSposnegRichards(x, Asym = Asym,
				K = K, Infl = Infl, M = M, RAsym = RAsym,
				Rk = 1, Ri = 1, RM = 1, modno = 11, pn.options = ",pnoptnm,"),
				data = userdata, ...)
			    }, silent = TRUE)",sep=""))))
			  dump <- rnassign()
			  print("--ASSESSING MODEL: richardsR5.lis  --")
			  richardsR5.lis <- rncheckfirst(richardsR5.lis,
			    existing = existing)
			  if (rncheck(richardsR5.lis, existing = existing) ==
			    TRUE)
			    richardsR5.lis <- eval(parse(text=sprintf("%s",paste("try({
			      nlsList(y ~ SSposnegRichards(x, Asym = Asym,
				K = K, Infl = Infl, M = M, RAsym = 1,
				Rk = 1, Ri = 1, RM = RM, modno = 5, pn.options = ",pnoptnm,"),
				data = userdata, ...)
			    }, silent = TRUE)",sep=""))))
			  dump <- rnassign()
			  submodel <- rankmod(richardsR11.lis, richardsR5.lis)
			  step4stat <- extraF(submodel, currentmodel, warn = F)
			  currentmodel <- tstmod(submodel, currentmodel)
			}
		    }
		    if (cnt == 5) {
			print("Step 5 of 6*********************************************************************")
			currentmodID1 <- FPCEnv$tempparam.select
			print("--ASSESSING MODEL: richardsR12.lis  --")
			richardsR12.lis <- rncheckfirst(richardsR12.lis,
			  existing = existing)
			if (rncheck(richardsR12.lis, existing = existing) ==
			  TRUE)
			  richardsR12.lis <- eval(parse(text=sprintf("%s",paste("try({
			    nlsList(y ~ SSposnegRichards(x, Asym = Asym,
			      K = K, Infl = Infl, M = M, RAsym = 1, Rk = 1,
			      Ri = 1, RM = 1, modno = 12, pn.options = ",pnoptnm,"), data = userdata, ...)
			  }, silent = TRUE)",sep=""))))
			dump <- rnassign()
			step5stat <- extraF(richardsR12.lis, currentmodel, warn = F)
			step5submod <- TRUE
			currentmodel <- tstmod(richardsR12.lis, currentmodel)
		    }
		    cnt <- cnt + 1
		    print("4 param")
		    print(cnt)
		} else {
		    if (cnt == 1) {
			print("Step 1 of a maximum of 6*********************************************************************")
			print("--ASSESSING MODEL: richardsR21.lis  --")
			richardsR21.lis <- rncheckfirst(richardsR21.lis,
			  existing = existing)
			if (rncheck(richardsR21.lis, existing = existing) ==
			  TRUE)
			  richardsR21.lis <- eval(parse(text=sprintf("%s",paste("try({
			    nlsList(y ~ SSposnegRichards(x, Asym = Asym,
			      K = K, Infl = Infl, M = 1, RAsym = RAsym,
			      Rk = Rk, Ri = Ri, RM = RM, modno = 21, pn.options = ",pnoptnm,"),
			      data = userdata, ...)
			  }, silent = TRUE)",sep=""))))
			currentmodel <- rnassign()
			if (class(currentmodel)[1] != "numeric")
			  step1submod <- TRUE
		    }
		    if (cnt == 2) {
			print("Step 2 of a maximum of 6*********************************************************************")
			print("--ASSESSING MODEL: richardsR22.lis  --")
			richardsR22.lis <- rncheckfirst(richardsR22.lis,
			  existing = existing)
			if (rncheck(richardsR22.lis, existing = existing) ==
			  TRUE)
			  richardsR22.lis <- eval(parse(text=sprintf("%s",paste("try({
			    nlsList(y ~ SSposnegRichards(x, Asym = Asym,
			      K = K, Infl = Infl, M = 1, RAsym = RAsym,
			      Rk = Rk, Ri = Ri, RM = 1, modno = 22, pn.options = ",pnoptnm,"),
			      data = userdata, ...)
			  }, silent = TRUE)",sep=""))))
			dump <- rnassign()
			print("--ASSESSING MODEL: richardsR27.lis  --")
			richardsR27.lis <- rncheckfirst(richardsR27.lis,
			  existing = existing)
			if (rncheck(richardsR27.lis, existing = existing) ==
			  TRUE)
			  richardsR27.lis <- eval(parse(text=sprintf("%s",paste("try({
			    nlsList(y ~ SSposnegRichards(x, Asym = Asym,
			      K = K, Infl = Infl, M = 1, RAsym = 1, Rk = Rk,
			      Ri = Ri, RM = RM, modno = 27, pn.options = ",pnoptnm,"), data = userdata, ...)
			  }, silent = TRUE)",sep=""))))
			dump <- rnassign()
			print("--ASSESSING MODEL: richardsR26.lis  --")
			richardsR26.lis <- rncheckfirst(richardsR26.lis,
			  existing = existing)
			if (rncheck(richardsR26.lis, existing = existing) ==
			  TRUE)
			  richardsR26.lis <- eval(parse(text=sprintf("%s",paste("try({
			    nlsList(y ~ SSposnegRichards(x, Asym = Asym,
			      K = K, Infl = Infl, M = 1, RAsym = RAsym,
			      Rk = 1, Ri = Ri, RM = RM, modno = 26, pn.options = ",pnoptnm,"),
			      data = userdata, ...)
			  }, silent = TRUE)",sep=""))))
			dump <- rnassign()
			print("--ASSESSING MODEL: richardsR28.lis  --")
			richardsR28.lis <- rncheckfirst(richardsR28.lis,
			  existing = existing)
			if (rncheck(richardsR28.lis, existing = existing) ==
			  TRUE)
			  richardsR28.lis <- eval(parse(text=sprintf("%s",paste("try({
			    nlsList(y ~ SSposnegRichards(x, Asym = Asym,
			      K = K, Infl = Infl, M = 1, RAsym = RAsym,
			      Rk = Rk, Ri = 1, RM = RM, modno = 28, pn.options = ",pnoptnm,"),
			      data = userdata, ...)
			  }, silent = TRUE)",sep=""))))
			dump <- rnassign()
			submodel <- rankmod(richardsR22.lis, richardsR27.lis,
			  richardsR26.lis, richardsR28.lis)
			step2stat <- extraF(submodel, currentmodel, warn = F)
			currentmodel <- tstmod(submodel, currentmodel)
		    }
		    if (cnt == 3) {
			print("Step 3 of a maximum of 6*********************************************************************")
			currentmodID3 <- FPCEnv$tempparam.select
			if (currentmodID3 == "NONE")
			  currentmodID3 = "22"
			if (currentmodID3 == "22") {
			  print("--ASSESSING MODEL: richardsR34.lis  --")
			  richardsR34.lis <- rncheckfirst(richardsR34.lis,
			    existing = existing)
			  if (rncheck(richardsR34.lis, existing = existing) ==
			    TRUE)
			    richardsR34.lis <- eval(parse(text=sprintf("%s",paste("try({
			      nlsList(y ~ SSposnegRichards(x, Asym = Asym,
				K = K, Infl = Infl, M = 1, RAsym = 1,
				Rk = Rk, Ri = Ri, RM = 1, modno = 34, pn.options = ",pnoptnm,"),
				data = userdata, ...)
			    }, silent = TRUE)",sep=""))))
			  dump <- rnassign()
			  print("--ASSESSING MODEL: richardsR33.lis  --")
			  richardsR33.lis <- rncheckfirst(richardsR33.lis,
			    existing = existing)
			  if (rncheck(richardsR33.lis, existing = existing) ==
			    TRUE)
			    richardsR33.lis <- eval(parse(text=sprintf("%s",paste("try({
			      nlsList(y ~ SSposnegRichards(x, Asym = Asym,
				K = K, Infl = Infl, M = 1, RAsym = RAsym,
				Rk = 1, Ri = Ri, RM = 1, modno = 33, pn.options = ",pnoptnm,"),
				data = userdata, ...)
			    }, silent = TRUE)",sep=""))))
			  dump <- rnassign()
			  print("--ASSESSING MODEL: richardsR35.lis  --")
			  richardsR35.lis <- rncheckfirst(richardsR35.lis,
			    existing = existing)
			  if (rncheck(richardsR35.lis, existing = existing) ==
			    TRUE)
			    richardsR35.lis <- eval(parse(text=sprintf("%s",paste("try({
			      nlsList(y ~ SSposnegRichards(x, Asym = Asym,
				K = K, Infl = Infl, M = 1, RAsym = RAsym,
				Rk = Rk, Ri = 1, RM = 1, modno = 35, pn.options = ",pnoptnm,"),
				data = userdata, ...)
			    }, silent = TRUE)",sep=""))))
			  dump <- rnassign()
			  submodel <- rankmod(richardsR34.lis, richardsR33.lis,
			    richardsR35.lis)
			  step3stat <- extraF(submodel, currentmodel, warn = F)
			  currentmodel <- tstmod(submodel, currentmodel)
			}
			if (currentmodID3 == "NONE")
			  currentmodID3 = "27"
			if (currentmodID3 == "27") {
			  print("--ASSESSING MODEL: richardsR34.lis  --")
			  richardsR34.lis <- rncheckfirst(richardsR34.lis,
			    existing = existing)
			  if (rncheck(richardsR34.lis, existing = existing) ==
			    TRUE)
			    richardsR34.lis <- eval(parse(text=sprintf("%s",paste("try({
			      nlsList(y ~ SSposnegRichards(x, Asym = Asym,
				K = K, Infl = Infl, M = 1, RAsym = 1,
				Rk = Rk, Ri = Ri, RM = 1, modno = 34, pn.options = ",pnoptnm,"),
				data = userdata, ...)
			    }, silent = TRUE)",sep=""))))
			  dump <- rnassign()
			  print("--ASSESSING MODEL: richardsR23.lis  --")
			  richardsR23.lis <- rncheckfirst(richardsR23.lis,
			    existing = existing)
			  if (rncheck(richardsR23.lis, existing = existing) ==
			    TRUE)
			    richardsR23.lis <- eval(parse(text=sprintf("%s",paste("try({
			      nlsList(y ~ SSposnegRichards(x, Asym = Asym,
				K = K, Infl = Infl, M = 1, RAsym = 1,
				Rk = 1, Ri = Ri, RM = RM, modno = 23, pn.options = ",pnoptnm,"),
				data = userdata, ...)
			    }, silent = TRUE)",sep=""))))
			  dump <- rnassign()
			  print("--ASSESSING MODEL: richardsR29.lis  --")
			  richardsR29.lis <- rncheckfirst(richardsR29.lis,
			    existing = existing)
			  if (rncheck(richardsR29.lis, existing = existing) ==
			    TRUE)
			    richardsR29.lis <- eval(parse(text=sprintf("%s",paste("try({
			      nlsList(y ~ SSposnegRichards(x, Asym = Asym,
				K = K, Infl = Infl, M = 1, RAsym = 1,
				Rk = Rk, Ri = 1, RM = RM, modno = 29, pn.options = ",pnoptnm,"),
				data = userdata, ...)
			    }, silent = TRUE)",sep=""))))
			  dump <- rnassign()
			  submodel <- rankmod(richardsR34.lis, richardsR23.lis,
			    richardsR29.lis)
			  step3stat <- extraF(submodel, currentmodel, warn = F)
			  currentmodel <- tstmod(submodel, currentmodel)
			}
			if (currentmodID3 == "NONE")
			  currentmodID3 = "26"
			if (currentmodID3 == "26") {
			  print("--ASSESSING MODEL: richardsR33.lis  --")
			  richardsR33.lis <- rncheckfirst(richardsR33.lis,
			    existing = existing)
			  if (rncheck(richardsR33.lis, existing = existing) ==
			    TRUE)
			    richardsR33.lis <- eval(parse(text=sprintf("%s",paste("try({
			      nlsList(y ~ SSposnegRichards(x, Asym = Asym,
				K = K, Infl = Infl, M = 1, RAsym = RAsym,
				Rk = 1, Ri = Ri, RM = 1, modno = 33, pn.options = ",pnoptnm,"),
				data = userdata, ...)
			    }, silent = TRUE)",sep=""))))
			  dump <- rnassign()
			  print("--ASSESSING MODEL: richardsR23.lis  --")
			  richardsR23.lis <- rncheckfirst(richardsR23.lis,
			    existing = existing)
			  if (rncheck(richardsR23.lis, existing = existing) ==
			    TRUE)
			    richardsR23.lis <- eval(parse(text=sprintf("%s",paste("try({
			      nlsList(y ~ SSposnegRichards(x, Asym = Asym,
				K = K, Infl = Infl, M = 1, RAsym = 1,
				Rk = 1, Ri = Ri, RM = RM, modno = 23, pn.options = ",pnoptnm,"),
				data = userdata, ...)
			    }, silent = TRUE)",sep=""))))
			  dump <- rnassign()
			  print("--ASSESSING MODEL: richardsR24.lis  --")
			  richardsR24.lis <- rncheckfirst(richardsR24.lis,
			    existing = existing)
			  if (rncheck(richardsR24.lis, existing = existing) ==
			    TRUE)
			    richardsR24.lis <- eval(parse(text=sprintf("%s",paste("try({
			      nlsList(y ~ SSposnegRichards(x, Asym = Asym,
				K = K, Infl = Infl, M = 1, RAsym = RAsym,
				Rk = 1, Ri = 1, RM = RM, modno = 24, pn.options = ",pnoptnm,"),
				data = userdata, ...)
			    }, silent = TRUE)",sep=""))))
			  dump <- rnassign()
			  submodel <- rankmod(richardsR33.lis, richardsR23.lis,
			    richardsR24.lis)
			  step3stat <- extraF(submodel, currentmodel, warn = F)
			  currentmodel <- tstmod(submodel, currentmodel)
			}
			if (currentmodID3 == "NONE")
			  currentmodID3 = "28"
			if (currentmodID3 == "28") {
			  print("--ASSESSING MODEL: richardsR35.lis  --")
			  richardsR35.lis <- rncheckfirst(richardsR35.lis,
			    existing = existing)
			  if (rncheck(richardsR35.lis, existing = existing) ==
			    TRUE)
			    richardsR35.lis <- eval(parse(text=sprintf("%s",paste("try({
			      nlsList(y ~ SSposnegRichards(x, Asym = Asym,
				K = K, Infl = Infl, M = 1, RAsym = RAsym,
				Rk = Rk, Ri = 1, RM = 1, modno = 35, pn.options = ",pnoptnm,"),
				data = userdata, ...)
			    }, silent = TRUE)",sep=""))))
			  dump <- rnassign()
			  print("--ASSESSING MODEL: richardsR29.lis  --")
			  richardsR29.lis <- rncheckfirst(richardsR29.lis,
			    existing = existing)
			  if (rncheck(richardsR29.lis, existing = existing) ==
			    TRUE)
			    richardsR29.lis <- eval(parse(text=sprintf("%s",paste("try({
			      nlsList(y ~ SSposnegRichards(x, Asym = Asym,
				K = K, Infl = Infl, M = 1, RAsym = 1,
				Rk = Rk, Ri = 1, RM = RM, modno = 29, pn.options = ",pnoptnm,"),
				data = userdata, ...)
			    }, silent = TRUE)",sep=""))))
			  dump <- rnassign()
			  print("--ASSESSING MODEL: richardsR24.lis  --")
			  richardsR24.lis <- rncheckfirst(richardsR24.lis,
			    existing = existing)
			  if (rncheck(richardsR24.lis, existing = existing) ==
			    TRUE)
			    richardsR24.lis <- eval(parse(text=sprintf("%s",paste("try({
			      nlsList(y ~ SSposnegRichards(x, Asym = Asym,
				K = K, Infl = Infl, M = 1, RAsym = RAsym,
				Rk = 1, Ri = 1, RM = RM, modno = 24, pn.options = ",pnoptnm,"),
				data = userdata, ...)
			    }, silent = TRUE)",sep=""))))
			  dump <- rnassign()
			  submodel <- rankmod(richardsR35.lis, richardsR29.lis,
			    richardsR24.lis)
			  step3stat <- extraF(submodel, currentmodel, warn = F)
			  currentmodel <- tstmod(submodel, currentmodel)
			}
		    }
		    if (cnt == 4) {
			print("Step 4 of a maximum of 6*********************************************************************")
			currentmodID2 <- FPCEnv$tempparam.select
			if (currentmodID3 == "NONE" & currentmodID2 ==
			  "NONE") {
			  currentmodID3 = "22"
			  currentmodID2 = "34"
			}
			if (currentmodID3 == "22" & currentmodID2 ==
			  "34") {
			  print("--ASSESSING MODEL: richardsR30.lis  --")
			  richardsR30.lis <- rncheckfirst(richardsR30.lis,
			    existing = existing)
			  if (rncheck(richardsR30.lis, existing = existing) ==
			    TRUE)
			    richardsR30.lis <- eval(parse(text=sprintf("%s",paste("try({
			      nlsList(y ~ SSposnegRichards(x, Asym = Asym,
				K = K, Infl = Infl, M = 1, RAsym = 1,
				Rk = 1, Ri = Ri, RM = 1, modno = 30, pn.options = ",pnoptnm,"),
				data = userdata, ...)
			    }, silent = TRUE)",sep=""))))
			  dump <- rnassign()
			  print("--ASSESSING MODEL: richardsR36.lis  --")
			  richardsR36.lis <- rncheckfirst(richardsR36.lis,
			    existing = existing)
			  if (rncheck(richardsR36.lis, existing = existing) ==
			    TRUE)
			    richardsR36.lis <- eval(parse(text=sprintf("%s",paste("try({
			      nlsList(y ~ SSposnegRichards(x, Asym = Asym,
				K = K, Infl = Infl, M = 1, RAsym = 1,
				Rk = Rk, Ri = 1, RM = 1, modno = 36, pn.options = ",pnoptnm,"),
				data = userdata, ...)
			    }, silent = TRUE)",sep=""))))
			  dump <- rnassign()
			  submodel <- rankmod(richardsR30.lis, richardsR36.lis)
			  step4stat <- extraF(submodel, currentmodel, warn = F)
			  currentmodel <- tstmod(submodel, currentmodel)
			}
			if (currentmodID3 == "NONE" & currentmodID2 ==
			  "NONE") {
			  currentmodID3 = "22"
			  currentmodID2 = "33"
			}
			if (currentmodID3 == "22" & currentmodID2 ==
			  "33") {
			  print("--ASSESSING MODEL: richardsR33.lis  --")
			  richardsR30.lis <- rncheckfirst(richardsR30.lis,
			    existing = existing)
			  if (rncheck(richardsR30.lis, existing = existing) ==
			    TRUE)
			    richardsR30.lis <- eval(parse(text=sprintf("%s",paste("try({
			      nlsList(y ~ SSposnegRichards(x, Asym = Asym,
				K = K, Infl = Infl, M = 1, RAsym = 1,
				Rk = 1, Ri = Ri, RM = 1, modno = 30, pn.options = ",pnoptnm,"),
				data = userdata, ...)
			    }, silent = TRUE)",sep=""))))
			  dump <- rnassign()
			  print("--ASSESSING MODEL: richardsR31.lis  --")
			  richardsR31.lis <- rncheckfirst(richardsR31.lis,
			    existing = existing)
			  if (rncheck(richardsR31.lis, existing = existing) ==
			    TRUE)
			    richardsR31.lis <- eval(parse(text=sprintf("%s",paste("try({
			      nlsList(y ~ SSposnegRichards(x, Asym = Asym,
				K = K, Infl = Infl, M = 1, RAsym = RAsym,
				Rk = 1, Ri = 1, RM = 1, modno = 31, pn.options = ",pnoptnm,"),
				data = userdata, ...)
			    }, silent = TRUE)",sep=""))))
			  dump <- rnassign()
			  submodel <- rankmod(richardsR30.lis, richardsR31.lis)
			  step4stat <- extraF(submodel, currentmodel, warn = F)
			  currentmodel <- tstmod(submodel, currentmodel)
			}
			if (currentmodID3 == "NONE" & currentmodID2 ==
			  "NONE") {
			  currentmodID3 = "22"
			  currentmodID2 = "35"
			}
			if (currentmodID3 == "22" & currentmodID2 ==
			  "35") {
			  print("--ASSESSING MODEL: richardsR36.lis  --")
			  richardsR36.lis <- rncheckfirst(richardsR36.lis,
			    existing = existing)
			  if (rncheck(richardsR36.lis, existing = existing) ==
			    TRUE)
			    richardsR36.lis <- eval(parse(text=sprintf("%s",paste("try({
			      nlsList(y ~ SSposnegRichards(x, Asym = Asym,
				K = K, Infl = Infl, M = 1, RAsym = 1,
				Rk = Rk, Ri = 1, RM = 1, modno = 36, pn.options = ",pnoptnm,"),
				data = userdata, ...)
			    }, silent = TRUE)",sep=""))))
			  dump <- rnassign()
			  print("--ASSESSING MODEL: richardsR31.lis  --")
			  richardsR31.lis <- rncheckfirst(richardsR31.lis,
			    existing = existing)
			  if (rncheck(richardsR31.lis, existing = existing) ==
			    TRUE)
			    richardsR31.lis <- eval(parse(text=sprintf("%s",paste("try({
			      nlsList(y ~ SSposnegRichards(x, Asym = Asym,
				K = K, Infl = Infl, M = 1, RAsym = RAsym,
				Rk = 1, Ri = 1, RM = 1, modno = 31, pn.options = ",pnoptnm,"),
				data = userdata, ...)
			    }, silent = TRUE)",sep=""))))
			  dump <- rnassign()
			  submodel <- rankmod(richardsR36.lis, richardsR31.lis)
			  step4stat <- extraF(submodel, currentmodel, warn = F)
			  currentmodel <- tstmod(submodel, currentmodel)
			}
			if (currentmodID3 == "NONE" & currentmodID2 ==
			  "NONE") {
			  currentmodID3 = "27"
			  currentmodID2 = "34"
			}
			if (currentmodID3 == "27" & currentmodID2 ==
			  "34") {
			  print("--ASSESSING MODEL: richardsR30.lis  --")
			  richardsR30.lis <- rncheckfirst(richardsR30.lis,
			    existing = existing)
			  if (rncheck(richardsR30.lis, existing = existing) ==
			    TRUE)
			    richardsR30.lis <- eval(parse(text=sprintf("%s",paste("try({
			      nlsList(y ~ SSposnegRichards(x, Asym = Asym,
				K = K, Infl = Infl, M = 1, RAsym = 1,
				Rk = 1, Ri = Ri, RM = 1, modno = 30, pn.options = ",pnoptnm,"),
				data = userdata, ...)
			    }, silent = TRUE)",sep=""))))
			  dump <- rnassign()
			  print("--ASSESSING MODEL: richardsR36.lis  --")
			  richardsR36.lis <- rncheckfirst(richardsR36.lis,
			    existing = existing)
			  if (rncheck(richardsR36.lis, existing = existing) ==
			    TRUE)
			    richardsR36.lis <- eval(parse(text=sprintf("%s",paste("try({
			      nlsList(y ~ SSposnegRichards(x, Asym = Asym,
				K = K, Infl = Infl, M = 1, RAsym = 1,
				Rk = Rk, Ri = 1, RM = 1, modno = 36, pn.options = ",pnoptnm,"),
				data = userdata, ...)
			    }, silent = TRUE)",sep=""))))
			  dump <- rnassign()
			  submodel <- rankmod(richardsR30.lis, richardsR36.lis)
			  step4stat <- extraF(submodel, currentmodel, warn = F)
			  currentmodel <- tstmod(submodel, currentmodel)
			}
			if (currentmodID3 == "NONE" & currentmodID2 ==
			  "NONE") {
			  currentmodID3 = "27"
			  currentmodID2 = "23"
			}
			if (currentmodID3 == "27" & currentmodID2 ==
			  "23") {
			  print("--ASSESSING MODEL: richardsR30.lis  --")
			  richardsR30.lis <- rncheckfirst(richardsR30.lis,
			    existing = existing)
			  if (rncheck(richardsR30.lis, existing = existing) ==
			    TRUE)
			    richardsR30.lis <- eval(parse(text=sprintf("%s",paste("try({
			      nlsList(y ~ SSposnegRichards(x, Asym = Asym,
				K = K, Infl = Infl, M = 1, RAsym = 1,
				Rk = 1, Ri = Ri, RM = 1, modno = 30, pn.options = ",pnoptnm,"),
				data = userdata, ...)
			    }, silent = TRUE)",sep=""))))
			  dump <- rnassign()
			  print("--ASSESSING MODEL: richardsR25.lis  --")
			  richardsR25.lis <- rncheckfirst(richardsR25.lis,
			    existing = existing)
			  if (rncheck(richardsR25.lis, existing = existing) ==
			    TRUE)
			    richardsR25.lis <- eval(parse(text=sprintf("%s",paste("try({
			      nlsList(y ~ SSposnegRichards(x, Asym = Asym,
				K = K, Infl = Infl, M = 1, RAsym = 1,
				Rk = 1, Ri = 1, RM = RM, modno = 25, pn.options = ",pnoptnm,"),
				data = userdata, ...)
			    }, silent = TRUE)",sep=""))))
			  dump <- rnassign()
			  submodel <- rankmod(richardsR30.lis, richardsR25.lis)
			  step4stat <- extraF(submodel, currentmodel, warn = F)
			  currentmodel <- tstmod(submodel, currentmodel)
			}
			if (currentmodID3 == "NONE" & currentmodID2 ==
			  "NONE") {
			  currentmodID3 = "27"
			  currentmodID2 = "29"
			}
			if (currentmodID3 == "27" & currentmodID2 ==
			  "29") {
			  print("--ASSESSING MODEL: richardsR36.lis  --")
			  richardsR36.lis <- rncheckfirst(richardsR36.lis,
			    existing = existing)
			  if (rncheck(richardsR36.lis, existing = existing) ==
			    TRUE)
			    richardsR36.lis <- eval(parse(text=sprintf("%s",paste("try({
			      nlsList(y ~ SSposnegRichards(x, Asym = Asym,
				K = K, Infl = Infl, M = 1, RAsym = 1,
				Rk = Rk, Ri = 1, RM = 1, modno = 36, pn.options = ",pnoptnm,"),
				data = userdata, ...)
			    }, silent = TRUE)",sep=""))))
			  dump <- rnassign()
			  print("--ASSESSING MODEL: richardsR25.lis  --")
			  richardsR25.lis <- rncheckfirst(richardsR25.lis,
			    existing = existing)
			  if (rncheck(richardsR25.lis, existing = existing) ==
			    TRUE)
			    richardsR25.lis <- eval(parse(text=sprintf("%s",paste("try({
			      nlsList(y ~ SSposnegRichards(x, Asym = Asym,
				K = K, Infl = Infl, M = 1, RAsym = 1,
				Rk = 1, Ri = 1, RM = RM, modno = 25, pn.options = ",pnoptnm,"),
				data = userdata, ...)
			    }, silent = TRUE)",sep=""))))
			  dump <- rnassign()
			  submodel <- rankmod(richardsR36.lis, richardsR25.lis)
			  step4stat <- extraF(submodel, currentmodel, warn = F)
			  currentmodel <- tstmod(submodel, currentmodel)
			}
			if (currentmodID3 == "NONE" & currentmodID2 ==
			  "NONE") {
			  currentmodID3 = "26"
			  currentmodID2 = "33"
			}
			if (currentmodID3 == "26" & currentmodID2 ==
			  "33") {
			  print("--ASSESSING MODEL: richardsR33.lis  --")
			  richardsR30.lis <- rncheckfirst(richardsR30.lis,
			    existing = existing)
			  if (rncheck(richardsR30.lis, existing = existing) ==
			    TRUE)
			    richardsR30.lis <- eval(parse(text=sprintf("%s",paste("try({
			      nlsList(y ~ SSposnegRichards(x, Asym = Asym,
				K = K, Infl = Infl, M = 1, RAsym = 1,
				Rk = 1, Ri = Ri, RM = 1, modno = 30, pn.options = ",pnoptnm,"),
				data = userdata, ...)
			    }, silent = TRUE)",sep=""))))
			  dump <- rnassign()
			  print("--ASSESSING MODEL: richardsR31.lis  --")
			  richardsR31.lis <- rncheckfirst(richardsR31.lis,
			    existing = existing)
			  if (rncheck(richardsR31.lis, existing = existing) ==
			    TRUE)
			    richardsR31.lis <- eval(parse(text=sprintf("%s",paste("try({
			      nlsList(y ~ SSposnegRichards(x, Asym = Asym,
				K = K, Infl = Infl, M = 1, RAsym = RAsym,
				Rk = 1, Ri = 1, RM = 1, modno = 31, pn.options = ",pnoptnm,"),
				data = userdata, ...)
			    }, silent = TRUE)",sep=""))))
			  dump <- rnassign()
			  submodel <- rankmod(richardsR30.lis, richardsR31.lis)
			  step4stat <- extraF(submodel, currentmodel, warn = F)
			  currentmodel <- tstmod(submodel, currentmodel)
			}
			if (currentmodID3 == "NONE" & currentmodID2 ==
			  "NONE") {
			  currentmodID3 = "26"
			  currentmodID2 = "23"
			}
			if (currentmodID3 == "26" & currentmodID2 ==
			  "23") {
			  print("--ASSESSING MODEL: richardsR30.lis  --")
			  richardsR30.lis <- rncheckfirst(richardsR30.lis,
			    existing = existing)
			  if (rncheck(richardsR30.lis, existing = existing) ==
			    TRUE)
			    richardsR30.lis <- eval(parse(text=sprintf("%s",paste("try({
			      nlsList(y ~ SSposnegRichards(x, Asym = Asym,
				K = K, Infl = Infl, M = 1, RAsym = 1,
				Rk = 1, Ri = Ri, RM = 1, modno = 30, pn.options = ",pnoptnm,"),
				data = userdata, ...)
			    }, silent = TRUE)",sep=""))))
			  dump <- rnassign()
			  print("--ASSESSING MODEL: richardsR25.lis  --")
			  richardsR25.lis <- rncheckfirst(richardsR25.lis,
			    existing = existing)
			  if (rncheck(richardsR25.lis, existing = existing) ==
			    TRUE)
			    richardsR25.lis <- eval(parse(text=sprintf("%s",paste("try({
			      nlsList(y ~ SSposnegRichards(x, Asym = Asym,
				K = K, Infl = Infl, M = 1, RAsym = 1,
				Rk = 1, Ri = 1, RM = RM, modno = 25, pn.options = ",pnoptnm,"),
				data = userdata, ...)
			    }, silent = TRUE)",sep=""))))
			  dump <- rnassign()
			  submodel <- rankmod(richardsR30.lis, richardsR25.lis)
			  step4stat <- extraF(submodel, currentmodel, warn = F)
			  currentmodel <- tstmod(submodel, currentmodel)
			}
			if (currentmodID3 == "NONE" & currentmodID2 ==
			  "NONE") {
			  currentmodID3 = "26"
			  currentmodID2 = "24"
			}
			if (currentmodID3 == "26" & currentmodID2 ==
			  "24") {
			  print("--ASSESSING MODEL: richardsR31.lis  --")
			  richardsR31.lis <- rncheckfirst(richardsR31.lis,
			    existing = existing)
			  if (rncheck(richardsR31.lis, existing = existing) ==
			    TRUE)
			    richardsR31.lis <- eval(parse(text=sprintf("%s",paste("try({
			      nlsList(y ~ SSposnegRichards(x, Asym = Asym,
				K = K, Infl = Infl, M = 1, RAsym = RAsym,
				Rk = 1, Ri = 1, RM = 1, modno = 31, pn.options = ",pnoptnm,"),
				data = userdata, ...)
			    }, silent = TRUE)",sep=""))))
			  dump <- rnassign()
			  print("--ASSESSING MODEL: richardsR25.lis  --")
			  richardsR25.lis <- rncheckfirst(richardsR25.lis,
			    existing = existing)
			  if (rncheck(richardsR25.lis, existing = existing) ==
			    TRUE)
			    richardsR25.lis <- eval(parse(text=sprintf("%s",paste("try({
			      nlsList(y ~ SSposnegRichards(x, Asym = Asym,
				K = K, Infl = Infl, M = 1, RAsym = 1,
				Rk = 1, Ri = 1, RM = RM, modno = 25, pn.options = ",pnoptnm,"),
				data = userdata, ...)
			    }, silent = TRUE)",sep=""))))
			  dump <- rnassign()
			  submodel <- rankmod(richardsR31.lis, richardsR25.lis)
			  step4stat <- extraF(submodel, currentmodel, warn = F)
			  currentmodel <- tstmod(submodel, currentmodel)
			}
			if (currentmodID3 == "NONE" & currentmodID2 ==
			  "NONE") {
			  currentmodID3 = "28"
			  currentmodID2 = "35"
			}
			if (currentmodID3 == "28" & currentmodID2 ==
			  "35") {
			  print("--ASSESSING MODEL: richardsR36.lis  --")
			  richardsR36.lis <- rncheckfirst(richardsR36.lis,
			    existing = existing)
			  if (rncheck(richardsR36.lis, existing = existing) ==
			    TRUE)
			    richardsR36.lis <- eval(parse(text=sprintf("%s",paste("try({
			      nlsList(y ~ SSposnegRichards(x, Asym = Asym,
				K = K, Infl = Infl, M = 1, RAsym = 1,
				Rk = Rk, Ri = 1, RM = 1, modno = 36, pn.options = ",pnoptnm,"),
				data = userdata, ...)
			    }, silent = TRUE)",sep=""))))
			  dump <- rnassign()
			  print("--ASSESSING MODEL: richardsR31.lis  --")
			  richardsR31.lis <- rncheckfirst(richardsR31.lis,
			    existing = existing)
			  if (rncheck(richardsR31.lis, existing = existing) ==
			    TRUE)
			    richardsR31.lis <- eval(parse(text=sprintf("%s",paste("try({
			      nlsList(y ~ SSposnegRichards(x, Asym = Asym,
				K = K, Infl = Infl, M = 1, RAsym = RAsym,
				Rk = 1, Ri = 1, RM = 1, modno = 31, pn.options = ",pnoptnm,"),
				data = userdata, ...)
			    }, silent = TRUE)",sep=""))))
			  dump <- rnassign()
			  submodel <- rankmod(richardsR36.lis, richardsR31.lis)
			  step4stat <- extraF(submodel, currentmodel, warn = F)
			  currentmodel <- tstmod(submodel, currentmodel)
			}
			if (currentmodID3 == "NONE" & currentmodID2 ==
			  "NONE") {
			  currentmodID3 = "28"
			  currentmodID2 = "29"
			}
			if (currentmodID3 == "28" & currentmodID2 ==
			  "29") {
			  print("--ASSESSING MODEL: richardsR36.lis  --")
			  richardsR36.lis <- rncheckfirst(richardsR36.lis,
			    existing = existing)
			  if (rncheck(richardsR36.lis, existing = existing) ==
			    TRUE)
			    richardsR36.lis <- eval(parse(text=sprintf("%s",paste("try({
			      nlsList(y ~ SSposnegRichards(x, Asym = Asym,
				K = K, Infl = Infl, M = 1, RAsym = 1,
				Rk = Rk, Ri = 1, RM = 1, modno = 36, pn.options = ",pnoptnm,"),
				data = userdata, ...)
			    }, silent = TRUE)",sep=""))))
			  dump <- rnassign()
			  print("--ASSESSING MODEL: richardsR25.lis  --")
			  richardsR25.lis <- rncheckfirst(richardsR25.lis,
			    existing = existing)
			  if (rncheck(richardsR25.lis, existing = existing) ==
			    TRUE)
			    richardsR25.lis <- eval(parse(text=sprintf("%s",paste("try({
			      nlsList(y ~ SSposnegRichards(x, Asym = Asym,
				K = K, Infl = Infl, M = 1, RAsym = 1,
				Rk = 1, Ri = 1, RM = RM, modno = 25, pn.options = ",pnoptnm,"),
				data = userdata, ...)
			    }, silent = TRUE)",sep=""))))
			  dump <- rnassign()
			  submodel <- rankmod(richardsR36.lis, richardsR25.lis)
			  step4stat <- extraF(submodel, currentmodel, warn = F)
			  currentmodel <- tstmod(submodel, currentmodel)
			}
			if (currentmodID3 == "NONE" & currentmodID2 ==
			  "NONE") {
			  currentmodID3 = "28"
			  currentmodID2 = "24"
			}
			if (currentmodID3 == "28" & currentmodID2 ==
			  "24") {
			  print("--ASSESSING MODEL: richardsR31.lis  --")
			  richardsR31.lis <- rncheckfirst(richardsR31.lis,
			    existing = existing)
			  if (rncheck(richardsR31.lis, existing = existing) ==
			    TRUE)
			    richardsR31.lis <- eval(parse(text=sprintf("%s",paste("try({
			      nlsList(y ~ SSposnegRichards(x, Asym = Asym,
				K = K, Infl = Infl, M = 1, RAsym = RAsym,
				Rk = 1, Ri = 1, RM = 1, modno = 31, pn.options = ",pnoptnm,"),
				data = userdata, ...)
			    }, silent = TRUE)",sep=""))))
			  dump <- rnassign()
			  print("--ASSESSING MODEL: richardsR25.lis  --")
			  richardsR25.lis <- rncheckfirst(richardsR25.lis,
			    existing = existing)
			  if (rncheck(richardsR25.lis, existing = existing) ==
			    TRUE)
			    richardsR25.lis <- eval(parse(text=sprintf("%s",paste("try({
			      nlsList(y ~ SSposnegRichards(x, Asym = Asym,
				K = K, Infl = Infl, M = 1, RAsym = 1,
				Rk = 1, Ri = 1, RM = RM, modno = 25, pn.options = ",pnoptnm,"),
				data = userdata, ...)
			    }, silent = TRUE)",sep=""))))
			  dump <- rnassign()
			  submodel <- rankmod(richardsR31.lis, richardsR25.lis)
			  step4stat <- extraF(submodel, currentmodel, warn = F)
			  currentmodel <- tstmod(submodel, currentmodel)
			}
		    }
		    if (cnt == 5) {
			print("Step 5 of 6*********************************************************************")
			currentmodID1 <- FPCEnv$tempparam.select
			print("--ASSESSING MODEL: richardsR32.lis  --")
			richardsR32.lis <- rncheckfirst(richardsR32.lis,
			  existing = existing)
			if (rncheck(richardsR32.lis, existing = existing) ==
			  TRUE)
			  richardsR32.lis <- eval(parse(text=sprintf("%s",paste("try({
			    nlsList(y ~ SSposnegRichards(x, Asym = Asym,
			      K = K, Infl = Infl, M = 1, RAsym = 1, Rk = 1,
			      Ri = 1, RM = 1, modno = 32, pn.options = ",pnoptnm,"), data = userdata, ...)
			  }, silent = TRUE)",sep=""))))
			dump <- rnassign()
			step5stat <- extraF(richardsR32.lis, currentmodel)
			step5submod <- TRUE
			currentmodel <- tstmod(richardsR32.lis, currentmodel)
		    }
		    cnt <- cnt + 1
		}
	    }
	    if (modelsig < 0.05) {
		mod1 <- "1."
		mod4 <- "12"
	    } else {
		mod1 <- "21"
		mod4 <- "32"
	    }
	    if (step5submod == FALSE) {
		step5submod <- NA
		step5stat <- NA
		mod4 <- NA
		mod5 <- NA
		step6stat <- NA
	    } else {
		print("Step 6 of 6*********************************************************************")
		if (class (mostreducedmod)[1] != "numeric"){
		mod5 <- mostreducednm
		step6stat <- extraF(mostreducedmod, currentmodel, warn = F)
		currentmodel <- tstmod(mostreducedmod, currentmodel)
		} else {
		step6stat <- NA
		mod5 <- NA
		}
	    }
	    options(warn = -1)
	    currentmodID3 <- as.numeric(currentmodID3)
	    currentmodID2 <- as.numeric(currentmodID2)
	    currentmodID1 <- as.numeric(currentmodID1)
	    options(warn = 0)
	    if(mod1 == "1.") mod1 <- "1"
	    modnames <- c(paste("richardsR", mod1, ".lis", sep = ""),
		paste("richardsR", currentmodID3, ".lis", sep = ""),
		paste("richardsR", currentmodID2, ".lis", sep = ""),
		paste("richardsR", currentmodID1, ".lis", sep = ""),
		paste("richardsR", mod4, ".lis", sep = ""), paste("richardsR",
		    mod5, ".lis", sep = ""))
	    modnames[modnames == "richardsRNA.lis"] <- ""
	    modnames[modnames == "richardsRNONE.lis"] <- ""
	    stepwisetable <- data.frame(`       Best Submodel at Step` = modnames)
	    testof <- rep("", 6)
	    xf <- rep(NA, 6)
	    dfn <- rep(NA, 6)
	    dfd <- rep(NA, 6)
	    pval <- rep(NA, 6)
	    RSSgen <- rep(NA, 6)
	    RSSsub <- rep(NA, 6)
	    countfit <- 0
	    for (i in 1:6) {
		if (i > 1) {
		    if (modnames[i] == "" & modnames[i - 1] == "") {
			testof[i] <- "No models converged at this step"
		    } else {
			testof[i] <- paste("|      ", modnames[i], " vs ",
			  modnames[i - 1], "      |", sep = "")
		    }
		    options(warn = -1)
		    assessfits <- try( eval(parse(paste("step",i, "stat", sep = ""))), silent = TRUE)
		    if (class(assessfits)[1] == "try-error") countfit <- countfit + 1
		    if (i == 6 & countfit == 0) stop("No models were successfully fitted. Aborting..... Please check your data or change argument options.")
		    currstat <- eval(parse(text = sprintf("%s", (paste("step",
			i, "stat", sep = "")))))
		    options(warn = 0)
		    xf[i] <- round(as.numeric(currstat[1]), 4)
		    dfn[i] <- as.numeric(currstat[2])
		    dfd[i] <- as.numeric(currstat[3])
		    pval[i] <- round(as.numeric(currstat[4]), 8)
		    RSSgen[i] <- as.numeric(currstat[5])
		    RSSsub[i] <- as.numeric(currstat[6])
		} else {
		    testof[i] <- "|     Reduced model    More complex model      |"
		}
	    }
	    modnames1 <- modnames
	    for (i in 2:6) {
		if (is.na(pval[i]) == FALSE) {
		    if (pval[i] > 0.05 & sign(xf[i]) != -1) {
			modnames1[i] <- modnames1[i]
		    } else {
		    	if (sign(xf[i]) != -1) {
				modnames1[i] <- modnames1[i - 1]			
		    	} else {
				modnames1[i] <- modnames1[i]		    
		  	}
		    }
		} else {
		 	modnames1[i] <- modnames1[i]
		}
	    }
	    for (i in 2:6) {
		if (modnames[i] == "" & modnames[i - 1] == "") {
		    testof[i] <- "No models converged at this step"
		} else {
		    if (modnames[i - 1] != modnames1[i - 1])
			testof[i] <- paste("|      ", modnames[i], " vs ",
			  modnames1[i - 1], "      |", sep = "")
		}
	    }
	    testof <- unlist(testof)
	    xf <- data.frame(sprintf("%.4f", as.numeric(unlist(xf))))
	    dfn <- data.frame(as.numeric(unlist(dfn)))
	    dfd <- data.frame(as.numeric(unlist(dfd)))
	    pval <- data.frame(sprintf("%.8f", as.numeric(unlist(pval))))
	    RSSgen <- data.frame(signif( as.numeric(unlist(RSSgen)),4))
	    RSSsub <- data.frame(signif( as.numeric(unlist(RSSsub)),4))
	    stepwisetable <- data.frame(`       Best Submodel at Step` = modnames1,
		Test = testof, `F-stat` = xf, df_n = dfn, df_d = dfd,
		P = pval, RSS_sub = RSSsub, RSS_gen = RSSgen)
	    names(stepwisetable) <- c("       Best Submodel at Step",
		"Test                 ", "F-stat", "df_n", "df_d", "P",
		"RSS_sub", "RSS_gen")
	    row.names(stepwisetable) <- c("Step 1", "Step 2", "Step 3",
		"Step 4", "Step 5", "Step 6")
	    stepwisetable[is.na(stepwisetable)] <- ""
	    print("###########  Minimal applicable model arrived at by stepwise reduction is saved as pn.bestmodel.lis  ################")
	    assign("pn.bestmodel.lis", currentmodel, envir = Envir)
	    print("writing output to environment:")
	    print(Envir)
	    get.mod(to.envir = Envir, write.mod = TRUE)
	    assign("userdata",userdata, envir = Envir)
	    options(warn=-1)
	    try(rm("model1",envir = FPCEnv),silent=T)
	    try(rm("tempparam.select",envir = FPCEnv),silent=T)
	    try(rm("tempmodnm",envir = FPCEnv),silent=T)
	    try(rm("legitmodel",envir = FPCEnv),silent=T)
	    options( warn = 0)
	    return(stepwisetable)
	    ##value<< A \code{\link{data.frame}} containing statistics produced by \code{\link{extraF}}
	    ## evaluations at each step, detailing the name of the general and best reduced model at each
	    ## step. The overall best model evaluated by the end of the function is saved in specified environment as
	    ## \eqn{pn.bestmodel.lis}
	    ## The naming convention for models is a concatenation of 'richardsR', the modno and '.lis'
	    ## (see \code{\link{SSposnegRichards}}).
	    ##seealso<< \code{\link{pn.mod.compare}}
	    ## \code{\link{extraF}}
	    ## \code{\link{SSposnegRichards}}
	    ## \code{\link{nlsList}}
	    ##note<< If object \eqn{pnmodelparams} does not exist, \code{\link{modpar}}
	    ## will be called automatically prior to model selection. During selection, text is output
    ## to the screen to inform the user of the progress of model selection
    ## (which model is being fitted)
}
, ex = function(){
#run model selection for posneg.data object (only first 3 group levels for example's sake)
subdata<-subset(posneg.data, as.numeric(row.names (posneg.data) ) < 40)
modseltable <- pn.modselect.step(subdata$age, subdata$mass,
    subdata$id, existing = FALSE)

#fit nlsList model initially and then run model selection
#for posneg.data object when at least one model is already fit
#note forcemod is set to 3 so that models 21-36 are evaluated
#(only first 4 group levels for example's sake)
subdata<-subset(posneg.data, as.numeric(row.names (posneg.data) ) < 40)
richardsR22.lis <- nlsList(mass ~ SSposnegRichards(age, Asym = Asym, K = K,
   Infl = Infl, M = 1, RAsym = RAsym, Rk = Rk, Ri = Ri, RM = 1 , modno = 22)
                        ,data = subdata)
modseltable <- pn.modselect.step(subdata$age, subdata$mass,
    subdata$id, forcemod = 3, existing = TRUE)

#run model selection ranked by residual standard error*sample size
#(only first 4 group levels for example's sake)
modseltable <- pn.modselect.step(subdata$age, subdata$mass,
    subdata$id, penaliz='1*(n)', existing = TRUE)  
}
)
