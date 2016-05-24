SSposnegRichards <- structure(function(x, Asym = NA, 
    K = NA, Infl = NA, M = NA, RAsym = NA, Rk = NA, Ri = NA, RM = NA, 
    modno, pn.options, Envir = ".GlobalEnv") {
    Envir1 <- try(FPCEnv$env,silent=T)
    env1ck <- try(is.environment(FPCEnv$env),silent=T)
    envck <- try(is.environment(Envir),silent=T)
    env.ck<-2
    if(envck == FALSE | class(envck) == "try-error") env.ck <- (env.ck - 1)
    if(env1ck == FALSE | class(env1ck) == "try-error") env.ck <- (env.ck - 1)
    if(env.ck == 2) {
    modselck<- try(get("mod.sel", envir = FPCEnv), silent =T)
   if(class(modselck)[1] != "try-error" & modselck == TRUE) {
     if(identical(Envir, Envir1) == FALSE & 
    	identical(Envir,.GlobalEnv) == TRUE){
    	Envir <- Envir1
    	}
    	}
    }
    if(env.ck == 1 & (envck == FALSE | class(envck) == "try-error")) Envir <- Envir1
    FPCEnv$env <- Envir
    pcall <- as.list((sys.call()))
    pnoptname <- as.character( pcall[["pn.options"]])
    nmpnmodelparams <- c("Asym", "K", "Infl", "M", "RAsym", "Rk",
   	"Ri", "RM", "first.y", "x.at.first.y", "last.y", "x.at.last.y",
   	"twocomponent.x","verbose","force4par","force.nonmonotonic")
    nmpnmodelparamsbounds <- c("Amin", "Amax", "Kmin", "Kmax", "Imin",
	"Imax", "Mmin", "Mmax", "RAmin", "RAmax", "Rkmin", "Rkmax",
	"Rimin", "Rimax", "RMmin", "RMmax") 
    try({
    parload <- try(get(pnoptname, envir = Envir), silent=T)
    if(class(parload) == "try-error") 
       stop("pn.options file not found in specified environment. Please check your call
       and the last call (if any) to modpar() to ensure they match")
    }, silent = TRUE)
    try({pnmodelparamsbounds <- parload[ names(parload) %in% nmpnmodelparamsbounds]
         pnmodelparams <- parload[ names(parload) %in% nmpnmodelparams] },silent = TRUE)
    fractM <- M
    fractRM <- RM
    if( !is.numeric(modno)) modno <- get(as.character(modno), envir = Envir)
    if (modno > 19) {
        fractM <- 1/pnmodelparams$M
        M <- pnmodelparams$M
     }
    if (modno < 17 | modno >= 18) {
    if (modno == 2 | modno == 22 | modno == 10 | modno == 30 | 
        modno == 11 | modno == 31 | modno == 12 | modno == 32 | 
        modno == 19 | modno == 13 | modno == 33 | modno == 14 | 
        modno == 34 | modno == 15 | modno == 35 | modno == 16 | 
        modno == 36 | modno == 20 | modno == 32) {
        fractRM <- 1/pnmodelparams$RM
        RM <- pnmodelparams$RM
      }
    if (modno == 3 | modno == 23 | modno == 5 | modno == 25 | 
        modno == 7 | modno == 27 | modno == 9 | modno == 29 | 
        modno == 10 | modno == 30 | modno == 12 | modno == 32 | 
        modno == 19 | modno == 14 | modno == 34 | modno == 16 | 
        modno == 36 | modno == 20 | modno == 32) {
        RAsym <- pnmodelparams$RAsym
    } else {
    }
    if (modno == 3 | modno == 23 | modno == 4 | modno == 24 | 
        modno == 5 | modno == 25 | modno == 6 | modno == 26 | 
        modno == 10 | modno == 30 | modno == 11 | modno == 31 | 
        modno == 12 | modno == 32 | modno == 19 | modno == 13 | 
        modno == 33 | modno == 20 | modno == 32) {
        Rk <- pnmodelparams$Rk
     }
    if (modno == 4 | modno == 24 | modno == 5 | modno == 25 | 
        modno == 8 | modno == 28 | modno == 9 | modno == 29 | 
        modno == 11 | modno == 31 | modno == 12 | modno == 32 | 
        modno == 19 | modno == 15 | modno == 35 | modno == 16 | 
        modno == 36 | modno == 20 | modno == 32) {
        Ri <- pnmodelparams$Ri
    }
    try({
    if (Re(as.complex(1 + M[1] * exp(-K[1] * (max(x) - Infl[1])))) < 0) {
     	    fractM <- round(1/M)
    } else {
            fractM <- 1/M
    }
    if (Re(as.complex(1 + RM[1] * exp(-Rk[1] * (max(x) - Ri[1])))) < 0) {
     	    fractRM <- round(1/RM)
    } else {
       	    fractRM <- 1/RM
    }
    },silent = TRUE)
    }
  if (!is.na(pnmodelparams$twocomponent.x )) {
        if (pnmodelparams$force4par == TRUE)
            stop("Cannot force a two component Richards model to have a single component.Either set force4par to FALSE
            	 or twocomponent.x to NA")
        c((Asym/Re(as.complex(1 + M * exp(-K * ((x[x <= pnmodelparams$twocomponent.x]) -
            Infl)))^(fractM))), (RAsym/Re(as.complex(1 + RM *
            exp(-Rk * ((x[x > pnmodelparams$twocomponent.x]) -
                Ri)))^(fractRM))))
   } else {
        if (modno >= 17 & modno < 18) {
            if (pnmodelparams$force4par == TRUE) {
                Asym/Re(as.complex(1 + exp(Infl - x)/M))
            } else {
            if(modno == 17.1 | modno == 17.3) RAsym <- (-Asym)
            if(modno == 17.2 | modno == 17.3) RM <- M
             (Asym/Re(as.complex(1 + exp(Infl - x)/M))) + 
                  (RAsym/Re(as.complex(1 + exp(Ri - x)/RM)))
            }
        } else {
            if (pnmodelparams$force4par == TRUE | (pnmodelparams$force.nonmonotonic == FALSE & modno == 12) |
            (pnmodelparams$force.nonmonotonic == FALSE & modno == 32) ) {
                (Asym/Re(as.complex(1 + M * exp(-K * ((x) - Infl)))^(fractM)))
            } else {
                (Asym/Re(as.complex(1 + M * exp(-K * ((x) - Infl)))^(fractM))) + 
                  (RAsym/Re(as.complex(1 + RM * exp(-Rk * ((x) - 
                    Ri)))^(fractRM)))
            }
        }
    }
}, initial = function(mCall, LHS, data, ...) {
    Envir1 <- try(FPCEnv$env,silent=T)
    env1ck <- try(is.environment(FPCEnv$env),silent=T)
    envck <- try(is.environment(Envir),silent=T)
    env.ck<-2
    if(envck == FALSE | class(envck) == "try-error") env.ck <- (env.ck - 1)
    if(env1ck == FALSE | class(env1ck) == "try-error") env.ck <- (env.ck - 1)
     if(env.ck == 2) {
        modselck<- try(get("mod.sel", envir = FPCEnv), silent =T)
       if(class(modselck)[1] != "try-error" & modselck == TRUE) {
         if(identical(Envir, Envir1) == FALSE & 
        	identical(Envir,.GlobalEnv) == TRUE){
        	Envir <- Envir1
        	}
        	}
    }
     if(env.ck == 1 & (envck == FALSE | class(envck) == "try-error")) Envir <- Envir1
    if(exists("Envir") == F) stop("Environment not previously specified - argument not
    successfully transfered: \n run modpar or manually assign Envir value to FPCEnv$env
    e.g. assign('env', FlexParamCurve:::FPCEnv, envir = FlexParamCurve:::FPCEnv")
    FPCEnv$env <- Envir
    xy <- sortedXyData(mCall[["x"]], LHS, data)
    substparams<-NA
    counterR<-0
    optnm <- mCall[["pn.options"]]
    optvarnm <- as.character(optnm)
    if(exists(optvarnm, envir = Envir) == FALSE) { 
    optvar<- FPCEnv$pnoptnm
     } else {
     optvar<- get(optvarnm, envir = Envir)
     }
     if(length(optvar$taper.ends) == 1){
    	taper.ends <- optvar$taper.ends
    	} else {
    	taper.ends <- 0.45
    					}
     if(length(optvar$modpar) == 1) {
     frommodpar <- TRUE
     				} else {
     frommodpar <- FALSE    				
     				    }
     if(length(optvar$width.bounds) == 1) {
     smallmult <- 1 * as.numeric(optvar$width.bounds)
     mediummult <- 2  * as.numeric(optvar$width.bounds)
     largemult <- 3  * as.numeric(optvar$width.bounds)
     					  } else {
     smallmult <- 1
     mediummult <- 2
     largemult <- 3
     					  }
     pntbndserr <- FALSE
     if(length(optvar$bounds.error) == 1) pntbndserr <- TRUE
  	    skel <- rep(list(1), 16)
	    pnmodelparams <- c(rep(NA, 12),rep(FALSE,4))
	    pnmodelparams <- relist(pnmodelparams, skel)	    
	    names(pnmodelparams) <- c("Asym", "K", "Infl", "M", "RAsym", "Rk",
		"Ri", "RM", "first.y", "x.at.first.y", "last.y", "x.at.last.y",
		"twocomponent.x","verbose","force4par", "force.nonmonotonic")  
	    skel <- rep(list(1), 16)
	    pnmodelparamsbounds <- c(rep(NA, 16))
	    pnmodelparamsbounds <- relist(pnmodelparamsbounds, skel)
	    names(pnmodelparamsbounds) <- c("Amin", "Amax", "Kmin", "Kmax", "Imin",
		"Imax", "Mmin", "Mmax", "RAmin", "RAmax", "Rkmin", "Rkmax",
		"Rimin", "Rimax", "RMmin", "RMmax")
	    userpar <- optvar[which(names(optvar) %in% names(pnmodelparams))]
	    userparbounds <- optvar[which(names(optvar) %in% names(pnmodelparamsbounds))]
	    parseval <- function(txt1,evtxt,txt2) {
	    callout<- parse(text=sprintf("%s",paste(txt1,as.character(evtxt),txt2,sep="")))
	    return(eval(callout))
    	        				  }  
	    options(warn=-1)
	    pnmodelparams[which(names(pnmodelparams) %in% names(userpar))] <- userpar
	    pnmodelparamsbounds[which(names(pnmodelparamsbounds) %in% names(userparbounds))] <- userparbounds	
	    options(warn=0)
	    parsprovided <- pnmodelparams[which(!is.na(pnmodelparams))] 
	    boundsprovided <- pnmodelparamsbounds[which(!is.na(pnmodelparamsbounds))]
	    splnms <- function(x) strsplit(names(x),"m")
	    try(if(names(parsprovided)[1] == "Asym") names(parsprovided)[1] <-"A",silent=TRUE)
	    try(if(names(parsprovided)[3] == "Infl") names(parsprovided)[3] <-"I",silent=TRUE)
	    matchbounds <-  names(pnmodelparamsbounds[which(sapply(splnms(boundsprovided), function(x) substring(x[[1]],1,2) ) %in% 
		sapply(names(parsprovided), function(x) substring(x[[1]],1,2)))])
	    "%w/o%" <- function(x, y) x[!x %in% y] #--  x without y
	    unmatchbounds <- (names(boundsprovided) %w/o% matchbounds) %w/o% c("first.y", "x.at.first.y", "last.y", "x.at.last.y",
		"twocomponent.x","verbose","force4par", "force.nonmonotonic")
	    modno <- mCall[["modno"]]
	    if(modno == 12 | modno == 32 | optvar$force4par == TRUE) unmatchbounds <- (names(boundsprovided) %w/o% 
	    	 matchbounds) %w/o% c("first.y", "x.at.first.y", "last.y", "x.at.last.y",
		"twocomponent.x","verbose","force4par", "force.nonmonotonic", "RAmin","RAmax","Rkmin","Rkmax","Rimin","Rimax","RMmin","RMmax")
	    .paramsestimated <- try( FPCEnv$.paramsestimated ,silent = TRUE)
	    if( class(.paramsestimated) == "try-error" ) .paramsestimated <- TRUE
	    if(modno != 18 & modno != 19 ){
	    if(length(unmatchbounds) > 0 & .paramsestimated == TRUE ) {
	    print("WARNING: In pn.options:  Some parameters specified are missing min/max bounds. Running modpar to populate these bounds.")
	    print("Note: this will also populate any parameters appropriate to the specified modno that are also missing")
       	    parseval("substparams <- modpar(xy$x, xy$y, first.y = pnmodelparams$first.y, x.at.first.y = pnmodelparams$x.at.first.y,
    		last.y = pnmodelparams$last.y,
    		x.at.last.y = pnmodelparams$x.at.last.y, twocomponent.x = pnmodelparams$twocomponent.x, verbose = pnmodelparams$verbose,
    		force4par = pnmodelparams$force4par, Envir = Envir, force.nonmonotonic = pnmodelparams$force.nonmonotonic
    		, pn.options =\"", as.character(optvarnm),"\")")
    	options(warn=-1)
    	pnmodelparams[which(names(pnmodelparams) %in% names(substparams))] <- substparams[which(names(substparams) %in% names(pnmodelparams))]
    	pnmodelparamsbounds[which(names(pnmodelparamsbounds) %in% names(substparams))] <- substparams[which(names(substparams) %in% names(pnmodelparamsbounds))]
   	options(warn=0)		  				     
  								}
  					} else {
  					}
    valexp <- sapply( c(unlist(pnmodelparams),unlist(pnmodelparamsbounds), 
        	unlist( optvar[names(optvar) %w/o% c(names(pnmodelparamsbounds),names(pnmodelparams)) ]) ), function(x) list(x))
    names( valexp[1:32] ) <- names( c(pnmodelparams[1:16],pnmodelparamsbounds) )
    valexp[c(1:12,17:32)] <- as.numeric( valexp[c(1:12,17:32)] )
    valexp[14:16] <- as.logical( valexp[14:16] )
    if(is.numeric(pnmodelparams$twocomponent.x)) {
        valexp[13] <- as.numeric( valexp[13] )
    } else {
  	valexp[13] <- as.logical( valexp[13] )     	
    }
    assign(optvarnm, valexp, envir = Envir) 
    modno <- mCall[["modno"]]
    if( !is.numeric(modno)) modno <- get(as.character(modno), envir = Envir)
    invars <- mCall
    modelparams <- valexp[1:16]
    if (is.na(modelparams[1]) == TRUE) {
        if (modno != 18 & modno != 19) {
            stop("ERROR: saved parameters empty or Asym missing: run function modpar before using selfStart functions: see ?modpar")
        }
    }
    if(!is.na(modelparams$twocomponent.x) & modelparams$force4par == TRUE) 
    	stop("Cannot force a two component Richards model to have a single component.Set force4par to FALSE")
    if(modelparams$force4par == TRUE) {
    if("RAsym" %in% names(invars)) stop("force4par is TRUE, therefore use models without RAsym, Rk, Ri or RM")
    if("Rk" %in% names(invars)) stop("force4par is TRUE, therefore use models without RAsym, Rk, Ri or RM")
    if("Ri" %in% names(invars)) stop("force4par is TRUE, therefore use models without RAsym, Rk, Ri or RM")
    if("RM" %in% names(invars)) stop("force4par is TRUE, therefore use models without RAsym, Rk, Ri or RM")
    }
    if (!"Asym" %in% names(invars)) stop(paste("Parameter Asym required for modno = ", modno, " but is absent from user-provided call", sep = ""))
    if (!"Infl" %in% names(invars)) stop(paste("Parameter Infl required for modno = ", modno, " but is absent from user-provided call", sep = ""))
    if (!"K" %in% names(invars) & (modno < 17 | modno >= 18) ) 
    	stop(paste("Parameter K required for modno = ", modno, " but is absent from user-provided call", sep = ""))
    if (modno <= 19) {
    	if (!"M" %in% names(invars)) stop(paste("Parameter M required for modno = ", modno, " but is absent from user-provided call", sep = ""))
    	if (.paramsestimated != FALSE & !"M" %in% names( pnmodelparams[which(!is.na(pnmodelparams))] )) 
    		stop(paste("Parameter M required for modno = ", modno, " but is absent from pn.options list in user-provided call", sep = ""))
    }
    if (modno != 2 & modno != 22 & modno != 10 & modno != 30 & 
        modno != 11 & modno != 31 & modno != 12 & modno != 32 & 
        modno != 19 & modno != 13 & modno != 33 & modno != 14 & 
        modno != 34 & modno != 15 & modno != 35 & modno != 16 & 
        modno != 36 & modno != 20 & modno != 32 & modno != 17.2 & modno != 17.3) {        
        if (!"RM" %in% names(invars)) stop(paste("Parameter RM required for modno = ", modno, " but is absent from user-provided call", sep = ""))
   	if (.paramsestimated != FALSE & !"RM" %in% names( pnmodelparams[which(!is.na(pnmodelparams))] )) 
    		stop(paste("Parameter RM required for modno = ", modno, " but is absent from pn.options list in user-provided call", sep = ""))
    }
    if (modno != 3 & modno != 23 & modno != 5 & modno != 25 & 
        modno != 7 & modno != 27 & modno != 9 & modno != 29 & 
        modno != 10 & modno != 30 & modno != 12 & modno != 32 & 
        modno != 19 & modno != 14 & modno != 34 & modno != 16 & 
        modno != 36 & modno != 20 & modno != 32 & modno != 17.1 & modno != 17.3) {
        if (!"RAsym" %in% names(invars)) stop(paste("Parameter RAsym required for modno = ", modno, " but is absent from user-provided call", sep = ""))
   	if (.paramsestimated != FALSE & !"RAsym" %in% names( pnmodelparams[which(!is.na(pnmodelparams))] )) 
    		stop(paste("Parameter RAsym required for modno = ", modno, " but is absent from pn.options list in user-provided call", sep = ""))
    }
    if (modno != 3 & modno != 23 & modno != 4 & modno != 24 & 
        modno != 5 & modno != 25 & modno != 6 & modno != 26 & 
        modno != 10 & modno != 30 & modno != 11 & modno != 31 & 
        modno != 12 & modno != 32 & modno != 19 & modno != 13 & 
        modno != 33 & modno != 20 & modno != 32 & modno != 17 &
        modno != 17.1 & modno != 17.2 & modno != 17.3) {
        if (!"Rk" %in% names(invars)) stop(paste("Parameter Rk required for modno = ", modno, " but is absent from user-provided call", sep = ""))
	if (.paramsestimated != FALSE & !"Rk" %in% names( pnmodelparams[which(!is.na(pnmodelparams))] )) 
    		stop(paste("Parameter Rk required for modno = ", modno, " but is absent from pn.options list in user-provided call", sep = ""))
    }
    if (modno != 4 & modno != 24 & modno != 5 & modno != 25 & 
        modno != 8 & modno != 28 & modno != 9 & modno != 29 & 
        modno != 11 & modno != 31 & modno != 12 & modno != 32 & 
        modno != 19 & modno != 15 & modno != 35 & modno != 16 & 
        modno != 36 & modno != 20 & modno != 32 ) {
        if (!"Ri" %in% names(invars)) stop(paste("Parameter Ri required for modno = ", modno, " but is absent from user-provided call", sep = ""))
	if (.paramsestimated != FALSE & !"Ri" %in% names( pnmodelparams[which(!is.na(pnmodelparams))] )) 
    		stop(paste("Parameter Ri required for modno = ", modno, " but is absent from pn.options list in user-provided call", sep = ""))    
    }
    if (is.na(modelparams$first.y) == FALSE) {
        if (is.na(modelparams$x.at.first.y))
        	modelparams$x.at.first.y <- 0
        modrngx <- diff(c( c(min(xy$x, na.rm = TRUE), max(xy$x, na.rm = TRUE)) ))
        xy <- rbind(xy, c(seq( (modelparams$x.at.first.y - (0.25*modrngx) ),
        	modelparams$x.at.first.y + (0.25*modrngx) ,length.out=10)
        	, rep(modelparams$first.y, 10) ))
        xy <- xy[order(xy$x), ]
    }
    if (is.na(modelparams$last.y) == FALSE) {
        if (is.na(modelparams$x.at.last.y)) 
            stop("x.at.last.y must be specified if last.y is not NA")
        if (is.na(modelparams$last.y)) 
            stop("last.y must be specified if x.at.last.y is not NA")
        modrngx <- diff(c( c(min(xy$x, na.rm = TRUE), max(xy$x, na.rm = TRUE)) ))
        xy <- rbind(xy, c(seq( (modelparams$x.at.last.y - (0.25*modrngx) ),
        	modelparams$x.at.last.y + (0.25*modrngx) ,length.out=10)
        	, rep(modelparams$last.y, 10) ))
        xy <- xy[order(xy$x), ]
    }
    roundparscale <- function(x) {
    options(warn=-1)
    if(x == 0 | length(x) == 0 | is.na(x) == TRUE) {
    	y = 1
    } else {
        y = (10^ceiling(log10(abs(x))))/10
        if(x < 0) y = y * -1
    }
    if(is.na(y) == TRUE) y = 1
    options(warn=0)
    return(y)
    				}
    sub <- subset(xy, xy$y == max(xy$y))
    tAsym <- mean(sub$x)
    SSposnegRichardsF <- function(x, Asym, K, Infl, M, RAsym, 
        Rk, Ri, RM) (Asym/Re(as.complex(1 + M * exp(-K * (x - 
        Infl)))^(1/M))) + (RAsym/Re(as.complex(1 + RM * exp(-Rk * 
        (x - Ri)))^(1/RM)))
    SSposnegRichardsFM <- function(x, Asym, K, Infl, M, RAsym, 
        Rk, Ri, RM) (Asym/Re(as.complex(1 + M * exp(-K * (x - 
        Infl)))^round(1/M))) + (RAsym/Re(as.complex(1 + RM * 
        exp(-Rk * (x - Ri)))^(1/RM)))
    SSposnegRichardsFRM <- function(x, Asym, K, Infl, M, RAsym, 
        Rk, Ri, RM) (Asym/Re(as.complex(1 + M * exp(-K * (x - 
        Infl)))^(1/M))) + (RAsym/Re(as.complex(1 + RM * exp(-Rk * 
        (x - Ri)))^round(1/RM)))
    SSposnegRichardsFMRM <- function(x, Asym, K, Infl, M, RAsym, 
        Rk, Ri, RM) (Asym/Re(as.complex(1 + M * exp(-K * (x - 
        Infl)))^round(1/M))) + (RAsym/Re(as.complex(1 + RM * 
        exp(-Rk * (x - Ri)))^round(1/RM)))
    SSposnegRichardsF17 <- function(x, Asym, Infl, M, RAsym, 
        Ri, RM) (Asym/Re(as.complex(1 + exp(Infl - x)/M))) + 
        (RAsym/Re(as.complex(1 + exp(Ri - x)/RM)))
    richards <- function(x, Asym, K, Infl, M) Asym/Re(as.complex(1 + 
        M * exp(-K * (x - Infl)))^(1/M))
    richards3 <- function(x, Asym, K, Infl) Asym/Re(as.complex(1 + 
        pnmodelparams$M * exp(-K * (x - Infl)))^(1/pnmodelparams$M))
    richards173 <- function(x, Asym, Infl, M) Asym/Re(as.complex(1 + 
        exp(Infl - x)/M))
    modelparamsbounds <- pnmodelparamsbounds
    if (is.na(modelparamsbounds[1]) == TRUE) {
        if (modno != 18 & modno != 19) {
            stop("ERROR: fix parameter bounds empty: run function modpar before using selfStart functions: see ?modpar")
        }
    }
    if (modno != 18 & modno != 19) {
    modelparams[1:8]<-lapply(modelparams[1:8],function(x) replace(x,is.na(x),1))
    modelparamsbounds<-lapply(modelparamsbounds,function(x) replace(x,is.na(x),1))
        Amax = modelparamsbounds$Amax
        Amin = modelparamsbounds$Amin
        Kmax = modelparamsbounds$Kmax
        Kmin = modelparamsbounds$Kmin
        Imax = modelparamsbounds$Imax
        Imin = modelparamsbounds$Imin
        Mmax = modelparamsbounds$Mmax
        Mmin = modelparamsbounds$Mmin
        RAmax = modelparamsbounds$RAmax
        RAmin = modelparamsbounds$RAmin
        Rkmax = modelparamsbounds$Rkmax
        Rkmin = modelparamsbounds$Rkmin
        Rimax = modelparamsbounds$Rimax
        Rimin = modelparamsbounds$Rimin
        RMmax = modelparamsbounds$RMmax
        RMmin = modelparamsbounds$RMmin
    }
    if (is.na(modelparams$twocomponent.x) == TRUE | modelparams$twocomponent.x == 
        TRUE) {
        movavg <- rep(NA, nrow(xy))
        windw <- (round(nrow(xy)/30) + 2)
        if (windw%%2) {
            subtr <- floor(windw/2)
            addtr <- floor(windw/2)
        } else {
            subtr <- floor(windw/2)
            addtr <- floor(windw/2) - 1
        }
        for (i in (subtr + 1):(nrow(xy) - addtr)) {
            movavg[i] <- mean(c(xy$y[(i - subtr):(i + addtr)]), 
                na.rm = TRUE)
        }
        signvalmov <- 1
        if (max(xy$x) != max(abs(xy$x))) 
            signvalmov <- -1
        movavg <- abs(movavg)
        xval <- xy$x[(subtr + 1):length(movavg)]
        xval <- abs(xval)
        movavg <- movavg[(subtr + 1):length(movavg)]/max(movavg, 
            na.rm = TRUE)
        tmp <- data.frame(xval, movavg)
        tmp <- subset(tmp, !is.na(tmp$movavg))
        xval <- tmp$xval
        movavg <- tmp$movavg
        slpe <- rep(NA, length(movavg))
        if(length(movavg) <2) {
        	        slpe = 1
        } else {
        for (i in 2:length(movavg)) {
            slpe[i] <- (movavg[i] - movavg[i - 1])/(xval[i] - 
                xval[i - 1])
        }
        }
         Intage <- NA
        stationry <- data.frame(xval, movavg)
        datalength <- nrow(stationry)
        stationry <- subset(stationry, slpe > -0.005 & slpe < 
            0.005)
        subid.1 <- as.numeric(row.names(stationry))
        teststart <- subid.1[1]
        peakage1<-min(xy$x, na.rm = TRUE)
        peakage2<-min(xy$x, na.rm = TRUE)
        if(nrow(stationry) > 0 ){
        if(teststart/datalength <= 0.2) {
            newstart <- ceiling(0.21*datalength)
            subid.2 <- subid.1[subid.1 >= newstart]
            stationry <- data.frame(xval, movavg)
            stationry <- stationry[subid.2,]
        }
                                 }
        if (length(stationry[, 1]) > 2) {
            savpeakmass2<-NA
            savpeakage2<-NA
            stdist <- rep(NA, length(stationry[, 1]))
            stdist[2:(length(stationry[, 1]))] <- stationry[(1:(length(stationry[, 
                1]) - 1)), 1] - stationry[2:(length(stationry[, 
                1])), 1]
            peakage1 <- subset(stationry[, 1], stdist == max(stdist, 
                na.rm = TRUE))[1]
            savpeakage1 <- peakage1
            peakmass1 <- subset(stationry[, 2], stdist == max(stdist, 
                na.rm = TRUE))[1]
            stationry <- subset(stationry, stationry[, 1] != 
                peakage1)
            cnting <- TRUE
            startcnt <- -1
            adjmass1 <- 0
            while (cnting == TRUE) {
                stdist <- rep(NA, length(stationry[, 1]))
                stdist <- stationry[1:(length(stationry[, 1]) - 
                  1), 1] - stationry[2:(length(stationry[, 1])), 
                  1]
                peakage2 <- subset(stationry[, 1], stdist == 
                  max(stdist, na.rm = TRUE))[1]
                peakmass2 <- subset(stationry[, 2], stdist == 
                  max(stdist, na.rm = TRUE))[1]
             if( length(peakage2) != 0 & length(peakage1) != 0 ) {
             if( !is.na(peakage2) & !is.na(peakage1) ) {
               	 if (abs(peakage2 - peakage1) > max(abs(xy$x)/8)) {
                 	 checkij <- subset(data.frame(xval, movavg), 
                	    xval >= peakage1 & xval <= peakage2)
                  if (nrow(checkij) > 1) {
                  	  cnting <- FALSE
                  } else {
                    	  cnting <- TRUE
                  }	
                 }
                 if (cnting == FALSE) {
                  if (startcnt == -1) {
                    savpeakmass2 <- peakmass2
                    savpeakage2 <- peakage2
                  } else {
                    tempp <- c(peakage1, peakmass1)
                    peakmass1 <- peakmass2
                    peakage1 <- peakage2
                    if (startcnt == 2) {
                      peakage2 <- tempp[1]
                      peakmass2 <- tempp[2]
                    } else {
                      peakage2 <- savpeakage2
                      peakmass2 <- savpeakmass2
                    }
                  }
                  if (startcnt != 0) {
                    if (peakmass1 > peakmass2) {
                      if (peakage1 > peakage2) {
                        temppk <- c(peakage1, peakmass1)
                        peakage1 <- peakage2
                        peakmass1 <- peakmass2
                        peakage2 <- temppk[1]
                        peakmass2 <- temppk[2]
                      } else {
                        cnting <- TRUE
                      }
                    } else {
                      checkij <- subset(data.frame(xval, movavg), 
                        xval >= peakage1 & xval <= peakage2)
                      checkij <- checkij[order(checkij$xval), 
                        ]
                      if (startcnt != 2) 
                        startcnt <- 0
                      for (ijcnt in 2:nrow(checkij)) {
                        if (checkij[ijcnt, 2] < checkij[(ijcnt - 
                          1), 2]) 
                          startcnt <- 1
                      }
                      if (startcnt == 0) {
                        cnting <- TRUE
                        stationry <- subset(stationry, stationry[, 
                          1] < peakage1)
                        if (length(stationry[, 1]) > 1) {
                          peakage2 <- max(stationry[, 1],na.rm = TRUE)
                          dataf <- data.frame(movavg, xval)
                          peakmass2 <- subset(dataf$movavg, dataf$xval == 
                            peakage2)[1]
                          peakmass1 <- subset(dataf$movavg, dataf$xval == 
                            peakage1)[1]
                          adjmass1 <- 1
                        }
                      }
                    }
                  }
                }
                if (nrow(stationry) < 3 & cnting == TRUE) {
                  if (startcnt == 0) {
                    stationry <- data.frame(xval, movavg)
                    stationry <- subset(stationry, slpe > -0.005 & 
                      slpe < 0.005)
                    stationry <- subset(stationry, stationry[, 
                      1] > savpeakage1)
                    if (nrow(stationry) > 1) {
                      peakage2 <- stationry[2, 1]
                      peakage1 <- savpeakage1
                      dataf <- data.frame(movavg, xval)
                      peakmass2 <- subset(dataf$movavg, dataf$xval == 
                        peakage2)[1]
                      peakmass1 <- subset(dataf$movavg, dataf$xval == 
                        peakage1)[1]
                      startcnt <- 2
                    }
                  } else {
                    Intage <- savpeakage1
                    cnting <- FALSE
                  }
                } else {
                  if (length(stationry[, 1]) < 3) {
                    Intage <- savpeakage1
                    cnting <- FALSE
                  }
                }
                if (adjmass1 == 0) {
                  stationry <- subset(stationry, stationry[, 
                    1] != peakage2)
                } else {
                  adjmass1 <- 0
                }
             } else { 
                    cnting <- FALSE
		    if(length(stationry[, 1])==0) {
        	    peakage1 <- subset(xy$x, xy$y == max(xy$y))[1]
        	    peakage2 <- max(xy$x)
        	    } else {
        	    peakage1 <- min(stationry[, 1], na.rm = TRUE)
        	    peakage2 <- max(stationry[, 1], na.rm = TRUE)
        	    }
             }
             }
            }
        } else {
        	    if(length(stationry[, 1])==0) {
        	    peakage1 <- subset(xy$x, xy$y == max(xy$y))[1]
        	    peakage2 <- max(xy$x)
        	    } else {
        	    peakage1 <- min(stationry[, 1], na.rm = TRUE)
        	    peakage2 <- max(stationry[, 1], na.rm = TRUE)
        	    }
        }
        if (peakage1 > peakage2) {
            peaksav <- c(peakage1, peakage2)
            peakage1 <- peaksav[2]
            peakage2 <- peaksav[1]
        }
	if (is.na(Intage) == TRUE) {
            check1 <- subset(data.frame(xval, movavg), xval >= 
                peakage1 & xval <= peakage2)
            check1 <- check1[order(check1$xval), ]
            cntj = 2
            cnting1 <- TRUE
            dwn <- FALSE
            turnpt <- NA
            while (cnting1 == TRUE) {
                if (cntj == nrow(check1) | nrow(check1) < 3) {
                  cnting1 <- FALSE
                  options(warn=-1)
                  try(Intage <- as.numeric(subset(check1[, 1], check1[, 
                    2] == max(check1[, 2], na.rm = TRUE))[1]), silent =TRUE)
                  options(warn=0)
                } else {
                  cntj = cntj + 1
                  if (check1[cntj - 1, 2] < (check1[(cntj - 2), 
                    2])) {
                    dwn <- TRUE
                    if (is.na(turnpt)) 
                      turnpt <- check1[(cntj - 1), 2]
                  }
                  if(!is.na(turnpt)) {
                  if (dwn == TRUE & check1[cntj, 2] > (turnpt + 
                    (max(movavg, na.rm = TRUE) * 0.05))) {
                    cnting1 <- FALSE
                    Intage <- check1[cntj, 1]
                  }
                  } else {
                  }
                }
            }
        }
 	mavrng <- diff(range(movavg))
        xrng <- diff(range(xval))
        plateu <- subset(data.frame(xval, movavg), 
        	movavg < (movavg[peakage2] + 0.15*mavrng) &
        	movavg > (movavg[peakage2] - 0.15*mavrng) )
        stplat <- sd(plateu[,2])
        plateu <- subset(data.frame(xval, movavg), 
	        	movavg < (movavg[peakage2] + stplat) &
	        	movavg > (movavg[peakage2] - stplat) )
        if(nrow(plateu) == 0) {
        availx <- data.frame(xval,movavg)
	closepeakage2 <- which( abs( availx[,1]-peakage2 ) == min( abs(availx[,1]-peakage2) ) )
	plateu <- data.frame(closepeakage2, movavg[closepeakage2])
	}
	try({
	if(Intage %in% plateu[,1] == FALSE){
	Intage <- plateu[ floor(nrow(plateu)/2),1]
	}
	},silent = TRUE)
	if(length(Intage) == 0) Intage <- NA
        if (is.na(Intage) | quantile(plateu[,1], probs = 0.90, na.rm = TRUE) > (peakage2 + 0.2*xrng) ) 
            Intage <- as.numeric(subset(xval, movavg == max(movavg, 
                na.rm = TRUE))[1])
        Intage <- Intage * signvalmov
        if (is.na(modelparams$twocomponent.x) == FALSE) {
            if (modelparams$twocomponent.x == TRUE) {
                modelparams$twocomponent.x <- Intage
                pnmodelparams$twocomponent.x <- Intage
                optvar$twocomponent.x <- Intage
                valexp <- sapply( c(unlist(pnmodelparams),unlist(pnmodelparamsbounds), 
        		unlist( optvar[names(optvar) %w/o% c(names(pnmodelparams),names(pnmodelparamsbounds)) ]) ), function(x) list(x))
   		names( valexp[1:32] ) <- names( c(pnmodelparams[1:16],pnmodelparamsbounds) )
		valexp[c(1:12,17:32)] <- as.numeric( valexp[c(1:12,17:32)] )
		valexp[14:16] <- as.logical( valexp[14:16] )
    		try(valexp$taper.ends <- as.numeric(valexp$taper.ends), silent = TRUE)
    		assign(optvarnm, valexp, envir = Envir )    		    
            }
        }
    } else {
        Intage <- modelparams$twocomponent.x
    }
    if (is.na(Intage) == TRUE) 
          Intage <- as.numeric(subset(xy$x, xy$y == max(xy$y, 
                na.rm = TRUE))[1])
    if (Intage == min(xy$x)) 
        Intage <- max(xy$x)
    if (modelparams$force4par == TRUE | (modelparams$force.nonmonotonic == FALSE & modno == 12) |
            (modelparams$force.nonmonotonic == FALSE & modno == 32) ) {
        xyE <- xy
        maxIval<- try(max(xyE$x, na.rm = TRUE) - 
        	(max(diff( range(xyE$x, na.rm = TRUE))*.05)), silent=TRUE)
        if(class(maxIval) == "try-error") maxIval <- max(xy$x, na.rm=TRUE)  
    } else {
        xyE <- subset(xy, xy$x <= Intage)
        maxIval<- try(max(xyE$x, na.rm = TRUE) - 
        	(max(diff( range(xyE$x, na.rm = TRUE))*.05)), silent=TRUE)
        if(class(maxIval) == "try-error") maxIval <- max(xy$x, na.rm=TRUE)      
        xyL <- data.frame(rep(NA, 1))
	xyL <- subset(xy, xy$x >= Intage)
        maxRival<- try(max(xyL$x, na.rm = TRUE) - 
        	(max(diff( range(xyL$x, na.rm = TRUE))*.05)), silent=TRUE)
  	if(class(maxRival) == "try-error") maxRival <- max(xy$x, na.rm=TRUE)        
    }
    if (nrow(xyE) < 5) {
        stop("ERROR: too few distinct input values to fit the positive Richards model, aborting")
    } else {
            moddata <- function(indata, absmin, Iage, taper.ends = 0) {
    	    savindata <- indata
            rng1 <- range(indata$x)
            rng2 <- range(indata$y)
            drng <- diff(rng2)
            indata <- indata[order(indata$x), ]
            if (max(indata$y[ ceiling( 0.8 * length(indata$y) ) : nrow(indata)], na.rm=TRUE)
            		< ((drng/2) + min(indata$y))) {
                starval <- (max(indata$y))
            } else {
                starval <- (min(indata$y))
            }
            Iagemod <- min(indata$x) - (diff(rng1) * abs(taper.ends))
            extrax <- seq(absmin, Iagemod, length.out = ceiling(length(indata$x)/20) + 
                1)
            extrax <- extrax[2:length(extrax - 1)]
            extrax <- sort(extrax)
            indata <- merge(indata, data.frame(x = extrax, y = rep(starval, 
                length(extrax))), all = TRUE)
            indata <- indata[order(indata$x), ]
            if (indata$y[nrow(indata)] < ((drng/2) + min(indata$y))) {
                endrval <- (min(indata$y))
            } else {
                endrval <- (max(indata$y))
            }
            Iagemod <- max(indata$x) + (diff(rng1) * abs(taper.ends))
            extrax <- seq(max(indata$x), Iagemod, length.out = ceiling(length(indata$x)/20) + 
                1)
            extrax <- extrax[2:length(extrax - 1)]
            extrax <- sort(extrax)
            if(taper.ends != 0) indata <- merge(indata, 
                data.frame(x = extrax, y = rep(endrval, 
                length(extrax))), all = TRUE)
            indata <- indata[order(indata$x), ]
            retdata <- list(indata)
            names(retdata) <-"indata"
            retdata[["savindata"]] <- savindata
            return(retdata)
        }
        trnsfrm <- function(indata, absmin, Iage, taper.ends) {
            savdata <- indata
            rng1 <- range(indata$x)
            rng2 <- range(indata$y)
            drng <- diff(rng2)
            orientx <- 0
            orienty <- 0
            dataorient <- 0
            indata <- indata[order(indata$x), ]
            indata <- moddata(indata, absmin, Iage, taper.ends = taper.ends)[["indata"]]
            savindata <- moddata(indata, absmin, Iage, taper.ends = taper.ends)[["savindata"]]
            names(indata) <- c("x","y")
            names(savindata) <- c("x","y")
            if ((rng1[1] <= 0 & rng1[2] <= 0)) {
                orientx <- 1
                indata$x <- (-indata$x)
                indata <- indata[order(indata$x), ]
            }
            if ((rng2[1] <= 0 & rng2[2] <= 0)) {
                orienty <- 1
                indata$y <- (-indata$y)
            }
            indata <- indata[order(indata$x), ]
            if (indata$y[nrow(indata)] < ((drng/2) + min(indata$y))) {
                if (indata$y[round(nrow(indata)/2)] < ((drng/2) + 
                  min(indata$y))) {
                  indata$y <- (abs(-(indata$y) + abs(max(indata$y))) + 
                    (min(indata$y)))
                  dataorient <- 1
                } else {
                  indata$x <- (abs(-(indata$x) + abs(max(indata$x))) + 
                    (min(indata$x)))
                  indata <- indata[order(indata$x), ]
                  dataorient <- 2
                }
            } else {
                if (indata$y[round(nrow(indata)/2)] < ((drng/2) + 
                  min(indata$y))) {
                  indata$y <- (abs(-(indata$y) + abs(max(indata$y))) + 
                    (min(indata$y)))   
		   indata$x <- (abs(-(indata$x) + abs(max(indata$x))) + 
                    (min(indata$x)))
                  indata <- indata[order(indata$x), ]
                  dataorient <- 3
                } else {
                }
            }
            rng <- range(indata$y)
            drng <- diff(rng)
            indata$prop <- (indata$y - rng[1L] + 0.05 * drng)/(1.1 * 
                drng)
            maxdat <- max(indata$x, na.rm = TRUE)
            transy <- (log(indata$prop/(1 - indata$prop)))
            Asym <- rep(NA, 31)
            K <- rep(NA, 31)
            Infl <- rep(NA, 31)
            M <- rep(NA, 31)
            indata$x <- indata$x/maxdat
            for (i in 1:31) {
                j <- c(-256, -96, -48, -24, -14, -7, -6, -5, 
                  -4, -3, -2, -1, -0.5, -0.25, -0.1, -1e-05, 
                  1e-05, 0.1, 0.25, 0.5, 1, 2, 3, 4, 5, 6, 12, 
                  24, 48, 96, 256)
                if (i == 1) 
                  irlogist <- as.vector(coef(lm(x ~ transy, data = indata)))
                try({
                  ir <- as.vector(coef(lm(x ~ transy, data = indata, 
                    weights = ((1 - indata$prop)^j[i]))))
                  M[i] <- exp(-abs(ir[1L] - irlogist[1]))
                  Infl[i] <- (ir[1L] + log(abs(M[i]))) * maxdat
                  K[i] <- (1/(ir[2L]))/maxdat
                  Asym[i] <- rng[2L]
                }, silent = TRUE)
            }
            indata$x <- indata$x * maxdat
            indata <- indata[order(indata$x), ]
            if (dataorient == 1 | dataorient == 3) {
                indata$y <- (abs(-(indata$y) + abs(max(indata$y))) + 
                  (min(indata$y)))
                if (dataorient != 3) {
                  if (orientx == 1 | orienty == 1) 
                    K <- (-K)
                  if (orientx == 0 & orienty == 0) 
                    K <- (-K)
                }
            }
            if (dataorient == 2 | dataorient == 3) {
                Infl <- max(abs(savdata$x)) - Infl
                if (dataorient != 3) {
                  if (orientx == 1 | orienty == 1) 
                    K <- -K
                  if (orientx == 0 & orienty == 0) 
                    K <- (-K)
                }
                indata$x <- (abs(-(indata$x) + abs(max(indata$x))) + 
                    (min(indata$x)))
                indata <- indata[order(indata$x), ]
            }
            if (orientx == 1) {
                indata$x <- (-indata$x)
                indata <- indata[order(indata$x), ]
                Infl <- (-Infl)
                K <- (-K)
            }
            if (orienty == 1) {
                indata$y <- (-indata$y)
                Asym <- (-Asym)
            }
            val <- data.frame(Asym, K, Infl, M)
            val <- na.omit(val)
            y0 <- Inf
            y1 <- Inf
            if(taper.ends != 0) {
       		evly <- function(val, y1, y0, i) {
                y1 <- sum((indata$y - richards(indata$x, val[i, 
                  1], val[i, 2], val[i, 3], val[i, 4]))^2)
                if (y1 < y0) {
                  y0 <- y1
                  indx <- i
                }
                return(c(y1, y0, indx))
            }            
            }else{
            	evly <- function(val, y1, y0, i) {
                y1 <- sum((savindata$y - richards(savindata$x, val[i, 
                  1], val[i, 2], val[i, 3], val[i, 4]))^2)
                if (y1 < y0) {
                  y0 <- y1
                  indx <- i
                }
                return(c(y1, y0, indx))
            }
            }
  		 if( max( val[ ,3] , na.rm = TRUE) < (min( indata$x, na.rm=TRUE)
  			+ (0.2* max(indata$x, na.rm=TRUE)) ) ) 
  			val[ ,3] <- val[ ,3] + (0.2* max(indata$x, na.rm=TRUE))
            for (i in 1:nrow(val)) {
                 try({
                  vald <- evly(val, y1, y0, i)
                  y1 <- vald[1]
                  y0 <- vald[2]
                  indx <- vald[3]
                }, silent = TRUE)
            }
    	    savy1 <- Inf
	    savy0 <- Inf
	    sav1y1 <- Inf
	    sav1y0 <- Inf          
            sav1indx <- indx
	    sav1y0 <- y0
	    y0 <- Inf
	    y1 <- Inf
	    sav1val<-val
            options(warn=-1)
            rngi <- 0.05*diff( c( min(indata$x, na.rm = TRUE), max(indata$x, na.rm = TRUE)) )
            rngyi <- diff( c( min(indata$y, na.rm = TRUE), max(indata$y, na.rm = TRUE)) )
            hway <- (min(indata$y,na.rm = TRUE) + (rngyi / 2))
            subdati <- subset(indata, indata$y >= (hway - 0.25*rngyi) 
            	& indata$y <= (hway + 0.25*rngyi) )
	    findi <- try(subdati$x[ ceiling(length(subdati$x) / 2) ] , silent = TRUE)
	    try(if( nrow(findi) == 0) findi <- indata$x[ which.max(indata$y[ indata$y <= 
	    	(min(indata$y,na.rm = TRUE) + (rngyi / 2)) ])], silent = TRUE)
	    if(class(findi)[1] == "try-error") findi <- indata$x[ which.max(indata$y[ indata$y <= 
	    	(min(indata$y,na.rm = TRUE) + (rngyi / 2)) ])]
	    try(
	    if( val[indx,3] < (findi - rngi) | val[indx,3] > (findi + rngi) ) {
	    val[,3] <- findi
	    	for (i in 1:nrow(val)) {
	                     try({
	                      vald <- evly(val, y1, y0, i)
	                      y1 <- vald[1]
	                      y0 <- vald[2]
               		      }, silent = TRUE)
	    	}
	    }
	    ,silent = TRUE)
	    savindx <- indx
	    savy0 <- y0
	    y0 <- Inf
	    y1 <- Inf
	    sav0val<-val
	    findk<-  function(a,k,i,m,desiredy,x) {
	    	k<- seq(-2*k,2*k,length.out=100)
	    	evalmod<-function(k) c(abs(desiredy-richards(x, a, k, i, m)) ,k)
		outmin<-data.frame(sapply( k, evalmod))
		if(which.min(outmin[1,])==1) {
			k<- seq(-200*k,-2*k,length.out=100)
			outmin<-data.frame(sapply( k, evalmod))
		}else{
			if(which.min(outmin[1,])==length(outmin[1,])) {
				k<- seq(2*k,200*k,length.out=100)
				outmin<-data.frame(sapply( k, evalmod))
			}
		}
	    return(outmin[ 2, which.min(outmin[1,])])
	    }
	    rngk <- 0.05*diff( c( min(indata$y, na.rm = TRUE), max(indata$y, na.rm = TRUE)) )
	    currenty <- richards(indata$x[1], val[indx,1], val[indx, 2], val[indx, 3], val[indx, 4]) 
	    try(
	    if( currenty < (indata$y[1] - rngk) | currenty > (indata$y[1] + rngk) ){
	    necessk <- findk( val[indx,1], val[indx,2], val[indx,3], val[indx,4], indata$y[1], indata$x[1] )
	    val[,2]<- necessk
	    	for (i in 1:nrow(val)) {
	          try({
	               vald <- evly(val, y1, y0, i)
	               y1 <- vald[1]
	               y0 <- vald[2]
	               indx <- vald[3]
	               }, silent = TRUE)
            	}
	    }
	    ,silent = TRUE)
	    if(savy0 < y0) {
	    indx <- savindx
	    val <- sav0val
	        if(sav1y0 < savy0) indx <- sav1indx
	        if(sav1y0 < savy0) val <- sav1val
	    } else {
	        if(sav1y0 < y0) indx <- sav1indx
	        if(sav1y0 < y0) val <- sav1val
	    }
            val <- c(val[indx,1], val[indx,2], 
            		 val[indx,3], val[indx,4])
            options(warn = 0)
            return(val)
        }
        val3 <- trnsfrm(xyE, min(xy$x), Intage, taper.ends = taper.ends)
        names(val3) <- c("Asym", "K", "Infl", "M")
        Asym <- val3[1]
        K <- val3[2]
        Infl <- val3[3]
        M <- val3[4]
        try(if(Asym > Amax) {
        if(pntbndserr == TRUE) stop("Upper bounds for Asym are too low: optimized Asym > Amax. If you wish bounds to be reestimated
        	automatically, specify bounds.error = FALSE in modpar") 
        cat("WARNING: Upper bounds for Asym are too low: optimized Asym > Amax. These will be adjusted. Alternatively
        	try runing modpar with option width.bounds set to more than 1 \n")
 	Amax <- as.numeric(Asym + (Asym * largemult))
      		}, silent = TRUE)
        try(if(Asym < Amin) {
        if(pntbndserr == TRUE) stop("Lower bounds for Asym are too high:optimized Asym < Amin. If you wish bounds to be reestimated
        	automatically, specify bounds.error = FALSE in modpar")         
        cat("WARNING: Lower bounds for Asym are too high: optimized Asym < Amin. These will be adjusted. Alternatively
        	try runing modpar with option width.bounds set to more than 1 \n")
        Amin <- as.numeric(Asym + (Asym * -largemult))
        		}, silent = TRUE)
        try(if(Infl > Imax) {
        if(pntbndserr == TRUE) stop("Upper bounds for Infl are too low: optimized Infl  > Imax. If you wish bounds to be reestimated
        	automatically, specify bounds.error = FALSE in modpar")         
        cat("WARNING: Upper bounds for Infl are too low: optimized Infl  > Imax. These will be adjusted. Alternatively
        	try runing modpar with option width.bounds set to more than 1 \n")
        Imax <- as.numeric(Infl + (Infl * largemult))
        		}, silent = TRUE)
        try(if(Infl < Imin) {
        if(pntbndserr == TRUE) stop("Lower bounds for Infl are too high: optimized Infl < Imin. If you wish bounds to be reestimated
         	automatically, specify bounds.error = FALSE in modpar")         
        cat("WARNING: Lower bounds for Infl are too high: optimized Infl < Imin. These will be adjusted. Alternatively
        	try runing modpar with option width.bounds set to more than 1 \n")
        Imin <- as.numeric(Infl + (Infl * -largemult))
        		}, silent = TRUE)
        if (modno >= 17 & modno < 18) {
            M <- ((Infl - mean(xyE$x)) * M) - (Infl - mean(xyE$x))
            names(M) <- "M"
        }
        val <- c(Asym, K, Infl, M)
        if (modno > 19 | (modno >= 17 & modno < 18)) {
            if (modno > 19) 
                val <- c(Asym, K, Infl)
            if (modno >= 17 & modno < 18) 
                val <- c(Asym, Infl, M)
            func21 <- function(val) {
                func23 <- function(x, K, Infl) {
                  x = x + (0 + (0+0i))
                  (1 + M * exp(-K * (x - Infl)))^(1/M)
                }
                func173 <- function(x, Infl, M) {
                  x = x + (0 + (0+0i))
                  (1 + exp(Infl - x)/M)
                }
                Asym = val[1]
                if (modno >= 17 & modno < 18) {
                  K = 1
                  Infl = val[2]
                  M = val[3]
                  names(M) <- "M"
                  P1 <- Re(func173(xyE$x, Infl, M))
                } else {
                  K = val[2]
                  Infl = val[3]
                  P1 <- Re(func23(xyE$x, K, Infl))
                }
                P1[is.na(P1)] <- 1e-290 * pnmodelparams$Asym
                P1[P1 == Inf] <- 1e-290 * pnmodelparams$Asym
                y1 <- Asym/P1
                evl <- sum((xyE$y - y1)^2)
   		if(!is.na(evl)){
                if (evl == Inf) {
                  evl <- 1e+290
                	} else {
                  evl <- 1e+290
                		}
               	}
                try(if (min(Im(as.complex(1 + M * exp(-K * ((0:max(xyE$x)) - 
                  Infl)))^(1/M))) < 0) {
                  evl <- 1e+200
                }, silent = T)
   		if(!is.na(evl)){
                if (evl == Inf) {
                  evl <- 1e+290
                	} else {
                  evl <- 1e+290
                		}
               	}
                return(evl)
            }
            oppar = 0
            is.na(oppar) <- TRUE
            func25 <- function(Infl) {
                repl <- Kmax
                if (M < 1e-20) {
                  op <- data.frame((max(xyE$x)) - Infl)
                  if (nrow(subset(op, op == 0)) > 0) 
                    repl <- (-(2.170611e-16 - 1e-06 + (log(abs(M)))/(Infl + 
                      1e-10)))
                }
                return(repl)
            }
            if (Infl >= Imax) Infl = Imax - ((Imax - min(xy$x, na.rm=TRUE)) * 0.95)
            if (Kmin < 1e-05 & Kmin >= 0 ) {
 	            Kmin = 1e-05
 	         } else {
            if (Kmin > (-1e-05) & Kmin < 0) Kmin = (-1e-05)
            				    }
            if (modno >= 17 & modno < 18) {
            	parscale1 <- sapply(c(Asym, Infl, M),roundparscale)
                try(oppar <- optim(c(Asym, Infl, M), func21, 
                  method = "L-BFGS-B", lower = c(Amin, Imin, 
                    Mmin), upper = c(Amax, Imax, Mmax), control = list(maxit = 2000, 
                    parscale = parscale1)), silent = TRUE)
            } else {
            	 parscale1 <- sapply(c(Asym, K, M),roundparscale)
                 try(oppar <- optim(c(Asym, K, Infl), func21, 
                  method = "L-BFGS-B", lower = c(Amin, Kmin, 
                    Imin), upper = c(Amax, func25(Infl), Imax), 
                  control = list(maxit = 2000, parscale = parscale1
                  )), silent = TRUE)
            }
            if (is.na(oppar[1]) == FALSE) {
                if (oppar$convergence < 52) {
                  val <- (c(oppar$par))
                } else {
                  if (!is.na(modelparams$twocomponent.x)) {
                    stop("ERROR: positive optimization failed,aborting")
                  } else {
                    print("Warning: positive optimization failed,using estimated parameters")
                    val <- c(Asym, K, Infl)
                    names(val)<-c("Asym", "K", "Infl")
                     if (modno >= 17 & modno < 18) {
                      	val <- c(Asym, Infl, M)
                      	names(val)<-c("Asym", "Infl", "M")
                     }
                  }
                }
            } else {
                if (!is.na(modelparams$twocomponent.x)) {
                  stop("ERROR: positive optimization failed,aborting")
                } else {
                  print("Warning: positive optimization failed,using estimated parameters")
   		  val <- c(Asym, K, Infl)
                  names(val)<-c("Asym", "K", "Infl")
                  if (modno >= 17 & modno < 18) {
                      	val <- c(Asym, Infl, M)
                      	names(val)<-c("Asym", "Infl", "M")
                  }
                }
            }
            if (modno >= 17 & modno < 18) {
                Asym = val[1]
                Infl = val[2]
                M = val[3]
                value <- c(Asym, Infl, M)
                names(value) <- c("Asym", "Infl", "M")
            } else {
                Asym = val[1]
                K = val[2]
                Infl = val[3]
                value <- c(Asym, K, Infl)
                names(value) <- c("Asym", "K", "Infl")
            }
            Kmax = K + (abs(K) * mediummult)
            Kmin = K - (abs(K) * mediummult)
 	    if (Kmin < 1e-05 & Kmin >= 0 ) {
 	            Kmin = 1e-05
 	         } else {
            if (Kmin > (-1e-05) & Kmin < 0) Kmin = (-1e-05)
            				    }
            Imax = Imax - (abs(Imax - Infl) * mediummult)
            Imin = Imin + (abs(Imin - Infl) * mediummult)
            if (modno >= 17 & modno < 18) {
            } else {
                while (abs(Imax * Kmax) > 700) Imax = Imax * 
                  0.9
                while (abs(Imin * Kmax) > 700) Imin = Imin * 
                  0.9
            }
            if( !is.na(Imax) ) {
       		if (Imax > maxIval) Imax <- maxIval
          	if(Imax < Infl) {
          	   maxIval <- Infl + ((maxIval - min(xy$x)) * 0.05)
          	   Imax <- maxIval  
          	   }
          	  		}
          	  if( !is.na(Imin) ) {
         	   if(Imin > Infl) Imin <- min(xy$x, na.rm=TRUE) 
         	   if(Imin > Infl) {
       		   Imin <- Infl - ((maxIval - min(xy$x)) * 0.05)
          	   			}
	    }
            Mmax = Mmax - (abs(Mmax - M) * largemult)
            Mmin = Mmin + (abs(Mmin - M) * largemult)
            pars <- 1
            is.na(pars) <- TRUE
            if (modno >= 17 & modno < 18) {
                try(pars <- as.vector(coef(nls(y ~ richards173(x, 
                  Asym, Infl, M), data = xyE, start = list(Asym = Asym, 
                  Infl = Infl, M = M), nls.control(maxiter = 1000, 
                  tol = 1), algorithm = "port", lower = list(Asym = Amin, 
                  Infl = Imin, M = Mmin), upper = list(Asym = Amax, 
                  I = Imax, M = Mmax)))), silent = TRUE)
            } else {
                try(pars <- as.vector(coef(nls(y ~ richards3(x, 
                  Asym, K, Infl), data = xyE, start = list(Asym = Asym, 
                  K = K, Infl = Infl), nls.control(maxiter = 1000, 
                  tol = 1), algorithm = "port", lower = list(Asym = Amin, 
                  K = Kmin, Infl = Imin), upper = list(Asym = Amax, 
                  K = Kmax, I = Imax)))), silent = TRUE)
            }
            if (is.na(pars[1]) == FALSE) {
                if (modno >= 17 & modno < 18) {
             if(!is.na(modelparams$twocomponent.x)) {
	                subxyp <- subset(xyE$x, !is.na(xyE$x) )
	                addRAsym <- richards(subxyp [length(subxyp )], pars[1L], pars[2L], pars[3L], pars[4L])
	                value <- c(pars[1L], K, pars[2L], pars[3L], 
	                    (Asym * 0.05)+addRAsym, K, pars[2L], tAsym)
	                } else {
	                value <- c(pars[1L], K, pars[2L], pars[3L], 
	                    (Asym * 0.05), K, pars[2L], tAsym)
            }           
                   names(value) <- c("Asym", "K", "Infl", "M", "RAsym", 
                  "Rk", "Ri", "RM")
                } else {
             if(!is.na(modelparams$twocomponent.x)) {
	                subxyp <- subset(xyE$x, !is.na(xyE$x) )
	                addRAsym <- richards(subxyp [length(subxyp )], pars[1L], pars[2L], pars[3L], pars[4L])
	                value <- c(pars[1L], pars[2L], pars[3L], M,
	                    (Asym * 0.05)+addRAsym, pars[2L], tAsym, M)
	                } else {
	                value <- c(pars[1L], pars[2L], pars[3L], M,
	                    (Asym * 0.05), pars[2L], tAsym, M)
            }   
                   names(value) <- c("Asym", "K", "Infl", "M", "RAsym", 
                  "Rk", "Ri", "RM")
                }
            } else {
                if (modelparams$verbose == TRUE) 
                  print("Warning: nls failed for positive Richards model, using optim parameters")
            if(!is.na(modelparams$twocomponent.x)) {
            subxyp <- subset(xyE$x, !is.na(xyE$x) )
            addRAsym <- richards(subxyp [length(subxyp )], Asym, K, Infl, M)
            value <- c(Asym, K, Infl, M, 
                (Asym * 0.05)+addRAsym, K, tAsym, M)
            } else {
            value <- c(Asym, K, Infl, M, 
                (Asym * 0.05), K, tAsym, M)
            }
                names(value) <- c("Asym", "K", "Infl", "M", "RAsym", 
                  "Rk", "Ri", "RM")
            }
        } else {
            func1 <- function(val) {
                Asym = val[1]
                K = val[2]
                Infl = val[3]
                M = val[4]
                func3 <- function(x, K, Infl, M) {
                  x = x + (0 + (0+0i))
                  (1 + M * exp(-K * (x - Infl)))^(1/M)
                }
                y1 <- Asym/Re(func3(xyE$x, K, Infl, M))
                evl <- sum((xyE$y - y1)^2)
                try(if (min(Im(as.complex(1 + M * exp(-K * ((max(xyE$x)) - 
                  Infl)))^(1/M))) < 0) {
                  evl <- 1e+200
                }, silent = T)
   		if(!is.na(evl)){
                if (evl == Inf) {
                  evl <- 1e+290
                	} else {
                  evl <- 1e+290
                		}
               	}
                return(evl)
            }
            oppar = 0
            is.na(oppar) <- TRUE
            func5 <- function(Infl, M) {
                repl <- Kmax
                if (M < 1e-20) {
                  op <- data.frame((max(xyE$x)) - Infl)
                  if (nrow(subset(op, op == 0)) > 0) 
                    repl <- (-(2.170611e-16 - 1e-06 + (log(abs(M)))/(Infl + 
                      1e-10)))
                }
                return(repl)
            }
            if (modno == 18 | modno == 19) {
            	Risav <- NA
                testpar <- get(optvarnm, envir = Envir )
                testpar <- testpar[ names(testpar) %in% names(pnmodelparamsbounds)]
                if (is.na(testpar[1])) {
                  Amax = Asym + (abs(Asym) * largemult)
                  Amin = Asym - (abs(Asym) * mediummult)
                  Kmax = K + (abs(K) * mediummult)
                  Kmin = K - (abs(K) * mediummult)
                  Imax = Infl + (abs(Infl) * largemult)
                  Imin = Infl + (abs(Infl) * -2*(largemult))
                  while (abs(Imax * Kmax) > 700) Imax = Imax * 
                    0.9
                  while (abs(Imin * Kmax) > 700) Imin = Imin * 
                    0.9
                  maxIval<- try(max(xyE$x, na.rm = TRUE) - 
                  	(max(diff( range(xyE$x, na.rm = TRUE))*.05)), silent=TRUE)
       		  if(class(maxIval) == "try-error") maxIval <- max(xy$x, na.rm=TRUE) 
            	 if( !is.na(Imax) ) {
        		if (Imax > maxIval) Imax <- maxIval
           	if(Imax < Infl) {
           	   maxIval <- Infl + ((maxIval - min(xy$x)) * 0.05)
           	   Imax <- maxIval  
           	   }
           	  		}
           	  if( !is.na(Imin) ) {
          	   if(Imin > Infl) Imin <- min(xy$x, na.rm=TRUE) 
          	   if(Imin > Infl) {
        		   Imin <- Infl - ((maxIval - min(xy$x)) * 0.05)
           	   			}
	    	}	
                  if (abs(M) < 0.1) {
                    Mmax = largemult
                    Mmin = -largemult
                  } else {
                    Mmax = M + abs(M)
                    Mmin = M - abs(M)
                  }
                  RAmax = Amax
                  RAmin = Amin
                  Rkmax = Kmax
                  Rkmin = Kmin
                  Rimax = Imax
                  Rimin = Imin
                  RMmax = Mmax
                  RMmin = Mmin
                  names(Amax) <- ("Amax")
 		  names(Amin) <- ("Amin")
                  names(Kmax) <- ("Kmax")
                  names(Kmin) <- ("Kmin")
 		  names(Imin) <- ("Imin")
                  names(Imax) <- ("Imax")
                  names(Mmax) <- ("Mmax")
 		  names(Mmin) <- ("Mmin")
                  names(RAmax) <- ("RAmax")
            	  names(RAmin) <- ("RAmin")
                  names(Rkmax) <- ("Rkmax")
  		  names(Rkmin) <- ("Rkmin")
 		  names(Rimin) <- ("Rimin")
                  names(Rimax) <- ("Rimax")
                  names(RMmax) <- ("RMmax")
 		  names(RMmin) <- ("RMmin")
                  maxRival<- try(max(xyL$x, na.rm = TRUE) - 
                  	(max(diff( range(xyL$x, na.rm = TRUE))*.05)), silent=TRUE)  
                  if(class(maxRival) == "try-error") maxRival <- max(xy$x, na.rm=TRUE)      
		  if( !is.na(Rimax) ) {if (Rimax > maxRival) Rimax <- maxRival}
                } else {
                  names(testpar) <- c("Amin", "Amax", "Kmin", "Kmax", "Imin",
       			 "Imax", "Mmin", "Mmax", "RAmin", "RAmax", "Rkmin", "Rkmax",
        		"Rimin", "Rimax", "RMmin", "RMmax")
                  Amax = testpar$Amax
                  Amin = testpar$Amin
                  Kmax = testpar$Kmax
                  Kmin = testpar$Kmin
                  Imax = testpar$Imax
                  Imin = testpar$Imin
                  Mmax = testpar$Mmax
                  Mmin = testpar$Mmin
                  RAmax = testpar$RAmax
                  RAmin = testpar$RAmin
                  Rkmax = testpar$Rkmax
                  Rkmin = testpar$Rkmin
                  Rimax = testpar$Rimax
                  Rimin = testpar$Rimin
                  RMmax = testpar$RMmax
                  RMmin = testpar$RMmin
                  maxIval<- try(max(xyE$x, na.rm = TRUE) - 
		       (max(diff( range(xyE$x, na.rm = TRUE))*.05)), silent=TRUE)
                  maxRival<- try(max(xyL$x, na.rm = TRUE) - 
		      (max(diff( range(xyL$x, na.rm = TRUE))*.05)), silent=TRUE)  
		  if(class(maxRival) == "try-error") maxRival <- max(xy$x, na.rm=TRUE)
             if( !is.na(Imax) ) {
       		if (Imax > maxIval) Imax <- maxIval
          	if(Imax < Infl) {
          	   maxIval <- Infl + ((maxIval - min(xy$x)) * 0.05)
          	   Imax <- maxIval  
          	   }
          	  		}
          	  if( !is.na(Imin) ) {
         	   if(Imin > Infl) Imin <- min(xy$x, na.rm=TRUE) 
         	   if(Imin > Infl) {
       		   Imin <- Infl - ((maxIval - min(xy$x)) * 0.05)
          	   			}
	    }	    
	     if( !is.na(Rimax) ) {
       		if (Rimax > maxRival) Rimax <- maxRival
          	if(Rimax < tAsym) {
          	   maxRival <- tAsym + ((maxRival - min(xy$x)) * 0.95)
          	   Rimax <- maxRival  
          	   }
          	  		}
          	  if( !is.na(Rimin) ) {
         	   if(Rimin > tAsym) Rimin <- min(xy$x, na.rm=TRUE) 
         	   if(Rimin > tAsym) {
       		   Rimin <- tAsym - ((maxIval - min(xy$x)) * 0.95)
          	   			}
             }	   			
	    }	    
	           if(is.na(Risav)) Risav <- tAsym
                skel <- rep(list(1), 16)
                if(!is.na(modelparams$twocomponent.x)) {
 		            subxyp <- subset(xyE$x, !is.na(xyE$x) )
		            addRAsym <- richards(subxyp [length(subxyp )], Asym, K, Infl, M)
		            } else {
		            addRAsym <- 0
            	}
                exportparams <- c(Asym, K, Infl, M, (Asym * 0.05) + addRAsym, 
                  K, Risav, M, modelparams$first.y, modelparams$x.at.first.y,
                  modelparams$last.y, 
                  modelparams$x.at.last.y, modelparams$twocomponent.x, 
                  modelparams$verbose, modelparams$force4par, modelparams$force.nonmonotonic)
                exportparams <- relist(exportparams, skel)
                names(exportparams) <- c("Asym", "K", "Infl", 
                  "M", "RAsym", "Rk", "Ri", "RM", "first.y", "x.at.first.y",
                  "last.y", "x.at.last.y", "twocomponent.x", 
                  "verbose", "force4par", "force.nonmonotonic")
                exportparams$twocomponent.x <- modelparams$twocomponent.x
                exportparams$verbose <- modelparams$verbose
                exportparams$force4par <- modelparams$force4par
                exportparams$force.nonmonotonic <- modelparams$force.nonmonotonic
                skel <- rep(list(1), 16)
                exportparamsbounds <- c(Amin, Amax, Kmin, Kmax, 
                  Imin, Imax, Mmin, Mmax, RAmin, RAmax, Rkmin, 
                  Rkmax, Rimin, Rimax, RMmin, RMmax)
                exportparamsbounds <- relist(exportparamsbounds, 
                  skel)
                names(exportparamsbounds) <- c("Amin", "Amax", 
                  "Kmin", "Kmax", "Imin", "Imax", "Mmin", "Mmax", 
                  "RAmin", "RAmax", "Rkmin", "Rkmax", "Rimin", 
                  "Rimax", "RMmin", "RMmax")
    		valexp <- sapply( c(unlist(exportparams),unlist(exportparamsbounds), 
            		unlist( optvar[names(optvar) %w/o% c(names(exportparams),names(exportparamsbounds)) ]) ), function(x) list(x))
       		names( valexp[1:32] ) <- names( c(exportparams[1:16],exportparamsbounds) )
        	valexp[c(1:12,17:32)] <- as.numeric( valexp[c(1:12,17:32)] )
        	valexp[14:16] <- as.logical( valexp[14:16] )
        	if(is.numeric(modelparams$twocomponent.x)) {
        	valexp[13] <- as.numeric( valexp[13] )
        	} else {
  		valexp[13] <- as.logical( valexp[13] )     	
        	}
  		try(valexp$taper.ends <- as.numeric(valexp$taper.ends), silent = TRUE)
 		assign(optvarnm, valexp, envir = Envir )
 		pnmodelparams<-valexp[1:16]
 	      }
            if (Infl >= Imax) Infl = Imax - ((Imax - min(xy$x, na.rm=TRUE)) * 0.95)
            if (Kmin < 1e-05 & Kmin >= 0 ) {
 	            Kmin = 1e-05
 	         } else {
            if (Kmin > (-1e-05) & Kmin < 0) Kmin = (-1e-05)
            				    }
            savminx <- min(xyE$x)
            savmaxx <- max(xyE$x)
            parscale1 <- sapply(c(Asym, K, Infl, M),roundparscale)
            try(oppar <- optim(c(Asym, K, Infl, M), func1, method = "L-BFGS-B", 
                lower = c(Amin, Kmin, Imin, Mmin), upper = c(Amax, 
                  func5(Infl, M), Imax, Mmax), control = list(maxit = 2000, 
                  parscale = parscale1)), silent = TRUE)
            if (is.na(oppar[1]) == TRUE) {
                xyE <- moddata(xyE, min(xy$x, Intage, taper.ends = taper.ends))[["indata"]]
                try(oppar <- optim(c(Asym, K, Infl, M), func1, 
                  method = "L-BFGS-B", lower = c(Amin, Kmin, 
                    Imin, Mmin), upper = c(Amax, func5(Infl, 
                    M), Imax, Mmax), control = list(maxit = 2000, 
                    parscale = parscale1)), silent = TRUE)     
            } else {
                if (oppar$convergence >= 52) {
                  xyE <- moddata(xyE, min(xy$x, Intage, taper.ends = taper.ends))[["indata"]]
                  try(oppar <- optim(c(Asym, K, Infl, M), func1, 
                    method = "L-BFGS-B", lower = c(Amin, Kmin, 
                      Imin, Mmin), upper = c(Amax, func5(Infl, 
                      M), Imax, Mmax), control = list(maxit = 2000, 
                      parscale = parscale1)), silent = TRUE)     
                } else {
                  evy1 <- sum((xyE$y - richards(xyE$x, as.numeric(oppar$par[1]), 
                    as.numeric(oppar$par[2]), as.numeric(oppar$par[3]), 
                    as.numeric(oppar$par[4])))^2)
                  evy2 <- sum((xyE$y - richards(xyE$x, Asym, 
                    K, Infl, M))^2)
                  if(is.na(evy1)) evy1 <- Inf
                  if(is.na(evy2)) evy2 <- Inf
                  if (evy1 > evy2 | !is.na(modelparams$twocomponent.x)) {
                    if (!is.na(modelparams$twocomponent.x)) 
                      savoppE <- oppar
                    xyE <- moddata(xyE, min(xy$x, Intage, taper.ends = taper.ends))[["indata"]]
                    try(oppar <- optim(c(Asym, K, Infl, M), func1, 
                      method = "L-BFGS-B", lower = c(Amin, Kmin, 
                        Imin, Mmin), upper = c(Amax, func5(Infl, 
                        M), Imax, Mmax), control = list(maxit = 2000, 
                        parscale = parscale1)), 
                      silent = TRUE)     
                    if (!is.na(modelparams$twocomponent.x)) {
                      if (class(oppar)[1] == "try-error") {
                        oppar <- savoppE
                      } else {
                        if (oppar$value > savoppE$value) 
                          oppar <- savoppE
                      }
                    }
                  }
                }
            }
            if (min(xyE$x) != savminx | max(xyE$x) != savmaxx) 
                xyE <- subset(xyE, xyE$x >= savminx & xyE$x <= 
                  savmaxx)
            xyE <- xyE[order(xyE$x), ]
            if (is.na(oppar[1]) == FALSE) {
                if (oppar$convergence < 52) {
                  val <- (c(oppar$par))
                } else {
                  if (!is.na(modelparams$twocomponent.x)) {
                    stop("ERROR: positive optimization failed,aborting")
                  } else {
                    print("Warning: positive optimization failed,using estimated parameters")
    		    val <- c(Asym, K, Infl, M)
                    names(val)<-c("Asym", "K", "Infl", "M")
                  }
                }
            } else {
                if (!is.na(modelparams$twocomponent.x)) {
                  stop("ERROR: positive optimization failed,aborting")
                } else {
                  print("Warning: positive optimization failed,using estimated parameters")
                  val <- c(Asym, K, Infl, M)
                  names(val)<-c("Asym", "K", "Infl", "M")
                }
            }
            Asym = val[1]
            K = val[2]
            Infl = val[3]
            M = val[4]
            pars = 0
            is.na(pars) <- TRUE
            Amax = Asym + diff(range(xy$y, na.rm = TRUE) * smallmult)
            Amin = Asym - diff(range(xy$y, na.rm = TRUE) * smallmult)
            Kmax = K + (abs(K) * smallmult)
            Kmin = K - (abs(K) * smallmult)
 	    if (Kmin < 1e-05 & Kmin >= 0 ) {
 	            Kmin = 1e-05
 	         } else {
            if (Kmin > (-1e-05) & Kmin < 0) Kmin = (-1e-05)
            				    }
            Imax = Infl + diff(range(xy$x, na.rm = TRUE)  * smallmult)
            Imin = Infl + diff(range(xy$x, na.rm = TRUE)  * -smallmult)
            while (abs(Imax * Kmax) > 700) Imax = Imax * 0.9
            while (abs(Imin * Kmax) > 700) Imin = Imin * 0.9
            	 if( !is.na(Imax) ) {
        		if (Imax > maxIval) Imax <- maxIval
           	if(Imax < Infl) {
           	   maxIval <- Infl + ((maxIval - min(xy$x)) * 0.05)
           	   Imax <- maxIval  
           	   }
           	  		}
           	  if( !is.na(Imin) ) {
          	   if(Imin > Infl) Imin <- min(xy$x, na.rm=TRUE) 
          	   if(Imin > Infl) {
        		   Imin <- Infl - ((maxIval - min(xy$x)) * 0.05)
           	   			}
	    	}	
              if (abs(M) < 0.1) {
                Mmax = mediummult
                Mmin = -mediummult
            } else {
                Mmax = M + abs(M * largemult)
                Mmin = M - abs(M * largemult)
            }
            try(pars <- as.vector(coef(nls(y ~ richards(x, Asym, 
                K, Infl, M), data = xyE, start = list(Asym = Asym, 
                K = K, Infl = Infl, M = M), nls.control(maxiter = 1000, 
                tol = 1), algorithm = "port", lower = list(Asym = Amin, 
                K = Kmin, Infl = Imin, M = Mmin), upper = list(Asym = Amax, 
                K = Kmax, I = Imax, M = Mmax)))), silent = TRUE)
            if (is.na(pars[1]) == FALSE) {
            if(!is.na(modelparams$twocomponent.x)) {
	                subxyp <- subset(xyE$x, !is.na(xyE$x) )
	                addRAsym <- richards(subxyp [length(subxyp )], pars[1L], pars[2L], pars[3L], pars[4L])
	                value <- c(pars[1L], pars[2L], pars[3L], pars[4L], 
	                    (Asym * 0.05)+addRAsym, pars[2L], tAsym, pars[4L])
	                } else {
	                value <- c(pars[1L], pars[2L], pars[3L], pars[4L], 
	                    (Asym * 0.05), pars[2L], tAsym, pars[4L])
            }
                names(value) <- c("Asym", "K", "Infl", "M", "RAsym", 
                  "Rk", "Ri", "RM")
            } else {
                if (modelparams$verbose == TRUE) 
                  print("Warning: nls failed for positive Richards model, using optim parameters (if no optimization occured, using estimated parameters)")
 	    if(!is.na(modelparams$twocomponent.x)) {
            subxyp <- subset(xyE$x, !is.na(xyE$x) )
            addRAsym <- richards(subxyp [length(subxyp )], Asym, K, Infl, M)
            value <- c(Asym, K, Infl, M, 
                (Asym * 0.05)+addRAsym, K, tAsym, M)
            } else {
            value <- c(Asym, K, Infl, M, 
                (Asym * 0.05), K, tAsym, M)
            }
                names(value) <- c("Asym", "K", "Infl", "M", "RAsym", 
                  "Rk", "Ri", "RM")
            }
        }
    }
    if (modelparams$force4par == FALSE) {
        xyL <- data.frame(rep(NA, 1))
        xyL <- subset(xy, xy$x >= Intage)
                maxRival<- try(max(xyL$x, na.rm = TRUE) - (max(diff( range(xyL$x, na.rm = TRUE))*.05)), silent=TRUE)
  		if(class(maxRival) == "try-error") maxRival = max(xy$x, na.rm=TRUE)
        if (is.na(modelparams$twocomponent.x) == FALSE) 
            xyL <- subset(xy, xy$x > Intage)
        if (nrow(xyL) < 3 | modno == 12 | modno == 20 | modno == 
            32 | modno == 19 | modelparams$force4par == TRUE) {
            if (nrow(xyL) < 3) {
                if (modelparams$verbose == TRUE) 
                  print("Warning: too few distinct input values to fit the negative Richards model, defaulting to RAsym= 5% of peak mass, Rk=k, Ri=age at peak, RM=M")
            counterR <- 0
            if(!is.na(modelparams$twocomponent.x)) {
            subxyp <- subset(xyE$x, !is.na(xyE$x) )
            addRAsym <- richards(subxyp [length(subxyp )], value[1L], value[2L], value[3L], value[4L])
            value <- c(value[1L], value[2L], value[3L], value[4L], 
                (Asym * 0.05)+addRAsym, value[2L], tAsym, value[4L])
            } else {
            value <- c(value[1L], value[2L], value[3L], value[4L], 
                (Asym * 0.05), value[2L], tAsym, value[4L])
            }
            names(value) <- c("Asym", "K", "Infl", "M", "RAsym", 
                "Rk", "Ri", "RM")
            }
            inputval <- 0
            is.na(inputval) <- TRUE
            if (modno == 1 | modno == 21 | modno == 18) {
                val1 <- data.frame(value[5L], value[6L], value[7L], 
                  value[8L])
                names(val1) <- c("RAsym", "Rk", "Ri", "RM")
                inputmin <- c(val1$RAmin, val1$Rkmin, val1$Rimin, 
                  val1$RMmin)
                inputmax <- c(val1$RAmax, val1$Rkmax, val1$Rimax, 
                  val1$RMmax)
            }
            if (modno == 2 | modno == 22) {
                val1 <- data.frame(value[5L], value[6L], value[7L])
                names(val1) <- c("RAsym", "Rk", "Ri")
                inputmin <- c(val1$RAmin, val1$Rkmin, val1$Rimin)
                inputmax <- c(val1$RAmax, val1$Rkmax, val1$Rimax)
            }
            if (modno == 3 | modno == 23 | modno == 17.1) {
                val1 <- data.frame(value[7L], value[8L])
                names(val1) <- c("Ri", "RM")
                inputmin <- c(val1$Rimin, val1$RMmin)
                inputmax <- c(val1$Rimax, val1$RMmax)
            }
            if (modno == 4 | modno == 24) {
                val1 <- data.frame(value[5L], value[8L])
                names(val1) <- c("RAsym", "RM")
                inputmin <- c(val1$RAmin, val1$RMmin)
                inputmax <- c(val1$RAmax, val1$RMmax)
            }
            if (modno == 5 | modno == 25) {
                val1 <- data.frame(value[8L])
                names(val1) <- c("RM")
                inputmin <- c(val1$RMmin)
                inputmax <- c(val1$RMmax)
            }
            if (modno == 6 | modno == 26 | modno == 17) {
                val1 <- data.frame(value[5L], value[7L], value[8L])
                names(val1) <- c("RAsym", "Ri", "RM")
                inputmin <- c(val1$RAmin, val1$Rimin, val1$RMmin)
                inputmax <- c(val1$RAmax, val1$Rimax, val1$RMmax)               
            }
            if (modno == 7 | modno == 27) {
                val1 <- data.frame(value[6L], value[7L], value[8L])
                names(val1) <- c("Rk", "Ri", "RM")
                inputmin <- c(val1$Rkmin, val1$Rimin, val1$RMmin)
                inputmax <- c(val1$Rkmin, val1$Rimax, val1$RMmax)
            }
            if (modno == 8 | modno == 28) {
                val1 <- data.frame(value[5L], value[6L], value[8L])
                names(val1) <- c("RAsym", "Rk", "RM")
                inputmin <- c(val1$RAmin, val1$Rkmin, val1$RMmin)
                inputmax <- c(val1$RAmax, val1$Rkmax, val1$RMmax)
            }
            if (modno == 9 | modno == 29) {
                val1 <- data.frame(value[6L], value[8L])
                names(val1) <- c("Rk", "RM")
                inputmin <- c(val1$Rkmin, val1$RMmin)
                inputmax <- c(val1$Rkmax, val1$RMmax)
            }
            if (modno == 10 | modno == 30 | modno == 17.3) {
                val1 <- data.frame(value[7L])
                names(val1) <- c("Ri")
                inputmin <- c(val1$Rimin)
                inputmax <- c(val1$Rimax)
            }
            if (modno == 11 | modno == 31) {
                val1 <- data.frame(value[5L])
                names(val1) <- c("RAsym")
                inputmin <- c(val1$RAmin)
                inputmax <- c(val1$RAmax)
            }
            if (modno == 12 | modno == 32 | modno == 19) {
            }
            if (modno == 13 | modno == 33 | modno == 17.2) {
                val1 <- data.frame(value[5L], value[7L])
                names(val1) <- c("RAsym", "Ri")
                inputmin <- c(val1$RAmin, val1$Rimin)
                inputmax <- c(val1$RAmax, val1$Rimax)
            }
            if (modno == 14 | modno == 34) {
                val1 <- data.frame(value[6L], value[7L])
                names(val1) <- c("Rk", "Ri")
                inputmin <- c(val1$Rkmin, val1$Rimin)
                inputmax <- c(val1$Rkmax, val1$Rimax)
            }
            if (modno == 15 | modno == 35) {
                val1 <- data.frame(value[5L], value[6L])
                names(val1) <- c("RAsym", "Rk")
                inputmin <- c(val1$RAmin, val1$Rkmin)
                inputmax <- c(val1$RAmax, val1$Rkmax)
            }
            if (modno == 16 | modno == 36) {
                val1 <- data.frame(value[6L])
                names(val1) <- c("Rk")
                inputmin <- c(val1$Rkmin)
                inputmax <- c(val1$Rkmax)
            }
            if (modno == 20) {
            }
        } else {
            if (is.na(modelparams$twocomponent.x) == TRUE) {
                signval <- 1
                if (max(abs(xy$y)) != max(xy$y)) 
                  signvalmax <- -1
                xyL$y <- xyL$y - (richards(xyL$x, value[1], value[2], 
                  value[3], value[4]) - (max(abs(xy$y)) * signval))
                xyL$y <- xyL$y - (min(abs(xyL$y)) * sign(min(xyL$y)))
                xyL <- xyL[order(xyL$x), ]
                signvalmin <- 1
                signvalmax <- 1
                if (max(abs(xyL$y)) != max(xyL$y)) 
                  signvalmax <- -1
                if (min(abs(xyL$y)) != min(xyL$y)) 
                  signvalmin <- -1
                if (xyL$y[nrow(xyL)] < ((diff(range(xyL$y))/2) + 
                  min(xyL$y))) {
                  xyL$y <- xyL$y - (max(abs(xyL$y)) * signvalmax)
                } else {
                  xyL$y <- xyL$y - (min(abs(xyL$y)) * signvalmin)
                }
            }
            val2 <- trnsfrm(xyL, min(xy$x), Intage, taper.ends = taper.ends)
            names(val2) <- c("RAsym", "Rk", "Ri", "RM")
            RAsym <- val2[1]
            Rk <- val2[2]
            Ri <- val2[3]
            RM <- val2[4]
            counterR <- 0
            try(if(RAsym > RAmax & modno != 18 & modno != 19 ) {
            if(pntbndserr == TRUE) stop("Upper bounds for RAsym are too low: optimized RAsym > RAmax. If you wish bounds to be reestimated
          	automatically, specify bounds.error = FALSE in modpar")         
            cat("WARNING: Upper bounds for RAsym are too low: optimized RAsym > RAmax. These will be adjusted. Alternatively
            	try runing modpar with option width.bounds set to more than 1 \n")
            RAmax <- as.numeric(RAsym + (RAsym * largemult))
            		}, silent = TRUE)
            try(if(RAsym < RAmin & modno != 18 & modno != 19) {
            if(pntbndserr == TRUE) stop("Lower bounds for RAsym are too high: optimized RAsym < RAmin. If you wish bounds to be reestimated
          	automatically, specify bounds.error = FALSE in modpar")         
            cat("WARNING: Lower bounds for RAsym are too high: optimized RAsym < RAmin. These will be adjusted. Alternatively
            	try runing modpar with option width.bounds set to more than 1 \n")
            RAmin <- as.numeric(RAsym + (RAsym * -largemult))
            		}, silent = TRUE)
            try(if(Ri > Rimax & modno != 18 & modno != 19) {
            if(pntbndserr == TRUE) stop("Upper bounds for RAsym are too low: optimized Ri  > Rimax. If you wish bounds to be reestimated
           	automatically, specify bounds.error = FALSE in modpar")         
            cat("WARNING: Upper bounds for Ri are too low: optimized Ri  > Rimax. These will be adjusted. Alternatively
            	try runing modpar with option width.bounds set to more than 1 \n")
            Rimax <- as.numeric(Ri + (Ri * largemult))
            		}, silent = TRUE)
            try(if(Ri < Rimin & modno != 18 & modno != 19) {
            if(pntbndserr == TRUE) stop("Lower bounds for RAsym are too high: optimized Ri < Rimin. If you wish bounds to be reestimated
           	automatically, specify bounds.error = FALSE in modpar")         
            cat("WARNING: Lower bounds for Ri are too high: optimized Ri < Rimin. These will be adjusted. Alternatively
            	try runing modpar with option width.bounds set to more than 1 \n")
            Rimin <- as.numeric(Ri + (Ri * -largemult))
        		}, silent = TRUE)
            if (modno == 18 | modno == 19) {
              if(modno == 18) {
                RAmax = RAsym + diff(range(xy$y, na.rm = TRUE) * largemult)
                RAmin = RAsym - diff(range(xy$y, na.rm = TRUE) * largemult)
                Rimax = Ri + diff(range(xy$x, na.rm = TRUE)  * largemult*2)
                Rimin = Ri + diff(range(xy$x, na.rm = TRUE)  * -largemult)
                while (abs(Rimax * Rkmax) > 700) Rimax = Rimax * 
                  0.9
                while (abs(Rimin * Rkmax) > 700) Rimin = Rimin * 
                  0.9
                if( !is.na(Rimax) ) {
 		        if (Rimax > maxRival) Rimax <- maxRival
 		        Risav <- tAsym
 		        if(Rimax < tAsym) Risav <- maxRival - ((Rimax - min(xy$x, na.rm=TRUE)) * 0.95)
 		        if(Rimin > tAsym) Risav <- min(xy$x, na.rm=TRUE) + ((Rimax - min(xy$x, na.rm=TRUE)) * 0.95)
       	     				  }
                if (abs(RM) < 0.1) {
                  RMmax = mediummult
                  RMmin = -mediummult
                } else {
                  RMmax = RM + abs(RM*mediummult)
                  RMmin = RM - abs(RM*mediummult)
                }
               } else {
               RAmax = Amax
	       RAmin = Amin
	       Rkmax = Kmax
	       Rkmin = Kmin
	       Rimax = Imax
               Rimin = Imin
               RMmax = Mmax
               RMmin = Mmin
               if( !is.na(Rimax) ) {if (Rimax > maxRival) Rimax <- maxRival}
               }
                skel <- rep(list(1), 16)
                exportparams <- c(Asym, K, Infl, M, RAsym, Rk, 
                  Ri, RM, modelparams$first.y, modelparams$x.at.first.y, 
                  modelparams$last.y, 
                  modelparams$x.at.last.y, modelparams$twocomponent.x, 
                  modelparams$verbose, modelparams$force4par,modelparams$force.nonmonotonic)
                exportparams <- relist(exportparams, skel)
                names(exportparams) <- c("Asym", "K", "Infl", 
                  "M", "RAsym", "Rk", "Ri", "RM", "first.y", "x.at.first.y",
                  "last.y", "x.at.last.y", "twocomponent.x", 
                  "verbose", "force4par", "force.nonmonotonic")
                exportparams$twocomponent.x <- modelparams$twocomponent.x
                exportparams$verbose <- modelparams$verbose
                exportparams$force4par <- modelparams$force4par
                exportparams$force.nonmonotonic <- modelparams$force.nonmonotonic
                skel <- rep(list(1), 16)
                exportparamsbounds <- c(Amin, Amax, Kmin, Kmax, 
                  Imin, Imax, Mmin, Mmax, RAmin, RAmax, Rkmin, 
                  Rkmax, Rimin, Rimax, RMmin, RMmax)
                exportparamsbounds <- relist(exportparamsbounds, 
                  skel)
                names(exportparamsbounds) <- c("Amin", "Amax", 
                  "Kmin", "Kmax", "Imin", "Imax", "Mmin", "Mmax", 
                  "RAmin", "RAmax", "Rkmin", "Rkmax", "Rimin", 
                  "Rimax", "RMmin", "RMmax")
     		valexp <- sapply( c(unlist(exportparams),unlist(exportparamsbounds), 
             		unlist( optvar[names(optvar) %w/o% c(names(exportparams),names(exportparamsbounds)) ]) ), function(x) list(x))
        		names( valexp[1:32] ) <- names( c(exportparams[1:16],exportparamsbounds) )
        	valexp[c(1:12,17:32)] <- as.numeric( valexp[c(1:12,17:32)] )
        	valexp[14:16] <- as.logical( valexp[14:16] )
        	if(is.numeric(modelparams$twocomponent.x)) {
        	valexp[13] <- as.numeric( valexp[13] )
        	} else {
  		valexp[13] <- as.logical( valexp[13] )     	
        	}
      		try(valexp$taper.ends <- as.numeric(valexp$taper.ends), silent = TRUE)
    		assign(optvarnm, valexp, envir = Envir )
    		pnmodelparams<-valexp[1:16]
		modelparams <- exportparams
                modelparamsbounds <- exportparamsbounds
                FPCEnv$.tmpposnegfile <- 0
            }
            if (modno == 1 | modno == 21 | modno == 18) {
                inputval <- c(RAsym, Rk, Ri, RM)
                inputmin <- c(RAmin, Rkmin, Rimin, RMmin)
                inputmax <- c(RAmax, Rkmax, Rimax, RMmax)
            }
            if (modno == 2 | modno == 22) {
                RM <- modelparams$RM
                inputval <- c(RAsym, Rk, Ri)
                inputmin <- c(RAmin, Rkmin, Rimin)
                inputmax <- c(RAmax, Rkmax, Rimax)
            }
            if (modno == 3 | modno == 23 | modno == 17.1) {
                Rk <- modelparams$Rk
                RAsym <- modelparams$RAsym
                inputval <- c(Ri, RM)
                inputmin <- c(Rimin, RMmin)
                inputmax <- c(Rimax, RMmax)
            }
            if (modno == 4 | modno == 24) {
                Rk <- modelparams$Rk
                Ri <- modelparams$Ri
                inputval <- c(RAsym, RM)
                inputmin <- c(RAmin, RMmin)
                inputmax <- c(RAmax, RMmax)
            }
            if (modno == 5 | modno == 25) {
                RAsym <- modelparams$RAsym
                Rk <- modelparams$Rk
                Ri <- modelparams$Ri
                inputval <- c(RM)
                inputmin <- c(RMmin)
                inputmax <- c(RMmax)
            }
            if (modno == 6 | modno == 26 | modno == 17) {
                Rk <- modelparams$Rk
                inputval <- c(RAsym, Ri, RM)
                inputmin <- c(RAmin, Rimin, RMmin)
                inputmax <- c(RAmax, Rimax, RMmax)
            }
            if (modno == 7 | modno == 27) {
                RAsym <- modelparams$RAsym
                inputval <- c(Rk, Ri, RM)
                inputmin <- c(Rkmin, Rimin, RMmin)
                inputmax <- c(Rkmax, Rimax, RMmax)
            }
            if (modno == 8 | modno == 28) {
                Ri <- modelparams$Ri
                inputval <- c(RAsym, Rk, RM)
                inputmin <- c(RAmin, Rkmin, RMmin)
                inputmax <- c(RAmax, Rkmax, RMmax)
            }
            if (modno == 9 | modno == 29) {
                RAsym <- modelparams$RAsym
                Ri <- modelparams$Ri
                inputval <- c(Rk, RM)
                inputmin <- c(Rkmin, RMmin)
                inputmax <- c(Rkmax, RMmax)
            }
            if (modno == 10 | modno == 30 | modno == 17.3) {
                Rk <- modelparams$Rk
                RAsym <- modelparams$RAsym
                RM <- modelparams$RM
                inputval <- c(Ri)
                inputmin <- c(Rimin)
                inputmax <- c(Rimax)
            }
            if (modno == 11 | modno == 31) {
                Rk <- modelparams$Rk
                Ri <- modelparams$Ri
                RM <- modelparams$RM
                inputval <- c(RAsym)
                inputmin <- c(RAmin)
                inputmax <- c(RAmax)
            }
            if (modno == 12 | modno == 32 | modno == 19) {
                print("error in model 12/32 code - please report")
            }
            if (modno == 13 | modno == 33 | modno == 17.2) {
                Rk <- modelparams$Rk
                RM <- modelparams$RM
                inputval <- c(RAsym, Ri)
                inputmin <- c(RAmin, Rimin)
                inputmax <- c(RAmax, Rimax)
            }
            if (modno == 14 | modno == 34) {
                RAsym <- modelparams$RAsym
                RM <- modelparams$RM
                inputval <- c(Rk, Ri)
                inputmin <- c(Rkmin, Rimin)
                inputmax <- c(Rkmax, Rimax)
            }
            if (modno == 15 | modno == 35) {
                Ri <- modelparams$Ri
                RM <- modelparams$RM
                inputval <- c(RAsym, Rk)
                inputmin <- c(RAmin, Rkmin)
                inputmax <- c(RAmax, Rkmax)
            }
            if (modno == 16 | modno == 36) {
                RAsym <- modelparams$RAsym
                Ri <- modelparams$Ri
                RM <- modelparams$RM
                inputval <- c(Rk)
                inputmin <- c(Rkmin)
                inputmax <- c(Rkmax)
            }
            if (modno == 20) {
                print("error in model 20 code - please report")
            }
            func2 <- function(val1) {
                tmpposnegfile <- FPCEnv$.tmpposnegfile
   		parload <- get(optvarnm, envir = Envir )
                pnmodelparams <- parload[ names(parload) %in% names(pnmodelparams)]
                value1 <- data.frame(RAsym = pnmodelparams$RAsym, 
                  Rk = pnmodelparams$Rk, Ri = pnmodelparams$Ri, 
                  RM = pnmodelparams$RM)
                if (tmpposnegfile == 0) {
                  value1$RAsym <- value1$RAsym * -1
                  value1$Ri <- max(xy$x) - value1$Ri
                }
                val1 <- (data.frame(t(val1)))
                if (length(val1$RAsym) == 1) 
                  value1$RAsym <- val1$RAsym
                if (length(val1$Rk) == 1) 
                  value1$Rk <- val1$Rk
                if (length(val1$Ri) == 1) 
                  value1$Ri <- val1$Ri
                if (length(val1$RM) == 1)
                  value1$RM <- val1$RM
                if (tmpposnegfile == 0 & length(val1$RM) == 
                  1) 
                  value1$RM <- 1
                RAsym <- value1$RAsym
                Rk <- value1$Rk
                Ri <- value1$Ri
                RM <- value1$RM
                if (modno >= 17 & modno < 18) {
                	if(modno == 17.3 | modno == 17.4) RM <- M
                  func4 <- function(x, Ri, RM) {
                    x = x + (0 + (0+0i))
                    (1 + exp(Ri - x)/RM)
                  }
                  P <- Re(func4(xyL$x, Ri, RM))
                } else {
                  func41 <- function(x, Rk, Ri, RM) {
                    x = x + (0 + (0+0i))
                    (1 + RM * exp(-Rk * (x - Ri)))^(1/RM)
                  }
                  P <- Re(func41(xyL$x, Rk, Ri, RM))
                }
                P[is.na(P)] <- 1e-290 * pnmodelparams$RAsym
                P[P == Inf] <- 1e-290 * pnmodelparams$RAsym
                if(modno == 17.2 | modno ==17.4) {
                y1 <- Asym/P
                }else{
                y1 <- RAsym/P
                }
                evl <- sum((xyL$y - y1)^2)
   		if(!is.na(evl)){
                if (evl == Inf) {
                  evl <- 1e+290
                	} else {
                  evl <- 1e+290
                		}
               	}
                options(warn = -1)
                try(if (min(Im(as.complex(1 + RM * exp(-Rk * 
                  ((max(xyL$x)) - Ri)))^(1/RM))) < 0) {
                  evl <- 1e+290
                }, silent = TRUE)
   		if(!is.na(evl)){
                if (evl == Inf) {
                  evl <- 1e+290
                	} else {
                  evl <- 1e+290
                		}
               	}
                options(warn = 0)
                return(evl)
            }
            if( is.na(inputval[1])) {
            	inputval<-val1
            	} else {
            	chinputval <- data.frame(t(inputval))
            	chinputmax <- data.frame(t(inputmax))
      	  	if(length(chinputval) == 1) chinputval <- data.frame(inputval)
          	if(length(chinputmax) == 1) chinputmax <- data.frame(inputmax)
            	names(chinputval) <- names(inputval)
            	chmx<-paste(substr(names(inputval),1,2),"max","")
            	chmx<-sub(" ","",chmx)
            	chmx<-sub(" ","",chmx)
            	names(chinputmax) <- chmx
            	if (length(chinputval$Ri) == 1) {
            		if (chinputval$Ri >= chinputmax$Rimax) 
            			chinputval$Ri <- chinputmax$Rimax - ((chinputmax$Rimax - min(xy$x, na.rm=TRUE)) * 0.95)
           	 }
            }
            oppar = 0
            is.na(oppar) <- TRUE
            savminx <- min(xyL$x)
            savmaxx <- max(xyL$x)
            xyL <- moddata(xyL, min(xy$x, Intage, taper.ends = taper.ends))[["indata"]]
            try(oppar <- optim(inputval, func2, method = "L-BFGS-B", 
                lower = inputmin, upper = inputmax, control = list(maxit = 2000)), 
                silent = TRUE)
            try(evy1L <- sum((xyL$y - richards(xyL$x, as.numeric(oppar$par[1]), 
                as.numeric(oppar$par[2]), as.numeric(oppar$par[3]), 
                as.numeric(oppar$par[4])))^2), silent = TRUE)
            try(evy2L <- sum((xyL$y - richards(xyL$x, inputval[1], 
                inputval[2], inputval[3], inputval[4]))^2), silent = TRUE)
            if (is.na(oppar[1]) == TRUE) {
                xyL <- moddata(xyL, min(xy$x, Intage, taper.ends = taper.ends))[["indata"]]
                try(oppar <- optim(inputval, func2, method = "L-BFGS-B", 
                  lower = inputmin, upper = inputmax, control = list(maxit = 2000)), 
                  silent = TRUE)
            } else {
                if (oppar$convergence >= 52) {
                  xyL <- moddata(xyL, min(xy$x, Intage, taper.ends = taper.ends))[["indata"]]
                  try(oppar <- optim(inputval, func2, method = "L-BFGS-B", 
                    lower = inputmin, upper = inputmax, control = list(maxit = 2000)), 
                    silent = TRUE)
                } else {
                  evy1L <- sum((xyL$y - richards(xyL$x, as.numeric(oppar$par[1]), 
                    as.numeric(oppar$par[2]), as.numeric(oppar$par[3]), 
                    as.numeric(oppar$par[4])))^2)
                  evy2L <- sum((xyL$y - richards(xyL$x, RAsym, 
                    Rk, Ri, RM))^2)
                  if(is.na(evy1L)) evy1L <- Inf
                  if(is.na(evy2L)) evy2L <- Inf
                  if (evy1L > evy2L | !is.na(modelparams$twocomponent.x)) {
                    if (!is.na(modelparams$twocomponent.x)) 
                      savoppL <- oppar
                    xyL <- moddata(xyL, min(xy$x, Intage, taper.ends = taper.ends))[["indata"]]
                    try(oppar <- optim(inputval, func2, method = "L-BFGS-B", 
                      lower = inputmin, upper = inputmax, control = list(maxit = 2000)), 
                      silent = TRUE)
                    if (!is.na(modelparams$twocomponent.x)) {
                      if (class(oppar)[1] == "try-error") {
                        oppar <- savoppL
                      } else {
                        if (oppar$value > savoppL$value) 
                          oppar <- savoppL
                      }
                    }
                  }
                }
            }
            if (min(xyL$x) != savminx | max(xyL$x) != savmaxx) 
                xyL <- subset(xyL, xyL$x >= savminx & xyL$x <= 
                  savmaxx)
            xyL <- xyL[order(xyL$x), ]
            if (is.na(oppar[1]) == FALSE) {
                if (oppar$convergence < 52) {
                  val1 <- data.frame(t(c(oppar$par)))
                  names(val1)<-names(oppar$par)
                } else {
                  if (modelparams$verbose == TRUE) 
                    print("Warning: negative optimization failed using default parameters (straight line)")
                   if(!is.na(modelparams$twocomponent.x)) {
		              subxyp <- subset(xyE$x, !is.na(xyE$x) )
		              addRAsym <- richards(subxyp [length(subxyp )], value[1L], value[2L], value[3L], value[4L])
		              val1 <- data.frame(RAsym=(Asym * 0.05)+addRAsym, Rk=value[2L], Ri= tAsym, RM=value[4L])
		              } else {
		              val1 <- data.frame(RAsym=(Asym * 0.05),Rk=value[2L], Ri= tAsym, RM= value[4L])
            	  }
                  inputval <- data.frame(t(inputval))
                  if (length(inputval$RAsym) == 0) 
                    val1[1] = NA
                  if (length(inputval$Rk) == 0) 
                    val1[2] = NA
                  if (length(inputval$Ri) == 0) 
                    val1[3] = NA
                  if (length(inputval$RM) == 0) 
                    val1[4] = NA
          	  val1<-as.numeric(val1)
                  names(val1) <- c("RAsym", "Rk", "Ri", "RM")
                  val1nms <- names(val1[!is.na(val1)])
                  val1 <- val1[!is.na(val1)]
                  val1 <- data.frame(t(val1))
                  names(val1) <- val1nms
                }
            } else {
                if (modelparams$verbose == TRUE) 
                  print("Warning: negative optimization failed using default parameters (straight line)")
                  if(!is.na(modelparams$twocomponent.x)) {
		              subxyp <- subset(xyE$x, !is.na(xyE$x) )
		              addRAsym <- richards(subxyp [length(subxyp )], Asym,value[2L],  tAsym, value[4L])
		              val1 <- data.frame(RAsym=(Asym * 0.05)+addRAsym, Rk=value[2L], Ri= tAsym, RM=value[4L])
		              } else {
		              val1 <- data.frame(RAsym=(Asym * 0.05),Rk=value[2L], Ri= tAsym, RM= value[4L])
            	  }
                inputval <- data.frame(t(inputval))
                if (length(inputval$RAsym) == 0) 
                  val1[1] = NA
                if (length(inputval$Rk) == 0) 
                  val1[2] = NA
                if (length(inputval$Ri) == 0) 
                  val1[3] = NA
                if (length(inputval$RM) == 0) 
                  val1[4] = NA
                val1<-as.numeric(val1)
                names(val1) <- c("RAsym", "Rk", "Ri", "RM")
                val1nms <- names(val1[!is.na(val1)])
		                  val1 <- val1[!is.na(val1)]
		                  val1 <- data.frame(t(val1))
                  names(val1) <- val1nms
            }
		parload <- get(optvarnm, envir = Envir )
                pnmodelparamsbounds <- parload[ names(parload) %in% names(pnmodelparamsbounds)]
        }
        if (modno == 12 | modno == 20 | modno == 32 | modno == 
            19 | modelparams$force4par == TRUE) {
            if (modno != 20 & modno != 32) {
                value <- value[1:4]
            } else {
                if (modno < 17 | modno >= 18) {
                value <- value[1:3]
                } else {
                  value <- data.frame(value[1], value[3], value[4])
                }
            }
        } else {
            if (length(val1$RAsym) == 1) {
                RAsym <- val1$RAsym
            } else {
                RAsym <- modelparams$RAsym
            }
            if (length(val1$Rk) == 1) {
                Rk <- val1$Rk
                } else {
                Rk <- modelparams$Rk
            }
            if (length(val1$Ri) == 1) {
                Ri <- val1$Ri
                } else {
                Ri <- modelparams$Ri
            }
            if (length(val1$RM) == 1) {
                RM <-   val1$RM
                } else {
                RM <- modelparams$RM
            }
            value <- data.frame(value[1L], value[2L], value[3L], 
                value[4L], RAsym, Rk, Ri, RM)
            names(value) <- c("Asym", "K", "Infl", "M", "RAsym", 
                "Rk", "Ri", "RM")
            Asym = value[1]
            K = value[2]
            Infl = value[3]
            M = value[4]
            RAsym = value[5]
            Rk = value[6]
            Ri = value[7]
            RM = value[8]
       	    value<-sapply(value,function(x) list(x))
            if (length(val1$RAsym) == 1) {
                RAsym <- val1$RAsym
            } else {
                value[5] <- NA
            }
            if (length(val1$Rk) == 1) {
                Rk <- val1$Rk
            } else {
                value[6] <- NA
            }
            if (length(val1$Ri) == 1) {
                Ri <- val1$Ri
            } else {
                value[7] <- NA
            }
            if (length(val1$RM) == 1) {
                RM <- val1$RM
            } else {
                value[8] <- NA
            }
            value<-unlist(value)
            if (modno > 19) 
                value[4] <- NA
            if (modno >= 17 & modno < 18) 
                value[2] <- NA
            savname <- names(value[!is.na(value)])
            value <- data.frame(t(value[!is.na(value)]))
            names(value) <- savname
            if (modno > 19) {
                posmin <- c(Amin, Kmin, Imin)
                posmax <- c(Amax, Kmax, Imax)
            } else {
                if (modno >= 17 & modno < 18) {
                  posmin <- c(Amin, Imin, Mmin)
                  posmax <- c(Amax, Imax, Mmax)
                } else {
                  posmin <- c(Amin, Kmin, Imin, Mmin)
                  posmax <- c(Amax, Kmax, Imax, Mmax)
                }
            }
            if (is.na(modelparams$twocomponent.x) == TRUE) {
                finalpars <- 0
                is.na(finalpars) <- TRUE
                richardsR <- function(Rparams) {
                  val2 <- data.frame(Asym = modelparams$Asym, 
                    K = modelparams$K, Infl = modelparams$Infl, 
                    M = modelparams$M, RAsym = modelparams$RAsym, 
                    Rk = modelparams$Rk, Ri = modelparams$Ri, 
                    RM = modelparams$RM)
                  val3 <- (data.frame(t(Rparams)))
                  if (length(val3$Asym) == 1) 
                    val2$Asym <- val3$Asym
                  if (length(val3$K) == 1) 
                    val2$K <- val3$K
                  if (length(val3$Infl) == 1) 
                    val2$Infl <- val3$Infl
                  if (length(val3$M) == 1) 
                    val2$M <- val3$M
                  if (length(val3$RAsym) == 1) 
                    val2$RAsym <- val3$RAsym
                  if (length(val3$Rk) == 1) 
                    val2$Rk <- val3$Rk
                  if (length(val3$Ri) == 1) 
                    val2$Ri <- val3$Ri
                  if (length(val3$RM) == 1) 
                    val2$RM <- val3$RM
                  Asym <- val2$Asym
                  K <- val2$K
                  Infl <- val2$Infl
                  M <- val2$M
                  RAsym <- val2$RAsym
                  Rk <- val2$Rk
                  Ri <- val2$Ri
                  RM <- val2$RM
                  if (is.na(exp(-K * (min(xy$x) - Infl))) == 
                    TRUE | (exp(-K * (min(xy$x) - Infl))) == 
                    Inf) {
                    K = modelparams$K
                    Infl = modelparams$Infl
                  }
                  if (is.na(exp(-K * (max(xy$x) - Infl))) == 
                    TRUE | (exp(-K * (min(xy$x) - Infl))) == 
                    Inf) {
                    K = modelparams$K
                    Infl = modelparams$Infl
                  }
                  if (is.na(exp(-Rk * (min(xy$x) - Ri))) == TRUE | 
                    (exp(-Rk * (min(xy$x) - Ri))) == Inf) {
                    Rk = modelparams$Rk
                    Ri = modelparams$Ri
                  }
                  if (is.na(exp(-Rk * (max(xy$x) - Ri))) == TRUE | 
                    (exp(-Rk * (min(xy$x) - Ri))) == Inf) {
                    Rk = modelparams$Rk
                    Ri = modelparams$Ri
                  }
                  options(warn = -1)
                  if(modno == 17.2 | modno == 17.4) RAsym <- Asym
                  if(modno == 17.3| modno == 17.4) RM <- M              
                  if (Re(as.complex(1 + M[1] * exp(-K[1] * (xy$x - 
                    Infl[1])))) < 0) {
                    if (Re(as.complex(1 + RM[1] * exp(-Rk[1] * (xy$x - 
                      Ri[1])))) < 0) {
                      if (modno >= 17 & modno < 18) {
                            y1 <- SSposnegRichardsF17(xy$x, Asym, 
                          Infl, M, RAsym, Ri, RM)
                        y1[is.na(y1)] <- 1e-290 * pnmodelparams$RAsym
                        y1[y1 == Inf] <- 1e-290 * pnmodelparams$RAsym
                        evl <- sum((xy$y - y1)^2)
   		if(!is.na(evl)){
                if (evl == Inf) {
                  evl <- 1e+290
                	} else {
                  evl <- 1e+290
                		}
               	}
                	try(if (min(Im(SSposnegRichardsF17((0:max(xy$x)), 
                          Asym, Infl, M, RAsym, Ri, RM)) < 0)) {
                          evl <- 1e+200
                        }, silent = TRUE)
                      } else {
                        y1 <- SSposnegRichardsFMRM(xy$x, Asym, 
                          K, Infl, M, RAsym, Rk, Ri, RM)
                        y1[is.na(y1)] <- 1e-290 * pnmodelparams$RAsym
                        y1[y1 == Inf] <- 1e-290 * pnmodelparams$RAsym
                        evl <- sum((xy$y - y1)^2)
   		if(!is.na(evl)){
                if (evl == Inf) {
                  evl <- 1e+290
                	} else {
                  evl <- 1e+290
                		}
               	}
                try(if (min(Im(SSposnegRichardsFMRM((0:max(xy$x)), 
                          Asym, K, Infl, M, RAsym, Rk, Ri, RM)) < 
                          0)) {
                          evl <- 1e+200
                        }, silent = TRUE)
                      }
                    } else {
                      if (modno >= 17 & modno < 18) {
                        y1 <- SSposnegRichardsF17(xy$x, Asym, 
                          Infl, M, RAsym, Ri, RM)
                        y1[is.na(y1)] <- 1e-290 * pnmodelparams$RAsym
                        y1[y1 == Inf] <- 1e-290 * pnmodelparams$RAsym
                        evl <- sum((xy$y - y1)^2)
   		if(!is.na(evl)){
                if (evl == Inf) {
                  evl <- 1e+290
                	} else {
                  evl <- 1e+290
                		}
               	}
                	try(if (min(Im(SSposnegRichardsF17((0:max(xy$x)), 
                          Asym, Infl, M, RAsym, Ri, RM)) < 0)) {
                          evl <- 1e+200
                        }, silent = TRUE)
                      } else {
                        y1 <- SSposnegRichardsFM(xy$x, Asym, 
                          K, Infl, M, RAsym, Rk, Ri, RM)
                        y1[is.na(y1)] <- 1e-290 * pnmodelparams$RAsym
                        y1[y1 == Inf] <- 1e-290 * pnmodelparams$RAsym
                        evl <- sum((xy$y - y1)^2)
   		if(!is.na(evl)){
                if (evl == Inf) {
                  evl <- 1e+290
                	} else {
                  evl <- 1e+290
                		}
               	}
                try(if (min(Im(SSposnegRichardsFM((0:max(xy$x)), 
                          Asym, K, Infl, M, RAsym, Rk, Ri, RM)) < 
                          0)) {
                          evl <- 1e+200
                        }, silent = TRUE)
                      }
                    }
                  } else {
                    if (Re(as.complex(1 + RM[1] * exp(-Rk[1] * (xy$x - 
                      Ri[1])))) < 0) {
                      if (modno >= 17 & modno < 18) {
                        y1 <- SSposnegRichardsF17(xy$x, Asym, 
                          Infl, M, RAsym, Ri, RM)
                        y1[is.na(y1)] <- 1e-290 * pnmodelparams$RAsym
                        y1[y1 == Inf] <- 1e-290 * pnmodelparams$RAsym
                        evl <- sum((xy$y - y1)^2)
   		if(!is.na(evl)){
                if (evl == Inf) {
                  evl <- 1e+290
                	} else {
                  evl <- 1e+290
                		}
               	}
                try(if (min(Im(SSposnegRichardsF17((0:max(xy$x)), 
                          Asym, Infl, M, RAsym, Ri, RM)) < 0)) {
                          evl <- 1e+200
                        }, silent = TRUE)
                      } else {
                        y1 <- SSposnegRichardsFRM(xy$x, Asym, 
                          K, Infl, M, RAsym, Rk, Ri, RM)
                        y1[is.na(y1)] <- 1e-290 * pnmodelparams$RAsym
                        y1[y1 == Inf] <- 1e-290 * pnmodelparams$RAsym
                        evl <- sum((xy$y - y1)^2)
   		if(!is.na(evl)){
                if (evl == Inf) {
                  evl <- 1e+290
                	} else {
                  evl <- 1e+290
                		}
               	}    
                try(if (min(Im(SSposnegRichardsFRM((0:max(xy$x)), 
                          Asym, K, Infl, M, RAsym, Rk, Ri, RM)) < 
                          0)) {
                          evl <- 1e+200
                        }, silent = TRUE)
                      }
                    } else {
                      if (modno >= 17 & modno < 18) {
                        y1 <- SSposnegRichardsF17(xy$x, Asym, 
                          Infl, M, RAsym, Ri, RM)
                        y1[is.na(y1)] <- 1e-290 * pnmodelparams$RAsym
                        y1[y1 == Inf] <- 1e-290 * pnmodelparams$RAsym
                        evl <- sum((xy$y - y1)^2)
    		if(!is.na(evl)){
                 if (evl == Inf) {
                   evl <- 1e+290
                 	} else {
                   evl <- 1e+290
                 		}
               	}
                	try(if (min(Im(SSposnegRichardsF17((0:max(xy$x)), 
                          Asym, Infl, M, RAsym, Ri, RM)) < 0)) {
                          evl <- 1e+200
                        }, silent = TRUE)
                      } else {
                        y1 <- SSposnegRichardsF(xy$x, Asym, K, 
                          Infl, M, RAsym, Rk, Ri, RM)
                        y1[is.na(y1)] <- 1e-290 * pnmodelparams$RAsym
                        y1[y1 == Inf] <- 1e-290 * pnmodelparams$RAsym
                        evl <- sum((xy$y - y1)^2)
   		if(!is.na(evl)){
                if (evl == Inf) {
                  evl <- 1e+290
                	} else {
                  evl <- 1e+290
                		}
               	}
                	try(if (min(Im(SSposnegRichardsF((0:max(xy$x)), 
                          Asym, K, Infl, M, RAsym, Rk, Ri, RM)) < 
                          0)) {
                          evl <- 1e+200
                        }, silent = TRUE)
                      }
                    }
                  }
   		if(!is.na(evl)){
                if (evl == Inf) {
                  evl <- 1e+290
                	} else {
                  evl <- 1e+290
                		}
               	} 
                options(warn = 0)
                  return(evl)
                }
                savname <- names(value)
                upbnds <- data.frame(value[1:length(value)])
                names(upbnds) <- savname
                dnbnds <- upbnds
                options(warn=-1)
                if ("Asym" %in% names(value)) {
                  upbnds$Asym = Asym + (abs(pnmodelparamsbounds$Amax - 
                    pnmodelparamsbounds$Amin) * 0.1)
                }
                if ("Asym" %in% names(value)) {
                  dnbnds$Asym = Asym - (abs(pnmodelparamsbounds$Amax - 
                    pnmodelparamsbounds$Amin) * 0.1)
                }
                if ("K" %in% names(value)) {
                  upbnds$K = K + (abs(pnmodelparamsbounds$Kmax - 
                    pnmodelparamsbounds$Kmin) * 0.25)
                }
                if ("K" %in% names(value)) {
                  dnbnds$K = K - (abs(pnmodelparamsbounds$Kmax - 
                    pnmodelparamsbounds$Kmin) * 0.25)
                }
                if ("Infl" %in% names(value)) {
                  upbnds$Infl = Infl + (abs(pnmodelparamsbounds$Imax - 
                    pnmodelparamsbounds$Imin) * 0.35)
                }
                if ("Infl" %in% names(value)) {
                  dnbnds$Infl = Infl - (abs(pnmodelparamsbounds$Imax - 
                    pnmodelparamsbounds$Imin) * 0.2)
                }
                if ("M" %in% names(value)) {
                  upbnds$M = M + abs(pnmodelparamsbounds$Mmax - 
                    pnmodelparamsbounds$Mmin) * 0.2
                }
                if ("M" %in% names(value)) {
                  dnbnds$M = M - abs(pnmodelparamsbounds$Mmax - 
                    pnmodelparamsbounds$Mmin) * 0.2
                }
                if ("RAsym" %in% names(value)) {
                  upbnds$RAsym = RAsym + (abs(pnmodelparamsbounds$RAmax - 
                    pnmodelparamsbounds$RAmin) * 0.5)
                }
                if ("RAsym" %in% names(value)) {
                  dnbnds$RAsym = RAsym - (abs(pnmodelparamsbounds$RAmax - 
                    pnmodelparamsbounds$RAmin) * 0.5)
                }
                if ("Rk" %in% names(value)) {
                  upbnds$Rk = Rk + (abs(pnmodelparamsbounds$Rkmax - 
                    pnmodelparamsbounds$Rkmin) * 0.25)
                }
                if ("Rk" %in% names(value)) {
                  dnbnds$Rk = Rk - (abs(pnmodelparamsbounds$Rkmax - 
                    pnmodelparamsbounds$Rkmin) * 0.25)
                }
                if ("Ri" %in% names(value)) {
                  upbnds$Ri = Ri + (abs(pnmodelparamsbounds$Rimax - 
                    pnmodelparamsbounds$Rimin) * 0.35)
                }
                if ("Ri" %in% names(value)) {
                  dnbnds$Ri = Ri - (abs(pnmodelparamsbounds$Rimax - 
                    pnmodelparamsbounds$Rimin) * 0.2)
                }
                if ("RM" %in% names(value)) {
                  upbnds$RM = RM + abs(pnmodelparamsbounds$RMmax - 
                    pnmodelparamsbounds$RMmin) * 0.2
                }
                if ("RM" %in% names(value)) {
                  dnbnds$RM = RM - abs(pnmodelparamsbounds$RMmax - 
                    pnmodelparamsbounds$RMmin) * 0.2
                }
                options(warn = 0)
                cnt <- 0
                finalpars <- 0
                is.na(finalpars) <- TRUE
                oppar1 <- 0
                is.na(oppar1) <- TRUE
                names(oppar1) <- "convergence"
                repoptm <- 1
                while (is.na(finalpars[1L]) == TRUE) {
                  if (cnt > 1) {
                    if ("Rk" %in% names(value)) {
                      if (cnt == 2) {
                      options(warn=-1)
                        value$Rk <- pnmodelparamsbounds$Rkmin
                        dnbnds$Rk <- pnmodelparamsbounds$Rkmin - 
                          abs(pnmodelparams$Rk) * 0.5
                        upbnds$Rk <- pnmodelparamsbounds$Rkmin + 
                          abs(pnmodelparams$Rk) * 0.5
                        repoptm <- 1
                        options(warn=0)
                      }
                    } else {
                      cnt = cnt + 1
                    }
                    if ("Rk" %in% names(value)) {
                      if (cnt == 3) {
                      options(warn=-1)
                        value$Rk <- pnmodelparamsbounds$Rkmax
                        dnbnds$Rk <- pnmodelparamsbounds$Rkmax - 
                          abs(pnmodelparams$Rk) * 0.5
                        upbnds$Rk <- pnmodelparamsbounds$Rkmax + 
                          abs(pnmodelparams$Rk) * 0.5
                        repoptm <- 1
                      options(warn=0)
                      }
                    } else {
                      cnt = cnt + 1
                    }
                    if (cnt == 4 & (modno < 17 | modno >= 18) & ("K" %in% names(value))) {
                     options(warn=-1)
                      if ("Rk" %in% names(value)) 
                        value$Rk <- pnmodelparams$Rk
                      if ("Rk" %in% names(value)) 
                        dnbnds$Rk <- pnmodelparamsbounds$Rkmin
                      if ("Rk" %in% names(value)) 
                        upbnds$Rk <- pnmodelparamsbounds$Rkmax
                      value$K <- pnmodelparamsbounds$Kmin
                      dnbnds$K <- pnmodelparamsbounds$Kmin - 
                        abs(pnmodelparams$K) * 0.5
                      upbnds$K <- pnmodelparamsbounds$Kmin + 
                        abs(pnmodelparams$K) * 0.5
                      repoptm <- 1
                       options(warn=0)
                    }
                    if (cnt == 5 & (modno < 17 | modno >= 18) & ("K" %in% names(value))) {
                     options(warn=-1)
                      value$K <- pnmodelparamsbounds$Kmax
                      dnbnds$K <- pnmodelparamsbounds$Kmax - 
                        abs(pnmodelparams$K) * 0.5
                      upbnds$K <- pnmodelparamsbounds$Kmax + 
                        abs(pnmodelparams$K) * 0.5
                      repoptm <- 1
                       options(warn=0)
                    }
                  }
                  options(warn=0)
                  dnbnds[dnbnds == 0] <- 1e-05
                  upbnds[upbnds == 0] <- 1e-05
                  value[value == 0] <- 1e-05
                  nmbndsav <- names(upbnds)
                  dnbnds <- unlist(dnbnds)
                  upbnds <- unlist(upbnds)
                  names(dnbnds) <- nmbndsav
                  names(upbnds) <- nmbndsav
                  oppar1 <- data.frame(52)
                  names(oppar1) <- c("convergence")
         	  if( is.na(inputval[1])) {
          	  	inputval<-val1
          	  	} else {
          	  	chinputval <- data.frame(t(inputval))
          	  	chinputmax <- data.frame(t(inputmax))
          	  	if(length(chinputval) == 1) chinputval <- data.frame(inputval)
          	  	if(length(chinputmax) == 1) chinputmax <- data.frame(inputmax)
 			names(chinputval) <- names(inputval)
            		chmx<-paste(substr(names(inputval),1,2),"max","")
            		chmx<-sub(" ","",chmx)
            		chmx<-sub(" ","",chmx)
            		names(chinputmax) <- chmx
            	  	if (length(chinputval$Ri) == 1) {
          	  		if (chinputval$Ri >= chinputmax$Rimax) 
          	  			chinputval$Ri <- chinputmax$Rimax - ((chinputmax$Rimax - min(xy$x, na.rm=TRUE)) * 0.95)
          	 	 }
          	  }
                  if (repoptm == 1) {
     		    SScurveeval <-function(value,xy){
     		    calltxt<-""
     		    nm <-names(value)
     		    value<-as.numeric(value)
     		    names(value)<-nm
     		    if("Asym" %in% names(value)) {
     		    calltxt <- "Asym=Asym"
     		    Asym <- as.numeric(value[names(value)=="Asym"])} else {
     		    Asym <- NA}
     		    if("K" %in% names(value)) {
     		    calltxt <- paste(calltxt,"K=K",sep=",")
     		    K <- as.numeric(value[names(value)=="K"])} else {
     		    K <- NA}
     		    if("Infl" %in% names(value)) {
     		    calltxt <- paste(calltxt,"Infl=Infl",sep=",")
     		    Infl <- as.numeric(value[names(value)=="Infl"])} else {
     		    Infl <- NA}
     		    if("M" %in% names(value)){
     		    calltxt <- paste(calltxt,"M=M",sep=",")
     		    M <- as.numeric(value[names(value)=="M"])} else {
     		    M <- NA}
     		    if("RAsym" %in% names(value)) {
     		    calltxt <- paste(calltxt,"RAsym=RAsym",sep=",")
     		    RAsym <- as.numeric(value[names(value)=="RAsym"])} else {
     		    RAsym <- NA}
     		    if("Rk" %in% names(value)) {
     		    calltxt <- paste(calltxt,"Rk=Rk",sep=",")
     		    Rk <- as.numeric(value[names(value)=="Rk"])} else {
     		    Rk <- NA}
     		    if("Ri" %in% names(value)) {
     		    calltxt <- paste(calltxt,"Ri=Ri",sep=",")
     		    Ri <- as.numeric(value[names(value)=="Ri"])} else {
     		    Ri <- NA}
     		    if("RM" %in% names(value)) {
     		    calltxt <- paste(calltxt, "RM=RM",sep=",")
     		    RM <- as.numeric(value[names(value)=="RM"])} else {
     		    RM <- NA}
     		    callfin <- paste("posnegRichards.eqn(xy$x,",calltxt,",modno = ",modno,",pn.options = ",paste("'",optvarnm,"'",sep=""),")",sep="")
     		    ynought <- sum(abs(xy$y-eval(parse(text=callfin)))^2)
     		    return(ynought)
     	            }
     		    if(counterR == 0) saveb4 <- value
		    counterR <- counterR + 1
		    if(counterR == 1) {
		    savparR <- list(Inf)
		    names(savparR) <- "1"
		    savyR <- list(Inf)
		    names(savyR) <- "1"
		    }
		     eval(parse(text=paste("savparR$'",counterR,"'<-value",sep="")))
		     eval(parse(text=paste("savyR$'",counterR,"'<-",SScurveeval(value,xy),sep="")))
		     parscaleR <- sapply(value,roundparscale)
                     try(oppar1 <- (optim(value, richardsR, method = "L-BFGS-B", 
                      lower = dnbnds, upper = upbnds, control = list(maxit = 1000,
                      parscale = parscaleR))), 
                      silent = TRUE)
                     counterR <- counterR + 1
		     eval(parse(text=paste("savparR$'",counterR,"'<-oppar1$par",sep="")))
		     eval(parse(text=paste("savyR$'",counterR,"'<-",SScurveeval( oppar1$par,xy),sep="")))
                  }
                  repoptm <- 0
                  cnt <- cnt + 1
                  if (oppar1$convergence < 52) {
                    finalpars <- data.frame(t(c(oppar1$par)))
                  } else {
                    if (cnt > 4) {
                      finalpars <- 1
                    } else {
                      is.na(finalpars) <- TRUE
                    }
                  }
                }
                if (cnt > 4) 
                  is.na(finalpars) <- TRUE
                options(warn = 0)
                if (is.na(finalpars[1] == TRUE)) {
                  if (modelparams$verbose == TRUE) 
                    print("Warning: simultaneous optimization of pre- and post- peak curves failed, using separately fitted parameters")
                  value<-NA
                  try(eval(parse(text=paste("value<-savparR$'",which.min(savyR),"'",sep=""))),silent=T)
                   if( class(value) == "try-error" ) value <- NA
                  if( is.na(value[1]) ) value <- saveb4
                } else {
                  finpars <- finalpars
                  if ("Rk" %in% names(finpars)) {
                    if (finpars$Rk < 1e-04 & finpars$Rk > -1e-04) {
                      finpars$Rk <- finalpars$Rk
                      finpars[1L] <- finalpars[1L] + (finalpars$Rk/(finalpars[1L]/finalpars$Rk))
                    }
                  }
                  value <- c(finpars)
                }
            }
        }
    }
    if (modelparams$force4par == TRUE | (modelparams$force.nonmonotonic == FALSE & modno == 12) |
            (modelparams$force.nonmonotonic == FALSE & modno == 32) ){
    	 if (modno != 20 & modno != 32) {
                   if (modno < 17 | modno >= 18) {
                         value <- value[1:4]
                   	 names(value) <- mCall[c("Asym", "K", "Infl", "M")]
                   } else {
                         value <- data.frame(value[1],value[3],value[4])
               	   	 names(value) <- mCall[c("Asym", "Infl", "M")]
                   }
   	 } else {
   	 		 value <- value[1:3]
              		 names(value) <- mCall[c("Asym", "K", "Infl")]
    	}
    value <- value[(names(value))!="NULL"]
    value <- value[!is.na(names(value))]
    } else {
        if (modno > 19) {
            vnames <- c("Asym", "K", "Infl")
        } else {
            vnames <- c("Asym", "K", "Infl", "M")
        }
        if (modno >= 17 & modno < 18) vnames <- c("Asym", "Infl", "M")
        if (modno == 1 | modno == 21 | modno == 18) {
            names(value) <- mCall[c(vnames, "RAsym", "Rk", "Ri", 
                "RM")]
        }
        if (modno == 2 | modno == 22) {
            names(value) <- mCall[c(vnames, "RAsym", "Rk", "Ri")]
        }
        if (modno == 3 | modno == 23 | modno == 17.1) {
            names(value) <- mCall[c(vnames, "Ri", "RM")]
        }
        if (modno == 4 | modno == 24) {
            names(value) <- mCall[c(vnames, "RAsym", "RM")]
        }
        if (modno == 5 | modno == 25) {
            names(value) <- mCall[c(vnames, "RM")]
        }
        if (modno == 6 | modno == 26 | modno == 17) {
            names(value) <- mCall[c(vnames, "RAsym", "Ri", "RM")]
        }
        if (modno == 7 | modno == 27) {
            names(value) <- mCall[c(vnames, "Rk", "Ri", "RM")]
        }
        if (modno == 8 | modno == 28) {
            names(value) <- mCall[c(vnames, "RAsym", "Rk", "RM")]
        }
        if (modno == 9 | modno == 29) {
            names(value) <- mCall[c(vnames, "Rk", "RM")]
        }
        if (modno == 10 | modno == 30 | modno == 17.3) {
            names(value) <- mCall[c(vnames, "Ri")]
        }
        if (modno == 11 | modno == 31) {
            names(value) <- mCall[c(vnames, "RAsym")]
        }
        if (modno == 12 | modno == 32 | modno == 19) {
            names(value) <- mCall[vnames]
        }
        if (modno == 13 | modno == 33 | modno == 17.2) {
            names(value) <- mCall[c(vnames, "RAsym", "Ri")]
        }
        if (modno == 14 | modno == 34) {
            names(value) <- mCall[c(vnames, "Rk", "Ri")]
        }
        if (modno == 15 | modno == 35) {
            names(value) <- mCall[c(vnames, "RAsym", "Rk")]
        }
        if (modno == 16 | modno == 36) {
            names(value) <- mCall[c(vnames, "Rk")]
        }
        if (modno == 20) {
            names(value) <- mCall[c("Asym", "K", "Infl")]
        }
    }
    options(warn = 0)
    if (modno == 18 | modno == 19 | frommodpar == TRUE) {
        Asym = as.numeric(value[1])
        K = as.numeric(value[2])
        Infl = as.numeric(value[3])
        M = as.numeric(value[4])
            Amax = Asym + diff(range(xy$y, na.rm = TRUE) * mediummult)
            Amin = Asym - diff(range(xy$y, na.rm = TRUE) * mediummult)
            Kmax = K + (abs(K) * smallmult)
            Kmin = K - (abs(K) * smallmult)
 	    if (Kmin < 1e-05 & Kmin >= 0 ) {
 	            Kmin = 1e-05
 	         } else {
            if (Kmin > (-1e-05) & Kmin < 0) Kmin = (-1e-05)
          				    }
            Imax = Infl + diff(range(xy$x, na.rm = TRUE)  * mediummult)
            Imin = Infl + diff(range(xy$x, na.rm = TRUE)  * -mediummult)
            while (abs(Imax * Kmax) > 700) Imax = Imax * 0.9
            while (abs(Imin * Kmax) > 700) Imin = Imin * 0.9
               	 if( !is.na(Imax) ) {
        		if (Imax > maxIval) Imax <- maxIval
           	if(Imax < Infl) {
           	   maxIval <- Infl + ((maxIval - min(xy$x)) * 0.05)
           	   Imax <- maxIval  
           	   }
           	  		}
           	  if( !is.na(Imin) ) {
          	   if(Imin > Infl) Imin <- min(xy$x, na.rm=TRUE) 
          	   if(Imin > Infl) {
        		   Imin <- Infl - ((maxIval - min(xy$x)) * 0.05)
           	   			}

	    	}
         if (abs(M) < 0.1) {
            Mmax = smallmult
            Mmin = -smallmult
        } else {
            Mmax = M + abs(M * smallmult)
            Mmin = M - abs(M * smallmult)
        }
        if (modno == 18 | modno == 1) {
            RAsym = as.numeric(value[5])
            Rk = as.numeric(value[6])
            Ri = as.numeric(value[7])
            RM = as.numeric(value[8])
            RAmax = RAsym + (abs(RAsym) * mediummult)
            RAmin = RAsym - (abs(RAsym) * mediummult)
            Rkmax = Rk + (abs(Rk) * smallmult)
            Rkmin = Rk - (abs(Rk) * smallmult)
            Rimax = Ri + (abs(Ri) * smallmult)
            Rimin = Ri + (abs(Ri) * -smallmult)
            while (abs(Rimax * Rkmax) > 700) Rimax = Rimax * 0.9
            while (abs(Rimin * Rkmax) > 700) Rimin = Rimin * 0.9     	 
            if( !is.na(Rimax) ) {
       		if (Rimax > maxRival) Rimax <- maxRival
          	if(Rimax < tAsym) {
          	   maxRival <- tAsym + ((maxRival - min(xy$x)) * 0.95)
          	   Rimax <- maxRival  
          	   }
          	  		}
          	  if( !is.na(Rimin) ) {
         	   if(Rimin > tAsym) Rimin <- min(xy$x, na.rm=TRUE) 
         	   if(Rimin > tAsym) {
       		   Rimin <- tAsym - ((maxIval - min(xy$x)) * 0.95)
          	   			}
	    }	
	    if (abs(RM) < 0.1) {
	                RMmax = smallmult
	                RMmin = -smallmult
	            } else {
	                RMmax = RM + abs(RM * smallmult)
	                RMmin = RM - abs(RM * smallmult)
        }
         } else {
            RAsym <- NA
            Rk <- NA
            Ri <- NA
            RM <- NA
            RAmax <- NA
            RAmin <- NA
            Rkmax <- NA
            Rkmin <- NA
            Rimax <- NA
            Rimin <- NA
            RMmax <- NA
            RMmin <- NA
        }
        skel <- rep(list(1), 16)
        exportparams <- c(Asym, K, Infl, M, RAsym, Rk, Ri, RM, 
            modelparams$first.y, modelparams$x.at.first.y, 
            modelparams$last.y, modelparams$x.at.last.y, 
            modelparams$twocomponent.x, modelparams$verbose, 
            modelparams$force4par, modelparams$force.nonmonotonic)
        exportparams <- relist(exportparams, skel)
        names(exportparams) <- c("Asym", "K", "Infl", "M", "RAsym", 
            "Rk", "Ri", "RM", "first.y", "x.at.first.y", "last.y", "x.at.last.y", 
            "twocomponent.x", "verbose", "force4par", "force.nonmonotonic")
        exportparams$twocomponent.x <- modelparams$twocomponent.x
        exportparams$verbose <- modelparams$verbose
        exportparams$force4par <- modelparams$force4par
        exportparams$force.nonmonotonic <- modelparams$force.nonmonotonic
        skel <- rep(list(1), 16)
        exportparamsbounds <- c(Amin, Amax, Kmin, Kmax, Imin, 
            Imax, Mmin, Mmax, RAmin, RAmax, Rkmin, Rkmax, Rimin, 
            Rimax, RMmin, RMmax)
        exportparamsbounds <- relist(exportparamsbounds, skel)
        names(exportparamsbounds) <- c("Amin", "Amax", "Kmin", 
            "Kmax", "Imin", "Imax", "Mmin", "Mmax", "RAmin", 
            "RAmax", "Rkmin", "Rkmax", "Rimin", "Rimax", "RMmin", 
            "RMmax")
     	valexp <- sapply( c(unlist(exportparams),unlist(exportparamsbounds), 
             	unlist( optvar[names(optvar) %w/o% c(names(exportparams),names(exportparamsbounds)) ]) ), function(x) list(x))
        	names( valexp[1:32] ) <- names( c(exportparams[1:16],exportparamsbounds) )
      	valexp[c(1:12,17:32)] <- as.numeric( valexp[c(1:12,17:32)] )
        	valexp[14:16] <- as.logical( valexp[14:16] )
        	if(is.numeric(modelparams$twocomponent.x)) {
        	valexp[13] <- as.numeric( valexp[13] )
        	} else {
  		valexp[13] <- as.logical( valexp[13] )     	
        	}
      	try(valexp$taper.ends <- as.numeric(valexp$taper.ends), silent = TRUE)
    	assign(optvarnm, valexp, envir = Envir) 
    }
    if (modelparams$verbose == TRUE) {
        print("Values from SSposnegRichards:")
        prnval <- unlist(value)
        names(prnval) <- names(value)
        print(prnval)
    } 
     value
}, pnames = c("Asym", "K", "Infl", "M", "RAsym", "Rk", "Ri", 
    "RM"), class = "selfStart") 
