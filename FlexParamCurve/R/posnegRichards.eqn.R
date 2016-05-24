posnegRichards.eqn <-
structure(function(x, Asym = NA,
    K = NA, Infl = NA, M = NA, RAsym = NA, Rk = NA, Ri = NA, RM = NA,
    modno, pn.options, Envir = .GlobalEnv) {
    Envir1 <- try(FPCEnv$env,silent=T)
    env1ck <- try(is.environment(FPCEnv$env),silent=T)
    envck <- try(is.environment(Envir),silent=T)
    env.ck<-2
    if(envck == FALSE | class(envck) == "try-error") env.ck <- (env.ck - 1)
    if(env1ck == FALSE | class(env1ck) == "try-error") env.ck <- (env.ck - 1)
    if(env.ck == 2) {
    if(identical(Envir, Envir1) == FALSE & 
    	identical(Envir,.GlobalEnv) == TRUE) Envir <- Envir1
    }
    if(env.ck == 1 & (envck == FALSE | class(envck) == "try-error")) Envir <- Envir1
    FPCEnv$env <- Envir
    params<-list(Asym = Asym, K = K, Infl = Infl, M = M, RAsym = RAsym,
        Rk = Rk, Ri = Ri, RM = RM, first.y = NA, x.at.first.y = NA, 
        last.y = NA, x.at.last.y = NA, twocomponent.x = NA, 
        verbose = NA, force4par = NA)
    pnmodelparams <- rep(list(NA),15)
    pnmodelparams[14]<-FALSE
    pnmodelparams[15]<-FALSE
    names(pnmodelparams)<-names(params)
    pnopt<- get(as.character( pn.options[1] ), envir = Envir)
    pnmodelparams[ names(pnmodelparams) %in% names(pnopt) ] <- pnopt[ names(pnopt) %in% names(pnmodelparams) ]
    pnmodelparams[ names(pnmodelparams) %in% names(params[!is.na(params)]) ] <- params[  names(pnmodelparams) %in% names(params[!is.na(params)])  ]   
    if(is.na(pnmodelparams[14])) pnmodelparams[14] <- FALSE
    if(is.na(pnmodelparams[15])) pnmodelparams[15] <- FALSE
    "%w/o%" <- function(x, y) x[!x %in% y] #--  x without y
    compareNA <- function(v1,v2) {
        same <- (v1 == v2)  |  (is.na(v1) & is.na(v2))
        same[is.na(same)] <- FALSE
        return(same)
    }
    adjRM<-FALSE
    parcheck<-rep("None",8)
    if(modno < 18 & modno >= 17) { 
    parcheck[1:2] <-c("Asym" , "Infl")
    } else {
    parcheck[1:3] <-c("Asym" , "K" , "Infl")
    }
    Asym <- pnmodelparams$Asym
    K <- pnmodelparams$K
    Infl <- pnmodelparams$Infl
    if (is.na(pnmodelparams$Asym)) stop(paste("Parameter Asym required for modno = ", modno,
    	" but is absent from user-provided call", sep = ""))
    if (is.na(pnmodelparams$Infl)) stop(paste("Parameter Infl required for modno = ", modno,
    	" but is absent from user-provided call", sep = ""))
    if (is.na(pnmodelparams$K) & (modno < 17 | modno >= 18)) stop(paste("Parameter K required for modno = ", modno,
    	" but is absent from user-provided call", sep = ""))
    fractM <- pnmodelparams$M
    fractRM <- pnmodelparams$RM
    if (modno > 19) {
        fractM <- 1/pnmodelparams$M
        M <- pnmodelparams$M
    } else {
        if (is.na(pnmodelparams$M)) stop(paste("Parameter M required for modno = ", modno,
        	" but is absent from user-provided call", sep = ""))
        parcheck[4] <- "M"
        fractM <- 1/pnmodelparams$M
        M <- pnmodelparams$M
    }
    if (modno == 2 | modno == 22 | modno == 10 | modno == 30 |
        modno == 11 | modno == 31 | modno == 12 | modno == 32 |
        modno == 19 | modno == 13 | modno == 33 | modno == 14 |
        modno == 34 | modno == 15 | modno == 35 | modno == 16 |
        modno == 36 | modno == 20 | modno == 32) {
        fractRM <- 1/pnmodelparams$RM
        RM <- pnmodelparams$RM
    } else {
    if( modno != 17.2 & modno != 17.3 & is.na(pnmodelparams$RM)) stop(paste("Parameter RM required for modno = ", modno,
    	" but is absent from user-provided call", sep = ""))
    	fractRM <- 1/pnmodelparams$RM
        RM <- pnmodelparams$RM
        if( modno != 17.2 & modno != 17.3) adjRM<-TRUE
    }
    if (modno == 3 | modno == 23 | modno == 5 | modno == 25 |
        modno == 7 | modno == 27 | modno == 9 | modno == 29 |
        modno == 10 | modno == 30 | modno == 12 | modno == 32 |
        modno == 19 | modno == 14 | modno == 34 | modno == 16 |
        modno == 36 | modno == 20 | modno == 32) {
        RAsym <- pnmodelparams$RAsym
    } else {
        if (modno != 17.1 & modno != 17.3 & is.na(pnmodelparams$RAsym)) stop(paste("Parameter RAsym required for modno = ", modno,
        	" but is absent from user-provided call", sep = ""))
        RAsym <- pnmodelparams$RAsym
        if (modno != 17.1 & modno != 17.3) parcheck[5] <- "RAsym"
    }
    if (modno == 3 | modno == 23 | modno == 4 | modno == 24 |
        modno == 5 | modno == 25 | modno == 6 | modno == 26 |
        modno == 10 | modno == 30 | modno == 11 | modno == 31 |
        modno == 12 | modno == 32 | modno == 19 | modno == 13 |
        modno == 33 | modno == 20 | modno == 32) {
        Rk <- pnmodelparams$Rk
    } else {
        if (modno != 17 & modno != 17.1 & modno != 17.2 & modno != 17.3 & is.na(pnmodelparams$Rk)) 
        	stop(paste("Parameter Rk required for modno = ", modno,
        	" but is absent from user-provided call", sep = ""))
        Rk <- pnmodelparams$Rk
        if (modno != 17 & modno != 17.1 & modno != 17.2 & modno != 17.3) parcheck[6] <- "Rk"
    }
    if (modno == 4 | modno == 24 | modno == 5 | modno == 25 |
        modno == 8 | modno == 28 | modno == 9 | modno == 29 |
        modno == 11 | modno == 31 | modno == 12 | modno == 32 |
        modno == 19 | modno == 15 | modno == 35 | modno == 16 |
        modno == 36 | modno == 20 | modno == 32) {
        Ri <- pnmodelparams$Ri
    } else {
        if (is.na(pnmodelparams$Ri)) stop(paste("Parameter Ri required for modno = ", modno,
        	" but is absent from user-provided call", sep = ""))
        Ri <- pnmodelparams$Ri
        parcheck[7] <- "Ri"
    }
    if(adjRM == TRUE)  parcheck[8] <- "RM"
      misspar<-"EMPTY"
      noparms = 8
      if(pnmodelparams$force4par == TRUE | modno ==12 | modno == 32) noparms = 4 
      for(cint in 1:noparms){
      if(names(params[cint]) == parcheck[cint]) {
      if(misspar == "EMPTY" ){
       if(compareNA(unlist(params[cint]),unlist(pnmodelparams[cint])) == FALSE) misspar <- names(params[cint])
			}else{
       if(compareNA(unlist(params[cint]),unlist(pnmodelparams[cint])) == FALSE) misspar <- paste(misspar,names(params[cint]),sep=" , ")
			  }	
			  			}
      }
      if(misspar != "EMPTY") print(paste("WARNING: Missing parameters in call (",misspar ,") have been substituted from the pn.options object: ",
    		pn.options[1], sep=""))
    if(modno < 17 | modno >= 18){
    try({
    if (Re(as.complex(1 + M[1] * exp(-K[1] * (max(x) - Infl[1])))) < 0) {
     	    fractM <- round(1/pnmodelparams$M)
    } else {
            fractM <- 1/pnmodelparams$M
    }
    if (Re(as.complex(1 + RM[1] * exp(-Rk[1] * (max(x) - Ri[1])))) < 0) {
     	    fractRM <- round(1/pnmodelparams$RM)
    } else {
       	    fractRM <- 1/pnmodelparams$RM
    }
    },silent = TRUE)
    }
    if(modno == 12 | modno == 32) {
    pnmodelparams$force4par <- TRUE
    pnmodelparams$twocomponent.x <- NA
    print("note: the model selected does not have any second curve parameters (i.e. is monotonic)")
    					}  
    if (!is.na(pnmodelparams$twocomponent.x)) {
        if (pnmodelparams$twocomponent.x == TRUE) 
            stop("Please run modpar or provide an estimate of value for x \nat which the 
    	    first curve ends and second curve begins: twocomponent.x=##")
        if (pnmodelparams$force4par == TRUE)
            stop("Cannot force a two component Richards model to have a single component.\nEither set force4par to FALSE
            	\n or twocomponent.x to NA")
        c(Asym/Re(as.complex(1 + M * exp(-K * (x[x <= pnmodelparams$twocomponent.x] -
            Infl)))^(fractM)), (RAsym/Re(as.complex(1 + RM *
            exp(-Rk * (x[x > pnmodelparams$twocomponent.x] -
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
            if (pnmodelparams$force4par == TRUE) {
                (Asym/Re(as.complex(1 + M * exp(-K * (x - Infl)))^(fractM)))
            } else {
                (Asym/Re(as.complex(1 + M * exp(-K * (x - Infl)))^(fractM))) +
                  (RAsym/Re(as.complex(1 + RM * exp(-Rk * (x -
                    Ri)))^(fractRM)))
            }
        }
    }
}, ex = function(){
    require(graphics)
    modpar(posneg.data$age, posneg.data$mass)
    y <- posnegRichards.eqn(10, 1000, 0.5, 25, 1, 100, 0.5, 125, 1, modno = 1)

    y <- posnegRichards.eqn(10 ,1000 ,0.5 ,25 ,1 ,100 ,0.5 ,125 ,1 ,modno = 12)

    plot(1:200 ,posnegRichards.eqn(1:200 ,1000 ,0.5 ,25 ,1 ,100 ,0.5 ,125 ,1 ,modno = 12),xlim=c(1, 200),
    xlab = "x", ylab = "y",pch = 1, cex = 0.7)
    }
)
