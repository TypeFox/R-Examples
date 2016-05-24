modpar <-
structure(function
		(x,
                     y,
                     pn.options = NA,
                     first.y = NA,
                     x.at.first.y = NA,
                     last.y = NA,
                     x.at.last.y = NA,
                     twocomponent.x = NA,
                     verbose = FALSE,
                     force8par = FALSE,
                     force4par = FALSE,
                     suppress.text = FALSE,
                     taper.ends = 0.45,
                     width.bounds = 1,
                     bounds.error = FALSE,
                     Envir =  .GlobalEnv,
                     force.nonmonotonic = FALSE,
                     ...
                     ) {
    if(!is.na(pn.options)) {
    if(is.character(pn.options) == FALSE) stop("character variable name required for pn.options") 
    			}
    ckfmtnum <- list(first.y[1], x.at.first.y[1], last.y[1], x.at.last.y[1])
    names(ckfmtnum) <- c("first.y", "x.at.first.y", "last.y", "x.at.last.y")
    fun1<-function(x){!is.na(x)}
    ts<-lapply(ckfmtnum,fun1)
    ckfmtnum <- ckfmtnum[unlist(ts)]
    if( length(ckfmtnum) > 0 ) {
    fun2<-function(x){!is.numeric(x)}
    tstfmtnum <- lapply(ckfmtnum,fun2)
    ckfmtnum <- ckfmtnum[unlist(tstfmtnum)]
    stop(paste("argument ", names(ckfmtnum[1]), " is not a number", sep ="") )
    }  
    options(warn = -1)
    if(exists("Envir", mode= "environment") == FALSE) stop("Envir must be a
    valid R environment (e.g. not a charater variable")
    savvalue<-NA
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
    FPCEnv$.paramsestimated <- FALSE
    if(!is.na(pn.options)) {
    pnopt<- as.character( pn.options[1] ) 
    } else {
    pnopt<- ".pntemplist"
    }
    pnoptnm <- pnopt
    if(!is.na(twocomponent.x) & force4par == TRUE) 
    	stop("Cannot force a two component Richards model to have a single component.Set force4par to FALSE")
    prntqut<-function(supp, tx) if(supp == FALSE) print(noquote(tx))
    prntqut( suppress.text, "modpar will attempt to parameterize your data using the following sequential procedures:")
    prntqut( suppress.text, "  (1) Extract parameter estimates for 8-parameter double-Richards curve in nls")
    prntqut( suppress.text, "  (2) Use getInitial to retrieve parameter estimates for 8-parameter double-Richards curve")
    prntqut( suppress.text, "  (3) Extract parameter estimates for 4-parameter Richards curve in ")
    prntqut( suppress.text, "  (4) Use getInitial to retrieve parameter estimates for 4-parameter Richards curve")
    prntqut( suppress.text, "if any approaches are successful, modpar will return these and terminate at that stage")
    prntqut( suppress.text, " ")
    detl <- TRUE
    if(verbose == TRUE) detl <- FALSE
    xy <- sortedXyData(x,y)
    skel <- rep(list(1), 16)
    initval <- c(rep(NA, 16))
    initval <- relist(initval, skel)
    names(initval) <- c("Asym", "K", "Infl", "M", "RAsym", "Rk",
        "Ri", "RM", "first.y", "x.at.first.y", "last.y", "x.at.last.y",
        "twocomponent.x","verbose","force4par", "force.nonmonotonic")
    initval$first.y <- first.y[1]
    initval$x.at.first.y <- x.at.first.y[1]
    initval$last.y <- last.y[1]
    initval$x.at.last.y <- x.at.last.y[1]
    initval$twocomponent.x <- twocomponent.x
    if(width.bounds != 1) initval$width.bounds <- width.bounds
    if(bounds.error == TRUE) initval$bounds.error <- bounds.error
    initval$verbose <- verbose
    initval$force4par <- force4par
    initval$force.nonmonotonic <- force.nonmonotonic
    initval$taper.ends <- taper.ends
    skel1 <- rep(list(1), 16)
    initval1 <- c(rep(NA, 16))
    initval1 <- relist(initval1, skel1)
    names(initval1) <- c("Amin", "Amax", "Kmin", "Kmax", "Imin",
        "Imax", "Mmin", "Mmax", "RAmin", "RAmax", "Rkmin", "Rkmax",
        "Rimin", "Rimax", "RMmin", "RMmax")
    formtassign <- function (x,y,type=0) {
    	options(warn=-1)
    	valexp <- sapply( c( unlist(x[1:16]),unlist(y) ) 
   		 , function(x) list(x))
    	names( valexp[1:32] ) <- names( c(x[1:16],y) )
    	valexp[c(1:12,17:32)] <- as.numeric( valexp[c(1:12,17:32)] )
    	valexp[14:16] <- as.logical( valexp[14:16] )
    	if(is.numeric(x$twocomponent.x)) {
	        valexp[13] <- as.numeric( valexp[13] )
	    } else {
	  	valexp[13] <- as.logical( valexp[13] )     	
   	}
   	if(type == 1) {
   		if(names( x[length(x)]) == "taper.ends") {
   		valexp[17] <- x$taper.ends
   		names(valexp[17]) <-"taper.ends"
   		return(valexp[1:17])
   		}else{
    		return(valexp[1:16])
    		}
    	} else {
    		if(type == 2) {
    			  if(names( x[length(x)]) == "taper.ends") {
			  valexp$taper.ends <- x$taper.ends
			  return(valexp[17:33])
			   }else{
			  return(valexp[17:32])
    			  }
    		}  else {
    			 if(names( x[length(x)]) == "taper.ends")
    			 	valexp$taper.ends <- x$taper.ends
			 return(valexp)
    		}
    		      }
         options(warn=0)
    		      			 }
    parseval <- function(txt1,evtxt,txt2) {
    callout<- parse(text=sprintf("%s",paste(txt1,as.character(evtxt),txt2,sep="")))
    options(warn =-1)
    outpt<- eval(callout)
    options(warn =0)
    return(outpt)
    					  }  					  
    valexp <- formtassign(initval, initval1)
    valexp$taper.ends <- taper.ends
    if(width.bounds != 1) valexp$width.bounds <- width.bounds
    if(bounds.error == TRUE) valexp$bounds.error <- bounds.error
    valexp$modpar <- TRUE
    assign(pnoptnm, valexp, envir = Envir)
    FPCEnv$pnoptnm <- valexp
    evlfit<-function(val1,val2,force8par,savoptions){
 	richards <- function(x, Asym, K, Infl, M) Asym/Re(as.complex(1 +
    	    M * exp(-K * (x - Infl)))^(1/M))
	 SSposnegRichardsF <- function(x, Asym, K, Infl, M, RAsym,
  	      Rk, Ri, RM) (Asym/Re(as.complex(1 + M * exp(-K * (x -
  	      Infl)))^(1/M))) + (RAsym/Re(as.complex(1 + RM * exp(-Rk *
  	      (x - Ri)))^(1/RM)))
	if(is.na(val1[5])){
		evl1<- sum((y-richards(x,as.numeric(val1[1]),as.numeric(val1[2]),
               	 as.numeric(val1[3]),as.numeric(val1[4])))^2)
		evl2<- sum((y-richards(x,as.numeric(val2[1]),as.numeric(val2[2]),
               	 as.numeric(val2[3]),as.numeric(val2[4])))^2)

	}else{
		evl1<- sum((y-SSposnegRichardsF(as.numeric(x),as.numeric(val1[1]),as.numeric(val1[2]),
               	 as.numeric(val1[3]),as.numeric(val1[4]),
		as.numeric(val1[5]),as.numeric(val1[6]),
		as.numeric(val1[6]),as.numeric(val1[8])))^2)
		evl2<- sum((y-SSposnegRichardsF(x,as.numeric(val2[1]),as.numeric(val2[2]),
               	 as.numeric(val2[3]),as.numeric(val2[4]),
		as.numeric(val2[5]),as.numeric(val2[6]),
		as.numeric(val2[6]),as.numeric(val2[8])))^2)
	}
	if(evl1<=evl2) {
	valfin<-val1
	assign(pnoptnm, savoptions, envir = Envir)
	FPCEnv$pnoptnm <- savoptions
	}else{
	valfin<-val2
	}
	if(length(val2)<5 & length(val1) >4 & force8par == TRUE) valfin <- val1
	assign(pnoptnm, savoptions, envir = Envir)
	FPCEnv$pnoptnm <- savoptions
	return(valfin)
	}
    value <- NA
    succ <- FALSE
    if(force4par == TRUE & is.na(twocomponent.x)) {
     if(detl == FALSE) prntqut( suppress.text, "Estimating parameter bounds....")
      value <- try(parseval("try(getInitial(y ~ SSposnegRichards(x, Asym = Asym,
    	             K = K, Infl = Infl, M = M, modno = 19,  pn.options =",  pnoptnm,"), data = xy), silent = detl)"),
    	             silent = detl)
      savvalue<-value
       value <-NA
      if(class(value[1]) == "try-error") stop ("Bounds unestimable")
      savoptions <- get(pnoptnm, envir = Envir)
    prntqut( suppress.text, "(3) Status of 4-parameter Richards curve nls fit:")
       value <- try(parseval("coef(nls(y ~ SSposnegRichards(x, Asym = Asym,
                    K = K, Infl = Infl, M = M, modno = 12, pn.options =",  pnoptnm,"), data = xy, ...))")
         		, silent = detl)
         if(is.na(value[1]) == FALSE & class(value)[1] != "try-error") {
         prntqut( suppress.text, "4 parameter nls fit successful")
         }
         if(is.na(value[1]) == TRUE | class(value)[1] == "try-error") {
          prntqut( suppress.text, "....4-parameter nls fit failed")
          prntqut( suppress.text, "(4) Status of 4-parameter Richards getInitial call:")
         value <- try(parseval("try(getInitial(y ~ SSposnegRichards(x, Asym = Asym,
	             K = K, Infl = Infl, M = M, modno = 12,  pn.options =",  pnoptnm,"), data = xy), silent = detl)"),
	             silent = detl)
            if(is.na(value[1]) == TRUE | class(value)[1] == "try-error") {
            stop("estimates not available for data provided. Please check data, call or provide estimates manually, see ?modpar")
            FPCEnv$.paramsestimated <- FALSE
         					} else {					
         					prntqut( suppress.text, "....4 parameter getInitial successful")
         					}
         }
  if(!is.na(value[1]) & class(value)[1] != "try-error" & !is.na(savvalue[1]) & class(savvalue)[1] != "try-error") 
  	value<-evlfit(savvalue,value,force8par,savoptions)
    }else{
    if(detl == FALSE) prntqut( suppress.text, "Estimating parameter bounds....")
    value <- try(parseval("try(getInitial(y ~ SSposnegRichards(x, Asym = Asym,
                K = K, Infl = Infl, M = M, RAsym = RAsym, Rk = Rk,
            Ri = Ri, RM = RM, modno = 18,  pn.options =",  pnoptnm,"), data = xy), silent = detl)"),
            silent = detl)
     if (is.na(value[1]) == TRUE | class(value)[1] == "try-error"){
        value <- try(parseval("try(getInitial(y ~ SSposnegRichards(x, Asym = Asym,
         	             K = K, Infl = Infl, M = M, modno = 19,  pn.options =",  pnoptnm,"), data = xy), silent = detl)"),
         	             silent = detl)
             if(class(value[1]) == "try-error") stop ("Bounds unestimable")
            bndsvals<-get(pnoptnm, envir = Envir) 
            initval[1:8] <- value
            initval[13]<-bndsvals[13]
            bndsvals<-bndsvals [17:32]
            bndsvals[9:16]<-bndsvals[1:8]
            valexp <- formtassign( initval , bndsvals)				  
	    assign(pnoptnm, valexp, envir = Envir)
            			} else {     
      if(class(value[1]) == "try-error") stop ("Bounds unestimable")
    savvalue<-value
      value <-NA
      savoptions <- get(pnoptnm, envir = Envir)    							}
    if(force4par == TRUE & !is.na(twocomponent.x)) stop("Cannot force a two component model to have 4 parameters")
    prntqut( suppress.text, "(1) Status of 8-parameter double-Richards curve fit in nls:")
     value <- try(parseval("try(coef(nls(y ~ SSposnegRichards(x, Asym = Asym,
        K = K, Infl = Infl, M = M, RAsym = RAsym, Rk = Rk, Ri = Ri,
        RM = RM, modno = 1, pn.options =",  pnoptnm,"), data = xy, ...)), silent = detl)"),
        silent=detl)
    if (is.na(value[1]) == TRUE | class(value)[1] == "try-error") {
        prntqut( suppress.text, "....8 parameter nls fit failed")
        prntqut( suppress.text, "(2) Status of 8-parameter double-Richards getInitial call")
        value <- try(parseval("try(getInitial(y ~ SSposnegRichards(x, Asym = Asym,
            K = K, Infl = Infl, M = M, RAsym = RAsym, Rk = Rk,
            Ri = Ri, RM = RM, modno = 1,  pn.options =",  pnoptnm,"), data = xy), silent = detl)"),
            silent = detl)
             if (is.na(value[1]) == FALSE & class(value)[1] != "try-error") {
             		succ <- TRUE
             		prntqut( suppress.text, "....8-parameter getInitial successful")
             		}
    } else {
        prntqut( suppress.text, "....8-parameter nls fit successful")
        succ <- TRUE
    }
    if(!is.na(value[1]) & class(value)[1] != "try-error" & !is.na(savvalue[1]) & class(savvalue)[1] != "try-error") 
    	value<-try(evlfit(savvalue,value,force8par,savoptions), silent = detl)
    if ((is.na(value[1]) == TRUE & is.na(twocomponent.x)) | (class(value)[1] == "try-error" & is.na(twocomponent.x)) ) {
        prntqut( suppress.text, "(3) Status of 4-parameter Richards curve nls fit:")
        prntqut( suppress.text, "if force8par==TRUE second curve parameters estimated as RAsym=Asym*0.05, Rk=K, Ri=Infl, RM=M")
         try({
            value <- try(parseval("coef(nls(y ~ SSposnegRichards(x, Asym = Asym,
                K = K, Infl = Infl, M = M, modno = 12, pn.options =",  pnoptnm,"), data = xy, ...))"),
                silent = detl)
            if (force8par == TRUE) {
                value <- c(value, -value[1] * 0.05, x-value[2],
                  value[3], value[4])
                names(value) <- c("Asym", "K", "Infl", "M", "RAsym",
                  "Rk", "Ri", "RM")
            }
        }, silent = detl)
        if(is.na(value[1]) == TRUE | class(value)[1] == "try-error") {
          prntqut( suppress.text, "....4-parameter nls fit failed")
          prntqut( suppress.text, "(4) Status of 4-parameter Richards getInitial call:")
          try({
 	    value <- try(parseval("getInitial(y ~ SSposnegRichards(x, Asym = Asym,
	             K = K, Infl = Infl, M = M, modno = 12,  pn.options = ", pnoptnm,"), data = xy)"), silent = detl)
	    if (force8par == TRUE) {
	        value <- c(value, -value[1] * 0.05, x - value[2],
	          value[3], value[4])
	        names(value) <- c("Asym", "K", "Infl", "M", "RAsym",
	        "Rk", "Ri", "RM")
	     }
          }, silent = detl)
          if(is.na(value[1]) == TRUE | class(value)[1] == "try-error") prntqut( suppress.text, "....4 parameter getInitial failed")
        } else {
        prntqut( suppress.text, "....4 parameter nls successful")
        }
        if(is.na(value[1]) == TRUE | class(value)[1] == "try-error")
            stop("**Estimates not available for data provided**. Please check data or provide estimates manually, see ?modpar")
    } else {
         if (succ != TRUE)
                  stop("**Estimates not available for data provided**. Please check data or provide estimates manually, see ?modpar")
    }
    }
    
    valexp<-get(pnoptnm, envir = Envir)
    assign(pnoptnm, valexp[1:length(valexp)-1], envir = Envir)
    FPCEnv$pnoptnm <- valexp[1:length(valexp)-1]
    valexp <- valexp[1:8]
    options(warn=-1)
    try(rm("1", envir = Envir), silent =T)
    try(rm("pnoptnm", envir = FPCEnv), silent =T)
    options(warn=0)
     return(valexp)
    }
, ex = function(){
        modpar(posneg.data$age,posneg.data$mass)

        modpar(posneg.data$age,posneg.data$mass)
        subdata<-subset(posneg.data, as.numeric(row.names (posneg.data) ) < 53)
        richardsR1.lis<-nlsList(mass~SSposnegRichards(age,Asym=Asym,K=K,
        	Infl=Infl,M=M,RAsym=RAsym,Rk=Rk,Ri=Ri,RM=RM,modno=1, pn.options = "myoptions")
                        ,data=subdata)

})
