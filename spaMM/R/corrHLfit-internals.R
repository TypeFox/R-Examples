dispFn <- function(x,xref=1) {log(x+xref)}
dispInv <- function(x,xref=1) {exp(x)-xref}
## for rho/trRho: rho=0 exactly is meaningful in adjacency model. rhoFn and rhoInv should be used only for other models.
rhoFn <- function(x,RHOMAX) {log(x/(RHOMAX-x))} ## rho should be constrained to <RHOMAX and trRho should diverge as rho approaches RHOMAX
rhoInv <- function(trRho,RHOMAX) {
  RHOMAX*exp(trRho)/(1+exp(trRho))
} 
nuFn <- function(nu,rho,NUMAX) {log(nu/(NUMAX-nu))} ## nu should be constrained to <NUMAX and trNu should diverge as nu approaches NUMAX
nuInv <- function(trNu,trRho,NUMAX) {NUMAX*exp(trNu)/(1+exp(trNu))}
## FR->FR rho/sqrt(nu) scaling => should be fixed nu and transformed rho to handle vectorial rho
## thus nu fns should be indep of nu and rho fns should be functions of nu ! :
if (FALSE) {
  futurerhoFn <- function(rho,nu) {
    rhosc <- rho/sqrt(nu)
    log(rhosc/(.spaMM.data$options$RHOMAX-rhosc))
  }
  futurerhoInv <- function(trRhosc,trNu) {
    nu <- .spaMM.data$options$NUMAX*exp(trNu)/(1+exp(trNu))
    rhosc <- .spaMM.data$options$RHOMAX*exp(trRhosc)/(1+exp(trRhosc))
    rhosc*sqrt(nu)
  }
}


makeLowerUpper <- function(canon.init, ## cf calls: ~ in user scale, must be a full list of relevant params
                           init.optim, ## ~in transformed scale
                           user.lower=list(),user.upper=list(),
                           corr.model="Matern",nbUnique,ranFix=list(),
                           control.dist=list(),
                           optim.scale, 
                           rhorange=NULL,
                           RHOMAX,NUMAX) {
  lower <- upper <- init.optim   ## init.optim not further used...
  if (corr.model %in% c("SAR_WWt","adjacency","ar1")) { ## adjacency model
    lower$rho <- user.lower$rho ## no transfo for adjacency model
    if (is.null(lower$rho)) lower$rho <- rhorange[1L] ## may remain NULL  
    upper$rho <- user.upper$rho ## no transfo again
    if (is.null(upper$rho)) upper$rho <- rhorange[2L]
  } else {
    if (corr.model=="AR1") {
      if ( ! is.null(canon.init$ARphi)) {
        ARphi <- user.lower$ARphi
        if (is.null(ARphi)) ARphi <- -0.999999
        lower$ARphi <- ARphi
        ARphi <- user.upper$ARphi
        if (is.null(ARphi)) ARphi <- 0.999999
        upper$ARphi <- ARphi
      }    
    } else { ## then Matern model....
      if (! is.null(canon.init$rho)) {
        rho <- user.lower$rho
        if (is.null(rho)) rho <- canon.init$rho/150
        if (optim.scale=="transformed") {
          lower$trRho <- rhoFn(rho,RHOMAX=RHOMAX)
        } else lower$rho <- rho
        rho <- user.upper$rho
        if (is.null(rho)) {
          if (inherits(nbUnique,"list")) nbUnique <- mean(unlist(nbUnique))
          rho <- canon.init$rho*2*nbUnique ## The following was a bit too low for experiments with nu=0.5 : 1/(maxrange/(2*nbUnique)) ## nb => unique rows !
          ## *modify* upper rho so that it does not exceed RHOMAX => /($RHOMAX+...)
          if (optim.scale=="transformed") rho <- 2*rho * RHOMAX/(RHOMAX+2*rho)
        }
        if (optim.scale=="transformed") {
          upper$trRho <- rhoFn(rho,RHOMAX=RHOMAX) 
        } else upper$rho <- rho
        rhoForNu <- canon.init$rho
      } else rhoForNu <- getPar(ranFix,"rho")
      if (! is.null(canon.init$nu)) {
        nu <- user.lower$nu
        if (is.null(nu)) nu <- canon.init$nu/100
        if (optim.scale=="transformed") {
          lower$trNu <- nuFn(nu,rhoForNu,NUMAX)
          #print(c(rhoForNu,nu,lower$trNu))
        } else lower$nu <-nu
        nu <- user.upper$nu
        if (is.null(nu)) {
          if ( ! is.null(dm <- control.dist$`dist.method`) && dm %in% c("Geodesic","Earth")) {
            nu <- 0.5
          } else {
            ## constructs upper nu from NUMAX => /(1+...)
            ## nu should not diverge otherwise it will diverge in Bessel_lnKnu, whatever the transformation used
            nu <- NUMAX * canon.init$nu/(1+canon.init$nu) 
          }
        }
        if (optim.scale=="transformed") {
          upper$trNu <- nuFn(nu,rhoForNu,NUMAX)
        } else upper$nu <- nu
        #print(c(rhoForNu,nu,upper$trNu))
      }
    } 
    ##### common to the different models except adjacency (because there are several places where NUgget+adjacency is not handled)
    if ( ! is.null(canon.init$Nugget)) {
      lower$Nugget <- 0
      upper$Nugget <- 0.999999
    }
  }
  if (! is.null(canon.init$phi)) {
    phi <- user.lower$phi
    if (is.null(phi)) phi <- canon.init$phi/1000
    lower$trPhi <- dispFn(phi)
    phi <- user.upper$phi
    if (is.null(phi)) phi <- canon.init$phi*1000
    ## if phi is badly initialized then it gets a default which may cause hard to catch problems in the bootstrap...
    upper$trPhi <- dispFn(phi)
  }
  if (! is.null(canon.init$lambda)) {
    lambda <- user.lower$lambda
    if (is.null(lambda)) lambda <- canon.init$lambda/1000
    lower$trLambda <- dispFn(lambda)
    lambda <- user.upper$lambda
    if (is.null(lambda)) lambda <- canon.init$lambda*1000
    upper$trLambda <- dispFn(lambda)
  }
  return(list(lower=lower,upper=upper))
}

checkDistMatrix <- function(distMatrix,data,coordinates) {
  classDistm <- class(distMatrix)
  if ( ! classDistm %in% c("matrix","dist")) {
    message(paste("(!) 'distMatrix' argument appears to be a '",classDistm,"',",sep=""))
    stop("not a 'matrix' or 'dist'. Check the input. I exit.")
  }
  ## chol() fails on distances matrices with repeated locations (which are pos SD)... but chol() not used by default
  ## the following code assumes that distMatrix deals only with unique locations, and checks this
  ## HENCE ******* distMatrix must refer to unique values of a grouping variable *********
  usernames <- rownames(distMatrix)
  checknames <- all(sapply(usernames,function(v) {v %in% rownames(data)})) ## 
  if (!checknames) {
    warning("The rownames of 'distMatrix' are not rownames of the 'data'. Further checking of 'distMatrix' is not possible.")
    nbUnique <- NA
  } else {
    uniqueGeo <- calcUniqueGeo(data=data[usernames,coordinates,drop=FALSE]) ## check that this corresponds to unique locations
    nbUnique <- nrow(uniqueGeo)
    if (nbUnique != nrow(distMatrix)) {
      stop("The dimension of 'distMatrix' does not match the number of levels of the grouping variable")
    } else { ## check order
      redondGeo <- data[,coordinates,drop=F]
      designRU <- apply(redondGeo,1,function(v) {which(apply(v==t(uniqueGeo),2,all))}) ## has no names
      ## eg 1 1 2 2 3 2 3 4 is valid for 8 obs, 4 unique locations
      designRU <- unique(as.vector(designRU)) ## should then be 1 2 3 4
      ## but if distMatrix in reverse order, the first row of redondGeo would match the 4th of uniqueGeo and then the following test is FALSE:
      if ( ! all (designRU==seq_len(length(designRU))) ) {
        stop("The rows of 'distMatrix' are not ordered as rows of the 'data'.")
      }
    } 
  }
  nbUnique ## if stop() did not occur
}

makeCheckGeoMatrices <- function(data,distMatrix=NULL,uniqueGeo=NULL,coordinates,dist.method="Euclidean") {
  isListData <- inherits(data,"list")
  if (is.null(distMatrix)) { 
    if ( is.null(uniqueGeo) ) { ## then construct it from the data ## this should be the routine case
      if (isListData) {
        uniqueGeo <- lapply(data,function(dd) {calcUniqueGeo(data=dd[,coordinates,drop=FALSE])})
        nbUnique <- lapply(uniqueGeo,nrow) 
      } else {
        uniqueGeo <- calcUniqueGeo(data=data[,coordinates,drop=FALSE])
        nbUnique <- nrow(uniqueGeo) 
      }
    } 
    ## (2): we need distMatrix *here* in all cases for the check
    if (isListData) {
      notnumericList <- lapply(uniqueGeo,function(uG) { ! unlist(lapply(uG,is.numeric))})
      if (any(unlist(notnumericList))) {
        firstwrong <- names(which(unlist(lapply(notnumericList,any))))[1L]
        stop(paste("In data list element '",firstwrong,"', variables ",paste(names(which(notnumericList[[firstwrong]])),collapse=" "),
                   " are not numeric, hence not suitable for computation of distance matrix.",sep=""))
      }
      distMatrix <- lapply(uniqueGeo,proxy::dist,method=dist.method)
    } else {
      notnumeric <- ! unlist(lapply(uniqueGeo,is.numeric))
      if (any(notnumeric)) stop(paste(paste(names(which(notnumeric)),collapse=" "),
                                      "are not numeric, hence not suitable for computation of distance matrix."))
      
      distMatrix <- proxy::dist(uniqueGeo,method=dist.method)
    }
  } else { ## there is a distMatrix, this is what will be used by HLCor
    if (isListData) {
      nbUnique <- lapply(seq_len(length(data)),function(dd) {checkDistMatrix(distMatrix,dd,coordinates)})
    } else nbUnique <- checkDistMatrix(distMatrix,data,coordinates)
    ## stops if problems, otherwise checkDistMatrix has no useful return value
  }
  return(list(nbUnique=nbUnique,uniqueGeo=uniqueGeo,distMatrix=distMatrix))
}




alternating <-function(init.optim,LowUp,anyOptim.args,maxIter,ranPars,HLCor.args,trace,Optimizer="L-BFGS-B",optimizers.args,maxcorners) {
  nam <- names(init.optim)
  if (any(c("trPhi","trLambda") %in% nam )) {
    mess <- pastefrom("Dispersion parameters non allowed in 'init.corrHLfit' with alternating algorithm.",prefix="(!) From ")
    stop(mess)
  }
  initcorr <- init.optim[nam %in% c("trRho","trNu","Nugget","ARphi")]
  HLfitLowUp <- LowUp
  HLfitLowUp$lower[c("trRho","trNu","Nugget","ARphi")] <- NULL
  HLfitLowUp$upper[c("trRho","trNu","Nugget","ARphi")] <- NULL
  corrLowUp <- LowUp
  corrLowUp$lower[c("trPhi","trLambda")] <- NULL
  corrLowUp$upper[c("trPhi","trLambda")] <- NULL
  anycorrOptim.args <- anyOptim.args
  iter <- 0
  conv <- 1
  currentLik <- -Inf
  while (iter < maxIter && conv > 1e-5 ) { ## if alternating: alternate HLCor and locoptim
    ranPars[names(initcorr)] <- initcorr
    attr(ranPars,"type")[names(initcorr)] <- "var"
    HLCor.args$ranPars <- ranPars ## variable locally
    oldLik <- currentLik
    if (is.character(trace$file)) {
      if(.spaMM.data$options$TRACE.UNLINK) unlink("HLCor.args.*.RData")
      zut <- paste(unlist(initcorr),collapse="")  
      save(HLCor.args,file=paste("HLCor.args.",zut,".RData",sep="")) ## for replicating the problem
    }
    givencorr <- do.call(HLCor,HLCor.args) ## optim disp and beta given corr param
    currentLik <- givencorr$APHLs$p_v ## iterations maximize p_v
    conv <- currentLik-oldLik
    anycorrOptim.args$ranPars$lambda <- givencorr$lambda
    anycorrOptim.args$ranPars$phi <- givencorr$phi
    anycorrOptim.args$etaFix <- list(beta=givencorr$fixef) ## LeeN01sm imply that eta should be fixed too 
    loclist <- list(initcorr,corrLowUp,anyObjfnCall.args=anycorrOptim.args,trace,Optimizer=Optimizer,
                    optimizers.args=optimizers.args,maxcorners=maxcorners,maximize=TRUE) 
    initcorr <- do.call(locoptim,loclist) 
    iter <- iter+1
  }
  optPars <- c(initcorr,givencorr$lambda,givencorr$phi)
  optPars
}

scalefnfn <- function(var,varMAX) {
  switch(var,
         trRho=function(x) rhoFn(x=x,RHOMAX=varMAX),
         trNu=function(x) nuFn(x,NUMAX=varMAX),
         trLambda=function(x) dispFn(x),
         function(x) x) 
}

invscalefnfn <- function(var,varMAX) {
  switch(var,
         trRho=function(x) rhoInv(x,RHOMAX=varMAX),
         trNu=function(x) nuInv(x,NUMAX=varMAX),
         trLambda=function(x) dispInv(x),
         function(x) x) 
}


## for diagnostic plots:
validRangefn <- function(st,varMAX) {switch(st,
                                            trRho=c(-Inf,RHOMAX=varMAX),
                                            trNu=c(-Inf,NUMAX=varMAX), ## 
                                            c(-Inf,Inf)
)}

labelfn <- function(st) {switch(st,
                                trRho=expression(rho),
                                trNu=expression(nu),
                                trLambda=expression(lambda),
                                st
)}

fullrho <- function(rho, coordinates,rho_mapping) {
  if ( is.null(rho_mapping) ) {
    if (length(rho)==1L ) rho <- rep(rho,length(coordinates))
    names(rho) <- coordinates
  } else {
    rho <- rho[rho_mapping]
    names(rho) <- names(rho_mapping)
  }  
  rho
}

# plots Krigobj$data$logLobj i.e. *not* the smoothing
SEMdiagnosticPlot2D <- function(Krigobj,MAX,smoothingOK,titlemain,titlesub, nextpoints, info, optrPar) {
  axesNames <- names(optrPar) 
  st <- axesNames[1L]
  ticks_x <- invscalefnfn(var=st,varMAX=MAX[[st]])(x=Krigobj$data[,st])
  tmp <- makeTicks(ticks_x, scalefn=scalefnfn(var=st,varMAX=MAX[[st]]),
                   axis=1L,logticks=TRUE,validRange=validRangefn(st=st,varMAX=MAX[[st]]))
  xlabs <- tmp$labels;xat <- tmp$at; #xph <- tmp$phantomat
  st <- axesNames[2L]
  ticks_y <- invscalefnfn(var=st,varMAX=MAX[[st]])(x=Krigobj$data[,st])
  tmp <- makeTicks(ticks_y, scalefn=scalefnfn(var=st,varMAX=MAX[[st]]),
                   axis=2L,logticks=TRUE,validRange=validRangefn(st=st,varMAX=MAX[[st]]))
  ylabs <- tmp$labels;yat <- tmp$at; #yph <- tmp$phantomat
  spaMMplot2D(Krigobj$data[,axesNames[1L]],Krigobj$data[,axesNames[2L]],Krigobj$data$logLobj,
              plot.title={
                title(main=titlemain,line=2,xlab=labelfn(axesNames[1L]),ylab=labelfn(axesNames[2L])) ## default cex.main=1.2, line ~1.7
                title(main=titlesub,line=0.8,cex.main=1.1) ## 
              },
              plot.axes={
                axis(1, at=xat, labels=xlabs)
                axis(2, at=yat, labels=ylabs)
              },
              decorations= {
                if ( ! is.null(nextpoints)) points(nextpoints,pch=15,cex=0.4)
                ## sampleNextPars() no longer returns the info attribute (easily put back in its code)
                if (smoothingOK && inherits(info,"volTriangulation")) apply(info$simplicesTable,1,function(v){
                  polygon(info$vertices[v,])
                });
                points(matrix(optrPar,nrow=1),pch="+",col="#008080")  ## col= colortools::complementary(spaMM.colors()[64])
              },
              map.asp=1/2 ## -(1-sqrt(5))/2 not wide enough
  )      
}

SEMdiagnosticPlot <- function(Krigobj,MAX,titlemain, optr) {
  nPredictors <- length(optr$par)
  intsqrt <- as.integer(sqrt(nPredictors))
  oldpar <- par(no.readonly=TRUE,mfrow=c(ceiling(nPredictors/intsqrt), intsqrt))
  for (st in names(optr$par)) {
    ticks_x <- invscalefnfn(var=st,varMAX=MAX[[st]])(x=Krigobj$data[,st])
    tmp <- makeTicks(ticks_x,scalefn=scalefnfn(st,varMAX=MAX[[st]]),axis=1,logticks=TRUE,
                     validRange=validRangefn(st,varMAX=MAX[[st]]))
    xlabs <- tmp$labels;xat <- tmp$at
    plot(Krigobj$data[,st],Krigobj$data$logLobj,
         axes=FALSE,frame=TRUE,xlab=labelfn(st),ylab="log(L)"
    )
    axis(1, at=xat, labels=xlabs)
    axis(2)
    points(optr$par[st],optr$value,pch=19,cex=2,col="red")  
  }
  title( titlemain, outer = TRUE )
  par(oldpar)                              
  invisible(NULL)
}

reformat_verbose <- function(verbose,For) {
  if (is.null(verbose)) verbose <- logical(0)
  if (is.list(verbose)) stop("is.list(verbose)")
  if (is.na(verbose["trace"])) verbose["trace"] <- FALSE
  if (is.na(verbose["SEM"])) verbose["SEM"] <- FALSE
  if (is.na(verbose["warn"])) verbose["warn"] <- switch(For, "HLCor" = TRUE, "corrHLfit" = FALSE)
  summstring <- paste(For,"Summary",sep="")
  verbose[summstring] <- verbose["summary"]
  if (is.na(verbose[summstring])) verbose[summstring] <- switch(For, "HLCor" = FALSE, "corrHLfit" = TRUE)
  verbose["summary"] <- FALSE ## this affects HLCor if reformat_verbose is called from corrHLfit, and HLfit if it is called from HLCor
  return(verbose)
}



calc_maxrange <- function(rho.size,distMatrix=NULL,uniqueGeo=NULL,rho_mapping,dist.method) {
  if (rho.size<2L) { ## can be 0 if no explicit rho in the input  
    if (inherits(distMatrix,"list")) {
      maxrange <- max(unlist(lapply(distMatrix,function(dd) max(c(-Inf,dd))))) ## les Inf to handle dist(0)...
      -min(unlist(lapply(distMatrix,function(dd) min(c( Inf,dd)))))
    } else maxrange <- max(distMatrix)-min(distMatrix)   
    if (maxrange<.Machine$double.eps) stop("Only one distance value: spatial model cannot be fitted.")
  } else { ## rho.size >1
    if (inherits(uniqueGeo,"list")) {
      maxrange <- lapply(unique(rho_mapping), function(idx) {
        ranges <- matrix(unlist(lapply(uniqueGeo,function(uu){
          if (nrow(uu)>1) {
            range(proxy::dist(uu[,rho_mapping==idx],method=dist.method))
          } else c(Inf,-Inf) ## encore des Inf to handle dist(0)...
        })),ncol=2)
        max(ranges[,2])-min(ranges[,1]) 
      })
    } else { ## single data set
      maxrange <- lapply(unique(rho_mapping), function(idx) {
        rng <- range(proxy::dist(uniqueGeo[,rho_mapping==idx],method=dist.method))
        rng[2]-rng[1] 
      })
    }  
    maxrange <- unlist(maxrange)
  }
  return(maxrange)
}

calc_init_dispPars <- function(init,init.optim,init.HLfit,ranFix) {
  ## does not modify init.HLfit, but keeps its original value. Also useful to keep ranFix for simple and safe coding
  if (is.null(getPar(ranFix,"lambda"))) { ## no ranFix$lambda: process init.optim
    init$lambda <- init.optim$lambda 
    if (!is.null(init$lambda)) {
      if (init$lambda<1e-4) init$lambda <- 1e-4
      init.optim$trLambda <- dispFn(init$lambda) 
      init.optim$lambda <- NULL
    }
  } else { ## ranFix$lambda present, do NOT put it in init.optim
    if (!is.null(init.optim$lambda)) stop("(!) Arguments 'ranFix$lambda' and 'init$lambda' conflict with each other.")  
  } 
  if (is.null(getPar(ranFix,"phi"))) {
    init$phi <- init.optim$phi 
    if (!is.null(init$phi)) {
      if (init$phi<1e-4) init$phi <- 1e-4
      init.optim$trPhi <- dispFn(init$phi)
      init.optim$phi <- NULL
    }
  } else {
    if (!is.null(init.optim$phi)) stop("(!) Arguments 'ranFix$phi' and 'init$phi' conflict with each other.")  
  } 
  return(list(init=init,init.optim=init.optim,init.HLfit=init.HLfit,ranFix=ranFix))
}

calc_inits_ARphi <- function(init,init.optim,init.HLfit,ranFix) {
  if (is.null(getPar(ranFix,"ARphi")) && (! is.numeric(init.HLfit$ARphi))) { 
    init$ARphi <- init.optim$ARphi 
    if (is.null(init$ARphi)) init$ARphi <- 0. 
    if (! is.null(init.HLfit$ARphi)) {
      init.HLfit$ARphi <- init$ARphi 
    } else {
      init.optim$ARphi <- init$ARphi
    }
  } 
  if ( (! is.null(init.HLfit$ARphi)) && (! is.numeric(init.HLfit$ARphi))) {
    init.HLfit$ARphi <- init.optim$ARphi
    init.optim$ARphi <- NULL
  }      
  return(list(init=init,init.optim=init.optim,init.HLfit=init.HLfit,ranFix=ranFix))
}

calc_inits_nu <- function(init,init.optim,init.HLfit,ranFix,control.dist,optim.scale,NUMAX) {
  if (is.null(getPar(ranFix,"nu")) && (! is.numeric(init.HLfit$nu))) { 
    init$nu <- init.optim$nu 
    if (is.null(init$nu)) {
      if ( ! is.null(dm <- control.dist$`dist.method`) && dm %in% c("Geodesic","Earth")) {
        init$nu <- 0.25  
      } else init$nu <- 0.5 
    }
    if (! is.null(init.HLfit$nu)) {
      init.HLfit$nu <- init$nu ## avant transformation
    } else {
      if (optim.scale=="transformed") {
        Fixrho <- getPar(ranFix,"rho")
        if (is.null(Fixrho)) { 
          init.optim$trNu <- nuFn(init$nu,init$rho,NUMAX=NUMAX) 
        } else init.optim$trNu <- nuFn(init$nu,Fixrho,NUMAX=NUMAX)
        init.optim$nu <- NULL
      } else init.optim$nu <- init$nu
    }
  } 
  if ( (! is.null(init.HLfit$nu)) && (! is.numeric(init.HLfit$nu))) {
    init.HLfit$nu <- init.optim$nu
    init.optim$nu <- NULL
  }
  return(list(init=init,init.optim=init.optim,init.HLfit=init.HLfit,ranFix=ranFix))
}

calc_inits_Auto_rho <- function(init,init.optim,init.HLfit,ranFix,rhorange) {
  if (is.null(getPar(ranFix,"rho")) && (! is.numeric(init.HLfit$rho))) { ## $rho NULL or NA
    if (is.null(init.optim$rho)) init.optim$rho <- mean(rhorange)
  }
  if ( (! is.null(init.HLfit$rho)) && (! is.numeric(init.HLfit$rho))) { ## init.HLfit$rho is NA
    init.HLfit$rho <- init.optim$rho
    init.optim$rho <- NULL 
  }
  return(list(init=init,init.optim=init.optim,init.HLfit=init.HLfit,ranFix=ranFix))
}


calc_inits_geostat_rho <- function(init,init.optim,init.HLfit,ranFix,maxrange,optim.scale,RHOMAX) {
  if (is.null(getPar(ranFix,"rho")) && (! is.numeric(init.HLfit$rho))) {
    init$rho <- init.optim$rho 
    if (is.null(init$rho)) {
      init$rho <- 30/(2*maxrange)
    } else if (any( narho <- is.na(init$rho))) init$rho[narho] <- 30/(2*maxrange[narho]) ## 05/2015 allows init.corrHLfit=list(rho=rep(NA,...
    if (! is.null(init.HLfit$rho)) { ## ! is.null, but ! is.numeric as tested above; but now it becomes numeric
      init.HLfit$rho <- init$rho ## avant transformation
    } else {
      if (optim.scale=="transformed") {
        init.optim$trRho <- rhoFn(init$rho,RHOMAX=RHOMAX) ## we're in Matern model here
        init.optim$rho <- NULL
      } else init.optim$rho <- init$rho
    }
  } 
  if ( (! is.null(init.HLfit$rho)) && (! is.numeric(init.HLfit$rho))) {
    init.HLfit$rho <- init.optim$rho
    init.optim$rho <- NULL
  }
  return(list(init=init,init.optim=init.optim,init.HLfit=init.HLfit,ranFix=ranFix))
}

calc_inits <- function(init.optim,init.HLfit,ranFix,corr.model,rhorange=NULL,maxrange=NULL,
                       optim.scale,control.dist,RHOMAX,NUMAX) { 
  inits <- as.list(match.call()[c("init.optim","init.HLfit","ranFix")])
  inits <- c(inits,list(init=list()))
  if (corr.model %in% c("Matern")) {
    arglist <- c(inits,list(maxrange=maxrange,optim.scale=optim.scale,RHOMAX=RHOMAX))
    inits <- do.call(calc_inits_geostat_rho,arglist)
    arglist <- c(inits,list(control.dist=control.dist, optim.scale=optim.scale, NUMAX=NUMAX))
    inits <- do.call(calc_inits_nu,arglist)
    # Nugget: remains NULL through all computations if init.optim$Nugget is NULL
    if (is.null(getPar(ranFix,"Nugget"))) { inits$init["Nugget"] <- init.optim$Nugget }  
  } else if ( corr.model  %in% c("SAR_WWt","adjacency","ar1") ) { 
    arglist <- c(inits,list(rhorange=rhorange))
    inits <- do.call(calc_inits_Auto_rho,arglist)
  } else if (corr.model=="AR1") {
    inits <- do.call(calc_inits_ARphi,inits)
  }  
  # phi, lambda
  inits <- do.call(calc_init_dispPars,inits)
  return(inits)
}


