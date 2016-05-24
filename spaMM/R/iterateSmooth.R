iterateSEMSmooth <- function(anyHLCor_obj_args, ## contains $processed 
                             LowUp, init.corrHLfit, #preprocess.formal.args, 
                            control.corrHLfit,verbose=interactive(),
                            MAX ## for diagnostic plots, list(rho, nu); where rho typically named num(1) but def'd as expanded. 
                            ) {
  eval_smoothtest <- function(Krigobj) { # perform LRT on the smoothing parameters...
    smoothtest <- as.list(attr(Krigobj,"HLCorcall"))
    smoothrho <- smoothtest$ranPars$rho
    ## there is trRho or rho whether smoothing was performed or not ## FR->FR how to ensure info is in only one place ???  
    ## if both trRho and rho, trRho is used if HLCor call
    if (is.null(smoothrho)) {
      RHOMAX <- attr(Krigobj,"optimInfo")$RHOMAX
      testedvalue <- rhoFn(rhoInv(smoothtest$ranPars$trRho,RHOMAX)*2,RHOMAX)
      if (any(is.nan(testedvalue))) { ## *2 exceeds max value => no real smoothing
        smoothtest <- FALSE
      } else {
        smoothtest <- eval(as.call(smoothtest))
        smoothtest <- Krigobj$APHLs$p_bv> (smoothtest$APHLs$p_bv+1.92) ## test of information about rho_smooth ## FR->FR should p_bv be used here ? 
      }
    } else {
      smoothtest$ranPars$rho <- smoothrho*2
      smoothtest <- eval(as.call(smoothtest))
      smoothtest <- Krigobj$APHLs$p_bv> (smoothtest$APHLs$p_bv+1.92) ## test of information about rho_smooth
    } 
  }
  
  write_diagnostics <- function(Krigobj,comment=NULL) {
    cat(paste("iter=",it," ",signif(optr$value,4), "+/-",signif(printRMSE,4), "; n_points=",nrow(Krigobj$data),
              "; smooth.lambda=",signif(Krigobj$lambda,4),sep=""))
    residform <- deparse(attr(Krigobj$resid.predictor,"oriFormula")) ## FR->FR maybe define safe extractor for predictor objects ?
    if (residform != "1" || (! is.null(comment) && verbose)) cat("\n   ") 
    if (residform != "1") cat(paste("residual variance predictor:",residform,", and",Krigobj$resid.family$link,"link"))
    if (residform != "1" && (! is.null(comment) && verbose)) cat("; ") 
    if ( ! is.null(comment) && verbose) cat(comment)
  }
  
  plot_diagnostics <- function(smoothRho=unlist(Krigobj$corrPars$rho),info=NULL,nextpoints=NULL) {
    if (length(lower)==2L) {
      zut <- signif(smoothRho,4)
      if (length(zut)>1) {
        titlesub <- bquote(paste(rho[f(rho)],"=",.(zut[1]),", ",rho[f(nu)],"=",.(zut[2]),
                                 "; max=",.(signif(optr$value,4)),"; predRMSE=",.(signif(printRMSE,4))))
      } else titlesub <- bquote(paste(rho[smooth],"=",.(zut),"; max=",.(signif(optr$value,4)),"; predRMSE=",.(signif(printRMSE,4))))
      SEMdiagnosticPlot2D(Krigobj, MAX=MAX, smoothingOK=smoothingOK, titlemain=titlemain, titlesub=titlesub, 
                          nextpoints=nextpoints, 
                          info=info, ## used only if (smoothingOK)
                          optrPar=optr$par)
      SEMdiagnosticPlot(Krigobj,MAX=MAX,"Raw profiles",optr) ## as the title says
      ## it would be nice to have contour lines on a SEMdiagnosticPlot2D -> spaMMplot2D but this requires a grid of values+ smoothing as in spaMM.filled.contour   
    } else { SEMdiagnosticPlot(Krigobj,MAX=MAX, titlemain=titlemain, optr) }
  } 

  pargrid <- sampleGridFromLowUp(LowUp,n=init.corrHLfit$nSmoothed) ## n may be NULL
  ## bits of codes needed whether PQL is run or not
  prevPredVars <- 0
  lower <- LowUp$lower ## list ! which elements may have length >1 !
  upper <- LowUp$upper ## list !
  if (interactive()) {
    predi <- getProcessed(anyHLCor_obj_args$processed,"predictor",from=1L)
    Xpv <- getProcessed(anyHLCor_obj_args$processed,"X.pv",from=1L)
  }  
  #
  ## PQL block removed from version 1.7.20
  ## now the SEM computations
  allsmooths <- list(initSmooth=control.corrHLfit$initSmooth, ## NULL by default
                     resid.model=control.corrHLfit[["smooth.resid.model"]] ## NULL by default -> default controlled by optimthroughSmooth 
  )  
  control.smooth <- allsmooths ## distinction between what goes in allsmooths and others is important ! nrepl will vary
  control.smooth$nrepl <- 20L ## number of points for which replicate estimates of likelihood are computed (modified later)
  anyHLCor_obj_args$`HLCor.obj.value` <- "logLapp"
  arglist <- list(pargrid=pargrid, 
                  anyHLCor.args=anyHLCor_obj_args, # contains $processed
                  control.smooth=control.smooth)
  precision <- control.corrHLfit$precision
  if (is.null(precision)) precision <- 0.02
  EIfac <- control.corrHLfit$EIfac
  if (is.null(EIfac)) EIfac <- 1.96
  #dit <- control.corrHLfit$dit ## NULL by default
  #if (is.null(dit)) dit <- 0 ## default: controls test predVar < prevPredVars[it-dit] for smoothingOK or not
  #prevPredVars <- rep(0,dit+1L)
  prevPredVars <- -Inf
  #
  it <- 1L
  optr <- do.call("optimthroughSmooth",arglist)  ## first SEM iteration, "it=1"   ############## CALL (with screen outputs)
  Krigobj <- optr$Krigobj
  predVar <- as.numeric(get_predVar(Krigobj,newdata=optr$par,variances=list(linPred=TRUE,disp=TRUE)))
  if (predVar<0) {
    predVar <- 2*precision ## affects sampleNextPars
    printRMSE <- NA
  } else {printRMSE <- sqrt(predVar)}
  if ( interactive() ) write_diagnostics(Krigobj=Krigobj)
  continue <- TRUE ## so there will always be a second iteration
  while ( continue ) { ## note that some SEM results have already been analyzed previous to the loop
    control.smooth <- allsmooths ## reinitialize optimthroughSmooth arguments with constant ones 
    prevPtls <- optr$forSmooth
    smoothRho <- unlist(Krigobj$corrPars$rho)
    smoothtest <- eval_smoothtest(Krigobj) ## tests whether some correlation structure has been detected
    tests <- c( predVar < prevPredVars[it], smoothtest )
    prevPredVars <- c(prevPredVars,predVar) ## after the test of predVar vs prevPredVars
    if ( smoothingOK <- all(tests) ){ 
      control.smooth$nrepl <- 0 ## only duplicate will be optr$par
      control.smooth$ranFix <- Krigobj$corrPars ## passing both nu and rho's
      sizes <- c(6,3) ## 6 point in simplex around max, 3 points by EI (+ nu bounds)
    } else {
      testMessages<-c( "high predvar",  "low LRT for scale parameters")
      comment <- paste(testMessages[!tests], collapse=" & ")
      sizes <- c(6,6) ## 6 point in simplex around max, 6 points by EI (+ nu bounds)
      ## To get rid of some possibly aberrant points that prevent good smoothing :
      prevPtls <- prevPtls[order(prevPtls$logLobj)[-c(1:2)],] ## FR->FR but aberrant points may not be the lowest... 
      control.smooth$ranFix <- Krigobj$corrPars["nu"] ## passing large original nu, always forcing smooth surfaces 
      control.smooth$nrepl <- ceiling(10/it - 0.0001)
    }
    nextpoints <- sampleNextPars(sizes=sizes,optr=optr,simplexExpand=0.95,D.resp=sqrt(predVar)/2) 
    if (smoothingOK) {
      info <- attr(nextpoints,"info") ## only used for diagnostic plot but removed by the following rbind:
      nextpoints <- rbind(nextpoints,optrpar=optr$par,optrpar=optr$par)
    } else {
      nextpoints <- rbind(nextpoints,optrpar=optr$par)
      ## need close pairs to estimate better the smoothing parameters
      ulower <- unlist(lower)
      uupper <- unlist(upper)
      epsilon <- (uupper-ulower)/1000
      ulower <- ulower+epsilon
      uupper <- uupper-epsilon ##useful for pmin, pmax 
      nearbypts <- sampleNearby(nextpoints,n=min(nrow(nextpoints),6),
                                stepsizes=(uupper-ulower)/(100*smoothRho),
                                margin=0.5) ## samples in a ring, not disc     
      ## FR->FR problem: nearbypts may extrapolate... particularly for small smoothRho. We correct:
      for (ii in seq_len(length(ulower))) {
        nearbypts[,ii] <- pmax(nearbypts[,ii],ulower[ii])
        nearbypts[,ii] <- pmin(nearbypts[,ii],uupper[ii])
      }
      nextpoints <- rbind(nextpoints,nearbypts)
    }
    ## and a bit of extrapolation
    #       if (it>1) {
    #         cS <- connectedSets(info$simplicesTable)
    #         outerpoints <- lapply(cS, function(v){
    #           v <- intersect(v,info$innerVertexIndices) ## only the really good points in the set
    #           pts <- info$vertices[v,,drop=FALSE]
    #           if (nrow(pts)>length(lower)+1) { ## more vertices than a simplex => can be redundant
    #             return(pts[unique(as.vector(convhulln(info$vertices[v,],"Pp"))),])
    #           } else return(pts) ## extrapolhull will handle special cases
    #         })
    #         extrap <- lapply(outerpoints,extrapolhull)
    #         extrap <- do.call(rbind,extrap)
    #         nextpoints <- rbind(nextpoints,extrap)
    #       }
    ##
    if (interactive() ) {
      titlemain <- bquote(paste(.(DEPARSE(predi)),", iter=",.(it)))
      if (nchar(eval(titlemain))>50) {
        titlemain <- bquote(paste(.(DEPARSE(nobarsNooffset(predi))),"+..., iter=",.(it)))
      }
      if (nchar(eval(titlemain))>50) {
        titlemain <- bquote(paste(.(substr(aschar,0,35)),"+... [length(",beta,")=",.(ncol(Xpv)),"], iter=",.(it)))
      }
      plot_diagnostics(smoothRho=smoothRho,info=info,nextpoints=nextpoints)
    }
    arglist <- list(pargrid=nextpoints,control.smooth=control.smooth,
                    anyHLCor.args=anyHLCor_obj_args, # contains $processed
                    prevPtls=prevPtls)
    it <- it+1L
    optr <- do.call("optimthroughSmooth",arglist) ## it>1    ############## CALL (with screen outputs)
    Krigobj <- optr$Krigobj
    predVar <- as.numeric(get_predVar(Krigobj,newdata=optr$par,variances=list(linPred=TRUE,disp=TRUE)))
    if (predVar<0) {
      predVar <- 2*precision ## affects sampleNextPars
      printRMSE <- NA
    } else {printRMSE <- sqrt(predVar)}
    if (interactive()) write_diagnostics(Krigobj=Krigobj,comment=comment)
    continue <- (it<3L || predVar > min(precision,prevPredVars[it])) 
  } ## end 'while' loop
  cat("\n") 
  if (interactive()) plot_diagnostics(info=NULL)
  attr(optr$value,"predVar") <- predVar
  return(optr)
}