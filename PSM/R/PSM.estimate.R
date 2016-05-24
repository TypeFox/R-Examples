`PSM.estimate` <-
function(Model,Data,Par,CI=FALSE,trace=0, control=NULL, fast=TRUE) {
  dimS <- length(Data)
  for(i in 1:dimS) {
    check <- ModelCheck(Model,Data[[i]],Par)
    if(!check$ok) {
      errmsg <- check$errmsg
      errmsg <- paste(errmsg, "- the error occured using data for individual", i)
      break
    }
  }
  if(!check$ok) stop(errmsg)
  Linear <- check$Linear

  if(trace>0)
    cat(ifelse( Linear, "* Linear model *\n", "* Non-linear model *\n"))
  
  # Check for fast option and Singular A
  if(Linear && fast) {
    #Get the ModelParameters
    tmp       <- Model$ModelPar(Par$Init)
    if( is.null(tmp$OMEGA) ) {
       # OMEGA IS NULL, but covariates can be present      
        tmpPhi <- Model$h( eta=NULL , theta=tmp$theta , covar=Data[[1]]$covar)                
    } else {
      # OMEGA is not null
      tmpDimEta <- dim(tmp$OMEGA)[1] 
      tmpPhi    <- Model$h( eta=rep(0, tmpDimEta) , theta=tmp$theta , covar=Data[[1]]$covar)
    }
    tmpMat    <- Model$Matrices(tmpPhi)
    matA  <- tmpMat$matA
    
    # Check for Model INPUT (U)
    if(is.null(Data[[1]][["U"]])) { #check if U exists.
      ModelHasInput <- FALSE
    } else if ( any(is.na(Data[[1]][["U"]])) ) { #check if it contains any NA
      ModelHasInput <- FALSE
    } else {
      ModelHasInput <- TRUE
    }
    tmpU     <- if( !ModelHasInput) { NA } else { Data[[1]][["U"]] }

    # Set Uk
    if(ModelHasInput) {
      Uk <- tmpU[,1,drop=FALSE]
    } else {Uk <- NA}
    
    tmpdimX <- nrow( Model$X0(Time=Data[[1]]$Time[1], phi=tmpPhi, U = Uk) )
    
    rankA <- qr(matA)$rank
    singA <- (rankA<tmpdimX)
    
    if( singA) {
      cat("Unable to use option \"fast\" on singular A matrix - Switching to fast=FALSE \n") 
      fast=FALSE 
    }
  }

  T0 <- proc.time()[3]

  if(!is.null(Par$LB)) {
    if(trace>1) cat("Using logit transformation of parameters \n")
    Par$Init <- logit(Par$Init,Par$LB,Par$UB) }

  # control
  if( is.null(control)) {
    # The user did not supply a control for the optimizer
    optimizer <- 'optim'
    control <- list(maxit=100, trace=trace, REPORT=1 )
    # put parameters on similar scale if bounds are missing
    if(is.null(Par$LB))
      control$parscale <- abs(Par$Init)+1e-3
  } else {
    # The _did_ supply control par
    optimizer <- control$optimizer
    if(is.null(optimizer))
      stop('Include choice of optimizer in control$optimizer')
    control <- control[-match('optimizer',names(control))]
  }

  switch( EXPR = optimizer,
         optim={
           if(trace>1)  cat( "Using Optimizer: \t optim\n")
           out <- optim(par = Par$Init, fn = APL.KF ,
                        gr = APL.KF.gr, method = "BFGS",
                        control = control, hessian = FALSE,
                        Model=Model, Pop.Data=Data, LB=Par$LB, UB=Par$UB,
                        GUIFlag=trace,fast=fast,Linear=Linear)
           NegLogL     <- out$value
           ParEstimate <- out$par
           #if(CI) Hess <- out$hessian
         },
         ucminf={
           if(trace>1)  cat( "Using Optimizer: \t ucminf\n")
           out <- ucminf(par = Par$Init, fn = APL.KF ,
                         gr = APL.KF.gr,
                         control = control, hessian = 0,
                         Model=Model, Pop.Data=Data, LB=Par$LB, UB=Par$UB,
                         GUIFlag=trace,fast=fast,Linear=Linear)
           NegLogL     <- out$value
           ParEstimate <- out$par
#         },
#         nlm={ 
#           if(trace>1)  cat( "Using Optimizer: \t nlm\n")
#           #typsize=Par$Init,stepmax=(.1*abs(Par$Init))+1e-3,          
#           out <- nlm(f=APL.KF, p=Par$Init, hessian=FALSE, print.level=trace,
#                      Model=Model, Pop.Data=Data, LB=Par$LB,
#                      UB=Par$UB, GUIFlag=trace,fast=fast,Linear=Linear,...)
#           NegLogL     <- out$minimum 
#           ParEstimate <- out$estimate
#           #if(CI) Hess <- out$hessian
         },
         stop( "Optimizer not recognized")      
         )
    
    
  if(!is.null(Par$LB)) {
    ParEstimate = invlogit(ParEstimate,Par$LB,Par$UB)
  }
  
  ci <- CI
  SD <- COR <- FALSE
  if(CI) {
    if(trace) cat( "Evaluating Hessian...\n")
    tmpfun <- function(THETA,Model,Data,THETAnames,GUIFlag) {
      names(THETA) <- THETAnames
      APL.KF(THETA,Model,Data,GUIFlag=GUIFlag)
    }
    Hess <- hessian(tmpfun,ParEstimate, method = "Richardson",
               Model=Model,Data=Data,GUIFlag=trace,THETAnames=names(Par$Init))
    COV <- solve(Hess)
    SD <- matrix(sqrt(diag(COV)),nrow=1)
    SDmat <- t(SD)%*%SD
    COR <- 1/SDmat*COV
    colnames(SD) <- colnames(COR) <- rownames(COR) <- names(Par$Init)

    ci <- matrix(c(
                   ParEstimate-1.96*sqrt(diag(solve(Hess))),
                   ParEstimate,
                   ParEstimate+1.96*sqrt(diag(solve(Hess)))),nrow=3,byrow=TRUE)
    rownames(ci) <- c("Lower CI95","MLE","Upper CI95")
    colnames(ci) <- names(Par$Init)
  } 
  
  totaltime = proc.time()[3]-T0
  
  if(trace>0) {
    minutes = floor(totaltime/60)
    tid <- paste("Runtime:  ",minutes,":",round(totaltime-60*minutes,2),"",sep="")
    print(tid,quote=FALSE)
  }

  list(NegLogL = NegLogL, THETA = ParEstimate, CI = ci, SD=SD, COR=COR, sec = totaltime, opt = out)

}

