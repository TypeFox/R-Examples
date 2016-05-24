nrm <-
function(reshOBJ,etastart=-0.1, ctrl=list())  
{
  call <- match.call()
  attr(call, "date") <- date()
  

  ## check
  RA <- qr(reshOBJ$Qmat)$rank
  if(RA < ncol(reshOBJ$Q))
  {
    stop(paste("Q matrix has only rank",RA," - but must have full column rank! \n"))  
  }
  
  
  #########################################
  ######### USER CONTROLS #################
  #########################################
  
  cont <- list(nodes=21, absrange=5, sigmaest=FALSE, exac=0.001, EMmax = 500, verbose=TRUE, NRmax=20, NRexac=0.01, Clist=NA, nonpar=FALSE, quads=NA)

  user_ctrlI <- match(names(ctrl),names(cont))
  if(any(is.na(user_ctrlI)))
      {
      notex <- names(ctrl)[is.na(user_ctrlI)]
      warning("The following options in ctrl do not exist: ", paste(notex,collapse=", "))
      ctrl       <- ctrl[!is.na(user_ctrlI)]
      user_ctrlI <- user_ctrlI[!is.na(user_ctrlI)]
      }
  
  cont[user_ctrlI] <- ctrl
  
  #-------------------------------------------------------------------
  
  ## generate starting values
  startOBJ <- startV_nrmMG(reshOBJ=reshOBJ,etastart=etastart,Clist=cont$Clist)
  ##generate quadrature nodes and weights
  if(all(is.na(cont$quads)) | ( !all(is.na(cont$quads)) & !cont$nonpar) )
  {
    quads <- quadIT(nodes=cont$nodes,absrange=cont$absrange,ngr=nlevels(reshOBJ$gr))
    wherefrom <- "automatically generated"
    if(!all(is.na(cont$quads))){warning("supplied quads are ignored because 'nonpar=FALSE'.")}
    
  } else {
    cquads(cont$quads) # check the quads
    quads <- cont$quads
    wherefrom <- "individual quads"
  }
  
  
  OLD     <- 0
  ZAEHL   <- 1
  PARS    <- startOBJ$ulstv
  NLev    <- nlevels(reshOBJ$gr)
  ESTlist <- list()
  mueERG  <- NA
  value0 <- 10000
  
  
  ##################################
  ## EM Algorithm  #################
  ##################################

  repeat
  {
    if(cont$verbose & (ZAEHL == 1 | NLev <= 1)){cat("\r Estep:",ZAEHL,"| Mstep:", ZAEHL -1,"\r")}
    
    # E - STEP
    erg_estep <-  Enrm(PARS,reshOBJ=reshOBJ,startOBJ=startOBJ,quads=quads,PREVinp=mueERG,nonpar=cont$nonpar)
        
    if(cont$nonpar) # nonpar estimation in any case (reference group is standardized 0,1)
    {
      quads <- quadIT(nodes=cont$nodes,absrange=cont$absrange,ngr=NLev, ergE=erg_estep)  
    }
    
    # M - STEP
    if(cont$verbose){cat("\r Estep:",ZAEHL,"| Mstep:", ZAEHL,"\r")}
    mPARS <- PARS
    
    #################################
    ## Newton Raphson Procedure  ####
    #################################
    
    for(i in 1:cont$NRmax)
      {
        fir1 <- de1nrm(mPARS,Estep=erg_estep,reshOBJ=reshOBJ,startOBJ=startOBJ,quads=quads)
        sec2 <- de2nrm(mPARS,Estep=erg_estep,reshOBJ=reshOBJ,startOBJ=startOBJ,quads=quads)
        
        newP <- mPARS - as.vector(fir1 %*% solve(sec2))
        
        if(sum(abs(newP - mPARS)) < cont$NRexac)
          {
          mPARS <- newP
          break 
          } else {
                  mPARS <- newP
                 }
      }


    
    ###### LIKELIHOOD BERECHNEN #######################
    value <- ZFnrm(mPARS,riqv_quer=erg_estep$riqv_querG,reshOBJ=reshOBJ,startOBJ=startOBJ,quads=quads)
    
          
    ###########################
    # FINISHED ? ##############
    ###########################
    
        if(abs(value0-value) <= cont$exac & ZAEHL > 2 | ZAEHL >= cont$EMmax)
          {
            if(cont$verbose){cat("\r Estep:",ZAEHL,"| Mstep:", ZAEHL," --> estimating EAP & co \r\n")}
            mueERG <- mueNRM(mPARS,reshOBJ=reshOBJ,startOBJ=startOBJ,quads=quads,sigmaest=cont$sigmaest,endest=TRUE)
            
            attr(quads,"wherefrom") <- wherefrom
            
            # convergence reached?
            if(abs(value0-value) <= cont$exac)
            {
              conv <- "convergence reached" 
            } else {
              conv <- "EM iteration limit reached! Increase EMmax."   
            }

            ESTlist[[1]]   <- mPARS
            ESTlist[[2]]   <- erg_estep
            ESTlist[[3]]   <- list(value=value, hessian=sec2)
            ESTlist[[4]]   <- ZAEHL
            ESTlist[[5]]   <- mueERG
            ESTlist[[6]]   <- quads
            ESTlist[[7]]   <- startOBJ
            ESTlist[[8]]   <- cont
            
            names(ESTlist) <- c("etapar","last_estep","last_mstep" ,"n_steps","erg_distr","QUAD","starting_values","ctrl")
            break
           }
    
    
    if(NLev > 1 & !cont$nonpar) # more than 1 group and NOT nonpar estimation
    {
      if(cont$verbose & ZAEHL >= 1){cat("\r Estep:",ZAEHL+1,"| Mstep:", ZAEHL,"\r")}  
      # estimate mu and sigma
      mueERG <- mueNRM(PARS,reshOBJ=reshOBJ,startOBJ=startOBJ,quads=quads,sigmaest=cont$sigmaest)
      # new quads
      quads <- quadIT(nodes=cont$nodes,absrange=cont$absrange,ngr=NLev,mu=mueERG$mean_est,sigma=mueERG$sig_est)
    }
    

    PARS <- mPARS
    ZAEHL <- ZAEHL + 1
    value0 <- value
  }  
  



  ## Person Parameters
  
  ##EAP_nrm <- PePNRM(quads=quads,mueERG=mueERG)
  ESTlist$EAPs <- mueERG$thetas
  ## centering
  
  if(all(!is.na(startOBJ$setC)))
  {
    
    etainclC <- vector(mode="numeric",length=ncol(reshOBJ$Qmat))
    etainclC[-startOBJ$setC$whichetas] <- ESTlist$etapar
    etainclC[startOBJ$setC$whichetas]  <- startOBJ$setC$whichconstant
    
    retrans1 <- as.vector(reshOBJ$Qmat %*% etainclC)
    
  } else 
      {
      retrans1 <- as.vector(reshOBJ$Qmat %*% ESTlist$etapar) 
      }
  
  
  ZLpar <- relist(retrans1, startOBJ$stwm1)
  
  ESTlist$ZLpar <- lapply(ZLpar,function(x)
          {
            names(x)  <- paste("Item",1:length(x),sep="")
            x
          }) 
  names(ESTlist$ZLpar) <- levels(reshOBJ$gr)
  
  #SE estimation ----------------
  
  if(all(!is.na(startOBJ$setC)))
  {
    
    comphess <- diag(reshOBJ$Qmat[,-startOBJ$setC$whichetas] %*% solve(ESTlist$last_mstep$hessian) %*% t(reshOBJ$Qmat[,-startOBJ$setC$whichetas]))
    comphesq <- sqrt(comphess*(-1))

    ## hier noch die constanten auf NA setzen!
    
    comphesq[rowSums(reshOBJ$Qmat[,startOBJ$setC$whichetas])  > 0] <- NA
    
    
  } else 
  {
  comphess <- diag(reshOBJ$Qmat %*% solve(ESTlist$last_mstep$hessian) %*% t(reshOBJ$Qmat))
  comphesq <- sqrt(comphess*(-1))
  #
  }
  
  notest <- which(rowSums(reshOBJ$Qmat) == 0)
  comphesq[notest] <- NA
  attr(comphesq,"listform")  <- relist(comphesq,startOBJ$stwm1)
  ESTlist$SE <- comphesq
  
  
  # also save the reshOBJ - and of course the queen
  ESTlist$reshOBJ   <- reshOBJ
  
  # -----------------------------
  # adding category information
  # -----------------------------
  ESTlist$Catinf <- Infnrm(ESTlist)
  
  
  cat(">>>>>>>",conv,"\n")
  attr(call,"convergence") <- conv
  
  
  ESTlist$call <- call
  
  class(ESTlist) <- "nrm"
  ## write out
  
  return(ESTlist)  
}
