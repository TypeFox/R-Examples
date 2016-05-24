index.batch <- function(
                        data,
                        labelnp,
                        ID,
                        Vcode = NULL,
                        olmeNB = NULL, #model fit,
                        subset = NULL, 
                        ## compute the index for a subset of patients when subset = a vector of patient ids.
                        ## by default (subset=NULL), compute the index for all patients.
                        ## The indecies of ID must agree to the indices of idat0$ID, returned from olmeNB 
                        qfun = "sum", # sum or max
                        IPRT = TRUE, #print control
                        i.se = TRUE, 
                        MC = FALSE, #for ar1 only
                        C = FALSE,
                        i.tol = 1E-75
                        ## If olmeNB  ==  NULL then this is required
                        ){
  if (nrow(data) !=length(labelnp)) stop("the length of labelnp does not agree with the number of rows in data")
  if (nrow(data) !=length(ID))      stop("the length of labelnp does not agree with the number of rows in data")
  if (C & olmeNB$RE=="NoN")         cat("\nC=TRUE option is not allowed for NoN. C=FALSE option will be used")
  if (C & i.tol < 1E-5) cat("i.tol=",i.tol,"could be too small for C=TRUE option!!")
  ## If the labelnp is TRUE/FALSE then the following code change them to 1/0
  ## If the labelnp is 1/0 then the following code keep them as 1/0
  labelnp[labelnp] <- 1; labelnp[labelnp] <- 0
  
  ftd <- formulaToDat(formula=olmeNB$formula,
                      data=data,ID=ID,
                      labelnp=labelnp,Vcode=Vcode)
  fulldata0 <- data.frame(ftd$dat)  # dat = (ID, Y, x1, x2, ...) numeric matrix
  lnp <- ftd$labelnp
  
  
  if (!olmeNB$AR) ftd$Vcode <- rep(1,length(fulldata0$ID))
  fulldat <- data.frame(ID=fulldata0$ID,Vcode=ftd$Vcode,CEL=fulldata0$CEL)
  ## If there is covariate, then add them in idat
  if (ncol( fulldata0 ) > 2) fulldat <- data.frame(fulldat,x=fulldata0[,3:ncol(fulldata0)])
  ## idat = dataframe(ID, Vcode, CEL, x.1, x.2, ...)

  if (is.null(subset)) fdID <- unique(fulldat$ID)  else fdID <- subset
  
  est <- as.data.frame(olmeNB$est)[,1] 
  Npat <- length(fdID)
  ## PP is output of this function
  PP <- matrix(NA, nrow=Npat, ncol=4)
  colnames(PP) <- c("p", "SE of logit(p_hat)", "low95", "up95")
  rownames(PP) <- unique(ftd$origID)

  if (i.se & olmeNB$RE == "NoN" )
    {
      if (is.null(olmeNB$bootstrap)) stop("if i.se == TRUE and olmeNB$RE == NoN then olmeNB$bootstrap must be supplied")
      PP.boot <- matrix(NA,nrow=Npat,ncol=length(olmeNB$bootstrap))
      rownames(PP.boot) <- unique(ftd$origID)
      colnames(PP.boot) <- paste("boot",1:length(olmeNB$bootstrap),sep="")
      colnames(PP)[2] <- "empirical SE(p)"
    }
  if (C & olmeNB$RE!="NoN"){
    
    IDnum <- getIDnum(Y=fulldata0$CEL,IDY=fulldata0$ID)
    tabID <- table(IDnum$mIDY)
    FirstScanOfPatient <- c(1,cumsum(table(IDnum$mIDY))+1)
    dif <- c(0,diff(ftd$Vcode))
    dif[FirstScanOfPatient] <- 0
    ## DT$ind contains the locations of the 1st repeated measures for each patient
    ## DT$diff  = 0 if its the first repeated measure and = 1 ifelse
    dif <- dif[-length(dif)]    # scan lag
    
    
    IDY <- ftd$dat[,1]
    Y <- ftd$dat[,2]

    if (ncol(ftd$dat) > 2){
      ZY <- cbind(1,ftd$dat[,3:ncol(ftd$dat),drop=FALSE])
    }else ZY <- matrix(1,nrow=nrow(ftd$dat),ncol=1)
    
    VcodeY <- ftd$Vcode
    lnpY <- ftd$labelnp
    
    reEst <- getEst(est=olmeNB$est,AR=olmeNB$AR)

    distRE <- ifelse(olmeNB$RE=="G",1,2) ## "G" <=> 1, "N" <=> 2
    ## point estimates
    typeSummary <- getSummaryType(qfun) 
    cpiall <- .Call("CPI_ALL",
                    as.double(Y), as.double(ZY), as.integer(IDnum$mIDY),as.double(reEst$alpha), 
                    as.double(reEst$theta),as.double(reEst$delta),as.double(reEst$beta),
                    as.double(distRE),as.integer(IDnum$maxni), as.integer(IDnum$NpatTot),
                    as.integer(lnpY), as.double(dif),as.integer(IPRT),as.double(typeSummary),as.double(i.tol)
                    )[[1]]
    PP[,1] <- cpiall
    ## Confidence interval
    if (i.se){
      if (IPRT) cat("\n compute SE...")
      for (ipat in 1 : IDnum$NpatTot)
        {
          if (IPRT) cat("\n computing SE of hat.p for patient",ipat,"...");
          probi <- cpiall[ipat]
          pick <- IDnum$mIDY %in% (ipat-1)
          Yi <- Y[pick]
          lnpYi <- lnpY[pick]
          ZYi <- ZY[pick,,drop=FALSE]
          diffi <- dif[pick]
          jac <- jacobian(func=logitCPIi, x=olmeNB$est[,1], Yi=Yi,ZYi=ZYi,distRE=distRE,
                          lnpYi=lnpYi,diffi=diffi,printing=FALSE,AR=olmeNB$AR,typeSummary=typeSummary,i.tol=i.tol)
          VarLogitHatP <- jac%*%olmeNB$V%*%t(jac)
          s <- sqrt(VarLogitHatP) ## In Logit scale
          ci <- pp.ci(c(probi,s))
          PP[ipat,2:4] <- c(s,ci)
        }
    }
    
  }else{
    if (IPRT ){
      if (olmeNB$RE != "NoN" & i.se==TRUE ) 
        { 
          cat("\nEstimated index (SE of logit(index)) [95%CILower,95%CIUpper]")
        }else{
          cat("\nEstimated index")
        }
    }
 
    for ( ipat in 1 : Npat ) ## 1 to N (the total number of patietns)
      {
        pick <- fulldat$ID==fdID[ipat]
        ## extract n_i by ncol(fulldat) observations, corresponding to the patient i
        idat <- fulldat[pick,,drop=FALSE] 
        if (nrow(idat)==0) next  #next if no data from this patient
        ilnp <- lnp[pick]
        if (! (1 %in% ilnp)) next  ## next if no new scans
        ##if (! (0 %in% ilnp)) next   ## next if no old scans !!!!!!CAUSUE problem for AR/semipara.

        xm <- NULL
        if ( ncol(idat)>3) xm <- as.matrix(idat[,-(1:3)])
        ## the component of ilnpTF = TRUE if the repeated measure is pre scans.
        ilnpTF <- (ilnp==0)       

        if (!olmeNB$AR)  ## independent model
          {
            ## total count on pre scans
            Y1 <- sum(idat$CEL[ilnpTF]) 
            ## if qfun = "sum" then Y2 = sum(idat$CEL[!ilnpTF]) else if qfun="max" then Y2= max(idat$CEL[!ilnpTF])
            Y2 <- eval(call(qfun, idat$CEL[!ilnpTF]))
            
            sn1 <- sum(ilnpTF)        ## number of pre scans
            sn2 <- sum(!ilnpTF)       ## number of new scans

            if (olmeNB$RE == "NoN"|i.se == FALSE){
              ## If nonparametric method is used for the random effects
              tem <- jCP(tpar=est, Y1=Y1, Y2=Y2, sn1=sn1, sn2=sn2, 
                         XM=xm, RE=olmeNB$RE, qfun=qfun, oth=olmeNB$gtb,
                         i.tol=i.tol)
              
              if(olmeNB$RE == "NoN" & i.se == TRUE){
                for (iboot in 1 : length(olmeNB$bootstrap))
                  {
                    olmeNB.boot <- olmeNB$bootstrap[[iboot]]
                    est.boot <- as.data.frame(olmeNB.boot$est)[,1]
                    PP.boot[ipat,iboot] <- jCP(tpar=est.boot, Y1=Y1, Y2=Y2, sn1=sn1, sn2=sn2, 
                                               XM=xm, RE=olmeNB.boot$RE, qfun=qfun, oth=olmeNB.boot$gtb,
                                               i.tol=i.tol)
                  }
              }
            }else{
              ## If distributional assumption was made for random effects
              tem <- CP.se(tpar=est, Y1=Y1,Y2=Y2, sn1=sn1, sn2=sn2,XM=xm,
                           V=olmeNB$V,RE=olmeNB$RE,qfun=qfun,i.tol=i.tol)
            }  
          }else{  ##AR(1) model
            
            ypre <- idat$CEL[ilnpTF]; ynew <- idat$CEL[!ilnpTF] ## vector of length # new, containing CEL
            stp <- c(0,diff(idat$Vcode))             
            Qf <- match.fun(qfun)
            newQ <- Qf(ynew, na.rm=TRUE)
            ## possible combinations for ynew under q()
            y2m <- getY2.1(newQ, sum(!is.na(ynew)), qfun)
            
            if (olmeNB$RE =="NoN"|i.se==FALSE) 
              {  
                tem <- jCP.ar1(tpar=est, ypre=ypre, ynew=ynew, y2m=y2m, stp=stp, XM=xm,   
                               RE=olmeNB$RE, MC=MC, qfun=qfun, oth=olmeNB$gi,i.tol=i.tol)
                
              }else{            
                tem <- CP.ar1.se(tpar=est,ypre=ypre,ynew=ynew,y2m=y2m, 
                                 XM=xm,stp=stp,RE=olmeNB$RE, V=olmeNB$V, MC=MC, qfun=qfun,i.tol=i.tol)
              }
          }
        
        if (i.se==FALSE){  
          tem1 <- rep(NA, 3)
        }else if (i.se==TRUE & olmeNB$RE == "NoN"){
          tem1 <- c(sd(PP.boot[ipat,]),quantile(PP.boot[ipat,],prob=c(0.025,0.975)))
        }else{
          tem1 <- pp.ci(tem)
        }
        PP[ipat,] <- c(tem, tem1)
        
        if (IPRT)
          {
            round.tem <- sprintf("%1.3f",tem)
            cat("\npatient", rownames(PP)[ipat], ":",round.tem[1])
            if (olmeNB$RE != "NoN" & i.se==TRUE) 
              { 
                round.tem1 <- sprintf("%1.3f",tem1) 
                cat(paste(" (",round.tem[2],")",
                    " [",round.tem1[1], ",", round.tem1[2], "]",sep=""))
              }
          }
      }
    
  }
  if (IPRT) cat("\n") 
  
  res <- list()
  res$condProbSummary <- PP
  res$para$CEL <- model.response(model.frame(formula=olmeNB$formula,data=data,na.action=function(z)z)) ## NA is kept
  res$para$labelnp <- labelnp ## original input which contains labelnp for corresponding NA
  res$para$ID <- ID ## original ID which contains ID for the corresponding NA
  res$para$qsum <- qfun
  
  if (i.se & olmeNB$RE == "NoN") res$bootstrap <- PP.boot
  class(res) <- "IndexBatch"
  return(res)
}

getEst <- function(est,AR)
  {
    if (is.matrix(est)) est <- est[,1]
    name <- rownames(est)
    alpha <- exp(est[1])
    theta <- exp(est[2])
    if (AR){
      Npara <- 3
      delta <- ilgt(est[3])
    }else{
      Npara <- 2
      delta <- 0
    }
    beta <- est[- (1:Npara)]
    return( list(alpha=alpha,theta=theta,delta=delta,beta=beta))
  }


logitCPIi <- function(est,Yi,ZYi,lnpYi,diffi,distRE,printing,AR,typeSummary,i.tol)
  {
    ## Return CPI for one patient in logit scale (unconstrained scale)
    EST<- getEst(est,AR)
    PP <- .Call("CPI_each",
                as.double(Yi),
                as.double(ZYi),  
                as.double(EST$alpha),
                as.double(EST$theta),
                as.double(EST$delta),
                as.double(EST$beta),  
                as.double(distRE), 
                as.integer(lnpYi), 
                as.double(diffi), 
                as.integer(printing),
                as.double(typeSummary),as.double(i.tol))
    return(lgt(PP))
  }




getSummaryType <- function(qfun){
  if (qfun == "sum"){ return(1);
                    }else if (qfun == "max"){ return(2);
                                            }else stop("qfun must be sum or max!")
}
getIDnum <- function(Y,IDY){
  IDnum <- list()
  IDnum$NpatTot <- length(unique(IDY))
  IDnum$maxni  <- max(c(tapply(Y,IDY,length)))
  IDnum$mIDY <- IDY - 1
  
  if (length(Y) == 0)
    {
      Y <- ZY <- 0;
      IDnum$mIDY <- -1000
    }  
  return(IDnum)
}



Pmax.non <-
  function(Y1=0,               #sum(ypre); ypre+ 
           Y2=1,               #max(ynew)
           u1=4.5,             #Ex(Ypre+)
           u2=c(1.5,1.5, 1.5), #Ex(Ynew), vector
           a=0.5, 
           gi=NULL)            #a vector (todo) allow a freq table [freq, valu]
  {
    pb=1/(gi*a+1)

    p=1
    if (Y2>0)
      { j2=j1=dnbinom(Y1, prob=pb, size=u1/a)
        for (i in 1:length(u2))
          { j2=pnbinom(Y2-1, prob=pb, size=u2[i]/a)*j2 } #pr(Ynew[i]<Y2)
        p=1-sum(j2)/sum(j1)
      }
    return(p)
  }



Psum.non <- function(Y1=0, Y2=1,     ## Y1=sum(y.pre), Y2=sum(y.new)
                     u1=1.5, u2=1.5, ## u1 = Ex(Y1); u2= Ex(Y2)
                     a=0.5, gi=NULL         ## a freqency table = [freq, g.value]
                     )
  {
    pb <- 1/(gi[,2]*a+1)
    p <- 1
    if ( Y2 > 0 )
      {
        ## j1 = Pr( Ypre = ypre | G = g[ib])*pi[ib], ib=1,...,Nb
        j1 <- dnbinom(Y1, prob=pb, size=u1/a)*gi[,1]
        ## j2 = a vector Pr( Yfol <= yfol - 1 | G = g[ib])* Pr( Ypre = ypre | G = g[ib])*pi[ib], ib=1,...,Nb
        ## Nb = The number of unique points for distribution of G
        ## pi is the proportion allocated to each unique point
        j2 <- pnbinom(Y2-1, prob=pb, size=u2/a)*j1 
        p <- 1-sum(j2)/sum(j1)
      }
    return(p)
  }


pp.ci <-
  function(x=c(0.01, 0.1), level=0.95, lg=TRUE)
  ## lg: logit transformation indicator
  ## x=c(phat, s.e.(phat)) if lg=F
  ## x=c(phat, s.e(logit(phat)) if lg=T
  ## level: confidence level
  ## x[1]: an estimate of conditional probability, p 
  ## x[2]: se (logit(phat)) if lg = TRUE
  {
    if (is.na(x[1])){ return(c(NA, NA))
                    }else if (x[1]==1) return(c(1,1))
                                        #lg=x[1]<0.2
    ll <- 0.5+level/2
    del <- qnorm(c(1-ll, ll))
    tem <- del*x[2]
    
    if (lg){
      return(ilgt(lgt(x[1])+tem))
    }else{
      tem1= x[1]+tem ##
      tem1[1] = max(tem1[1], 0)
      tem1[2] = min(tem1[2], 1)
      return(tem1)
    }
  }





ppv.ci <- function(x = rbind(c(1, 0), c(0.1, 0.11)), level=0.95, lg=T)
                                        #see pp.ci
                                        # x is a n by 2 matrix
  {   
    pp=cbind(x[,1], x[,1])
    ss=cbind(x[,2], x[,2])

    ll=0.5+level/2
    del=qnorm(c(1-ll, ll))
    tem=ss%*%diag(del)
    if (lg) 
      { r1 = pp*exp(tem)
        res=r1/(1-pp+r1)
      }else{ tem1= pp+tem 
             tem1[tem1[,1]<0,1] = 0
             tem1[tem1[,2]>1, 2] = 1
             res=tem1
           } #res[is.na(x[,1]), 1]=res[is.na(x[,1]), 2]= NA
                                        #res[x[,1]==1,2] = res[x[,1]==1,2] =1
    return(res)
  }

## pmarg.gauss.fun <-
##   function(y=1, u=1.5, a=0.5, th=3, RE="G", gau.n=32)
##   {   ainv=1/a 
##       uVa=u/a

##       if (RE=="G")  
##         { sc=th 
##           sh=1/th
##           GAU <- gauss.quad.prob(n=gau.n, RE="gamma", alpha=sh, beta=sc)
##         }
##       if (RE=="N")  
##         { tem=log(th+1)
##           u.ln=-tem/2
##           s.ln=sqrt(tem)
##           GAU <- gauss.quad.prob(n=gau.n,RE="normal", mu=u.ln, sigma=s.ln)
##           GAU$nodes=exp(GAU$nodes)   
##         }

##       p=GAU$nodes/(GAU$nodes+ainv)
      
##       val=sum(pnbinom(y, size=uVa, prob=1-p)*GAU$weights)
##       return(val)
##     }

