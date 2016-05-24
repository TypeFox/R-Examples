lmeNB <- function(formula,data,ID,p.ini=NULL,IPRT=FALSE,AR=FALSE,
                  RE=c("G","N","NoN"),deps=1e-03,Vcode=NULL,C=FALSE,
                  i.tol=1e-7,
                  o.tol=sqrt(.Machine$double.eps),##1.e-3,
                  labelnp,
                  maxit=100,semi.boot=0,
                  u.low=0)##,replaceI=TRUE) #YZ: Jan 9, 2015
  {
    ## A Wrapper function to call
    ## fitParaIND
    ## fitSemiIND
    ## fitParaAR1
    ## fitSemiAR1
    if (i.tol < 1E-5 & C) cat("\ni.tol might be too small for C=TRUE option!!")
    if (length(RE)>1) RE <- RE[1]
    if (length(ID) != nrow(data)) stop("The length of Vcode must be the same as the nrow of data.")
    
    if (is.null(Vcode)){
      if (AR) stop("AR model requires Vcode object")
    }else{
      if (length(Vcode) != nrow(data)) stop("The length of Vcode must be the same as the nrow of data.")
    }
    
    if (C & RE != "NoN"){
      ## Adjust ID and Vcode so that they have the same length as the cleaned data
      ftd <- formulaToDat(formula=formula,data=data,ID=ID,Vcode=Vcode)
      dat <- ftd$dat
      Vcode <- ftd$Vcode
      ## dat = (ID, Y, x1, x2, ...) numeric matrix
      DT <- getDT(dat) ## list

      ## DT$dif is not used if AR=FALSE, the vector DT$dif has to be a vector of length sum_i ni
      ## but can hold any double values
      if (AR){
        DT$dif <- c(0,diff(Vcode))   
        DT$dif[DT$ind] <- 0
      }else{
        DT$dif <- rep(1,nrow(dat))
      }
      ## DT$ind contains the locations of the 1st repeated measures for each patient
      ## DT$diff  = 0 if its the first repeated measure and = 1 ifelse
      DT$dif <- DT$dif[1:DT$totN]    # scan lag


      if (AR){
        Npara <- 4## AR: log(alpha), log(theta), logit(delta), beta0,..
        Npara.name <- c("log(alpha)", "log(theta)", "logit(delta)","Intercept")
      }else{
        Npara <- 3
        Npara.name <- c("log(alpha)", "log(theta)","Intercept")
      }
      
      if (is.null(p.ini))
        {
          p.ini <- rep(0, Npara + DT$cn)
          p.ini[Npara] <- mean(DT$y) ## beta0 intercept  
        }else{
          if (length(p.ini) != Npara + DT$cn)
            stop("The length of p.ini does not agree to ",Npara,"+",DT$cn,"(=# covariates)")
          }
  
      if (IPRT){
       
        cat(paste("\n\n estimates:",paste(Npara.name,collapse=""),
                  "beta1,... and the negative of the log-likelihood",sep=""))
      }
      if (RE == "G"){
        REint <- 1
      }else if (RE == "N") REint <- 2
      
      tt <- optim(p.ini,  ## c(log(a), log(th), lgt(dt), b0, b1, ...)
                  llk.C, ##ar1.lk likelihood function for gamma/log-normal RE model 
                  hessian=TRUE,  
                  control=list(reltol=o.tol),
                  REint=REint,
                  DT=DT,#input data (output from getDT; lag)
                  IPRT=IPRT,AR=AR,absTol=i.tol
                  )
      nlk <- tt$value #neg likelihood
      
      vcm <- solve(tt$hessian)
      if (is.matrix(vcm)) colnames(vcm) <- rownames(vcm) <- c(Npara.name, DT$xnames)
      p.est <- cbind(tt$p, sqrt(diag(vcm)))
      row.names(p.est) <- c(Npara.name, DT$xnames)
      re <- list(opt=tt, nlk=nlk, V=vcm, est=p.est, RE=RE,##idat=data.frame(dat),
                 Vcode=Vcode,AR=AR,formula=formula)
      class(re) <- "LinearMixedEffectNBFreq"
      return(re)

      
    }else{
      
      if (!AR){
        
        if (RE=="G" || RE == "N"){
          return( fitParaIND(formula=formula,data=data,ID=ID,p.ini=p.ini,IPRT=IPRT,RE=RE))
        }else if (RE=="NoN"){

          re <- lmeNB.sp.ind(formula=formula,  data=data,       
                             ID=ID,  labelnp=labelnp,
                             p.ini=p.ini,  ## initial values for the parameters (log(a), log(th), b0, b1, ...)
                             IPRT=IPRT,      ## printing control
                             deps=deps,   ## stop iteration when max(|p.new - p.old|) < deps
                             maxit = maxit,
                             semi.boot=semi.boot, u.low=u.low)##,replaceI=replaceI) #YZ: Jan 9, 2015
          return(re)

        }else{
          stop("RE must be G, N or semipara!!")
        }
        
      }else{
        
        if (RE=="G" || RE == "N"){
          return(fitParaAR1(formula=formula, data=data, ID=ID, Vcode=Vcode,
                             p.ini=p.ini, IPRT = IPRT, RE = RE, 
                             i.tol = i.tol, o.tol = o.tol
                             ))
        }else if (RE=="NoN"){
          return(fitSemiAR1(formula=formula, data=data, ID=ID, Vcode=Vcode,
                              p.ini = p.ini, IPRT = IPRT, deps = deps, maxit=maxit))
        }else{
          stop("RE must be G, N or NoN!!")
        }
        
      }
    }
  }





print.LinearMixedEffectNBFreq <- function(x,...)
  {

    cat("\n -----Negative binomial mixed effect regression----- ")
    if (x$RE == "G") RE <- "gamma"
    else if (x$RE == "N") RE <- "log-normal"
    else if (x$RE == "NoN") RE <- "semi-parametric"

    depn <- ifelse(x$AR,"AR(1)","independent")
    cat("\n ---Random effect is from ",RE," distribution---")
    cat("\n The ",depn," correlation structure for the conditional dist'n of response given random effects.")

    cat("\n Formula: ");print(x$formula)
    cat("\n Estimates: \n")
    crit <- qnorm(0.975)
    if (ncol(x$est) ==4)
      {
        printForm <- data.frame(cbind(sprintf("%1.3f",x$est[,1]),sprintf("%1.3f",x$est[,2]),
                                     sprintf("%1.3f",x$est[,1]-crit*x$est[,2]),
                                     sprintf("%1.3f",x$est[,1]+crit*x$est[,2])))
        rownames(printForm) <-rownames(x$est)
        colnames(printForm) <-c("Value","Std.Error","lower CI","upper CI")
      }else{
        ## When RE="NoN" and semi.boot=0 then the estimated variance and 95% CIs are not returned.
        printForm <- x$est
        colnames(printForm) <- "Value"
      }
    print(printForm) 
    cat("\n Estimated covariance matrix: \n")
    print(x$V)
    cat("------------------------")
    if (x$RE != "NoN") cat("\n Log-likelihood",-x$nlk)
    cat("\n")
  }


formulaToDat <- function(formula, data, ID,labelnp=NULL,Vcode=NULL) 
{
  
  ## datamatrix exclude the covariate values corresponding to the NA CELs or NA selected covariates
  datamatrix <- model.matrix(object=formula,data=data)
  covNames <- colnames(datamatrix)
  if ("(Intercept)" %in% colnames(datamatrix) )
    {
      datamatrix <- datamatrix[,-1,drop=FALSE]
    }else{
      stop("A model without an intercept term is not accepted!!")
    }
  ## y exclude the covariate values corresponding to the NA CELs or NA selected covariates
  y <- model.response(model.frame(formula=formula,data=data,na.action=na.omit))

  ## Update ID so that length(ID) = nrow(datamatrix) 
  formulaNonMiss <- update(formula,    ~ . + NonMiss)
  dataNonMiss <- data; dataNonMiss$NonMiss <- 1:nrow(data)
  datamatNonMiss <- model.matrix(object=formulaNonMiss,data=dataNonMiss)
  pick.nonMiss <-  1:nrow(data) %in% datamatNonMiss[,colnames(datamatNonMiss)=="NonMiss"] 
  origID <- ID[pick.nonMiss]
  ## Update Vcode so that length(Vcode) = nrow(datamatrix)
  ## Vcode is passed only in the AR1 model funcitons
  if (!is.null(Vcode)) nMissVcode <- Vcode[pick.nonMiss] else nMissVcode <- NULL
  ## Update labelnp so that length(labelnp) = nrow(datamatrix)
  if (!is.null(labelnp)) nMisslabelnp <- labelnp[pick.nonMiss] else nMisslabelnp <- NULL
  
  upID <- as.character(origID)
  ## If ID is a character vector of length sum ni,
  ## it is modified to an integer vector, indicating the first appearing patient
  ## as 1, the second one as 2, and so on..
  uniID <- unique(upID)
  numID <- rep(NA,length(upID))
  for (i in 1 : length(uniID))
    {
      numID[upID == uniID[i]] <- i

      s <- sum(upID == uniID[i])
      for (j in 1 : s ){
        ## upID contains the original patientID-ivec where ivec = 1,...,ni
        ## upID must be a character string as the renamed object is a character string
        upID[upID == uniID[i]] <- paste(uniID[i],1:s,sep="--")
      }
    }
  dat <- cbind(ID=numID,CEL=y,datamatrix)
  rownames(dat) <- upID
  
  
  ##if (! is.null(labelnp) | !is.null(Vcode)){
  dat <-list(dat=dat,labelnp=nMisslabelnp,origID=origID,Vcode=nMissVcode)
  ##}

  return(dat) #dat=(ID, Y, x1, x2, ...)
}

getDT <- function(dat)             #dat=(ID, Y, x1, x2, ...)
{
  dat <- dat[order(dat[,1]),]
  ID <- dat[,1]
  ydat <- dat[,2]
  ncov <- ncol(dat)-2
   
  if (ncov>0){
    xdat <- as.matrix(dat[,3:ncol(dat)])
    xnames <- colnames(dat)[-(1:2)] 
    ##print(xnames)
  }else{
    xdat <- xnames <- NULL
  } 

  Ysum <- tapply(ydat, ID, sum) #total lesion count for each patient

  Ni <- table(ID) # number of scans for each patient
  totN <- sum(Ni) # total number scans of cohort
  IND <- c(0,cumsum(Ni))+1 #location of the 1st row for each patient
  N <- length(Ni) #number of patients
  
  return(list(id=ID, y=ydat, x=xdat, cn=ncov, ys=Ysum, np=N, totN=totN, ni=Ni,  ind=IND, xnames=xnames))
}






lmeNB.sp.ind <- function(formula,     ## an object of class "formula"
                         ## (or one that can be coerced to that class): a symbolic description of the model to be fitted.
                         data,        ## a data frame, list or environment (or object coercible by as.data.frame to a data frame)
                         ## containing the variables in the model.
                         ID,          ## a vector of length n*ni containing patient IDs of observations in data
                         labelnp,
                         p.ini=NULL,  ## initial values for the parameters (log(a), log(th), b0, b1, ...)
                         IPRT=TRUE,      ## printing control
                         deps=1E-4,   ## stop iteration when max(p.new - p.old) < deps
                         maxit = 1000,
                         semi.boot=500,
                         u.low=0)##,replaceI=TRUE)#YZ: Jan 9, 2015
  {
    uniID <- unique(ID)
    Npat <- length(uniID)
    re <- fitSemiIND(formula=formula,data=data,ID=ID,p.ini=p.ini,
                     IPRT=IPRT,deps=deps,maxit=maxit, 
                     u.low=u.low)##,replaceI=replaceI)#YZ: Jan 9, 2015


    if (semi.boot >0){
      par <- NULL
      re.boot <- list()
      
      Nfollowup <- tapply(labelnp,ID,function(x) sum(x==1) )
      tab.Nfollowup <- table(Nfollowup)
      uni.Nfollowup <- as.numeric(names(tab.Nfollowup))
      ## cat("tab.Nfollowup",tab.Nfollowup)
      ## cat("uni.Nfollowup",uni.Nfollowup)
      for (iboot in 1 : semi.boot)
        {
          if (IPRT & iboot%%100 == 0)cat("\n",iboot,"bootstraps are done...")
          ## To account for the varying follow-up times of the patients, the bootstrap
          ## sampling is stratified according to the follow-up time
          bootSamp <- NULL
          for (iNfup in 1 : dim(tab.Nfollowup))
            {
              Nfup.pick <- uni.Nfollowup[iNfup]
              freq.pick <- tab.Nfollowup[iNfup]
              bootSamp <- c(bootSamp,sample(uniID[Nfollowup==Nfup.pick],freq.pick,replace=TRUE))
            }
          ## bootSamp <- sample(uniID,Npat,replace=TRUE)
          boot.data <- boot.ID <- NULL
          for (ipatboot in 1 : Npat)
            {
              pick <- ID == bootSamp[ipatboot]
              boot.ID <- c(boot.ID,rep(ipatboot,sum(pick)))
              boot.data <- rbind(boot.data,data[pick,])
            }
          boot.data <- data.frame(boot.data)
          re.boot[[iboot]] <- fitSemiIND(formula=formula,data=boot.data,ID=boot.ID,
                                         IPRT=IPRT,deps=deps,maxit=maxit,u.low=u.low)##,replaceI=replaceI)
          par <- rbind(par,c(re.boot[[iboot]]$est))
        }
      re$V <- var(par)
      ## need to check this one!
      CI <- t(apply(par,2,quantile,prob=c(0.025,0.975)))
      re$est <- cbind(re$est,SE=sqrt(diag(re$V)),CI)
      rownames(re$V) <- colnames(re$V) <- rownames(re.boot[[iboot]]$est)
      re$bootstrap <- re.boot
    }
    return(re)
  }



fitSemiIND <-  function(formula,     ## an object of class "formula"
                        ## (or one that can be coerced to that class): a symbolic description of the model to be fitted.
                        data,        ## a data frame, list or environment (or object coercible by as.data.frame to a data frame)
                        ## containing the variables in the model.
                        ID,          ## a vector of length n*ni containing patient IDs of observations in data
                        p.ini=NULL,  ## initial values for the parameters (log(a), log(th), b0, b1, ...)
                        IPRT=TRUE,      ## printing control
                        deps=1.e-4,   ## stop iteration when max(p.new - p.old) < deps
                        maxit = 100,
			u.low = 0
                        ## YZ: lower limit for u when calculating weights (Jan 9, 2015)
                        ##,replaceI=TRUE
                        )
{
  dat <- formulaToDat(formula=formula,data=data,ID=ID)  # dat = (ID, Y, x1, x2, ...) numeric matrix
  DT <- getDT(dat$dat)  
  #DT=list(id, y, x, cn=ncov, ys=Ysum, np=N, totN=totN, ni=Ni, ind=IND, xnames=xnames))
  if (is.null(p.ini))
  {
    p.ini=rep(0, 3+DT$cn) 
    p.ini[3]=mean(DT$y)   
  }
  p.new <- p.ini
  Wh <- paraSave <- NULL; ##weightAll <- list()

  #initialize diffPara, ahat
  diffPara <- 1    #max(p.new-p.old)
  ahat <- 0
  counter <- 0 
  while(diffPara > deps & maxit > counter)
    {
      p.ini <- p.new           #reset the initial values at each iteration
      ## the initial values for the numerical method for the generalized least square are set to
      ## the inputted p.ini when counter = 0
      ## the solution to the least square (weight matrix is identity) when counter >= 1
      if (counter<2) b.ini <- p.ini[-(1:2)] #YZ: Jan 9, 2015
      ## ------------------------------------------------------- ##
      ## Main Step 1 : estimate b using Generalized Least Square
      ## ------------------------------------------------------- ##
      if (DT$cn>0) #with covariates
        { temB <- optim(b.ini,## initial values of coefficients b0,b1,...
                        estb, ## weighted residuals function
                        hessian=FALSE,dat=DT, ## output of getDT
                        Wh=Wh, ## list of N elements, containing ni by ni weight matrix for each patient
                        PRT=FALSE)
          bhat <- temB$par
        }else{         #without covariates
          tem <- optimize(estb, lower=-5, upper=5, dat=DT, Wh=Wh, PRT=F)
          bhat <- tem$min
        }
      ## ------------------------------------------------------------------- ##
      ## Main Step 2: estimate g_i i=1,...,N, g_i = w_i*(y_i+/mu_i+)+(1-w_i)
      ## ------------------------------------------------------------------- ##
      
      ## Substep 1: the estimation of mu_i+
      ##         hat.mu_i+ = exp( X_i^T bhat)
                                        # compute exp(b0+b1*x1+ ...)
      uh <- rep(exp(bhat[1]), DT$totN)
      if (DT$cn>0) {
        tem <- exp(DT$x%*%bhat[-1])
        uh <- tem*uh
      }
                                        #uhat=u.i+
      uhat <- as.vector(tapply(uh, DT$id, sum))
      
      ## Substep 2: the estimation of w_i = sqrt(var(G_i)/var(Y_i+/mu_i+))
      ## weights for gi
      if (ahat == 0) #first iteration, set wi=1
        {
          wi <- rep(1, DT$np)
        }else{
          gsig <- that+(1+ahat*(that+1))/uhat
          wi <- sqrt(that/gsig)   
        }
      
      gh0 <- DT$ys/uhat ## An estimate of y_i+/mu_i+ for each i
      
      gh <- wi*gh0 + 1 - wi ## An estimate of w_i*(y_i+/mu_i+)+(1-w_i) for each i
      
      if (IPRT) { ##print(summary(wi)) ; print(summary(gh0));
        cat("\n iteration:",counter)
        ## cat("\n The distribution of hat.g \n")
        ## print(summary(gh))
      }
      
      ## normalized gh so that the average is 1
      gh <- gh/mean(gh)
      gh <- as.vector(gh)
      
      ## frequency table for gh
      tem <- sort(round(gh,6))
      gh1 <- unique(tem)
      ghw <- table(tem)/DT$np
      ## the number of unique values of gh by 1 vector, containing the proportion of observations with that value
      gtb <- cbind(ghw, gh1) 
      rownames(gtb)=NULL
      ## ------------------------------------------------ ## 
      ## Main Step 3: estimate a using profile likelihood ##
      ## ------------------------------------------------ ##
      tt <- optimize(lk.a.fun, lower=-5, upper=5, dat=DT, Iprt=FALSE, gs=gtb,
                     uhat=uh ## a vector of length n*sn containing hat.mu_ij 
                     )
      ## lk.a.fun is a function of log(a)
      ## exponentiate back 
      ahat <- exp(tt$min)
      ## -------------------------------------- ## 
      ## Main Step 4: moment estimate of var(G) ##
      ## -------------------------------------- ## 
      that <- estSg3(DT, ui=uhat, uij=uh) ## uh is a vector of length length(Y)
      
      ## weights for GLS
      ## YZ: Jan 9, 2015
      Wh <- tapply(uh, DT$id, getWhs, th=that, a=ahat, u.low=u.low)
      
      p.new <- c(tt$min, log(that), bhat)
      paraSave <- rbind(paraSave,p.new)
      ##weightAll[[counter + 1]] <- Wh
      if (IPRT) {
        cat("\n log(alpha), log(theta), beta0, beta1, ...\n");
        print(p.new)
      }
      ## YZ (Jan 9, 2015): perhaps exclude log(theta) from comparison
      diffPara <- max(abs(p.new-p.ini))
      counter <- counter +1
      
      ## If p.new and p.ini contains Inf at the same spot, then
      ## diffPara = NaN
      whereInf <- (abs(p.new)==Inf)* (abs(p.new)==Inf)
      if ( sum(whereInf))
        {
          diffPara <- max(abs(p.new[whereInf]-p.ini[whereInf]))
        }
    }

  ##names(weightAll) <-  rownames(paraSave) <- paste("sim",1:nrow(paraSave),sep="")
  if (counter == maxit) warning("The maximum number of iterations occur!")
  p.new <- matrix(p.new,ncol=1)
  colnames(paraSave)  <- rownames(p.new) <- c("log_a", "log_var(G)", "(Intercept)", DT$xnames)
  
  re <- list(opt = tt, 
             diffPara=diffPara, ##YZ: 9 Jan 2015
             V = NULL, est = p.new, gtb = gtb,counter = counter,
             gi = as.vector(gh), RE = "NoN", 
             AR = FALSE,formula=formula,
             paraAll=paraSave
             ##,weightAll=weightAll
             )
  class(re) <- "LinearMixedEffectNBFreq"
  return(re)
}




estb <- function(b, # (b0, b1, ...)
                 dat, #data see fitParaIND
                 Wh=NULL, ## inverse of Var(Yi)
                 ## a list of N elements.
                 ## each element is ni by ni matrix of weights
                 ## If Wh=NULL then the weight matrix are all regarded as identity
                 PRT=TRUE)
{
  ## GLS compute the weighted sum of the residuals
  mu.i <- exp(b[1]) #mu.i =exp(b0)
  ## If the number of covariate is nonzero:
  if (dat$cn>0) 
  {
    residual <- exp(dat$x%*%b[-1])
    mu.i <- residual*mu.i #mu.i = exp(b0+b1*x1 + ...) n*sn by 1 vector
  }
  residual <- dat$y-mu.i
  ## If weights are not provided then  W_i = diag(dt$vn) identity
  ## so the weighted sum of residuals is
  ## sum_i (y_i-X_i^T beta)^T (y_i-X_i^T beta) = sum(residual^2)
  if (is.null(Wh)){ 
    val <- sum(residual^2)
  }else{
    val <- 0
    for (ipat in 1 : dat$np)
    {
      ## pick is a vector of length ni containing the locations of
      ## repeated measures of patient ipat
      weightMat <- Wh[[ipat]]$weightMat
      pick <- dat$ind[ipat]:(dat$ind[ipat+1]-1)
      val <- val+t(residual[pick])%*%weightMat%*%residual[pick]
    }
  }
  if (PRT) print(val)
  return(val)
}

estSg3 <-
function(
         dat,  ## data from getDT
         ui,   ## u.i+
         uij   ## u.ij, vector
         )
{ 
  ssy <- sum(dat$ys^2)-sum(dat$y^2) ## SS
  ssu <- sum(ui^2)-sum(uij^2)
  that <- ssy/ssu
  return(max(that-1,0)) 
}





lk.a.fun <- function(para=-0.5, #log(a)
                     dat, Iprt=FALSE, 
                     gs, ## gi in freqency table form
                     uhat    ## hat u.ij, same length as Y.ij
                     ) 
{
  ## psudo profile likelihood as a function of alpha
  ## Li(a; hatb, hatg, y)
  ## = 1/N sum_l=1^N prod_j=1^ni Pr(Y_ij=y_ij|G_l=g_l,hatb,a)
  ## = sum_l=1^N* prod_j=1^ni Pr(Y_ij=y_ij|G_l=g*_l,hatb,a)*prop(g*_l)
  ## = sum_l=1^N* prod_j=1^ni choose(y_ij+r_ij-1,y_ij) (1-p(g*_l,a))^r_ij p(g*_l,a)^{y_ij}*prop(g*_l)
  ## = [prod_j=1^ni choose(y_ij+r_ij-1,y_ij)] sum_l=1^N* (1-p(g*_l,a))^r_i+ p(g*_l,a)^{y_i+}*prop(g*_l)
  ## 
  ## log Li
  ## = [sum_j=1^ni log choose(y_ij+r_ij-1,y_ij)] + log [sum_l=1^N* (1-p(g*_l,a))^r_i+ p(g*_l,a)^{y_i+}*prop(g*_l)]
  ## 
  ## g*_l represents the unique value of the approximations
  ## prop(g*_l) represents the proportion of approximated g_l's that have the value g*_l
  ## N* represents the number of unique g*
  ## p = p(g*_l,a) = 1/(g*_l*a + 1)
  
  if (Iprt) cat(para)
  
  a <- exp(para)
  
  ainv <- 1/a  
  
  th2 <- uhat/a ## r_ij = exp(X_ij^T hat.beta)/alpha

  if (!all(is.finite(th2))) return(Inf) ## alpha is not in the acceptable range
  ## -log [sum_j=1^ni log choose(y_ij+r_ij-1,y_ij)] = sum( - lgamma(dat$y+th2) + lgamma(th2) + lgamma(dat$y))
  ##  but the last term is not necessary because it does not depend on a
  i.pos = dat$y>0 #YZ: Jan 9, 2015
  th2[th2<1.e-100] = 1.e-100
  nllk <- sum( - lgamma((dat$y+th2)[i.pos]) + lgamma(th2[i.pos]))
  us <- tapply(th2, dat$id, sum)
  
  p <- gs[,2]/(gs[,2]+ainv)
  ## please note that p = 1 - p in this manuscript
  ## gs[,2] is a vector of length # unique g_i containing the values of unique g_i
  
  q <- 1-p
  for (i in 1:dat$np)
    {
      tem <- sum(
                 p^dat$ys[i]*q^us[i]*gs[,1]
                 ## gs[,1] is a vector of length # unique g_i
                 ## containing the proportion of g_i
                 )
      nllk <- nllk-log(tem)
    }
  
  if (Iprt) cat(" nllk=", nllk, "\n")
  return(nllk)
}

getWhs <-
  function(
           u=c(1.5, 1.35, 1.2, 1.1, 1, 0.9), ##mu.ij, vector
           th=exp(1.3), ## var(G)
           a=exp(-0.5),  ## alpha
           u.low=0 ## YZ: Jan 9, 2015
           ##,replaceI=TRUE ## Yumi: Jan 13, 2015
           )
{
  ## Independent model has:
  ## Var(Y_ij) = mu_ij^2*var(G) + mu_ij*(1+(var(G)+1)*a)
  ## Cov(Y_ij,Y_ij') = mu_ij * mu_ij'*Var(G)
  pick.toosmall <- u<u.low
  u[pick.toosmall] <- u.low #YZ: Jan 9, 2015
  w1 <- u%*%t(u)*th
  if (length(u)==1){
    w2 <- (1+(th+1)*a)*u ## Modified by Yumi 2013 Mar8
  }else{
    w2 <- (1+(th+1)*a)*diag(u)
  }
  w12 <- w1 + w2
  ##cat("\nrcond",rcond(w12),"det",det(w12))
  if (rcond(w12) < 1E-6)
    {
      ## if Var(Y) is not invertible, then the weights of weighted least squares
      ## are replaced by the identity matrix
      ##if (replaceI) W <- diag(length(u)) else
      stop("the weight matrix is singular at weighted least square stage!! u.low is too small.")
    }else{
      W <- solve(w12)
    }
  ##return(W) ## inverse of Var(Yi)
  return(list(weightMatrix.i=W,out=pick.toosmall)) ## inverse of Var(Yi)
}







## mle.a1.fun <- function(dat, p.ini=NULL, IPRT=TRUE, deps=1e-4)##,replaceI=TRUE)
##   { 
##     DT=getDT(dat)  
##                                         #DT=(id, y, x, cn=ncov, ys=Ysum, np=N, totN=totN, ni=Ni, ind=IND, xnames=xnames))
    
##     if (is.null(p.ini))
##       { p.ini=rep(0, 3+DT$cn) 
##         p.ini[3]=mean(DT$y)
##       }
##     p.new=p.ini
##     W=NULL

##   dth=1 #max(p.new-p.ini)
##   while(dth > deps)
##   { p.ini=p.new
##     b.ini=p.ini[-(1:2)]
    
##     if (DT$cn>0) 
##     { temB=optim(b.ini, estb, hessian=T, dat=DT, Wh=W, PRT=F)
##       bhat=temB$par }
##     else
##     { tem=optimize(estb, lower=-5, upper=5, dat=DT, Wh=W, PRT=F)
##       bhat=tem$min
##     }

##     uh=rep(exp(bhat[1]), DT$totN)
##     if (DT$cn>0) { tem=exp(DT$x%*%bhat[-1])
##     uh=tem*uh }

##     uhat=as.vector(tapply(uh, DT$id, sum))
##     gh0=DT$ys/uhat

##     #for cases of gh0=0, uniformly between 0 to 1/2 of the smallest non zero value
##     gh_1=min(gh0[gh0>0])/2

##     ll=sum(gh0==0)
 
##     gh=gh0
##     ll1=ll%%10
##     ll2=ll%/%10
##     ll3=rep(1:0, c(ll1, 10-ll1))+ll2
##     gh[gh0==0] = rep(seq(0, gh_1, length=10), ll3)
   
##     gh=gh/mean(gh)
  
##     gh=as.vector(gh)

##     #frequency table
##     tem=sort(round(gh,6))
##     gh1=unique(tem)
##     ghw=table(tem)/DT$np
##     gtb=cbind(ghw, gh1)
##     rownames(gtb)=NULL

##     tt=optimize(lk.a.fun, lower=-5, upper=5, dat=DT, Iprt=F, gs=gtb, uhat=uh)
##     ahat=exp(tt$min)

##      that=estSg3(DT, uhat, uh)
##      W=tapply(uh, DT$id, getWhs, th=that, a=ahat)##,replaceI=replaceI)

##     p.new=c(tt$min, log(that), bhat)
##     if (IPRT) print(p.new)
##     dth=max(abs(p.new-p.ini))
##   }  

##   if (DT$cn>0) vcm=solve(temB$hessian)
##   names(p.new)=c("log_a", "log_th", "log_u0", DT$xnames)
##   return(list(opt=tt, vcm=vcm, est=p.new, gi=as.vector(gh), gtb=gtb, RE="NoN"))
## }




## ===========  predicting the values of RE  ============= 

## dbb <- function(x, N, u, v) {
##   beta(x+u, N-x+v)/beta(u,v)*choose(N,x)
## }

dbb <- function(x,N,al,bet,logTF=FALSE)
  {
    ## beta-binomial density
    temp = lbeta(x+al, N-x+bet)-lbeta(al,bet) + lchoose(N,x);
    if (logTF){
      return (temp);
    }else return (exp(temp));
  }

get.prob <- function(gY,alpha)1/(gY*alpha + 1)
densYijGivenYij_1AndGY <- function(Yik,Yik_1,
                                   sizeYik,sizeYik_1,
                                   delta,prob=NULL,gY,alpha,cdf=FALSE)
  {
    if (Yik <0) return(0)
    ## Yik | Yi(k-1), GY
    ## under AR(1) model
    if (is.null(prob)) prob <- get.prob(gY=gY,alpha=alpha)
    support.max <- min(Yik,Yik_1)
    fYij.Yik_1 <- 0
    sizeNB <- sizeYik - delta*sizeYik_1
    if (sizeNB <= 0) return(0)
    
    for (t in 0 : support.max )
      {
        if (cdf){
          temp <- pnbinom(Yik - t,size=sizeNB,prob=prob)
        }else{
          temp <- dnbinom(Yik - t,size=sizeNB,prob=prob)
        }
        fYij.Yik_1 <- fYij.Yik_1 +
          dbb(x = t,N = Yik_1,
              al = sizeYik_1*delta,
              bet = sizeYik_1*(1-delta),logTF=FALSE)*temp
        
      }
    return(fYij.Yik_1)
  }
dens_Yi.gY <- function(Yi,sizeYi,prob,ARY=FALSE,delta=NULL,Vcode=1:length(Yi),
                       logTF=FALSE)
  {
    if (ARY)
      {
        out <- dnbinom(Yi[1],size=sizeYi[1],prob=prob,log=TRUE)
        for (ik in 2 : length(Yi))
          {
            Yik <- Yi[ik]
            sizeYik <- sizeYi[ik]
            Yik_1 <- Yi[ik-1]
            sizeYik_1 <- sizeYi[ik-1]

            fYij.Yik_1 <- densYijGivenYij_1AndGY(Yik=Yik,Yik_1=Yik_1,
                                                 sizeYik=sizeYik,
                                                 sizeYik_1=sizeYik_1,
                                                 delta=delta,prob=prob)
            out <- out + log(fYij.Yik_1)
          }
      }else{
        out <- sum(dnbinom(Yi,size=sizeYi,prob=prob,log=TRUE))
      }
    if (logTF) return(out) else return(exp(out))
  }



int.denRE <- function(gs,Yis,Xis,bet,alpha,theta,AR=FALSE,RE="G",delta=NULL,Vcode=NULL)
  {
    temp <- rep(NA,length(gs))
    for (ig in 1 : length(gs))
      {
        g <- gs[ig]
        size <- exp(Xis%*%bet)/alpha
        prob <- 1/(g*alpha+1)
        if (prob == 0){
          temp[ig] <- 0
        }else{
          temp2 <- dens_Yi.gY(Yi=Yis,sizeYi=size,prob=prob,
                              ARY=AR,delta=delta,Vcode=Vcode)
          if (RE=="N"){
            dRE <- dlnorm(g,meanlog=-log(theta+1)/2,sdlog=sqrt(log(theta+1)))
          }else if (RE=="G"){
            dRE <- dgamma(g,shape=1/theta, scale=theta)
          }
          temp[ig] <- temp2*dRE
        } 
      }
    ##cat("\n",temp)
    return(temp)
  }
int.numRE <- function(gs,Yis,Xis,bet,alpha,theta,delta,Vcode,AR,RE, expG)
  {
    rr <- rep(NA,length(gs))
    for (ig in 1 : length(gs))
      {
        g <- g.temp <- gs[ig]
        if (expG) g.temp <- log(g)
        dens <- int.denRE(gs=g,Yis=Yis,Xis=Xis,bet=bet,RE=RE,
                          alpha=alpha,theta=theta,AR=AR,delta=delta,Vcode=Vcode)
        if (dens == 0) rr[ig] <- 0 else rr[ig] <- g.temp*dens
      }
    ##cat("\n\n",gs)
    ##cat("\n",rr)
    return(rr)
  }


RElmeNB <- function(theta,alpha,betas,delta,
                    formula,ID,Vcode=NULL,data,
                    AR,RE,rel.tol=.Machine$double.eps^0.8,
                    expG=FALSE)
  {
    ## Compute E(G|Yi)
    ## expG = TRUE then parametrization is prob=1/(exp(G)*alpha + 1 )
    ## else if expG = FALSE the parametrization is prob=1/(G*alpha + 1)
    ## par <- fM$par
    ## theta <- par$theta#getByName(obj=fM$par,name="theta")
    ## alpha <- par$alpha#getByName(obj=fM$par,name="alpha")
    ## betas <- par$betaY#getByName(obj=fM$par,name="betaY",all=TRUE)
    ## delta <- par$delta#getByName(obj=fM$par,name="delta")    
    ## temp <-  fM$data$CEL
    ## Y <- temp$resp
    ## ID <- temp$ID
    ## Z <- temp$Z
    ## Vcode <- temp$Vcode
    ftd <- formulaToDat(formula=formula,data=data,ID=ID,Vcode=Vcode)
    ID <- ftd$dat[,1]
    Y <- ftd$dat[,2]
    if (ncol(ftd$dat)>2)
      {
        X <- cbind(rep(1,nrow(ftd$dat)), ftd$dat[,-(1:2)])
      }else if (ncol(ftd$dat)==2){
        X <- matrix(rep(1,nrow(ftd$dat)),ncol=1)
      }else{
        stop("!!")
      }
    rownames(X) <- rownames(ftd$dat)
    Vcode <- ftd$Vcode
    
    uniID <- unique(ID)
    if (is.vector(betas)) betas <- matrix(betas,ncol=1)
    ##if (delta == -Inf) AR <- FALSE else AR <- TRUE
    gss <- uniIDorig <-  rep(NA,length(uniID))
    for (ipat in 1 : length(uniID) )
      {
        uniIDorig[ipat] <-
          unlist(strsplit(rownames(X)[ID==uniID[ipat]][1],"--"))[1]
        ipick <- ID==uniID[ipat]
        Yis <- Y[ipick]
        Xis <- X[ipick,,drop=FALSE]
        Vcodei <- Vcode[ipick]
        num <- integrate(f=int.numRE,lower=0,upper=Inf,
                         Yis=Yis,Xis=Xis,bet=betas,alpha=alpha,theta=theta,
                         AR=AR,delta=delta,Vcode=Vcodei,RE=RE,
                         rel.tol = rel.tol,expG=expG
                         )$value
        den <- integrate(f=int.denRE,lower=0,upper=Inf,
                         Yis=Yis,Xis=Xis,bet=betas,alpha=alpha,theta=theta,
                         AR=AR,delta=delta,Vcode=Vcodei,RE=RE,
                         rel.tol = rel.tol)$value
        
        gss[ipat] <- num/den
        ##cat("\n ipat",ipat,"g",gs[ipat]," num",num," den",den,"\n")
        ##print(Yis)
        ##print(Xis)
      }
    names(gss) <- uniIDorig
    return(gss)
  }




## lmeNB.sp.algo <-  function(formula,     ## an object of class "formula"
##                         ## (or one that can be coerced to that class): a symbolic description of the model to be fitted.
##                         data,        ## a data frame, list or environment (or object coercible by as.data.frame to a data frame)
##                         ## containing the variables in the model.
##                         ID,          ## a vector of length n*ni containing patient IDs of observations in data
##                         p.ini=NULL,  ## initial values for the parameters (log(a), log(th), b0, b1, ...)
##                         IPRT=TRUE,      ## printing control
##                         deps=1.e-4,   ## stop iteration when max(p.new - p.old) < deps
##                         maxit = 100
##                         )
## {
##   dat <- formulaToDat(formula=formula,data=data,ID=ID)  # dat = (ID, Y, x1, x2, ...) numeric matrix
##   DT <- getDT(dat$dat)  
##   #DT=list(id, y, x, cn=ncov, ys=Ysum, np=N, totN=totN, ni=Ni, ind=IND, xnames=xnames))
##   if (is.null(p.ini))
##   {
##     p.ini=rep(0, 3+DT$cn) 
##     p.ini[3]=mean(DT$y)   
##   }
##   p.new <- p.ini
##   W <- paraSave <- NULL

##   #initialize diffPara, ahat
##   diffPara <- 1    #max(p.new-p.old)
##   ahat <- 0
##   uhat <- rep(mean(tapply(DT$ys,DT$id,sum)),DT$np)
##   counter <- 0
##   while(diffPara > deps & maxit > counter)
##     {
##       p.ini <- p.new         
##       b.ini <- p.ini[-(1:2)]
##       ## ==================================================================== ## 
##       ## Main Step 1: estimate g_i i=1,...,N, g_i = w_i*(y_i+/mu_i+)+(1-w_i)  ##
##       ## ==================================================================== ##
##       ## the estimation of w_i = sqrt(var(G_i)/var(Y_i+/mu_i+))
##       ## weights for gi
##       if (ahat == 0) #first iteration, set wi=1
##         {
##           wi <- rep(1, DT$np)
##         }else{
##           gsig <- varG+(1+ahat*(varG+1))/mui.plus
##           wi <- sqrt(varG/gsig)   
##         }
##       gh0 <- DT$ys/mui.plus ## An estimate of y_i+/mu_i+ for each i
##       gh <- wi*gh0+1-wi ## An estimate of w_i*(y_i+/mu_i+)+(1-w_i) for each i
      
##       ## normalized gh so the average is 1
##       gh <- gh/mean(gh)
##       gh <- as.vector(gh)
      
##       ## frequency table for gh
##       tem <- sort(round(gh,6))
##       gh1 <- unique(tem)
##       ghw <- table(tem)/DT$np
##       ## the number of unique values of gh by 1 vector, containing the proportion of observations with varG value
##       gtb <- cbind(ghw, gh1) 
##       rownames(gtb)=NULL
##       ## moment estimate of var(G) 
##       varG <- estSg3(DT, ui=mui.plus, uij=muij) ## uh is a vector of length length(Y)

##       ## ========================================= ##
##       ## Main Step 2 : estimate log.alpha and beta ##
##       ## ========================================= ##
##       temB <- optim(c(log.alpha,b.ini), 
##                     psuedoNllk,
##                     DT=DT, Iprt=FALSE, 
##                     gs=gs## gi in freqency table form
##                     )
##       ahat <- exp(temB$par[1])
##       bhat <- temB$par[-1]
##       temp <- get.mu(bhat,DT)

##       muij <- temp$muij
##       mui.plus <- temp$mui.plus
       
##       p.new <- c(ahat, varG, bhat)
##       paraSave <- rbind(paraSave,p.new)
##       if (IPRT) {
##         cat("\n iteration:",counter)
##         cat("\n alpha, theta, beta0, beta1, ...\n");
##         print(p.new)
##       }
##       diffPara <- max(abs(p.new-p.ini))
##       counter <- counter +1
##     }
##   vcm <- NULL
  
##   rownames(paraSave) <- paste("sim",1:nrow(paraSave),sep="")
##   if (counter == maxit) warning("The maximum number of iterations occur!")
##   p.new <- matrix(p.new,ncol=1)
##   colnames(paraSave)  <- rownames(p.new) <- c("log_a", "log_var(G)", "(Intercept)", DT$xnames)
##   re <- list(opt = tt, V = vcm, est = p.new, gtb = gtb,counter=counter,
##               gi = as.vector(gh), RE = "NoN", 
##               ##idat = data.frame(dat),
##               cor="ind",formula=formula)
##   class(re) <- "LinearMixedEffectNBFreq"
##   return(re)
## }


## psuedoNllk <- function(para, 
##                        DT, Iprt=FALSE, 
##                        gs ## gi in freqency table form
##                        ) 
## {
##   ## para = c(logAlpha,beta0,beta1,...)
##   ## psudo profile likelihood as a function of alpha
##   ## Li(a; hatb, hatg, y)
##   ## = 1/N sum_l=1^N prod_j=1^ni Pr(Y_ij=y_ij|G_l=g_l,hatb,a)
##   ## = sum_l=1^N* prod_j=1^ni Pr(Y_ij=y_ij|G_l=g*_l,hatb,a)*prop(g*_l)
##   ## = sum_l=1^N* prod_j=1^ni choose(y_ij+r_ij-1,y_ij) (1-p(g*_l,a))^r_ij p(g*_l,a)^{y_ij}*prop(g*_l)
##   ## = [prod_j=1^ni choose(y_ij+r_ij-1,y_ij)] sum_l=1^N* (1-p(g*_l,a))^r_i+ p(g*_l,a)^{y_i+}*prop(g*_l)
##   ## 
##   ## log Li
##   ## = [sum_j=1^ni log choose(y_ij+r_ij-1,y_ij)] + log [sum_l=1^N* (1-p(g*_l,a))^r_i+ p(g*_l,a)^{y_i+}*prop(g*_l)]
##   ## 
##   ## g*_l represents the unique value of the approximations
##   ## prop(g*_l) represents the proportion of approximated g_l's that have the value g*_l
##   ## N* represents the number of unique g*
##   ## p = p(g*_l,a) = 1/(g*_l*a + 1)
  
##   if (Iprt) cat(para)
  
##   alpha <- exp(para[1])
##   betas <- para[-1]
##   muij <- get.mu(bhat=betas,DT=DT)$muij
##   ainv <- 1/alpha  
  
##   th2 <- muij/alpha ## r_ij = exp(X_ij^T hat.beta)/alpha

##   if (!all(is.finite(th2)) | any(th2==0)) return(Inf) ## alpha is not in the acceptable range

##   ## -log [sum_j=1^ni log choose(y_ij+r_ij-1,y_ij)] = sum( - lgamma(dat$y+th2) + lgamma(th2) + lgamma(dat$y))
##   ##  but the last term is not necessary because it does not depend on a
##   nllk <- sum( - lgamma(DT$y+th2) + lgamma(th2))
##   us <- tapply(th2, DT$id, sum)
  
##   p <- gs[,2]/(gs[,2]+ainv)
##   ## please note that p = 1 - p in this manuscript
##   ## gs[,2] is a vector of length # unique g_i containing the values of unique g_i
  
##   q <- 1-p
##   for (i in 1:DT$np)
##     {
##       tem <- sum(
##                  p^DT$ys[i]*q^us[i]*gs[,1]
##                  ## gs[,1] is a vector of length # unique g_i
##                  ## containing the proportion of g_i
##                  )
##       nllk <- nllk-log(tem)
##     }
  
##   if (Iprt) cat(" nllk=", nllk, "\n")
##   return(nllk)
## }
