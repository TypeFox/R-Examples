ibr.fit <- function(x,y,criterion="gcv",df=1.5,Kmin=1,Kmax=1e+06,smoother="k",kernel="g",rank=NULL,control.par=list(),cv.options=list()) {
  cl <- match.call() 
  crit <-c("aic","aicc","gcv","bic","gmdl","rmse","map")
  crite <- pmatch(criterion,crit)
  if (any(is.na(crite))) stop(paste("parameter criterion must be in",paste(crit,collapse=", "))) else criterion <- crit[crite]
  if ((any(crite>5))&&(length(crite)>2)) stop("RMSE or MAP must be used without any other criterion")
  iterautre <- NULL
  choixkautre <- NULL
  ##  criterion <- match.arg(criterion,crit)
  smoothertable <- c("k","tps","ds","lrtps","lrds")
  smoother <- match.arg(smoother,smoothertable)
  lowrank <- ifelse(substr(smoother,0,2)=="lr",TRUE,FALSE)
  if (!is.matrix(x)) {
    x <- data.matrix(x)
    warning("x is coerced to matrix by data.matrix function\n")
  }
  n <- nrow(x)
  p <- ncol(x)
  if (length(y)!=n) stop("numbers of observations in x and y do not match\n")
  if (smoother=="k") contr.sp <- list(bandwidth=NULL,iter=NULL,really.big=FALSE,dftobwitmax=1000,exhaustive=FALSE,m=NULL,s=NULL,dftotal=FALSE,accuracy=0.01,dfmaxi=2*n/3,fraction=c(100, 200, 500, 1000, 5000,10^4,5e+04,1e+05,5e+05,1e+06),scale=FALSE,criterion="strict",aggregfun=function(x) {10^(floor(log10(x[2]))+2)})
  else  contr.sp <- list(bandwidth=NULL,iter=NULL,really.big=FALSE,dftobwitmax=1000,exhaustive=FALSE,m=NULL,s=NULL,dftotal=FALSE,accuracy=0.01,dfmaxi=2*n/3,fraction=c(100, 200, 500, 1000, 5000,10^4,5e+04,1e+05,5e+05,1e+06),scale=TRUE,criterion="strict",aggregfun=function(x) {10^(floor(log10(x[2]))+2)})
  contr.sp[(names(control.par))] <- control.par
  if (!(is.logical(contr.sp$dftotal))) stop("contr.sp$dftotal must be logical\n")
  if ((!is.null(contr.sp$bandwidth))&(!is.numeric(contr.sp$bandwidth) || (contr.sp$bandwidth<0))) stop("invalid bandwidth\n")
  if ((!is.null(contr.sp$iter))&(!is.numeric(contr.sp$iter) || (contr.sp$iter<0) || (floor(contr.sp$iter)!=contr.sp$iter))) stop("invalid number of iterations\n")
  iter <- contr.sp$iter
   if ((contr.sp$dfmaxi<=0)|(contr.sp$dfmaxi>n)) stop("invalid dfmaxi\n")
  crit <-c("recalc","strict","aggregation")
  crite <- pmatch(contr.sp$criterion,crit)
  if (is.na(crite)) stop(paste("control.par$criterion must be in",crit))
  if (length(criterion)==2&&(all(criterion%in%c("map","rmse")))&&(!contr.sp$exhaustive)) stop("when RMSE and MAP are both selected, exhaustive search (in control.par list) must be chosen")
  if ((smoother=="tps")|(smoother=="lrtps")) {
    contr.sp$s <- 0
    if (is.null(contr.sp$m)) contr.sp$m <- floor(p/2)+1 else {
      if (contr.sp$m<=(p/2)) stop("order of thin plate splines is invalid (need to be greater than p/2)\n")
    }
    if (!is.numeric(contr.sp$m) || (contr.sp$m<0) || (floor(contr.sp$m)!=contr.sp$m)) stop("invalid spline order\n")
    ddlmin <- choose(contr.sp$m+p-1,contr.sp$m-1)
    if (ddlmin>=n) stop(paste("ddl min is equal to",ddlmin,"and the number observations is",n))
  }
  if (smoother=="k") {
    contr.sp$m <- NULL
}
  smoothobject <- NULL
  if ((smoother=="ds")|(smoother=="lrds")) {
    if (is.null(contr.sp$m)) contr.sp$m <- 2 ## default penalty order 2
    if (is.null(contr.sp$s)) contr.sp$s <- (p-1)/2 ## default pseudo cubic
    if ((!is.numeric(contr.sp$m))|(!is.numeric(contr.sp$s))) stop("contr.par$m or contr.par$s is not numeric...\n")
    contr.sp$m <- round(contr.sp$m)     ## m is integer
    contr.sp$s <- round(contr.sp$s*2)/2 ## s is in halfs
    if (contr.sp$m< 1) contr.sp$m <- 1  ## m > 0
    ## check that -p/2 < s < p/2...
    if (contr.sp$s >= p/2) { 
      contr.sp$s <- (p-1)/2
      warning("contr.par$s value reduced")
    } 
    if (contr.sp$s <= -p/2) { 
      contr.sp$s <- -(p-1)/2
      warning("contr.par$s value increased")
    }
    
    ## m + s > p/2 for continuity...
    if ((contr.sp$m+contr.sp$s)<=p/2) {
      contr.sp$s <- 1/2 + p/2 - contr.sp$m
      if (contr.sp$s>=p/2) stop("No suitable contr.par$s try increasing contr.par$m")
      warning("contr.par$s value modified to give continuous function")
  }
        ddlmin <- choose(contr.sp$m+p-1,contr.sp$m-1)
    if (ddlmin>=n) stop(paste("ddl min is equal to",ddlmin,"and the number observations is",n))
  }
  
  if (lowrank) {
      if (is.null(rank)) stop("rank argument for lowrank splines must be chosen...")
      if (rank>n) stop(paste("rank argument must be less than",n))
        bs <- substr(smoother,3,4)
        listvarx <- colnames(x)
    } 
  moy <- NULL
  ec <- NULL
  if (contr.sp$scale) {
      if (smoother=="k") warning("when using kernel smoother, you do not need to scale\n")
      if (lowrank) ec <- apply(x,2, sd)*sqrt((n-1)/n) else  ec <- apply(x,2, sd)
      x <- scale(x,scale=ec)
      moy <- attr(x,"scaled:center")
  }
  if (all(criterion%in%c("rmse","map"))) {
    cv <- list(bwchange=FALSE,ntest=floor(nrow(x)/10),ntrain=NULL,Kfold=FALSE,type="random",seed=NULL,npermut=20)
    cv[(names(cv.options))] <- cv.options
    if (!all(sapply(cv[1],is.logical))) stop("invalid cv$bwchange or cv$Kfold: must be logical\n")
    if (!all(sapply(cv[c(2,3,6,8)], FUN=function(x) is.numeric(x)||is.null(x)))) stop("invalid cv parameters: must be numeric or NULL\n")
    if (any(names(cv.options)=="ntrain")) cv$ntest <- NULL
  } else cv <- NULL
  if (!(lowrank)&&((n>1000)&(! contr.sp$really.big))) stop("number of observations is greater than 1000, set control.par$really.big to TRUE if you really want to do the requested calculations (but computational time -eigen decomposition- could be prohibitive)\n")
  if (smoother=="k") {
    m <- NULL
    if (!is.null(contr.sp$bandwidth)) {
      if (length(contr.sp$bandwidth)==1) {
        bandwidth <- rep(contr.sp$bandwidth,p)
      } else {
        if (length(contr.sp$bandwidth)!=p)  stop(paste("the length of bandwidth vector have to be",p,"or 1\n"))
        bandwidth <- contr.sp$bandwidth
      }
      listeA <- calcA(X = x, bx = bandwidth, kernelx = kernel)
    } else {
      if (kernel=="g") {
        if (df<=1) stop("degree of freedom should be greater than 1\n")   
        departbw <- apply(x,2,FUN=function(z) 3*abs(diff(range(z))))
        calculus <- .C("gaustotal",
                      as.double(departbw),
       if (contr.sp$dftotal) as.double(rep(1e-10,p)) else
                       as.double(1e-10),
                 as.double(x),as.integer(n),as.integer(p),
                 as.double(.Machine$double.eps^0.25),
                 maxit=as.integer(contr.sp$dftobwitmax),
                 as.double(df),as.integer(contr.sp$dftotal),
                 A=double(n^2),Ddemi=double(n),df=double(1),
                 bandwidth=double(p),PACKAGE="ibr")
        if (calculus$maxit==-1) warning("failed to find appropriate bandwidth. Try to increase dftobwitmax argument of control.par list ?\n")
        listeA <- list(A=matrix(calculus$A,n,n),Ddemi=calculus$Ddemi,df=calculus$df)
        bandwidth <- calculus$bandwidth
      rm(calculus)
      } else {
        if (contr.sp$dftotal) {
          dfobjectif <- df
          dfmini <- 1.05
          repeat {
          ddlminimum <- departnoyau(dfmini,x,kernel,contr.sp$dftobwitmax,n,p,dfobjectif)
          if (ddlminimum<0) break else dfmini <- dfmini-0.02
        }
          if (abs(ddlminimum-dfobjectif)>0.1) {
            ddlmaximum <- 1
            dfmaxi <- 1.2
            repeat {
              ddlmaximum <- departnoyau(dfmaxi,x,kernel,contr.sp$dftobwitmax,n,p,dfobjectif)
              if (ddlmaximum>0) break else dfmaxi <- dfmaxi+0.2
            }
            if (abs(ddlmaximum-dfobjectif)>0.1) {
              res <- uniroot(departnoyau,c(dfmini,dfmaxi),tol=contr.sp$accuracy,x=x,kernel=kernel,dftobwitmax=contr.sp$dftobwitmax,n=n,p=p,dfobjectif)
              df <- res$root
            } else df <- dfmaxi
          } else df <- dfmini
        }
        bandwidth <- bwchoice(x,df,kernel,contr.sp$dftobwitmax)
        listeA <- calcA(X=x,bx=bandwidth,kernelx=kernel)
      }
    }
    dfstart <- listeA$df
    listeA.eig <- eigen(listeA$A,symmetric=TRUE)
    tPADmdemiY <- t(listeA.eig$vectors*(1/listeA$Ddemi))%*%y
    eigenvaluesA <- listeA.eig$values
    if (any(zapsmall(eigenvaluesA-1,digits=9)==0)) {
      ddlmini <-  sum(zapsmall(eigenvaluesA-1,digits=9)==0)
    } else ddlmini <- NA
    if (any(zapsmall(eigenvaluesA,digits=9)==0)) {
      index0 <-  which(zapsmall(eigenvaluesA,digits=9)==0)[1]
    } else index0 <- NA
    DdemiPA <- (listeA$Ddemi*listeA.eig$vectors)
    rm(listeA.eig)
    if (is.null(iter)) {
      if (Kmax<=Kmin) stop("Kmax should be greater than Kmin\n")
      if (any(c(Kmin,Kmax)<=0)) stop("Kmin and Kmax should be greater than 0")
      if (all(criterion%in%c("rmse","map"))) {
        if (cv$bwchange) {
          if (is.null(df)) stop("df needs to be set\n")
          bx <- NULL } else bx <- bandwidth
        if (contr.sp$exhaustive) {
          choixkautre <- iterchoiceAcve(x,y,bx,df,kernel,ddlmini,cv$ntest,cv$ntrain,cv$Kfold,cv$type,cv$npermut,cv$seed,Kmin,Kmax)
          iterautre <- (Kmin:Kmax)[unlist(lapply(choixkautre,FUN=which.min))]
          iter <- switch(contr.sp$criterion,strict=iterautre[1],
                 aggregation=contr.sp$aggregfun(iterautre),
                 recalc= { Kmax2 <- contr.sp$aggregfun(iterautre) ; iterautre[1] <- (Kmin:Kmax2)[which.min(choixkautre[[1]])] ; iterautre[1]})
          if (contr.sp$criterion=="aggregation") {
            choixk <- NA
          } else {
            choixk <- switch(criterion[1],aic=choixkautre$aic[iter],aicc=choixkautre$aicc[iter],gcv=choixkautre$gcv[iter],bic=choixkautre$bic[iter],gmdl=choixkautre$gmdl[iter])
          }
#          iter <- iterautre[1]
#          allcrit <- switch(criterion,rmse=choixk$rmse,map=choixk$map)
#          iter <- (Kmin:Kmax)[which.min(allcrit)]
        } else {
          prov <- iterchoiceAcv(x,y,bx,df,kernel,ddlmini,cv$ntest,cv$ntrain,cv$Kfold,cv$type,cv$npermut,cv$seed,Kmin,Kmax,criterion,contr.sp$fraction)
          iter <- prov$iter
          choixk <- prov$objective
         if (criterion=="rmse") choixk <- sqrt(choixk)
        }
      } else {
        if (contr.sp$exhaustive) {
          choixkautre <- iterchoiceAe(y,Kmin:Kmax,eigenvaluesA,tPADmdemiY,DdemiPA,ddlmini,contr.sp$dfmaxi)
          iterautre <- (Kmin:Kmax)[unlist(lapply(choixkautre[c("aic","aicc","gcv","bic","gmdl")],FUN=which.min))]
          iter <- switch(contr.sp$criterion,strict=iterautre[1],
                 aggregation=contr.sp$aggregfun(iterautre),
                 recalc= { Kmax2 <- contr.sp$aggregfun(iterautre) ; iterautre[1] <- (Kmin:Kmax2)[which.min(choixkautre[[1]])] ; iterautre[1] })
#          iter <- iterautre[1]
          if (contr.sp$criterion=="aggregation") {
            choixk <- NA
          } else {
            choixk <- switch(criterion[1],aic=choixkautre$aic[iter],aicc=choixkautre$aicc[iter],gcv=choixkautre$gcv[iter],bic=choixkautre$bic[iter],gmdl=choixkautre$gmdl[iter])
          }
#          allcrit <- switch(criterion[1],aic=choixk$aic,aicc=choixk$aicc,gcv=choixk$gcv,bic=choixk$bic,gmdl=choixk$gmdl)
#          iter <- (Kmin:Kmax)[which.min(allcrit)]
        } else {
          prov <- iterchoiceA(n,Kmin,Kmax,eigenvaluesA,tPADmdemiY,DdemiPA,ddlmini,contr.sp$dfmaxi,y,criterion[1],contr.sp$fraction)
          iter <- prov$iter
          choixk <- prov$objective
          iterautre <- iter
          choixkautre <- prov$objective
          if (length(criterion)>1) {
            for (i in 2:length(criterion)) {
              prov <- iterchoiceA(n,Kmin,Kmax,eigenvaluesA,tPADmdemiY,DdemiPA,ddlmini,contr.sp$dfmaxi,y,criterion[i],contr.sp$fraction)
              iterautre <- c(iterautre,prov$iter)
              choixkautre <- c(choixkautre,prov$objective)
            }
            iter <- switch(contr.sp$criterion,strict=iterautre[1],
                           aggregation={
                   choixk <- NA
                   contr.sp$aggregfun(iterautre) }, recalc= {
                     Kmax2 <- contr.sp$aggregfun(iterautre)
                     prov <- iterchoiceA(n,Kmin,Kmax2,eigenvaluesA,tPADmdemiY,DdemiPA,ddlmini,contr.sp$dfmaxi,y,criterion[1],contr.sp$fraction)
                     choixk <- prov$objective
                     iterautre[1] <- prov$iter
                     prov$iter })
          }
        }
      }
      if ((((Kmax-iter)/(Kmax-Kmin)<1e-5)|(Kmax-iter)<3)&(!contr.sp$exhaustive)) warning(paste("Number of iterations is chosen close to the boundary of grid search: ",Kmax,".\nIncrease the maximum number of iterations or use contr.sp$exhaustive search\n",sep=""))
      if ((iter==max(Kmax))&(contr.sp$exhaustive)) warning(paste("Number of iterations is chosen at the boundary of grid search: ",Kmax,".\nIncrease the maximum number of iterations\n",sep=""))
    } else {
      criterion="user"
      choixk <- NULL
  }
    beta <- betaA(n,eigenvaluesA,tPADmdemiY,DdemiPA,ddlmini,k=iter,index0)
    listefit <- fittedA(n,eigenvaluesA,tPADmdemiY,DdemiPA,ddlmini,k=iter)
    residuals <- y-listefit$fit
  }
  if ((smoother=="tps")|(smoother=="ds")|lowrank) {
    bandwidth <- contr.sp$bandwidth
    if (length(df)>1) stop("only one df is possible with Splines\n")
    ddlmini <- choose(contr.sp$m+p-1,contr.sp$m-1)
    if (lowrank)     {
        if (is.null(bandwidth)) {
            lambda <- lambdachoicelr(x,ddlmini*df,m=contr.sp$m,contr.sp$s,rank,itermax=contr.sp$dftobwitmax,bs,listvarx)
            bandwidth <- lambda
        } else {
            lambda <- bandwidth
        }
        S2 <- lrsmoother(x,bs,listvarx,lambda=lambda,m=contr.sp$m,s=contr.sp$s,rank)
        eigenvaluesS1 <- S2$values
        dfstart <- sum(eigenvaluesS1)
        if (any(zapsmall(eigenvaluesS1,digits=9)==0)) {
            index0 <-  which(zapsmall(eigenvaluesS1,digits=9)==0)[1]
        } else index0 <- rank+1
        eigenvaluesS1 <- c(eigenvaluesS1,rep(0,n-rank))
        eigenvaluesS1[eigenvaluesS1<0] <- 0
        U <- S2$vectors
        Rm1U <- S2$Rm1U
        tUy <- c(as.vector(crossprod(U, y)),rep(0,n-rank))
        smoothobject <- S2$smoothobject
        rm(S2)
    } else {
        if (is.null(bandwidth)) {
            lambda <- lambdachoice(x,ddlmini*df,m=contr.sp$m,contr.sp$s,itermax=contr.sp$dftobwitmax,smoother)
            bandwidth <- lambda
        } else {
            lambda <- bandwidth
        }
        S1 <- dssmoother(x, y,lambda=lambda,m=contr.sp$m,s=contr.sp$s)
        vp1.S1 <- eigen(S1$H,symmetric=TRUE)
        U <- vp1.S1$vect
        eigenvaluesS1 <- vp1.S1$values
        rm(vp1.S1)
        tUy <- as.vector(crossprod(U, y))
        if (any(zapsmall(eigenvaluesS1,digits=9)==0)) {
            index0 <-  which(zapsmall(eigenvaluesS1,digits=9)==0)[1]
        } else index0 <- NA
        dfstart <- sum(eigenvaluesS1)
    }
    if (is.null(iter)) {
      if (Kmax<=Kmin) stop("Kmax should be greater than Kmin\n")
      if (any(c(Kmin,Kmax)<=0)) stop("Kmin and Kmax should be greater than 0")
      if (all(criterion%in%c("rmse","map"))) {
        if (cv$bwchange) {
          if (is.null(df)) stop("df needs to be set\n")
          lambda <- NULL }
        if (contr.sp$exhaustive) {
            if (lowrank)     {
                choixkautre <- iterchoiceS1lrcve(x,y,lambda,rank,bs,listvarx,df,ddlmini,cv$ntest,cv$ntrain,cv$Kfold,cv$type,cv$npermut,cv$seed,Kmin,Kmax,contr.sp$m,contr.sp$s)
            } else {
                choixkautre <- iterchoiceS1cve(x,y,lambda,df,ddlmini,cv$ntest,cv$ntrain,cv$Kfold,cv$type,cv$npermut,cv$seed,Kmin,Kmax,contr.sp$m,contr.sp$s)
            }
          iterautre <- (Kmin:Kmax)[unlist(lapply(choixkautre,FUN=which.min))]
          iter <- switch(contr.sp$criterion,strict=iterautre[1],
                 aggregation=contr.sp$aggregfun(iterautre),
                 recalc= { Kmax2 <- contr.sp$aggregfun(iterautre) ; iterautre[1] <- (Kmin:Kmax2)[which.min(choixkautre[[1]])] ; iterautre[1] })
#          iter <- iterautre[1]
          if (contr.sp$criterion=="aggregation") {
            choixk <- NA
          } else {
            choixk <- switch(criterion[1],aic=choixkautre$aic[iter],aicc=choixkautre$aicc[iter],gcv=choixkautre$gcv[iter],bic=choixkautre$bic[iter],gmdl=choixkautre$gmdl[iter])
          }
#          iter <- iterautre[1]
#          choixk <- switch(criterion[1],aic=choixkautre$aic[iter],aicc=choixkautre$aicc[iter],gcv=choixkautre$gcv[iter],bic=choixkautre$bic[iter],gmdl=choixkautre$gmdl[iter])
        } else {
            if (lowrank)     {
                 prov <- iterchoiceS1lrcv(x,y,lambda,rank,bs,listvarx,df,ddlmini,cv$ntest,cv$ntrain,cv$Kfold,cv$type,cv$npermut,cv$seed,Kmin,Kmax,criterion,contr.sp$m,contr.sp$s,contr.sp$fraction)
             } else {
                 prov <- iterchoiceS1cv(x,y,lambda,df,ddlmini,cv$ntest,cv$ntrain,cv$Kfold,cv$type,cv$npermut,cv$seed,Kmin,Kmax,criterion,contr.sp$m,contr.sp$s,contr.sp$fraction)
             }
          iter <- prov$iter
          choixk <- prov$objective
          if (criterion=="rmse") choixk <- sqrt(choixk)
        }
      } else {
        if (contr.sp$exhaustive) {
          choixkautre <- iterchoiceS1e(y,Kmin:Kmax,tUy,eigenvaluesS1,ddlmini,contr.sp$dfmaxi)
          iterautre <- (Kmin:Kmax)[unlist(lapply(choixkautre[c("aic","aicc","gcv","bic","gmdl")],FUN=which.min))]
          iter <- switch(contr.sp$criterion,strict=iterautre[1],
                 aggregation=contr.sp$aggregfun(iterautre),
                 recalc= { Kmax2 <- contr.sp$aggregfun(iterautre) ; iterautre[1] <- (Kmin:Kmax2)[which.min(choixkautre[[1]])] ; iterautre[1]})
          if (contr.sp$criterion=="aggregation") {
            choixk <- NA
          } else {
            choixk <- switch(criterion[1],aic=choixkautre$aic[iter],aicc=choixkautre$aicc[iter],gcv=choixkautre$gcv[iter],bic=choixkautre$bic[iter],gmdl=choixkautre$gmdl[iter])
          }
#          allcrit <- switch(criterion[1],aic=choixkautre$aic,aicc=choixkautre$aicc,gcv=choixkautre$gcv,bic=choixkautre$bic,gmdl=choixkautre$gmdl)
#          iter <- (Kmin:Kmax)[which.min(allcrit)]
        } else {
          prov <- iterchoiceS1(n,Kmin,Kmax,tUy,eigenvaluesS1,ddlmini,contr.sp$dfmaxi,y,criterion[1],contr.sp$fraction)
          iter <- prov$iter
          choixk <- prov$objective
          iterautre <- iter
          choixkautre <- prov$objective
          if (length(criterion)>1) {
            for (i in 2:length(criterion)) {
              prov <- iterchoiceS1(n,Kmin,Kmax,tUy,eigenvaluesS1,ddlmini,contr.sp$dfmaxi,y,criterion[i],contr.sp$fraction)
              iterautre <- c(iterautre,prov$iter)
              choixkautre <- c(choixkautre,prov$objective)
            }
            iter <- switch(contr.sp$criterion,strict=iterautre[1],
                           aggregation={
                   choixk <- NA
                   contr.sp$aggregfun(iterautre) }, recalc= {
                     Kmax2 <- contr.sp$aggregfun(iterautre)
              prov <- iterchoiceS1(n,Kmin,Kmax2,tUy,eigenvaluesS1,ddlmini,contr.sp$dfmaxi,y,criterion[1],contr.sp$fraction)
                     choixk <- prov$objective
                     iterautre[1] <- prov$iter
                     prov$iter }) 
           }
        }
      }
      if ((((Kmax-iter)/(Kmax-Kmin)<1e-5)|(Kmax-iter)<3)&(!contr.sp$exhaustive)) warning(paste("Number of iterations is chosen close to the boundary of grid search: ",Kmax,".\n  Increase the maximum number of iterations or use contr.sp$exhaustive search\n",sep=""))
      if ((iter==max(Kmax))&(contr.sp$exhaustive)) warning(paste("Number of iterations is chosen at the boundary of grid search: ",Kmax,".\n  Increase the maximum number of iterations\n",sep=""))
    } else {
      criterion="user"
      choixk <- NULL
    }
    if (lowrank) {
        beta <- betaS1lr(n,U,tUy,eigenvaluesS1,ddlmini,iter,lambda,rank,Rm1U,index0)
        listefit <- fittedS1lr(n,U,tUy,eigenvaluesS1,ddlmini,iter,rank)
      if (listefit$trace>=0.99*rank) warning(paste("Final rank is chosen equal to the lowrank of spline: ",rank,".\n  Increase the rank argument\n",sep=""))
    } else {
        listebeta <- betaS1(n,U,tUy,eigenvaluesS1,ddlmini,iter,lambda,S1$Sgu,S1$Qgu,index0)
        beta <- list(d=listebeta$dgub,c=listebeta$cgub)
        listefit <- fittedS1(n,U,tUy,eigenvaluesS1,ddlmini,iter)
    }
    residuals <- y- listefit$fit
}
  res <- list(beta=beta,residuals=residuals,fitted=listefit$fit,iter=iter,initialdf=dfstart,finaldf=listefit$trace,bandwidth=bandwidth,call=cl,parcall=list(p=p,m=contr.sp$m,s=contr.sp$s,scaled=contr.sp$scale,mean=moy,sd=ec,critmethod=contr.sp$criterion,rank=rank,criterion=criterion,smoother = smoother, kernel = kernel,smoothobject=smoothobject),criteria=choixk,alliter=iterautre,allcriteria=choixkautre)    
  return(res)
}
