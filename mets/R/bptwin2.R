bptwin2 <- function(formula, data, id, zyg, DZ, group=NULL,
                   num=NULL,
                   weight=NULL,
                   biweight=function(x) 1/min(x),
                   strata=NULL,
                   messages=1,
                   control=list(trace=0),
                   type="ace",
                   eqmean=TRUE,
                   pairsonly=FALSE,
                   samecens=TRUE,
                   allmarg=samecens&!is.null(weight),
                   stderr=TRUE,                  
                   robustvar=TRUE,                   
                   p, indiv=FALSE,
                   constrain,
                   bound=FALSE,
                   varlink,
                   ...) {

###{{{ setup

    OSon <- FALSE
    idx2 <- NULL
    
  mycall <- match.call()
  formulaId <- unlist(Specials(formula,"cluster"))
  formulaStrata <- unlist(Specials(formula,"strata"))
  formulaSt <- paste("~.-cluster(",formulaId,")-strata(",paste(formulaStrata,collapse="+"),")")
  formula <- update(formula,formulaSt)
  if (!is.null(formulaId)) {
    id <- formulaId
    mycall$id <- id
  }
  if (!is.null(formulaStrata)) strata <- formulaStrata
  mycall$formula <- formula
 
  if (!is.null(strata)) {
    dd <- split(data,interaction(data[,strata]))
    nn <- unlist(lapply(dd,nrow))
    dd[which(nn==0)] <- NULL
    if (length(dd)>1) {
      fit <- lapply(seq(length(dd)),function(i) {
        if (messages>0) message("Strata '",names(dd)[i],"'")
        mycall$data <- dd[[i]]
        eval(mycall)
      })
      res <- list(model=fit)
      res$strata <- names(res$model) <- names(dd)
      class(res) <- c("twinlm.strata","biprobit")
      res$coef <- unlist(lapply(res$model,coef))
      res$vcov <- blockdiag(lapply(res$model,vcov.biprobit))
      res$N <- length(dd)
      res$idx <- seq(length(coef(res$model[[1]])))
      rownames(res$vcov) <- colnames(res$vcov) <- names(res$coef)
      return(res)
    }
  }

##################################################
### No strata
  if (is.null(control$method)) {
      if (!samecens & !is.null(weight)) {
        control$method <- "bhhh"
      } else {
        if (requireNamespace("ucminf",quietly=TRUE)) {
          control$method <- "gradient"
        } else control$method <- "nlminb"
      }
  }
  if (length(grep("flex",tolower(type)))>0) { type <- "u"; eqmean <- FALSE }

  yvar <- paste(deparse(formula[[2]]),collapse="")
  data <- data[order(data[,id]),]
  idtab <- table(data[,id])
  if (sum(idtab>2)) stop("More than two individuals with the same id ")
  
  if (pairsonly) {
    data <- data[as.character(data[,id])%in%names(idtab)[idtab==2],]
    idtab <- table(data[,id])
  }
  if (is.logical(data[,yvar])) data[,yvar] <- data[,yvar]*1
  if (is.factor(data[,yvar])) data[,yvar] <- as.numeric(data[,yvar])-1  

  if (missing(DZ)) {
    DZ <- levels(as.factor(data[,zyg]))[1]
    message("Using '",DZ,"' as DZ",sep="")
  }

  idx1 <- which(data[,zyg]%in%DZ) ## DZ
  if (length(idx1)==0) stop("No DZ twins found")
  idx0 <- which(!(data[,zyg]%in%DZ)) ## MZ
  if (length(idx1)==0) stop("No MZ twins found")
  zyg2 <- rep(1,nrow(data)); zyg2[idx0] <- 0;
  data[,zyg] <- zyg2 ## MZ=0, DZ=1

  if (!is.null(group)) data[,group] <- as.factor(data[,group])
 

  ff <- paste(as.character(formula)[3],"+",
              paste(c(id,zyg,weight,num),collapse="+"))
  ff <- paste("~",yvar,"+",ff)
  formula0 <- as.formula(ff)
  opt <- options(na.action="na.pass")
  Data <- model.matrix(formula0,data)
  options(opt)
  rnames1 <- setdiff(colnames(Data),c(yvar,id,weight,zyg,num))
  nx <- length(rnames1) 
  if (nx==0) stop("Zero design not allowed")
  

  bidx0 <- seq(nx)
  midx0 <- bidx0; midx1 <- midx0+nx
  dS0. <- rbind(rep(1,4),rep(1,4),rep(1,4)) ## MZ
  dS1. <- rbind(c(1,.5,.5,1),rep(1,4),c(1,.25,.25,1)) ## DZ
  dS2.  <- rbind(c(1,0,0,1),rep(1,4),c(1,0,0,1),c(0,1,1,0))
  
  ##mytr <- function(x) x; dmytr <- function(x) 1
  ##mytr <- function(x) x^2; dmytr <- function(x) 2*x
  ##mytr <- function(z) 1/(1+exp(-z)); dmytr <- function(z) exp(-z)/(1+exp(-z))^2
  ACDU <- sapply(c("a","c","d","e","u"),function(x) length(grep(x,tolower(type)))>0)
  ACDU <- c(ACDU)

  if (missing(varlink) || !is.null(varlink)) {
      mytr <- exp; dmytr <- exp; myinvtr <- log
      trname <- "exp"; invtrname <- "log"
  } else {
      mytr <- myinvtr <- identity; dmytr <- function(x) rep(1,length(x))
      trname <- ""; invtrname <- ""
  }
  
  dmytr2 <- function(z) 4*exp(2*z)/(exp(2*z)+1)^2
  mytr2 <- tanh;  myinvtr2 <- atanh
  trname2 <- "tanh"; invtrname2 <- "atanh"

  if (!is.null(group)) {
    mytr <- function(x) c(exp(x[-length(x)]),mytr2(x[length(x)]))
    myinvtr <- function(z) c(log(z[-length(z)]),myinvtr2(z[length(z)]))
    dmytr <- function(x) c(exp(x[-length(x)]),dmytr2(x[length(x)]))
  }

  if (ACDU["u"]) {
    ##      datanh <- function(r) 1/(1-r^2)
    dmytr <- dmytr2
    mytr <- mytr2;  myinvtr <- myinvtr2
    trname <- trname2; invtrname <- invtrname2
    dS0 <- rbind(c(0,1,1,0))
    vidx0 <- 1
    vidx1 <- 2
    vidx2 <- 3
    dS2 <- dS1 <- dS0
    nvar <- length(vidx0)+length(vidx1)
    if (OSon) nvar <- nvar+length(vidx2)
  } else {
    nvar <- sum(ACDU[1:3])
    vidx0 <- vidx1 <- seq(nvar); vidx2 <- seq(nvar+1)
    if (OSon) nvar <- nvar+1
    dS0 <- dS0.[ACDU[1:3],,drop=FALSE]
    dS1 <- dS1.[ACDU[1:3],,drop=FALSE]
    dS2 <- dS2.[which(c(ACDU[1:3],TRUE)),,drop=FALSE]
  }  
  if (eqmean) {
    bidx2 <- bidx1 <- bidx0
  } else {
    bidx1 <- bidx0+nx
    bidx2 <- bidx1+nx    
    if (OSon) nx <- 3*nx else nx <- 2*nx;
  }
  
  vidx0 <- vidx0+nx; vidx1 <- vidx1+nx; vidx2 <- vidx2+nx
  vidx <- nx+seq_len(nvar)
  midx <- seq_len(nx)
  plen <- nx+nvar

  Am <- matrix(c(1,.5,.5,1),ncol=2)
  Dm <- matrix(c(1,.25,.25,1),ncol=2)
  Vm <- matrix(c(1,0,0,1),ncol=2)
  Rm <- matrix(c(0,1,1,0),ncol=2)

##################################################

  ## system.time(Wide <- reshape(as.data.frame(Data),idvar=c(id,zyg),timevar=time,direction="wide"))
  ##  system.time(Wide <- as.data.frame(fast.reshape(Data,id=c(id),sep=".")))
  Wide <- as.data.frame(fast.reshape(Data,id=c(id,zyg),sep=".",idcombine=FALSE,labelnum=TRUE))
  yidx <- paste(yvar,1:2,sep=".")
  rmidx <- c(id,yidx,zyg)
    
  W0 <- W1 <- W2 <- NULL
  if (!is.null(weight)) {
    widx <- paste(weight,1:2,sep=".")
    rmidx <- c(rmidx,widx)
    W0 <- as.matrix(Wide[which(Wide[,zyg]==0),widx,drop=FALSE])
    W1 <- as.matrix(Wide[which(Wide[,zyg]==1),widx,drop=FALSE])
    W2 <- as.matrix(Wide[which(Wide[,zyg]==2),widx,drop=FALSE])
  }
  XX <- as.matrix(Wide[,setdiff(colnames(Wide),rmidx)])
  XX[is.na(XX)] <- 0

  Y0 <- as.matrix(Wide[which(Wide[,zyg]==0),yidx,drop=FALSE])
  Y1 <- as.matrix(Wide[which(Wide[,zyg]==1),yidx,drop=FALSE])
  Y2 <- as.matrix(Wide[which(Wide[,zyg]==2),yidx,drop=FALSE])
  XX0 <- XX[which(Wide[,zyg]==0),,drop=FALSE]
  XX1 <- XX[which(Wide[,zyg]==1),,drop=FALSE]
  XX2 <- XX[which(Wide[,zyg]==2),,drop=FALSE]
  
##################################################

###}}} setup

###{{{ Mean/Var function

  ##Marginals etc.
  MyData0 <- ExMarg(Y0,XX0,W0,dS0,eqmarg=TRUE,allmarg=allmarg)
  MyData1 <- ExMarg(Y1,XX1,W1,dS1,eqmarg=TRUE,allmarg=allmarg)
  MyData2 <- ExMarg(Y2,XX2,W2,dS2,eqmarg=TRUE,allmarg=allmarg)

  N <- cbind(length(idx0),length(idx1),length(idx2)); 
  N <- cbind(N,
             2*nrow(MyData0$Y0)+if (!pairsonly) NROW(MyData0$Y0_marg) else 0, 
             2*nrow(MyData1$Y0)+if (!pairsonly) NROW(MyData1$Y0_marg) else 0,
             2*nrow(MyData2$Y0)+if (!pairsonly) NROW(MyData2$Y0_marg) else 0,
             NROW(MyData0$Y0),NROW(MyData1$Y0),NROW(MyData2$Y0))

  
  colnames(N) <- c("Total.MZ","Total.DZ","Total.OS","Complete.MZ","Complete.DZ","Complete.OS","Complete pairs.MZ","Complete pairs.DZ","Complete pairs.OS")
  rownames(N) <- rep("",nrow(N))
  if (!OSon) N <- N[,-c(3,6,9),drop=FALSE]
  
  if (samecens & !is.null(weight)) {
    MyData0$W0 <- cbind(apply(MyData0$W0,1,biweight))
    if (!is.null(MyData0$Y0_marg))
      MyData0$W0_marg <- cbind(apply(MyData0$W0_marg,1,biweight))
  }
  if (samecens & !is.null(weight)) {
    MyData1$W0 <- cbind(apply(MyData1$W0,1,biweight))
    if (!is.null(MyData1$Y0_marg))
      MyData1$W0_marg <- cbind(apply(MyData1$W0_marg,1,biweight))
  }
  if (samecens & !is.null(weight)) {
    MyData2$W0 <- cbind(apply(MyData2$W0,1,biweight))
    if (!is.null(MyData2$Y0_marg))
      MyData2$W0_marg <- cbind(apply(MyData2$W0_marg,1,biweight))
  }

  rm(Y0,XX0,W0,Y1,XX1,W1,Y2,XX2,W2)
  

  Sigma <- function(p0) {
    Sigma2 <- NULL
    p0[vidx] <- mytr(p0[vidx])    
    if (ACDU["u"]) {
      pos0 <- ifelse(OSon, plen-2, plen-1)
      Sigma0 <- diag(2) + p0[pos0]*Rm
      Sigma1 <- diag(2) + p0[pos0+1]*Rm
      if (OSon) Sigma2 <- diag(2) + p0[pos0+2]*Rm
    } else {
      ii <- ACDU; ii[4:5] <- FALSE
      pv <- ACDU*1;  pv[ii] <- p0[vidx]
      Sigma0 <- Vm*pv["e"] + pv["a"] + pv["c"] + pv["d"]
      Sigma1 <- Vm*pv["e"] + pv["a"]*Am + pv["c"] + pv["d"]*Dm
      Sigma2 <- Vm*pv["e"] + pv["c"] + (pv["a"] + pv["d"])*Vm +
        pv["os"]*(pv["a"] + pv["d"])*Rm
      if (OSon) {
        dS2 <- dS2.
        dS2[c(1,3),2:3] <- pv["os"]
        dS2[4,2:3] <- pv["a"]+pv["d"]
        dS2 <- dS2[which(c(ACDU[1:3],TRUE)),]
      }
    }
    return(list(Sigma0=Sigma0,Sigma1=Sigma1,Sigma2=Sigma2,dS2=dS2))
  }

  ## p0 <- op$par
  ## ff <- function(p) as.vector(Sigma(p)$Sigma2)
  ## numDeriv::jacobian(ff,p0)
  ## Sigma(p0)$dS2
  ## dmytr(p0[vidx])
  ## Sigma(p0)$dS2[1,]*dmytr(p0[vidx])[1]
  ## Sigma(p0)$dS2[2,]*dmytr(p0[vidx])[2]
  ## Sigma(p0)$dS2[3,]*dmytr(p0[vidx])[3]

###}}} Mean/Var function
  
###{{{ U  

  p0 <- rep(-1,plen); ##p0[vidx] <- 0
  if (!missing(varlink) && is.null(varlink)) p0 <- rep(0.5,plen)
  if (OSon) p0[length(p0)] <- 0.3
  if (type=="u")
    p0[vidx] <- 0.3
  if (!is.null(control$start)) {
    p0 <- control$start
    control$start <- NULL
  } else {
    X <- rbind(MyData0$XX0[,midx0,drop=FALSE],MyData0$XX0[,midx1,drop=FALSE])
    Y <- rbind(MyData0$Y0[,1,drop=FALSE],MyData0$Y0[,2,drop=FALSE])
    g <- suppressWarnings(glm(Y~-1+X,family=binomial(probit)))
    p0[midx] <- coef(g)
  }


  U <- function(p,indiv=FALSE) {
    b0 <- cbind(p[bidx0])
    b1 <- cbind(p[bidx1])
    b2 <- cbind(p[bidx2])
    b00 <- b0; b11 <- b1; b22 <- b2
    if (bound) p[vidx] <- min(p[vidx],20)
    S <- Sigma(p)
    lambda <- eigen(S$Sigma0)$values
    if (any(lambda<1e-12 | lambda>1e9)) stop("Variance matrix out of bounds")
    
    Mu0 <- with(MyData0, cbind(XX0[,midx0,drop=FALSE]%*%b00,
                               XX0[,midx1,drop=FALSE]%*%b00))
    U0 <- with(MyData0, .Call("biprobit0",
                             Mu0,
                             S$Sigma0,dS0,Y0,XX0,W0,!is.null(W0),samecens))

    if (!is.null(MyData0$Y0_marg) &&!pairsonly) {
      mum <- with(MyData0, XX0_marg%*%b00)
      dSmarg <- dS0[,1,drop=FALSE]
       U_marg <- with(MyData0, .Call("uniprobit",
                                   mum,XX0_marg,
                                   S$Sigma0[1,1],t(dSmarg),Y0_marg,
                                   W0_marg,!is.null(W0_marg),TRUE))
      U0$score <- rbind(U0$score,U_marg$score)
      U0$loglik <- c(U0$loglik,U_marg$loglik)
    }

    Mu1 <- with(MyData1, cbind(XX0[,midx0,drop=FALSE]%*%b11,
                               XX0[,midx1,drop=FALSE]%*%b11))

    U1 <- with(MyData1, .Call("biprobit0",
                             Mu1,
                             S$Sigma1,dS1,Y0,XX0,W0,!is.null(W0),samecens))
    if (!is.null(MyData1$Y0_marg) &&!pairsonly) {
      mum <- with(MyData1, XX0_marg%*%b11)
      dSmarg <- dS1[,1,drop=FALSE]
      U_marg <- with(MyData1, .Call("uniprobit",
                                    mum,XX0_marg,
                                    S$Sigma1[1,1],t(dSmarg),Y0_marg,
                                    W0_marg,!is.null(W0_marg),TRUE))
      U1$score <- rbind(U1$score,U_marg$score)
      U1$loglik <- c(U1$loglik,U_marg$loglik)
    }

    U2 <- val2 <- NULL
    if (OSon) {
      Mu2 <- with(MyData2, cbind(XX0[,midx0,drop=FALSE]%*%b22,
                                 XX0[,midx1,drop=FALSE]%*%b22))
      U2 <- with(MyData2, .Call("biprobit0",
                                Mu2,
                                S$Sigma2,S$dS2,Y0,XX0,W0,!is.null(W0),samecens))
      if (!is.null(MyData2$Y0_marg) &&!pairsonly) {
        mum <- with(MyData2, XX0_marg%*%b22)
        dSmarg <- S$dS2[,1,drop=FALSE]
        U_marg <- with(MyData2, .Call("uniprobit",
                                      mum,XX0_marg,
                                      S$Sigma2[1,1],t(dSmarg),Y0_marg,
                                      W0_marg,!is.null(W0_marg),TRUE))
        U2$score <- rbind(U2$score,U_marg$score)
        U2$loglik <- c(U2$loglik,U_marg$loglik)
      }
    }

    
    if (indiv) {
      ll0 <- U0$loglik
      ll1 <- U1$loglik
      val0 <- U0$score[MyData0$id,,drop=FALSE]
      val1 <- U1$score[MyData1$id,,drop=FALSE]
      
      N0 <- length(MyData0$id)
      idxs0 <- seq_len(N0)
      if (length(MyData0$margidx)>0) {
        for (i in seq_len(N0)) {
          idx0 <- which((MyData0$idmarg)==(MyData0$id[i]))+N0
          idxs0 <- c(idxs0,idx0)
          val0[i,] <- val0[i,]+colSums(U0$score[idx0,,drop=FALSE])
        }
        val0 <- rbind(val0, U0$score[-idxs0,,drop=FALSE])
        ll0 <- c(ll0,ll0[-idxs0])
      }
      
      N1 <- length(MyData1$id)
      idxs1 <- seq_len(N1)
      if (length(MyData1$margidx)>0) {
        for (i in seq_len(N1)) {
          idx1 <- which((MyData1$idmarg)==(MyData1$id[i]))+N1
          idxs1 <- c(idxs1,idx1)
          val1[i,] <- val1[i,]+colSums(U1$score[idx1,,drop=FALSE])
        }
        val1 <- rbind(val1, U1$score[-idxs1,,drop=FALSE])
        ll1 <- c(ll1,ll1[-idxs1])        
      }

      if (OSon) {
        ll2 <- U2$loglik
        val2 <- U2$score[MyData2$id,,drop=FALSE]
        N2 <- length(MyData2$id)
        idxs2 <- seq_len(N2)
        if (length(MyData2$margidx)>0) {
          for (i in seq_len(N2)) {
            idx2 <- which((MyData2$idmarg)==(MyData2$id[i]))+N2
            idxs2 <- c(idxs2,idx2)
            val2[i,] <- val2[i,]+colSums(U2$score[idx2,,drop=FALSE])
          }
          val2 <- rbind(val2, U2$score[-idxs2,,drop=FALSE])
          ll2 <- c(ll2,ll2[-idxs2])        
        }
      } 
      
      val <- matrix(0,ncol=plen,nrow=nrow(val0)+nrow(val1) + NROW(val2))
      val[seq_len(nrow(val0)),c(bidx0,vidx0)] <- val0
      val[nrow(val0)+seq_len(nrow(val1)),c(bidx1,vidx1)] <- val1
      if (OSon) {
        val[nrow(val0)+nrow(val1)+seq_len(nrow(val2)),c(bidx2,vidx2)] <- val2
      }
      trp <- dmytr(p[vidx])
      for (i in seq(length(vidx))) {
        val[,vidx[i]] <- val[,vidx[i]]*trp[i]
      }
      attributes(val)$logLik <- c(U0$loglik,U1$loglik,U2$loglik)
      return(val)
      
    }
    val <- numeric(plen)
    val[c(bidx0,vidx0)] <- colSums(U0$score)
    val[c(bidx1,vidx1)] <- val[c(bidx1,vidx1)]+colSums(U1$score)
    if (OSon) val[c(bidx2,vidx2)] <- val[c(bidx2,vidx2)]+colSums(U2$score)
    val[vidx] <- val[vidx]*dmytr(p[vidx])
    attributes(val)$logLik <- sum(U0$loglik)+sum(U1$loglik)+sum(U2$loglik)
    return(val)
  }

###}}} U

###{{{ optim

  if (!missing(p)) return(U(p,indiv=indiv))


  f <- function(p) crossprod(U(p))[1]
  f0 <- function(p) -sum(attributes(U(p))$logLik)
  g0 <- function(p) -as.numeric(U(p))
  h0 <- function(p) crossprod(U(p,indiv=TRUE))

  
  if (!missing(constrain)) {
    freeidx <- is.na(constrain)
    f <- function(p) {      
      p1 <- constrain; p1[freeidx] <- p
      res <- U(p1)[freeidx]
      crossprod(res)[1]
    }
    f0 <- function(p) {
      p1 <- constrain; p1[freeidx] <- p
      -sum(attributes(U(p1))$logLik)
    }
    g0 <- function(p) {
      p1 <- constrain; p1[freeidx] <- p
      -as.numeric(U(p1)[freeidx])
    }
    p0 <- p0[is.na(constrain)]    
  }


  ## Derivatives, Sanity check 
  ## ff <- function(p) attributes(U(p,indiv=FALSE))$logLik
  ## pp <- c(0,-.1,.1,0.5)
  ## numDeriv::grad(ff,pp)
  ## U(pp,indiv=FALSE)

  
  controlstd <- list(hessian=0)
  controlstd[names(control)] <- control
  control <- controlstd
  
  nlminbopt <- intersect(names(control),c("eval.max","iter.max","trace","abs.tol","rel.tol","x.tol","step.min"))
  ucminfopt <- intersect(names(control),c("trace","grtol","xtol","stepmax","maxeval","grad","gradstep","invhessian.lt"))
  optimopt <- names(control) 

  op <- switch(tolower(control$method),
               nlminb=nlminb(p0,f0,gradient=g0,control=control[nlminbopt]),
               optim=optim(p0,fn=f0,gr=g0,control=control[ucminfopt]),
               ucminf=,
               quasi=,
               gradient=ucminf::ucminf(p0,fn=f0,gr=g0,control=control[ucminfopt],hessian=0),
               ## ,
               ## bhhh={
               ##   controlnr <- list(stabil=FALSE,
               ##                     gamma=0.1,
               ##                     gamma2=1,
               ##                     ngamma=5,
               ##                     iter.max=200,
               ##                     epsilon=1e-12,
               ##                     tol=1e-9,
               ##                     trace=1,
               ##                     stabil=FALSE)
               ##   controlnr[names(control)] <- control
               ##   lava:::NR(start=p0,NULL,g0, h0,control=controlnr)
##               },
               ##                 op <- switch(mycontrol$method,
               ##                              ucminf=ucminf(p0,f,control=mycontrol[ucminfopt],hessian=F),
               ##                optim=optim(p0,f,control=mycontrol[ucminfopt],...),
                 nlminb(p0,f,control=control[nlminbopt]))

  if (stderr) {
    UU <- U(op$par,indiv=TRUE)    
    I <- -numDeriv::jacobian(U,op$par)
    tol <- 1e-15
    V <- Inverse(I,tol)
    sqrteig <- attributes(V)$sqrteig
    J <- NULL
    if (robustvar) {
      J <- crossprod(UU)
      V <- V%*%J%*%V
    }
    if (any(sqrteig<tol)) warning("Near-singular covariance matrix (pseudo-inverse used)")
  } else {
    UU <- matrix(NA,ncol=length(op$par),nrow=1)
    I <- J <- V <- matrix(NA,ncol=length(op$par),nrow=length(op$par))
  }

###}}} optim

###{{{ return

  suppressWarnings(cc <- cbind(op$par,sqrt(diag(V))))
  cc <- cbind(cc,cc[,1]/cc[,2],2*(pnorm(abs(cc[,1]/cc[,2]),lower.tail=FALSE)))
  colnames(cc) <- c("Estimate","Std.Err","Z","p-value")
  vnames1 <- NULL
  trnam <- " "

  if (!eqmean) {
    rnam <- rnames1
    rnames1 <- c(paste(rnam,"MZ",sep=trnam),paste(rnam,"DZ",sep=trnam))
    if (OSon) rnames1 <- c(rnames1,paste(rnam,"OS",sep=trnam))

  }

  if (ACDU["u"]) {
    groups <- c("MZ","DZ"); if (OSon) groups <- c(groups,"OS")
    rnames <- c(rnames1,paste(c("atanh(rho)","atanh(rho)"),groups,sep=trnam))
  } else {
      rnames0 <- c("var(A)","var(C)","var(D)")[ACDU[1:3]]
      if (invtrname!="")
          rnames0 <- paste(invtrname,"(",rnames0,")",sep="")
      rnames <- c(rnames1,rnames0)
      if (OSon) rnames <- c(rnames,"atanh(rho(G[OS]))")
  }
  if (!missing(constrain)) rnames <- rnames[freeidx]
  rownames(cc) <- rnames
  rownames(V) <- colnames(V) <- rnames
  S <- Sigma(op$par)

  npar <- list(intercept=attributes(terms(formula))$intercept,
               pred=nrow(attributes(terms(formula))$factor)-1,
               var=sum(ACDU[-4]),
               ACDU=ACDU[-4]*1)
  
  npar[unlist(lapply(npar,length))==0] <- 0
  
  val <- list(coef=cc,vcov=V,score=UU,logLik=attributes(UU)$logLik,opt=op,
              Sigma0=S$Sigma0, Sigma1=S$Sigma1, Sigma2=S$Sigma2,
              dS0=dS0, dS1=dS1, dS2=dS2,
              N=N,
              midx0=midx0, midx1=midx1,
              vidx0=vidx0, vidx1=vidx1, vidx2=vidx2,
              eqmean=eqmean, I=I,J=J, robustvar=robustvar,
              transform=list(tr=mytr, invtr=myinvtr, dtr=dmytr,
                name=trname, invname=invtrname),
              SigmaFun=Sigma,
              npar=npar,
              OS=OSon
              )
  class(val) <- c("bptwin","biprobit")
  return(val)
}

###}}} return


