
##' @export
biprobit.vector <- function(x,id,X=NULL,Z=NULL,
                            weights=NULL,
                            biweight=function(x) { u <- min(x); ifelse(u==0,u,1/u) },
                            ## biweight=function(x) { u <- min(x); ifelse(u==0,0,1/max(u,lava.options()$min.weight)) },
                            eqmarg=TRUE,add=NULL,control=list(),p=NULL,
                            cells,
                            ...) {

    if (is.factor(x)) x <- factor(x,labels=c(0,1))
    else x <- factor(x*1,levels=c(0,1))
    x <- factor(x,levels=c(0,1),labels=c(0,1))
    namX <- colnames(X)
    namZ <- colnames(Z)
    if (is.null(namX)) namX <- "(Intercept)"
    if (is.null(namZ)) namZ <- "(Intercept)"
    namZ <- paste("r:",namZ,sep="")

    resh <- function(x,...,nam,onecol=0,df=1) {
        x <- c(list(x),list(...))
        res <- mis <- c()
        for (i in seq_along(x)) {
            if (length(x[[i]])==0) {
                res <- c(res,list(NULL))
            } else {                
                 if (i%in%df) {
                    xx <- fast.reshape(x[[i]],id=id)
                } else {
                    xx <- as.matrix(fast.reshape(x[[i]],id=id))
                }
                if (i%in%onecol) {
                    xx <- xx[,seq(ncol(x[[i]])),drop=FALSE]
                    colnames(xx) <- colnames(x[[i]])
                }
                idx <- which(is.na(xx),arr.ind=TRUE)                
                if (length(idx)>0) mis <- c(mis,idx[,1])
                res <- c(res,list(xx))
            }
        }
        mis <- unique(mis)
        if (length(mis)>0) {
            for (i in seq_along(res)) {
                if (length(res[[i]])>0) res[[i]] <- res[[i]][-mis,,drop=FALSE]
            }
        }
        if (!missing(nam)) names(res) <- nam
        return(res)
    }
    DD <- resh(data.frame(x),X,Z,weights,nam=c("y","x","z","w"),onecol=3)
    Y00 <- matrix(c(0,0, 1,0, 0,1, 1,1),ncol=2,byrow=TRUE)

    if (is.null(X) && is.null(Z)) {
        pos <- factor(interaction(DD$y))
        ipos <- unique(as.numeric(pos))
        Tab <- rbind(as.vector(table(DD$y))); colnames(Tab) <- c("00","10","01","11")
        Y0 <- Y00
        NN <- as.vector(Tab)
        midx1 <- 1; midx2 <- 2
        X0 <- matrix(1,ncol=2,nrow=4)        
        Z0 <- matrix(1,ncol=1,nrow=4)
        namX <- "(Intercept)"
        namZ <- "r:(Intercept)"
    } else {
        pos2 <- fast.pattern(cbind(DD$x,DD$z))
        pos <- interaction(interaction(DD$y),pos2$group)
        XZ0 <- pos2$pattern; colnames(XZ0) <- c(colnames(DD$x),colnames(DD$z))
        NN2 <- unlist(by(DD$y,pos,nrow))            
        NN <- rep(0,4*nrow(XZ0))
        ipos <- unique(as.numeric(pos))
        NN[ipos] <- NN2[ipos]

        ## tt <- with(DD, by(y,as.list(as.data.frame(cbind(x,z))),FUN=function(x) as.vector(table(x[,1:2])),simplify=FALSE))
        ## XZ0 <- do.call("expand.grid",lapply(attributes(tt)$dimnames,as.numeric))
        ## rem <- which(unlist(lapply(tt,is.null)))
        ## if (length(rem)>0) {
        ##     tt <- tt[-rem]
        ##     XZ0 <- XZ0[-rem,,drop=FALSE]
        ## }
        ## NN <- Reduce("c",tt);
        ## Reduce("cbind",tt)
        suppressWarnings(Tab <- cbind(matrix(NN,ncol=4,byrow=TRUE),XZ0)); colnames(Tab)[1:4] <- c("00","10","01","11")
        
        Y0 <- matrix(rep(t(Y00),length(NN)/4),ncol=2,byrow=TRUE)
        if (is.null(X)) X0 <- matrix(1,ncol=2,nrow=nrow(Y0)) else {
            X0 <- as.matrix(XZ0[rep(seq(nrow(XZ0)),each=4),seq(ncol(X)*2),drop=FALSE])
        }
        if (is.null(Z)) Z0 <- matrix(1,ncol=1,nrow=nrow(Y0)) else {
            Z1 <- as.matrix(XZ0[,ncol(XZ0)-ncol(Z)+seq(ncol(Z)),drop=FALSE])
            Z0 <- Z1[rep(seq(nrow(XZ0)),each=4),,drop=FALSE]
            colnames(Z0) <- colnames(Z1)
        }
    }

    nx <- ncol(X0)/2
    midx1 <- seq(nx)
    midx2 <- midx1+nx
    midx <- seq(2*nx)
    blen <- ifelse(eqmarg,nx,2*nx)
    zlen <- ncol(Z0)      
    plen <- blen+zlen
    
    MyData <- ExMarg(Y0,X0,W0=NULL,NULL,
                     midx1=1,midx2=2,eqmarg=eqmarg,allmarg=FALSE,Z0)
 
    datanh <- function(r) 1/(1-r^2)
    dtanh <- function(z) 4*exp(2*z)/(exp(2*z)+1)^2
    vartr <- tanh
    dvartr <- dtanh; varitr <- atanh
    trname <- "tanh"; itrname <- "atanh"    
    varcompname <- "Tetrachoric correlation"    
    ##msg <- "Variance of latent residual sterm = 1 (standard probit link)"
    msg <- NULL

    model <- list(tr=vartr,name=trname,inv=itrname,invname=itrname,deriv=dvartr,varcompname=varcompname,eqmarg=eqmarg,blen=blen,zlen=zlen,...)
    SigmaFun <- function(p,Z=MyData$Z0,cor=TRUE,...) {
        if (!cor) {
            r <- vartr(p[1])
            Sigma <- matrix(c(1,r,r,1),2)
            attributes(Sigma)$dvartr <- dvartr
            return(Sigma)
        }
        val <- Z%*%p
        dr <- apply(Z,2,function(x) x*dvartr(val))
        structure(list(rho=vartr(val),lp=val,drho=dr),dvartr=dvartr,vartr=vartr)
    }
    

    if (length(weights)<2) {
        w0 <- NN
    } else {        
        w1 <- apply(DD$w,1,biweight)
        w2 <- unlist(by(w1,pos,sum))
        w0 <- rep(0,4*nrow(Tab))
        w0[ipos] <- w2[ipos]
    }

    if (!missing(cells)) {
        idx <- unlist(apply(MyData$Y0,1,function(x) all(x==cells)))
        w0[!idx] <- 0
    }
    if (!is.null(control$start)) {
        p0 <- control$start
        control$start <- NULL
    } else {
        p0 <- rep(0,plen)
        events <- Tab[,2]+Tab[,3]+2*Tab[,4]
        totals <- rowSums(Tab[,1:4,drop=FALSE])
        if (is.null(X)) xx1 <- rep(1,length(totals)) else  xx1 <- Tab[,midx1+4,drop=FALSE]
        b1 <- glm(cbind(events,totals)~-1+xx1,family=binomial("probit"))
        p0[midx1] <- coef(b1)
        if (!eqmarg) {
            if (is.null(X)) xx2 <- rep(1,length(totals)) else  xx2 <- Tab[,midx2+4,drop=FALSE]
            b2 <- glm(cbind(events,totals)~-1+xx2,family=binomial("probit"))
            p0[midx2] <- coef(b2)
        }     
    }
   
    U <- function(p,w0) {
        val <- Ubiprobit(p,SigmaFun,eqmarg,nx,MyData,indiv=TRUE)
        logl <- w0*as.vector(attributes(val)$logLik)
        score <- apply(val,2,function(x) w0*x)
        return(structure(score,logLik=logl))
    }


    f0 <- function(p) -sum(attributes(U(p,w0))$logLik)
    g0 <- function(p) -as.numeric(colSums(U(p,w0)))   
    if (!is.null(p)) op <- list(par=p) else {
        suppressWarnings(op <- nlminb(p0,f0,gradient=g0,control=control))
    }
    
    iI <- Inverse(numDeriv::jacobian(g0,op$par))
    V <- iI
    UU <- U(op$par,w0)
    logLik <- sum(attributes(UU)$logLik)
    UU <- U(op$par,1)    
    if (length(weights)>1) {
        UU <- apply(UU[pos,],2,function(x) w1*x)
    } else {
        UU <- UU[pos,]
    }
    meat <- crossprod(UU)
    if (length(weights>1)) V <- iI%*%meat%*%iI
    
    mycall <- match.call()    
    cc <- cbind(op$par,sqrt(diag(V)))
    cc <- cbind(cc,cc[,1]/cc[,2],2*(pnorm(abs(cc[,1]/cc[,2]),lower.tail=FALSE)))
    colnames(cc) <- c("Estimate","Std.Err","Z","p-value")
    if (!eqmarg)
        rownames(cc) <- c(paste(namX,rep(c(1,2),each=length(namX)),sep="."),
                          namZ)
    else
    rownames(cc) <- c(namX,namZ)
    rownames(V) <- colnames(V) <- rownames(cc)

    npar <- list(intercept=1,
                 pred=blen-2+eqmarg)
    if (!eqmarg) npar <- lapply(npar,function(x) x*2)
    npar$var <- zlen
    N <- with(MyData, c(pairs=sum(NN)))
    
    val <- c(list(coef=cc,
                  N=N,
                  vcov=V,
                  bread=iI,
                  score=rbind(UU),
                  logLik=logLik,
                  opt=op,
                  call=mycall,
                  model=model,
                  msg=msg,
                  table=Tab,
                  npar=npar,
                  SigmaFun=SigmaFun),add)
    class(val) <- "biprobit"
    return(val)
}


##' Bivariate Probit model
##' 
##' @export
##' @aliases biprobit biprobit.vector biprobit.time
##' @param x formula (or vector)
##' @param data data.frame
##' @param id The name of the column in the dataset containing the cluster id-variable.
##' @param rho Formula specifying the regression model for the dependence parameter
##' @param num Optional name of order variable
##' @param strata Strata
##' @param eqmarg If TRUE same marginals are assumed (exchangeable)
##' @param indep Independence
##' @param weights Weights
##' @param biweight Function defining the bivariate weight in each cluster
##' @param samecens Same censoring
##' @param randomeffect If TRUE a random effect model is used (otherwise correlation parameter is estimated allowing for both negative and positive dependence)
##' @param vcov Type of standard errors to be calculated
##' @param pairs.only Include complete pairs only?
##' @param allmarg Should all marginal terms be included
##' @param control Control argument parsed on to the optimization routine. Starting values may be parsed as '\code{start}'.
##' @param messages Control amount of messages shown 
##' @param constrain Vector of parameter constraints (NA where free). Use this to set an offset.
##' @param table Type of estimation procedure
##' @param p Parameter vector p in which to evaluate log-Likelihood and score function
##' @param ... Optional arguments
##' @examples
##' data(prt)
##' prt0 <- subset(prt,country=="Denmark")
##' a <- biprobit(cancer~1+zyg, ~1+zyg, data=prt0, id="id")
##' b <- biprobit(cancer~1+zyg, ~1+zyg, data=prt0, id="id",pairs.only=TRUE)
##' predict(b,newdata=Expand(prt,zyg=c("MZ")))
##' predict(b,newdata=Expand(prt,zyg=c("MZ","DZ")))
##' 
##' \donttest{ ## Reduce Ex.Timings
##' m <- lvm(c(y1,y2)~x)
##' covariance(m,y1~y2) <- "r"
##' constrain(m,r~x+a+b) <- function(x) tanh(x[2]+x[3]*x[1])
##' distribution(m,~x) <- uniform.lvm(a=-1,b=1)
##' ordinal(m) <- ~y1+y2
##' d <- sim(m,1000,p=c(a=0,b=-1)); d <- d[order(d$x),]
##' dd <- fast.reshape(d)
##' 
##' a <- biprobit(y~1+x,rho=~1+x,data=dd,id="id")
##' summary(a, mean.contrast=c(1,.5), cor.contrast=c(1,.5))
##' with(predict(a,data.frame(x=seq(-1,1,by=.1))), plot(p00~x,type="l"))
##' 
##' pp <- predict(a,data.frame(x=seq(-1,1,by=.1)),which=c(1))
##' plot(pp[,1]~pp$x, type="l", xlab="x", ylab="Concordance", lwd=2, xaxs="i")
##' confband(pp$x,pp[,2],pp[,3],polygon=TRUE,lty=0,col=Col(1))
##' 
##' pp <- predict(a,data.frame(x=seq(-1,1,by=.1)),which=c(9)) ## rho
##' plot(pp[,1]~pp$x, type="l", xlab="x", ylab="Correlation", lwd=2, xaxs="i")
##' confband(pp$x,pp[,2],pp[,3],polygon=TRUE,lty=0,col=Col(1))
##' with(pp, lines(x,tanh(-x),lwd=2,lty=2))
##' 
##' xp <- seq(-1,1,length.out=6); delta <- mean(diff(xp))
##' a2 <- biprobit(y~1+x,rho=~1+I(cut(x,breaks=xp)),data=dd,id="id")
##' pp2 <- predict(a2,data.frame(x=xp[-1]-delta/2),which=c(9)) ## rho
##' confband(pp2$x,pp2[,2],pp2[,3],center=pp2[,1])
##' 
##' 
##' }
##' 
##' ## Time
##' \dontrun{
##'     a <- biprobit.time(cancer~1, rho=~1+zyg, id="id", data=prt, eqmarg=TRUE,
##'                        cens.formula=Surv(time,status==0)~1,
##'                        breaks=seq(75,100,by=3),fix.censweights=TRUE)
##' 
##'     a <- biprobit.time2(cancer~1+zyg, rho=~1+zyg, id="id", data=prt0, eqmarg=TRUE,
##'                        cens.formula=Surv(time,status==0)~zyg,
##'                        breaks=100)
##' 
##'     a1 <- biprobit.time2(cancer~1, rho=~1, id="id", data=subset(prt0,zyg=="MZ"), eqmarg=TRUE,
##'                        cens.formula=Surv(time,status==0)~1,
##'                        breaks=100,pairs.only=TRUE)
##' 
##'     a2 <- biprobit.time2(cancer~1, rho=~1, id="id", data=subset(prt0,zyg=="DZ"), eqmarg=TRUE,
##'                         cens.formula=Surv(time,status==0)~1,
##'                         breaks=100,pairs.only=TRUE)
##' 
##'     prt0$trunc <- prt0$time*runif(nrow(prt0))*rbinom(nrow(prt0),1,0.5)
##'     a3 <- biprobit.time(cancer~1, rho=~1, id="id", data=subset(prt0,zyg=="DZ"), eqmarg=TRUE,
##'                         cens.formula=Surv(trunc,time,status==0)~1,
##'                         breaks=100,pairs.only=TRUE)
#
##' 
##'     plot(a,which=3,ylim=c(0,0.1))
##' }
biprobit <- function(x, data, id, rho=~1, num=NULL, strata=NULL, eqmarg=TRUE,
                             indep=FALSE, weights=NULL, 
                             biweight,
                             samecens=TRUE, randomeffect=FALSE, vcov="robust",
                             pairs.only=FALSE,                             
                             allmarg=samecens&!is.null(weights),
                             control=list(trace=0),
                             messages=1, constrain=NULL,                     
                             table=pairs.only,
                             p=NULL,...) {

  mycall <- match.call()
  if (missing(biweight)) {
      biweight <- mycall$biweight <- function(x) { u=min(x); ifelse(u==0,0,1/min(u)) }
  }
  formulaId <- unlist(Specials(x,"cluster"))
  if (is.null(formulaId)) {
    formulaId <- unlist(Specials(x,"id"))
  }
  ##  formulaOffset <- unlist(Specials(x,"offset"))
  formulaStrata <- unlist(Specials(x,"strata"))
  formulaSt <- paste("~.-cluster(",formulaId,")",
                     "-id(",formulaId,")",
                     "-strata(",paste(formulaStrata,collapse="+"),")",sep="")
  formula <- update(x,formulaSt)

  if (!is.null(formulaId)) {
    id <- formulaId
    mycall$id <- id
  }
  if (!is.null(formulaStrata)) strata <- formulaStrata
  mycall$x <- formula

  if (!is.null(strata)) {
    dd <- split(data,interaction(data[,strata]))
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
  
  if (missing(id)) {    
    if (!is.null(weights)) {
      weights <- data[,weights]
      return(glm(formula,data=data,family=binomial(probit),weights=weights,...))
    }
    return(glm(formula,data=data,family=binomial(probit),...))    
  }

  yx <- getoutcome(formula)

  if (pairs.only) {
      X <- Z <- NULL
      zf <- getoutcome(rho); if (length(attr(zf,"x"))>0) Z <- model.matrix(rho,data);
      if (table && NCOL(Z)<10 && length(unique(sample(Z,min(1000,length(Z)))))<10) { ## Not quantitative?
          if (!is.null(weights)) weights <- data[,weights]
          if (length(attr(yx,"x")>0)) X <- model.matrix(x,data);
          return(biprobit.vector(data[,yx],X=X,Z=Z,id=data[,id],weights,biweight=biweight,eqmarg=eqmarg,add=list(formula=formula,rho.formula=rho),control=control,p=p,...))
      }
  }
  
  mycall <- match.call()
  DD <- procdatabiprobit(formula,data,id,num=num,weights=weights,pairs.only=pairs.only,rho,...)
  rnames1 <- DD$rnames1

  nx <- length(rnames1)
  ## if (nx==0) stop("Zero design not allowed")  
  midx1 <- seq(nx)
  midx2 <- midx1+nx
  midx <- seq(2*nx)
  blen <- ifelse(eqmarg,nx,2*nx)
  zlen <- ncol(DD$Z0)      
  plen <- blen+zlen
  
  datanh <- function(r) 1/(1-r^2)
  dtanh <- function(z) 4*exp(2*z)/(exp(2*z)+1)^2
  vartr <- tanh
  dvartr <- dtanh; varitr <- atanh
  trname <- "tanh"; itrname <- "atanh"    
  Sigma1 <- diag(2)  
  Sigma2 <- matrix(c(0,1,1,0),2,2)
  dS0 <- rbind(c(0,1,1,0))
  
  varcompname <- "Tetrachoric correlation"
  msg <- NULL
  if (randomeffect) msg <- "Variance of latent residual term = 1 (standard probit link)"
  if (randomeffect) {
      dS0 <- rbind(rep(1,4))
      vartr <- dvartr <- exp; inv <- log
      trname <- "exp"; itrname <- "log"
      Sigma2 <- 1
      varcompname <- NULL
  }
  model <- list(tr=vartr,name=trname,inv=itrname,invname=itrname,deriv=dvartr,varcompname=varcompname,dS=dS0,eqmarg=eqmarg,randomeffect=randomeffect,blen=blen,zlen=zlen)

  MyData <- with(DD,ExMarg(Y0,XX0,W0,dS0,midx1,midx2,eqmarg=eqmarg,allmarg=allmarg,Z0))
  if (samecens & !is.null(weights)) {
      MyData$W0 <- cbind(apply(MyData$W0,1,biweight))
      if (!is.null(MyData$Y0_marg)) {
          MyData$W0_marg <- cbind(apply(MyData$W0_marg,1,biweight))
      }
  }
  
  SigmaFun <- function(p,Z=MyData$Z0,cor=!randomeffect,...) {
      if (!cor) {
          r <- vartr(p[1])
          Sigma <- matrix(c(1,r,r,1),2)
          if (indep) Sigma <- diag(2)
          attributes(Sigma)$dvartr <- dvartr
          return(Sigma)
      }
      val <- Z%*%p
      dr <- apply(Z,2,function(x) x*dvartr(val))
      structure(list(rho=vartr(val),lp=val,drho=dr),dvartr=dvartr,vartr=vartr)
  }

  U <- function(p,indiv=FALSE) {
      gamma <- p[seq(zlen)+blen]
      ##if (bound) gamma <- min(gamma,20)
      Sigma <- SigmaFun(gamma)
      if (randomeffect) {
          lambda <- eigen(Sigma)$values
          if (any(lambda<1e-12 | lambda>1e9)) stop("Variance matrix out of bounds")
      }
      Mu_marg <- NULL
      if (eqmarg) {
          B <- cbind(p[midx1])
          Mu <- with(MyData,
                     cbind(XX0[,midx1,drop=FALSE]%*%B,XX0[,midx2,drop=FALSE]%*%B))     
          if (!is.null(MyData$Y0_marg)) 
              Mu_marg <- with(MyData, XX0_marg%*%B)
      } else {
          B1 <- cbind(p[midx1])
          B2 <- cbind(p[midx2])
          Mu <- with(MyData,
                     cbind(XX0[,midx1,drop=FALSE]%*%B1,XX0[,midx2,drop=FALSE]%*%B2))
          if (!is.null(MyData$Y0_marg))
              Mu_marg <- with(MyData, rbind(X0_marg1%*%B1,X0_marg2%*%B2))
      }

      if (randomeffect) {
          U <- with(MyData, .Call("biprobit2",
                                  Mu,XX0,
                                  Sigma,dS0*attributes(Sigma)$dvartr(p[plen]),Y0,W0,
                                  !is.null(W0),TRUE,eqmarg,FALSE))
      } else {
          U <- with(MyData, .Call("biprobit2",
                                  Mu,XX0,
                                  Sigma$rho,Sigma$drho,
                                  Y0,W0,
                                  !is.null(W0),TRUE,eqmarg,TRUE))
      }
      
      if (!is.null(MyData$Y0_marg)) {
          if (randomeffect) {
              U_marg <- with(MyData, .Call("uniprobit",
                                           Mu_marg,XX0_marg,
                                           Sigma[1,1],dS0_marg*attributes(Sigma)$dvartr(p[plen]),Y0_marg,
                                           W0_marg,!is.null(W0_marg),TRUE))
          } else {
              U_marg0 <- matrix(0,length(MyData$Y0_marg),ncol=plen)
              U_marg <- with(MyData, .Call("uniprobit",
                                           Mu_marg,XX0_marg,
                                           1,matrix(ncol=0,nrow=0),Y0_marg,
                                           W0_marg,!is.null(W0_marg),TRUE))
              U_marg0[,seq(blen)] <-  U_marg[[1]]
              U_marg[[1]] <- U_marg0
          }
          U$score <- rbind(U$score,U_marg$score)
          U$loglik <- c(U$loglik,U_marg$loglik)
      }
      
      if (indiv) {
          val <- U$score
          if (!is.null(MyData$idmarg) && !pairs.only) {
              val <- with(MyData, cluster.index(c(id,idmarg),mat=U$score))
          }          
          ## val <- U$score[MyData$id,,drop=FALSE]
          ## N <- length(MyData$id)
          ## idxs <- seq_len(N)
          ## for (i in seq_len(N)) {
          ##     idx <- which((MyData$idmarg)==(MyData$id[i]))+N
          ##     idxs <- c(idxs,idx)
          ##     val[i,] <- val[i,]+colSums(U$score[idx,,drop=FALSE])
          ## }      
          ## val <- rbind(val, U$score[-idxs,,drop=FALSE])
          attributes(val)$logLik <- U$loglik
          return(val)
      }
      val <- colSums(U$score)
      attributes(val)$logLik <- sum(U$loglik)
      return(val)
  }  
  
  p0 <- rep(0,plen)  
  if (!is.null(control$start)) {
    p0 <- control$start
    control$start <- NULL
  } else {
      g <- suppressWarnings(glm(formula,data,family=binomial(probit)))
      p0[midx1] <- coef(g)
      if (!eqmarg) p0[midx2] <- coef(g)    
  }

  f <- function(p) crossprod(U(p))[1]
  f0 <- function(p) -sum(attributes(U(p))$logLik)
  g0 <- function(p) -as.numeric(U(p))
  h0 <- function(p) crossprod(U(p,indiv=TRUE))

  if (!is.null(constrain)) {
      if (length(constrain)!=length(p0)) stop("Wrong length of constraints (should be NA at positions not to be fixed)")
      fix <- which(!is.na(constrain))
      free <- which(is.na(constrain))
      p0 <- p0[free]
      U0 <- U
      U <- function(p,indiv=FALSE) {
          p1 <- constrain
          p1[free] <- p
          res <- U0(p1,indiv)
          if (is.matrix(res)) {              
              return(structure(res[,free,drop=FALSE],logLik=attributes(res)$logLik))
          } 
          return(structure(res[free],logLik=attributes(res)$logLik))
      }
  }
  
  if (is.null(control$method)) {
    ##    control$method <- ifelse(samecens & !is.null(weights), "bhhh","quasi")
    control$method <- "quasi"
  }
  control$method <- tolower(control$method)
  
  if (is.null(p)) {
      if (control$method=="score") {
          control$method <- NULL
          op <- nlminb(p0,f,control=control,...)
      } else if (control$method=="quasi") {
          control$method <- NULL
          suppressWarnings(op <- nlminb(p0,f0,gradient=g0,control=control))
          ## }
  ## else if (control$method=="bhhh") {
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
  ##   op <- lava:::NR(start=p0,NULL,g0, h0,control=controlnr)
      } else {
          control$method <- NULL
          op <- nlminb(p0,f0,control=control)
      }
  } else op <- list(par=p)

  UU <- U(op$par,indiv=TRUE)
  J <- crossprod(UU)
  ##  iJ <- Inverse(J)
  iI <- Inverse(-numDeriv::jacobian(U,op$par))
  V <- switch(vcov,
              robust=,
              sandwich=iI%*%J%*%iI,##iJ%*%I%*%iJ,
              score=,
              outer=Inverse(J),
              hessian=iI              
              )
  
  cc <- cbind(op$par,sqrt(diag(V)))
  cc <- cbind(cc,cc[,1]/cc[,2],2*(pnorm(abs(cc[,1]/cc[,2]),lower.tail=FALSE)))
  if (!is.null(constrain)) {
      cc0 <- matrix(NA,nrow=length(constrain),ncol=4)
      cc0[free,] <- cc
      cc0[fix,1] <- constrain[fix]
      cc <- cc0
      V0 <- matrix(0,nrow=length(constrain),ncol=length(constrain))
      V0[free,free] <- V
      V <- V0
  }
  colnames(cc) <- c("Estimate","Std.Err","Z","p-value")

  p1 <- "("; p2 <- ")"
  if (itrname=="log") rhonam <- "U" else {
      rhonam <- DD$znames
      p1 <- p2 <- ""; itrname <- "r:" 
  }
  if (!eqmarg)
    rownames(cc) <- c(paste(rnames1,rep(c(1,2),each=length(rnames1)),sep="."),
                      paste(itrname,p1,rhonam,p2,sep=""))
  else
    rownames(cc) <- c(rnames1,paste(itrname,p1,rhonam,p2,sep=""))
  rownames(V) <- colnames(V) <- rownames(cc)

  npar <- list(intercept=attributes(terms(formula))$intercept,
              pred=nrow(attributes(terms(formula))$factor)-1)
  if (!eqmarg) npar <- lapply(npar,function(x) x*2)
  npar$var <- 1##nrow(cc)-sum(unlist(npar))
  N <- with(MyData, c(n=nrow(XX0)*2+length(margidx), pairs=nrow(XX0)))
  val <- list(coef=cc,N=N,vcov=V,bread=iI,score=UU,logLik=attributes(UU)$logLik,opt=op, call=mycall, model=model,msg=msg,npar=npar,
              SigmaFun=SigmaFun,rho.formula=rho,formula=formula,constrain=constrain)
  class(val) <- "biprobit"
  return(val)
}



procdatabiprobit <- function(formula,data,id,num=NULL,weights=NULL,pairs.only=FALSE,rho=~1,...) {

    data <- data[order(data[,id]),]
    idtab <- table(data[,id])
    if (pairs.only) {
        data <- data[which(as.character(data[,id])%in%names(idtab)[idtab==2]),]
        idtab <- table(data[,id])
    }

    ff <- paste(as.character(formula)[3],"+",
                paste(c(id,num),collapse="+"))
    yvar <- paste(deparse(formula[[2]]),collapse="")
    if (!is.null(weights))
        ff <- paste(weights,"+",ff)
    ff <- paste("~",yvar,"+",ff)

    if (is.logical(data[,yvar])) data[,yvar] <- data[,yvar]*1
    if (is.factor(data[,yvar])) data[,yvar] <- as.numeric(data[,yvar])-1  

    formula0 <- as.formula(ff)  
    opt <- options(na.action="na.pass")
    Data <- model.matrix(formula0,data)
    options(opt)
    rnames1 <- setdiff(colnames(Data),c(yvar,num,id,weights))
    X0 <- as.matrix(Data[,rnames1])

    ex <- 1+!is.null(num)
    rho <- update(rho,paste("~.+",
                           paste(c(id,num),collapse="+")))
    Z0 <- model.matrix(rho,data);
    znames <- setdiff(colnames(Z0),c(id,num))
    znames1 <- paste(znames,1,sep="")
    Z0 <- as.matrix(subset(fast.reshape(Z0,id=id),select=znames1))
    colnames(Z0) <- znames1
    Wide <- fast.reshape(as.data.frame(Data),id=id,num=num,sep=".",labelnum=TRUE)
    
    W0 <- NULL
    yidx <- paste(yvar,1:2,sep=".")
    rmidx <- c(id,yidx)
    if (!is.null(weights)) {
        W <- cbind(data[,weights])
        widx <- paste(weights,1:2,sep=".")
        W0 <- as.matrix(Wide[,widx])
        rmidx <- c(rmidx,widx)
    }   
    Y0 <- as.matrix(Wide[,yidx])
    XX0 <- as.matrix(Wide[,setdiff(colnames(Wide),rmidx)])
    XX0[is.na(XX0)] <- 0

    list(Y0=Y0,XX0=XX0,W0=W0,Z0=Z0,znames=znames,rnames1=rnames1)
}



Ubiprobit <- function(p,Rho,eqmarg,nx,MyData,indiv=FALSE) {
    midx1 <- seq(nx)
    midx2 <- midx1+nx
    midx <- seq(2*nx)
    blen <- ifelse(eqmarg,nx,2*nx)
    zlen <- length(p)-blen
    plen <- blen+zlen
    gamma <- p[blen+seq(zlen)]
    Sigma <- Rho(gamma)
    ## lambda <- eigen(Sigma)$values
    ## if (any(lambda<1e-12 | lambda>1e9)) stop("Variance matrix out of bounds")
    Mu_marg <- NULL
    if (eqmarg) {
        B <- cbind(p[midx1])
        Mu <- with(MyData,
                   cbind(XX0[,midx1,drop=FALSE]%*%B,XX0[,midx2,drop=FALSE]%*%B))     
        ##      Mu <- with(MyData, matrix(X0%*%B,ncol=2,byrow=TRUE))
        if (!is.null(MyData$Y0_marg))
            Mu_marg <- with(MyData, XX0_marg%*%B)
    } else {
        B1 <- cbind(p[midx1])
        B2 <- cbind(p[midx2])
        Mu <- with(MyData,
                   cbind(XX0[,midx1,drop=FALSE]%*%B1,XX0[,midx2,drop=FALSE]%*%B2))
        if (!is.null(MyData$Y0_marg))
            Mu_marg <- with(MyData, rbind(X0_marg1%*%B1,X0_marg2%*%B2))
    }

    U <- with(MyData, .Call("biprobit2",
                            Mu,XX0,
                            Sigma$rho,Sigma$drho,Y0,W0,
                            !is.null(W0),TRUE,eqmarg,TRUE))
    
    if (!is.null(MyData$Y0_marg)) {
        U_marg0 <- matrix(0,length(MyData$Y0_marg),ncol=plen)
        U_marg <- with(MyData, .Call("uniprobit",
                                     Mu_marg,XX0_marg,
                                     1,matrix(ncol=0,nrow=0),Y0_marg,
                                     W0_marg,!is.null(W0_marg),TRUE))
        U_marg0[,seq(blen)] <-  U_marg[[1]]
        U_marg[[1]] <- U_marg0
        
        U$score <- rbind(U$score,U_marg$score)
        U$loglik <- c(U$loglik,U_marg$loglik)
    }

    if (indiv) {
        val <- with(MyData, cluster.index(c(id,idmarg),mat=U$score))
        attributes(val)$logLik <- U$loglik
        return(val)
    }
    val <- colSums(U$score)
    attributes(val)$logLik <- sum(U$loglik)
    return(val)
}


