##' Dose response calculation for binomial regression models
##'
##' @title Dose response calculation for binomial regression models
##' @param model Model object or vector of parameter estimates
##' @param intercept Index of intercept parameters
##' @param slope Index of intercept parameters
##' @param prob Index of mixture parameters (only relevant for
##' \code{zibreg} models)
##' @param x Optional weights
##' length(x)=length(intercept)+length(slope)+length(prob)
##' @param level Probability at which level to calculate dose
##' @param ci.level Level of confidence limits
##' @param vcov Optional estimate of variance matrix of parameter
##' estimates
##' @param family Optional distributional family argument
##' @param EB Optional ratio of treatment effect and adverse effects
##' used to find optimal dose (regret-function argument)
##' @author Klaus K. Holst
##' @export
PD <- function(model,intercept=1,slope=2,prob=NULL,x,level=0.5,
               ci.level=0.95,vcov,family, EB=NULL) {
    if (is.vector(model)) {
        beta <- model
        if (missing(vcov)) stop("vcov argument needed")
        if (missing(family)) stop("family argument needed")
    } else beta <- coef(model)
    if (missing(vcov)) vcov <- stats::vcov(model)
    if (missing(family)) family <- stats::family(model)
  N <- length(intercept)+length(slope)+length(prob)
  if (length(intercept)<length(beta)) {
    B.intercept <- rep(0,length(beta));
    if (!missing(x)) {
      if (length(x)!=N) stop("x should be of same length as the total length of 'intercept','slope','prob'")
      B.intercept[intercept] <- x[seq_len(length(intercept))]
    } else B.intercept[intercept] <- 1
  } else {
    B.intercept <- intercept
  }
  if (length(slope)<length(beta)) {
    B.slope <- rep(0,length(beta));
    if (!missing(x))
      B.slope[slope] <- x[length(intercept)+seq_len(length(slope))]
    else
      B.slope[slope] <- 1
  } else {
    B.slope <- slope
  }
  if (!is.null(prob)) {
    if (length(prob)<length(beta)) {
      B.prob <- rep(0,length(beta));
      if (!missing(x))
        B.prob[prob] <- x[length(intercept)+length(slope)+seq_len(length(prob))]
      else
        B.prob[prob] <- 1
    } else {
      B.prob <- prob
    }
  }
  if (is.null(prob)) B.prob <- NULL
  B <- rbind(B.intercept,B.slope,B.prob)
  S <- B%*%vcov%*%t(B)
  b <- as.vector(B%*%beta)

  f <- function(b) {
    mylevel <- level
    if (!is.null(EB)) {
      if (is.null(prob)) stop("Index of mixture-probability parameters needed")
      pi0 <- family$linkinv(b[3])
      mylevel <- 1-(1-pi0)/pi0*(EB)/(1-EB)
    }
    return(structure((family$linkfun(mylevel)-b[1])/b[2],level=mylevel))
  }

  xx <- f(b)
  Dxx <- -1/b[2]*rbind(1,xx)
  if (!is.null(EB))
    Dxx <- numDeriv::grad(f,b)
  se <- diag(t(Dxx)%*%S%*%Dxx)^0.5
  res <- cbind(Estimate=xx,"Std.Err"=se)
  alpha <- 1-ci.level
  alpha.str <- paste(c(alpha/2,1-alpha/2)*100,"",sep="%")
  res <- cbind(res,res[,1]-qnorm(1-alpha/2)*res[,2],res[,1]+qnorm(1-alpha/2)*res[,2])
  colnames(res)[3:4] <- alpha.str
  rownames(res) <- paste0(round(1000*attributes(xx)$level)/10,"%")
  structure(res,b=b)
}


TN.zibreg <- function(object,data=model.frame(object),p=coef(object),intercept=1,slope=2,alpha=0.95,...) {
    pp <- predict(object,link=FALSE,p=p,newdata=data)
    X <- attributes(pp)$grad$beta
    Z <- attributes(pp)$grad$gamma
    db1 <- db2 <- matrix(0,nrow(X),ncol(X))
    db1[,intercept] <- X[,intercept]
    db2[,slope[1]] <- 1; db2[,slope[-1]] <- X[,slope[-1]]
    b1 <- as.vector(db1%*%p[object$beta.idx])
    b2 <- as.vector(db2%*%p[object$beta.idx])
    ginv <-  object$family$linkinv
    dginv <- object$family$mu.eta ## D[linkinv]
    g <- object$family$linkfun
    dg <- function(x) 1/dginv(g(x)) ## Dh^-1 = 1/(h'(h^-1(x)))
    pi0 <- ginv(pp[,2])
    A2 <- dginv(pp[,2])
    dpi0 <- rbind(apply(Z,2,function(z) A2*z))
    h <- function(pi0) (alpha+pi0-1)/(alpha*pi0)
    dh <- function(pi0) (1-alpha)/(alpha*pi0^2)
    lev <- h(pi0)
    eta <- g(lev)
    detad2 <- rbind(apply(dpi0,2,function(z) dg(lev)*dh(pi0)*z))
    val <- (eta-b1)/b2
    dvald1 <- -(db1+db2*val)/b2
    return(structure(val,grad=cbind(dvald1,detad2/b2),varnames="theta"))
  ##  structure(g(coef(object)),grad=grad(g,coef(object)))
}



##' Regression model for binomial data with unkown group of immortals (zero-inflated binomial regression)
##'
##' @title Regression model for binomial data with unkown group of immortals
##' @param formula Formula specifying
##' @param formula.p Formula for model of disease prevalence
##' @param data data frame
##' @param family Distribution family (see the help page \code{family})
##' @param offset Optional offset
##' @param start Optional starting values
##' @param var Type of variance (robust, expected, hessian, outer)
##' @param ... Additional arguments to lower level functions
##' @author Klaus K. Holst
##' @export
##' @examples
##'
##' ## Simulation
##' n <- 2e3
##' x <- runif(n,0,20)
##' age <- runif(n,10,30)
##' z0 <- rnorm(n,mean=-1+0.05*age)
##' z <- cut(z0,breaks=c(-Inf,-1,0,1,Inf))
##' p0 <- lava:::expit(model.matrix(~z+age) %*% c(-.4, -.4, 0.2, 2, -0.05))
##' y <- (runif(n)<lava:::tigol(-1+0.25*x-0*age))*1
##' u <- runif(n)<p0
##' y[u==0] <- 0
##' d <- data.frame(y=y,x=x,u=u*1,z=z,age=age)
##' head(d)
##'
##' ## Estimation
##' e0 <- zibreg(y~x*z,~1+z+age,data=d)
##' e <- zibreg(y~x,~1+z+age,data=d)
##' compare(e,e0)
##' e
##' PD(e0,intercept=c(1,3),slope=c(2,6))
##'
##' B <- rbind(c(1,0,0,0,20),
##'            c(1,1,0,0,20),
##'            c(1,0,1,0,20),
##'            c(1,0,0,1,20))
##' prev <- summary(e,pr.contrast=B)$prevalence
##'
##' x <- seq(0,100,length.out=100)
##' newdata <- expand.grid(x=x,age=20,z=levels(d$z))
##' fit <- predict(e,newdata=newdata)
##' plot(0,0,type="n",xlim=c(0,101),ylim=c(0,1),xlab="x",ylab="Probability(Event)")
##' count <- 0
##' for (i in levels(newdata$z)) {
##'   count <- count+1
##'   lines(x,fit[which(newdata$z==i)],col="darkblue",lty=count)
##' }
##' abline(h=prev[3:4,1],lty=3:4,col="gray")
##' abline(h=prev[3:4,2],lty=3:4,col="lightgray")
##' abline(h=prev[3:4,3],lty=3:4,col="lightgray")
##' legend("topleft",levels(d$z),col="darkblue",lty=seq_len(length(levels(d$z))))
zibreg <- function(formula,formula.p=~1,data,family=stats::binomial(),offset=NULL,start,var="hessian",...) {
  md <- cbind(model.frame(formula,data),model.frame(formula.p,data))
  y <- md[,1]
  X <- model.matrix(formula,data)
  Z <- model.matrix(formula.p,data)
  beta.idx <- seq(ncol(X)); gamma.idx <- seq(ncol(Z))+ncol(X)
  if (missing(start)) start <- rep(0,ncol(X)+ncol(Z))
  op <- nlminb(start,function(x)
               -zibreg_logL(x[beta.idx],x[gamma.idx],y,X,Z),
               gradient=function(x)
               -zibreg_score(x[beta.idx],x[gamma.idx],y,X,Z),...)
  beta <- op$par[beta.idx]; gamma <- op$par[gamma.idx]
  cc <- c(beta,gamma)
  names(cc) <- c(colnames(X),paste0("pr:",colnames(Z)))
  bread <- Inverse(zibreg_information(beta,gamma,y,X,Z,offset,type="hessian",...))
  if (tolower(var[1])%in%c("robust","sandwich")) {
      meat <- zibreg_information(beta,gamma,y,X,Z,offset,family,type="outer",...)
      V <- bread%*%meat%*%bread
  } else {
      V <- bread
  }
  colnames(V) <- rownames(V) <- names(cc)
  res <- list(coef=cc,opt=op,beta=beta,gamma=gamma,
              beta.idx=beta.idx,gamma.idx=gamma.idx,bread=bread,
              formula=formula,formula.p=formula.p, y=y, X=X, Z=Z, offset=offset, vcov=V, model.frame=md,family=family)
  class(res) <- "zibreg"
  res$fitted.values <- predict(res)
  return(res)
}

##' @export
vcov.zibreg <- function(object,...) object$vcov

##' @export
coef.zibreg <- function(object,...) object$coef

##' @export
family.zibreg <- function(object,...) object$family

##' @export
predict.zibreg <- function(object,p=coef(object),gamma,newdata,link=TRUE,subdist=FALSE,...) {
  newf <- as.formula(paste("~",as.character(object$formula)[3]))
  if (missing(newdata)) {
    X <- object$X; Z <- object$Z
  } else {
    X <- model.matrix(newf,newdata)
    Z <- model.matrix(object$formula.p,newdata)
  }
  if (length(p)==length(object$beta)+length(object$gamma)) {
    gamma <- p[object$gamma.idx]
    p <- p[object$beta.idx]
  }
  g <- object$family$linkfun
  ginv <- object$family$linkinv
  dginv <- object$family$mu.eta ## D[linkinv]
  Xbeta <- as.vector(X%*%p)
  Zgamma <- as.vector(Z%*%gamma)
  if (!link) {
    res <- cbind(beta=Xbeta,gamma=Zgamma)
    return(structure(res,grad=list(beta=X,gamma=Z)))
  }
  Pred <- ginv(Xbeta)
  p0 <- ginv(Zgamma)
  A1 <- dginv(Xbeta)
  A2 <- dginv(Zgamma)
  if (subdist) {
    dgamma <- apply(Z,2,function(z) A2*z)
    dbeta <- apply(X,2,function(x) A1*x)
    res <- cbind(subdist=Pred,pr=p0)
    return(structure(res,grad=list(subdist=dbeta,pr=dgamma)))
  }
  Pred <- p0*Pred
  A1 <- p0*A1
  A2 <- Pred*dginv(Zgamma)
  dgamma <- apply(Z,2,function(z) A2*z)
  dbeta <- apply(X,2,function(x) A1*x)
  attributes(Pred)$grad <- cbind(dbeta,p0*dgamma)
  return(Pred)
}

##' @export
residuals.zibreg <- function(object,newdata,...) {
  if (missing(newdata)) {
    y <- object$y
  } else {
    y <- model.frame(object$formula,newdata)[,1]
  }
  y-predict(object,newdata=newdata,...)
}

##' @export
summary.zibreg <- function(object,level=0.95,pr.contrast,...) {
  alpha <- 1-level
  alpha.str <- paste(c(alpha/2,1-alpha/2)*100,"",sep="%")
  cc <- cbind(coef(object),diag(vcov(object))^0.5)
  pval <- 2*(pnorm(abs(cc[,1]/cc[,2]),lower.tail=FALSE))
  qq <- qnorm(1-alpha/2)
  cc <- cbind(cc[,1],cc[,1]-qq*cc[,2],cc[,1]+qq*cc[,2],pval)
  colnames(cc) <- c("Estimate",alpha.str,"P-value")
  pr.names <- unlist(lapply(rownames(cc)[object$gamma.idx],
                            function(x) substr(x,4,nchar(x))))
  if (missing(pr.contrast)) {
    withIntercept <- pr.names[1]=="(Intercept)"
    pr.contrast <- diag(length(object$gamma.idx))
    if (withIntercept) pr.contrast[,1] <- 1
  }
  pr.cc <- cbind(pr.contrast%*%cc[object$gamma.idx,1],
                 diag((pr.contrast)%*%vcov(object)[object$gamma.idx,object$gamma.idx]%*%t(pr.contrast))^0.5)
  pr.cc <- object$family$linkinv(cbind(pr.cc[,1],pr.cc[,1]-qq*pr.cc[,2],pr.cc[,1]+qq*pr.cc[,2]))
  colnames(pr.cc) <- colnames(cc)[1:3]
  ## B <- cbind(0,cbind(0,pr.contrast))
  ## print(compare(object,contrast=B))
  pr.rnames <- c()
  for (i in seq_len(nrow(pr.contrast))) {
    Bidx <- which(pr.contrast[i,]!=0)
    Bval <- pr.contrast[i,Bidx]; Bval[Bval==1] <- ""
    pr.rnames <- c(pr.rnames,
                   paste0(paste0(Bval,paste0("{",pr.names[Bidx],"}"),collapse=" + ")))
  }
  rownames(pr.cc) <- pr.rnames

  return(structure(list(coef=cc, prevalence=pr.cc),class="summary.zibreg"))
}


##' @export
print.summary.zibreg <- function(x,...) {
  print(x$coef,...)
  cat("\nPrevalence probabilities:\n")
  print(x$prevalence,...)
}

##' @export
print.zibreg <- function(x,...) {
  print(summary(x,...))
}

##' @export
logLik.zibreg <- function(object,beta=object$beta,gamma=object$gamma,data,offset=object$offset,indiv=FALSE,...) {
  if (!missing(data)) {
    y <- model.frame(object$formula,data)[,1]
    X <- model.matrix(object$formula,data)
    Z <- model.matrix(object$formula.p,data)
    return(zibreg_logL(beta,gamma,y,X,Z,offset,object$family,indiv=indiv,...))
  }
  zibreg_logL(beta,gamma,object$y,object$X,object$Z,offset,object$family,indiv=indiv,...)
}
zibreg_logL <- function(beta,gamma,y,X,Z,offset=NULL,family=stats::binomial(),indiv=FALSE,...) {
  g <- family$linkfun
  ginv <- family$linkinv
  dginv <- family$mu.eta ## D[linkinv]
  n <- nrow(X)
  Xbeta <- as.vector(X%*%beta)
  Zgamma <- as.vector(Z%*%gamma)
  p0 <- ginv(Zgamma)
  if (!is.null(offset)) Xbeta <- Xbeta+offset
  Pr <- p0*ginv(Xbeta)
  loglik <- y*log(Pr)+(1-y)*log(1-Pr)
  if (indiv) return(loglik)
  loglik <- sum(loglik)
  structure(loglik,nobs=n,df=length(beta)+length(gamma),class="logLik")
}

##' @export
score.zibreg <- function(x,beta=x$beta,gamma=x$gamma,data,offset=x$offset,indiv=FALSE,...) {
  if (!missing(data)) {
    y <- model.frame(x$formula,data)[,1]
    X <- model.matrix(x$formula,data)
    Z <- model.matrix(x$formula.p,data)
    s <- zibreg_score(beta,gamma,y,X,Z,offset,x$family,indiv=indiv,...)
  } else {
    s <- zibreg_score(beta,gamma,x$y,x$X,x$Z,offset,x$family,indiv=indiv,...)
  }
  if (indiv) colnames(s) <- names(x$coef) else names(s) <- names(x$coef)
  return(s)
}

zibreg_score <- function(beta,gamma,y,X,Z,offset=NULL,family=stats::binomial(),indiv=FALSE,...) {
  g <- family$linkfun
  ginv <- family$linkinv
  dginv <- family$mu.eta ## D[linkinv]
  n <- nrow(X)
  Xbeta <- as.vector(X%*%beta)
  Zgamma <- as.vector(Z%*%gamma)
  p0 <- ginv(Zgamma)
  if (!is.null(offset)) Xbeta <- Xbeta+offset
  Pr <- p0*ginv(Xbeta)
  A0 <- (y/Pr  - (1-y)/(1-Pr))
  A1 <- A0*p0*dginv(Xbeta)
  A2 <- A0*ginv(Xbeta)*dginv(Zgamma)
  dbeta <- apply(X,2,function(x) A1*x)
  dgamma <- apply(Z,2,function(z) A2*z)
  ss <- cbind(dbeta,dgamma)
  if (indiv) return(ss)
  colSums(ss)
}

##' @export
information.zibreg <- function(x,beta=x$beta,gamma=x$gamma,data,offset=x$offset,type=c("robust","outer","obs"),...) {
  if (!missing(data)) {
    y <- model.frame(x$formula,data)[,1]
    X <- model.matrix(x$formula,data)
    Z <- model.matrix(x$formula.p,data)
    I <- zibreg_information(beta,gamma,y,X,Z,offset,x$family,type=type,...)
  } else {
    I <- zibreg_information(beta,gamma,x$y,x$X,x$Z,offset,x$family,type=type,...)
  }
  colnames(I) <- rownames(I) <- names(x$coef)
  return(I)
}

zibreg_information <- function(beta,gamma,y,X,Z,offset=NULL,family=stats::binomial(),type=c("outer","obs","robust"),...) {
  if (tolower(type[1])%in%c("obs","hessian")) {
    beta.idx <- seq(ncol(X)); gamma.idx <- seq(ncol(Z))+ncol(X)
    I <- -numDeriv::jacobian(function(x)
                   zibreg_score(x[beta.idx],x[gamma.idx],y,X,Z,offset,family,...),c(beta,gamma))
    return(I)
  }
  if (tolower(type[1])%in%c("robust","sandwich")) {
    I <- zibreg_information(beta,gamma,y,X,Z,offset,family,type="obs")
    J <- zibreg_information(beta,gamma,y,X,Z,offset,family,type="outer")
    return(J%*%Inverse(I)%*%J)
  }
  S <- zibreg_score(beta,gamma,y,X,Z,offset,family,indiv=TRUE,...)
  crossprod(S)
}
