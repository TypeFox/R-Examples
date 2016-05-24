testModels <- function(model, null.model, method=c("D1","D2","D3"), use=c("wald","likelihood"), df.com=NULL){
# model comparison and hypothesis tests for k-dimensional estimands

  # *** 
  # general errors
  if(!"list"%in%class(model) & !"list"%in%class(null.model) & !is.null(null.model)) stop("The 'model' and 'null.model' arguments must be lists of fitted statistical models.")
  if(!"list"%in%class(model) & is.null(null.model)) stop("The 'model' argument must be a list of fitted statistical models.")
  if(!missing(method) & length(method)>1) stop("Only one 'method' may be supplied.")
  if(!missing(use) & length(use)>1) stop("Only one of 'wald' or 'likelihood' may be supplied.")
  method <- match.arg(method)
  use <- match.arg(use)

  # ***
  # warnings
  if(!is.null(df.com) & method!="D1") warning("Complete-data degrees of freedom are not available for use with '",method,"', and thus were ignored.")
  if(use=="likelihood" & method!="D2") warning("The 'likelihood' option is not available with method '",method,"', and thus was ignored.")

  # ***
  # select extraction methods
  cls <- class(model[[1]])

  # default method (lm)
  coef.method <- vc.method <- lr.method <- "default"
  if(cls[1]=="lm") vc.method <- lr.method <- "lm"

  # merMod (lme4)
  if(any(grepl("merMod",cls)) & coef.method=="default"){
    if(!requireNamespace("lme4", quietly=TRUE)) stop("The 'lme4' package must be installed to handle 'merMod' class objects.")
    coef.method <- vc.method <- lr.method <- "lmer"
  }
  
  # lme (nlme)
  if(any(grepl("^.?lme$",cls)) & coef.method=="default"){
    if(!requireNamespace("nlme", quietly=TRUE)) stop("The 'nlme' package must be installed to handle 'lme' class objects.")
    coef.method <- vc.method <- lr.method <- "nlme"
  }
  
  # ***
  #!
  if(method=="D1"){

    fe <- switch(coef.method,
      lmer=.getCOEF.lmer(model,null.model),
      nlme=.getCOEF.nlme(model,null.model),
      default=.getCOEF.default(model,null.model)
    )

    m <- length(model)
    Qhat <- fe$Qhat
    Uhat <- fe$Uhat
    if(is.null(dim(Qhat))) dim(Qhat) <- c(1,m)
    if(is.null(dim(Uhat))) dim(Uhat) <- c(1,1,m)
    k <- dim(Qhat)[1]
    
    Qbar <- apply(Qhat,1,mean)
    Ubar <- apply(Uhat,c(1,2),mean)
    B <- cov(t(Qhat))
    r <- (1+m^(-1))*sum(diag(B%*%solve(Ubar)))/k
    Ttilde <- (1 + r)*Ubar
    
    # D1 (Li, Raghunathan and Rubin, 1991)
    val <- t(Qbar) %*% solve(Ttilde) %*% Qbar / k
    t <- k*(m-1)
  
    if(!is.null(df.com)){
      a <- r*t/(t-2)
      vstar <- ( (df.com+1) / (df.com+3) ) * df.com
      v <- 4 + ( (vstar-4*(1+a))^(-1) + (t-4)^(-1) * ((a^2*(vstar-2*(1+a))) / 
           ((1+a)^2*(vstar-4*(1+a)))) )^(-1)
    } else {
      if (t>4){ 
        v <- 4 + (t-4) * (1 + (1 - 2*t^(-1)) * (r^(-1)))^2
      }else{
        v <- t * (1 + k^(-1)) * ((1 + r^(-1))^2) / 2
      }
    }
    p <- 1-pf(val, k, v)

    out <- matrix(c(val,k,v,p,r),ncol=5)
    colnames(out) <- c("F.value","df1","df2","p.value","RIV")
    
    out <- list(
      call=match.call(),
      test=out,
      m=m,
      adj.df=!is.null(df.com),
      df.com=df.com,
      method="D1",
      use="wald",
      reml=FALSE
    )
  }

  # ***
  #!
  if(method=="D2"){

    if(use=="wald"){

      reml=FALSE

      fe <- switch(coef.method,
        lmer=.getCOEF.lmer(model,null.model),
        nlme=.getCOEF.nlme(model,null.model),
        default=.getCOEF.default(model,null.model)
      )

      m <- length(model)
      Qhat <- fe$Qhat
      Uhat <- fe$Uhat
      if(is.null(dim(Qhat))) dim(Qhat) <- c(1,m)
      if(is.null(dim(Uhat))) dim(Uhat) <- c(1,1,m)
      k <- dim(Qhat)[1]
      dW <- sapply(1:m, function(z) t(Qhat[,z]) %*% solve(Uhat[,,z]) %*% Qhat[,z])

    }

    # TODO: likelihood test for (single) model fit (with null.model=NULL)
    if(use=="likelihood"){

      # check for REML and refit
      reml1 <- sapply(model, .is.REML, lr.method=lr.method)
      reml0 <- sapply(null.model, .is.REML, lr.method=lr.method)
      reml <- ( any(reml0) | any(reml1) )
      if(reml){
        model[reml1] <- lapply(model[reml1], .update.ML, lr.method=lr.method)
        null.model[reml0] <- lapply(null.model[reml0], .update.ML, lr.method=lr.method)
      }

      dW <- switch(coef.method,
        lmer=.getLR.lmer(model,null.model),
        nlme=.getLR.nlme(model,null.model),
        default=.getLR.default(model,null.model)
      )

      m <- length(model)
      k <- attr(dW,"df")
      if(is.null(k)) stop("Degrees of freedom for the model comparison could not be detected.")
  
    }

    # D2 (Li, Meng et al., 1991)
    dWbar <- mean(dW)
    r <- (1+m^(-1)) * var(sqrt(dW))
    val <- (dWbar/k - (m+1)/(m-1) * r) / (1+r)

    v <- k^(-3/m) * (m-1) * (1+r^(-1))^2
    p <- 1-pf(val, k, v)
   
    out <- matrix(c(val,k,v,p,r),ncol=5)
    colnames(out) <- c("F.value","df1","df2","p.value","RIV")
    
    out <- list(
      call=match.call(),
      test=out,
      k=k,
      m=m,
      adj.df=FALSE,
      df.com=NULL,
      method="D2",
      use=use,
      reml=reml
    )
  }

  # ***
  #!
  if(method=="D3"){

    # error checking
    if(!lr.method%in%c("lm","lmer","nlme")) stop("The 'D3' method is currently not supported for models of class '",cls,"'.")
    if(!grepl("^lme$",cls[1]) & lr.method=="nlme") stop("The 'D3' method is currently only supported for linear mixed-effects models.")
    if(!grepl("^l?merMod$",cls[1]) & lr.method=="lmer") stop("The 'D3' method is currently only supported for linear mixed-effects models.")

    # check for REML and refit
    reml1 <- sapply(model, .is.REML, lr.method=lr.method)
    reml0 <- sapply(null.model, .is.REML, lr.method=lr.method)
    reml <- ( any(reml0) | any(reml1) )
    if(reml){
      model[reml1] <- lapply(model[reml1], .update.ML, lr.method=lr.method)
      null.model[reml0] <- lapply(null.model[reml0], .update.ML, lr.method=lr.method)
    }

    # LR at fit-specific estimates
    dL <- switch( lr.method, lmer=.getLR.lmer(model,null.model), nlme=.getLR.nlme(model,null.model),
                  lm=.getLR.default(model,null.model))

    fe0 <- switch( coef.method, lmer=.getCOEF.lmer(null.model), nlme=.getCOEF.nlme(null.model),
                   default=.getCOEF.default(null.model) )
    vc0 <- switch( vc.method, lmer=.getVC.lmer(null.model), nlme=.getVC.nlme(null.model),
                   lm=.getVC.lm(null.model,ML=TRUE) )
    fe1 <- switch( coef.method, lmer=.getCOEF.lmer(model), nlme=.getCOEF.nlme(model),
                  default=.getCOEF.default(model) )
    vc1 <- switch( vc.method, lmer=.getVC.lmer(model), nlme=.getVC.nlme(model),
                   lm=.getVC.lm(model,ML=TRUE) )

    dLbar <- mean(dL)
    m <- length(model)
    k <- attr(dL,"df")

    # LR at average estimates
    switch( lr.method, 

      lmer={

        if(length(vc0$vlist)>2) stop("The 'D3' method is only supported for models of class 'merMod' with a single cluster variable.")

        Q0 <- fe0$Qhat
        Q1 <- fe1$Qhat
        if(is.null(dim(Q0))) dim(Q0) <- c(1,m)
        if(is.null(dim(Q1))) dim(Q1) <- c(1,m)
        V0 <- lapply(vc0$vlist, function(z) unname(apply(z,1:2,mean)) )
        V1 <- lapply(vc1$vlist, function(z) unname(apply(z,1:2,mean)) )
        psi0bar <- list(beta=rowMeans(Q0),D=V0[[1]],sigma2=V0[[2]][1,1])
        psi1bar <- list(beta=rowMeans(Q1),D=V1[[1]],sigma2=V1[[2]][1,1])
        dLt <- .getLR.lmer(model,null.model,psi=psi1bar,null.psi=psi0bar)

      },

      nlme={

        if(length(vc0$vlist)>2) stop("The 'D3' method is only supported for models of class 'lme' with a single cluster variable.")

        Q0 <- fe0$Qhat
        Q1 <- fe1$Qhat
        if(is.null(dim(Q0))) dim(Q0) <- c(1,m)
        if(is.null(dim(Q1))) dim(Q1) <- c(1,m)
        V0 <- lapply(vc0$vlist, function(z) unname(apply(z,1:2,mean)) )
        V1 <- lapply(vc1$vlist, function(z) unname(apply(z,1:2,mean)) )
        psi0bar <- list(beta=rowMeans(Q0),D=V0[[1]],sigma2=V0[[2]][1,1])
        psi1bar <- list(beta=rowMeans(Q1),D=V1[[1]],sigma2=V1[[2]][1,1])
        dLt <- .getLR.nlme(model,null.model,psi=psi1bar,null.psi=psi0bar)

      },

      lm={

        Q0 <- fe0$Qhat
        Q1 <- fe1$Qhat
        if(is.null(dim(Q0))) dim(Q0) <- c(1,m)
        if(is.null(dim(Q1))) dim(Q1) <- c(1,m)
        V0 <- lapply(vc0$vlist, function(z) unname(apply(z,1:2,mean)) )
        V1 <- lapply(vc1$vlist, function(z) unname(apply(z,1:2,mean)) )
        psi0bar <- list(beta=rowMeans(Q0),sigma2=V0[[1]][1,1])
        psi1bar <- list(beta=rowMeans(Q1),sigma2=V1[[1]][1,1])
        dLt <- .getLR.default(model,null.model,psi=psi1bar,null.psi=psi0bar)

      }

    )
  
    # D3 (Meng & Rubin, 1992)
    dLtilde <- mean(dLt)
    r <- (m+1) * (k*(m-1))^(-1) * (dLbar-dLtilde)
    val <- dLtilde / (k*(1+r))
  
    t <- k*(m-1)
    if( t>4 ){
      v <- 4 + (t-4) * (1 + (1-2*t^(-1)) * r^(-1))^2
    }else{
      v <- t * (1+k^(-1)) * (1+r^(-1))^2 / 2
    }
  
    p <- 1- pf(val, k, v)

    out <- matrix(c(val,k,v,p,r),ncol=5)
    colnames(out) <- c("F.value","df1","df2","p.value","RIV")

    out <- list(
      call=match.call(),
      test=out,
      m=m,
      adj.df=FALSE,
      df.com=NULL,
      method="D3",
      use="likelihood",
      reml=reml
    )
  }

  class(out) <- "mitml.testModels"
  out
}

.is.REML <- function(x, lr.method){

  reml <- FALSE
  if(lr.method=="lmer") reml <- lme4::isREML(x)
  if(lr.method=="nlmer") reml <- x$method=="REML"
  reml

}

.update.ML <- function(x, lr.method){

  if(lr.method=="lmer") x <- update(x, REML=FALSE)
  if(lr.method=="nlme") x <- update(x, data=x$data, method="ML")
  x

}
