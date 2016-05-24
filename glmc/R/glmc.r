glmc <- function (formula,
                  family = gaussian, data,
                  na.action, start = NULL, etastart, mustart, offset,
                  control = glmc.control(...),
                  model = TRUE, glm.method= "glm.fit",optim.method="Nelder-Mead",
                  emplik.method="Owen",
                  optim.hessian=FALSE,
                  x = FALSE, y = TRUE, 
                  Amat=NULL, confn=NULL,...) 

{
#
#
    call <- match.call()
    structure(1,CALL=match.call())
    if (is.character(family)) 
        family <- get(family, mode = "function", envir = parent.frame())
    if (is.function(family)) 
        family <- family()
    if (is.null(family$family)) {
        print(family)
        stop("`family' not recognized")
    }
    
    if (missing(data)){ 
        data <- environment(formula)
    }else{
      if(is.null(Amat) & is.null(confn) & !is.na(charmatch("constraints",names(data)))){
        Amat <- as.matrix(data)[,c(grep("constraints",names(data)))]
        if(is.vector(Amat)){Amat <- matrix(Amat,ncol=1)}
      }
      if(!is.null(Amat) & is.character(Amat)){
        Amat <- data[Amat[1]]
      }
    }
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "na.action", 
        "etastart", "mustart", "offset"), names(mf), 0)
    mf <- mf[c(1, m)]

    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    switch(glm.method, model.frame = return(mf), glm.fit = 1, glm.fit.null = 1, 
        stop("invalid `glm.method': ", glm.method))
    mt <- attr(mf, "terms")
    Y <- model.response(mf, "numeric")
    X <- if (!is.empty.model(mt)) 
        model.matrix(mt, mf, contrasts)
    else matrix(, NROW(Y), 0)
    

    weights <- NULL

    offset <- model.offset(mf)
    if (!is.null(weights) && any(weights < 0)) 
        stop("Negative wts not allowed")
    if (!is.null(offset) && length(offset) != NROW(Y)) 
        stop("Number of offsets is ", length(offset), ", should equal ", 
            NROW(Y), " (number of observations)")
    mustart <- model.extract(mf, "mustart")
    etastart <- model.extract(mf, "etastart")
    
     if(length(control$optim.control$parscale)==1){
     control$optim.control$parscale <- NULL  }
     if(length(control$optim.control$ndeps)==1){
     control$optim.control$ndeps <- NULL  }
     
    if(is.null(weights))
       weights=rep(1,nrow(as.matrix(Y)))

    if(is.null(Amat) & is.null(confn)){
    
      nra <- 0
      fit <- glm.fit(x = X, y = Y, weights = weights, start = start, 
        etastart = etastart, mustart = mustart, offset = offset, 
        family = family, control = control$glm.control,
        intercept = attr(mt, "intercept") > 0)
      
    if (any(offset) && attr(mt, "intercept") > 0) {
        fit$null.deviance <- glm.fit(x = X[, "(Intercept)", drop = FALSE], 
            y = Y, offset = offset, family = family, weights=weights,
            control = control$glm.control, intercept = TRUE)$deviance
      }

      if (model) 
        fit$model <- mf
        fit$na.action <- attr(mf, "na.action")


      if (x) 
        fit$x <- X
      if (!y) 
        fit$y <- NULL
        fit$hessian<-NULL
      fit <- c(fit, list(call = call, formula = formula, terms = mt, 
        data = data, offset = offset, control = control$glm.control, glm.method = glm.method, 
        contrasts = attr(X, "contrasts"), xlevels = glmc.getXlevels(mt,mf)))
      class(fit) <- c("glm", "lm")
    
    }else{
      
     if(!is.null(Amat)){
      
     initmu<-rep(0,ncol(Amat))
     if(emplik.method=="Owen"){
     initwt<-el.owen(Amat,initmu,control=control$weights.control)}else{
     initwt<-el.CSW(Amat,initmu,control=control$weights.control)}
      
      lam<-initwt$lambda     

      wf<-(initwt$wts/nrow(as.matrix(Y))) 
     
 
      final.weights<-wf*weights


     fit<-glm.fit(x=X,y=Y,weights=final.weights,start=start,etastart=etastart,mustart=mustart,offset=offset,family=family,control=control$glm.control,intercept=attr(mt,"intercept")>0)

      t<- fit$coefficient


      if(fit$converged!=TRUE||fit$boundary!=FALSE){
      
      v2 <- optim(t, glmc.inner2,
                method=optim.method,
                control=control$optim.control,
                hessian=optim.hessian,
                rt=1, Y=Y, X=X, family=family, Amat=Amat, confn=NULL, weights=weights,controlwt=control$weights.control,emplikmethod=emplik.method,...)

      t <- v2$par
      rt <- 0
      wf <- glmc.inner2(t,rt,Y,X,family,Amat,confn=NULL,weights,controlwt=control$weights.control,emplikmethod=emplik.method)
     
      final.weights <- (wf$el$wts/nrow(as.matrix(Y)))*weights
      
      fit <- glm.fit(x = X, y = Y, weights = final.weights, start = start, 
                     etastart = etastart, mustart = mustart, 
                     offset = offset, 
                     family = family, 
                     control = control$glm.control, intercept = attr(mt, "intercept") > 0)
    }

}else{

      fit1 <- glm.fit(x = X, y = Y, weights = weights, start = start, 
            etastart = etastart, mustart = mustart, offset = offset, 
            family = family, control = control$glm.control, intercept = attr(mt, 
                "intercept") > 0)
      tstart <- fit1$coefficients
      

v2 <- optim(tstart, glmc.inner2,
                method=optim.method,
                control=control$optim.control,
                hessian=optim.hessian,
                rt=1, Y=Y, X=X, family=family, Amat=NULL, confn=confn, weights=weights,controlwt=control$weights.control,emplikmethod=emplik.method,...)

      t <- v2$par
      rt <- 0
      wf <- glmc.inner2(t,rt,Y,X,family,Amat=NULL,confn=confn,weights,controlwt=control$weights.control,emplikmethod=emplik.method)

      final.weights <- (wf$el$wts/nrow(as.matrix(Y)))*weights
      Amat=wf$Amat
      fit <- glm.fit(x = X, y = Y, weights = final.weights, start = start, 
                     etastart = etastart, mustart = mustart, 
                     offset = offset, 
                     family = family, 
                     control = control$glm.control, intercept = attr(mt, "intercept") > 0)
          }
      
      
      if (any(offset) && attr(mt, "intercept") > 0) {
        fit$null.deviance <- glm.fit(x = X[, "(Intercept)", drop = FALSE], 
            y = Y, weights = NULL, offset = offset, family = family, 
            control = control$glm.control, intercept = TRUE)$deviance
      }
 
      if (model) 
          fit$model <- mf
      fit$na.action <- attr(mf, "na.action")
      if (x) 
        fit$x <- X
      if (!y) 
        fit$y <- NULL
      
       eta<-X%*%t
       mu<-family$linkinv(eta)

       scale<-final.weights/weights

      if(optim.hessian==TRUE){
      fit$hessian=(-1)*solve(v2$hessian)}else{
      #Z<-t(t(X)*diag(as.array(family$mu.eta(eta)/family$variance(mu))))
      Z<-X*drop(family$mu.eta(eta)/family$variance(mu))
      G<-t(Z*c(family$variance(mu)*final.weights))%*%Z
      G2<-t(Z*c(((Y-mu)*final.weights)^2))%*%Z
      Tv<-t(Z*c((Y-mu)*final.weights))%*%(Amat*c(scale))
      H<-t(Amat*c(scale))%*%(Amat*c(scale))
      hessian<-solve(G)%*%(G2-Tv%*%solve(H)%*%t(Tv))%*%solve(G)

      fit$hessian<-hessian
      }
 
      fit$final.weights<-final.weights
      fit$prior.weights<-weights

      fit <- c(fit, list(call = call, formula = formula, terms = mt, 
        data = data, offset = offset, control = control, glm.method = glm.method,emplik.method=emplik.method,optim.method=optim.method, 
        contrasts = attr(X, "contrasts"), xlevels = glmc.getXlevels(mt, 
            mf)))
      class(fit) <- c("glmc","glm", "lm")
    }
    return(fit)
}

glmc.inner2 <- function(t,rt,Y,X,family,Amat,confn,weights,controlwt,emplikmethod,...){
    ncx <- ncol(as.matrix(X))
    nr <- nrow(as.matrix(Y))
    if(!is.null(Amat)){
        nra <- nrow(Amat)
        if(nra!=nr)
         {stop("Number of rows in the constraint matrix ", nra, ", should equal ", NROW(Y), " (number of observations)")}
     nc <- ncol(Amat)
    }else{
    
    confn1 <- function(t) confn(t,Y,X,...)
#    print(confn1)
    Amat <- confn1(t)
    nra <- nrow(Amat)
    if(nra!=nr)
       {stop("Number of rows in the constraint matrix ", nra, ", should equal ", 
            NROW(Y), " (number of observations)")}
     nc <- ncol(Amat)
    }

    if(is.null(weights)){
     weights <- rep(1,nr)
    }else{
     weights <- weights
    }

    eta <- X%*%t
    mu<-family$linkinv(eta)

    link <- family$mu.eta(mu)*(Y-mu)*weights/family$variance(mu)

    mat <- matrix(0,nrow=nr,ncol=(nc+ncx))

    mat[,(1:ncx)] <- X*as.numeric(link)
   
    if(!is.matrix(Amat)){
     mat[,((ncx+1):(nc+ncx))] <- as.matrix(Amat)
    }else{
    mat[,((ncx+1):(nc+ncx))] <- Amat
    }
    
    
    mean <- rep(0,ncol(mat))
    if(emplikmethod=="Owen"){
    el <- el.owen(mat,mean,control=controlwt)}else{
    el <- el.CSW(mat,mean,control=controlwt)}

    constrerr<-abs(sum(el$wts)/nr-1)
    wtsum=sum(el$wts/nr)
     

    if(rt==1)
     {if(constrerr>=10^(-3)){
       return((-1)*nr*log(sum(nr/el$wts)))}else{
      return(el$"-2LLR"/(-2))}
    }
    else
     {list(el=el,Amat=Amat)}
}

el.owen <- function (x, mu, lam, control) 
{
    x <- as.matrix(x)
    n <- nrow(x)
    p <- ncol(x)
    mu <- as.vector(mu)


    maxit<-control$maxit
    gradtol<-control$gradtol
    svdtol<-control$svdtol
    itertrace<-control$itertrace

    if (length(mu) != p) 
        stop("mu must have same dimension as observation vectors.")
    if (n <= p) 
        stop("Need more observations than length(mu) in el.test().")
    z <- t(t(x) - mu)
    TINY <- sqrt(.Machine$double.xmin)
    scale <- mean(abs(z)) + TINY
    z <- z/scale
    if (!missing(lam)) {
        lam <- as.vector(lam)
        lam <- lam * scale
        if (logelr(z, rep(0, p), lam) > 0) 
            lam <- rep(0, p)
    }
    if (missing(lam)) 
        lam <- rep(0, p)
    if (svdtol < TINY) 
        svdtol <- TINY
    if (gradtol < TINY) 
        gradtol <- TINY
    nwts <- c(3^-c(0:3), rep(0, 12))
    gwts <- 2^(-c(0:(length(nwts) - 1)))
    gwts <- (gwts^2 - nwts^2)^0.5
    gwts[12:16] <- gwts[12:16] * 10^-c(1:5)
    nits <- 0
    gsize <- gradtol + 1
    while (nits < maxit && gsize > gradtol) {
        arg <- 1 + z %*% lam
        wts1 <- as.vector(llogp(arg, 1/(n)))
        wts2 <- as.vector(-llogpp(arg, 1/(n)))^0.5
        grad <- as.matrix(-z * wts1)
        grad <- as.vector(rowsum(grad, rep(1, nrow(grad))))
        gsize <- mean(abs(grad))
        hess <- z * wts2
        svdh <- svd(hess)
        if (min(svdh$d) < max(svdh$d) * svdtol) 
            svdh$d <- svdh$d + max(svdh$d) * svdtol
        nstep <- svdh$v %*% (t(svdh$u)/svdh$d)
        nstep <- as.vector(nstep %*% matrix(wts1/wts2, n, 1))
        gstep <- -grad
        if (sum(nstep^2) < sum(gstep^2)) 
            gstep <- gstep * (sum(nstep^2)^0.5/sum(gstep^2)^0.5)
        ologelr <- -sum(llog(arg, 1/(n)))
        ninner <- 0
        for (i in 1:length(nwts)) {
            nlogelr <- logelr(z, rep(0, p), lam + nwts[i] * nstep + 
                gwts[i] * gstep)
            if (nlogelr < ologelr) {
                lam <- lam + nwts[i] * nstep + gwts[i] * gstep
                ninner <- i
                break
            }
        }
        nits <- nits + 1
        if (ninner == 0) 
            nits <- maxit
        if (itertrace) 
            print(c(lam, nlogelr, gsize, ninner))
    }
    list("-2LLR" = -2 * nlogelr, Pval = 1 - pchisq(-2 * nlogelr, 
        df = p), lambda = lam/scale, grad = grad * scale, hess = t(hess) %*% 
        hess * scale^2, wts = wts1, nits = nits)
}

el.CSW<-function(x, mu, lam,control){
    x <- as.matrix(x)
    n <- nrow(x)
    p <- ncol(x)
    mu <- as.vector(mu)
    if (length(mu) != p) 
        stop("mu must have same dimension as observation vectors.")
    if (n <= p) 
        stop("Need more observations than length(mu) in el.test().")

    maxit<-control$maxit
    gradtol<-control$gradtol
    svdtol<-control$svdtol
    itertrace<-control$itertrace
    
    TINY <- sqrt(.Machine$double.xmin)

    if (gradtol < TINY) 
        gradtol <- TINY

    if (!missing(lam)) {
        lam <- as.vector(lam)}else
    {lam<-as.vector(rep(0,p))}
        
    k<-0
    nits=0
    gamma=1

    norm=1
    lold=0
    while(norm>=gradtol && nits<=maxit){
     divisor=1/(1+x%*%lam)
     grad=t(x)%*%divisor
     hess=(-1)*solve(t(x*c(divisor^2))%*%x)%*%grad

     delta=gamma*hess
     divisorn=1+x%*%(lam-delta)
     ninner=0
     while(any(divisorn<=0)==TRUE||(sum(log(divisorn))<(-1)*sum(log(divisor)))){
      gamma=gamma/2
      delta=gamma*hess
      divisorn=1+x%*%(lam-delta)
      ninner=ninner+1

      }
     lam=lam-delta
     nits=nits+1
     gamma=sqrt(1/nits)
     norm=sqrt(t(hess)%*%hess)
     
     nlogelr=logelr(x,mu,lam)
     
     wts1=1/divisorn 
     if (itertrace)
           print(c(lam, nlogelr, norm, ninner))
     }

list("-2LLR" = -2 * nlogelr, Pval = 1 - pchisq(-2 * nlogelr, 
        df = p), lambda = lam, grad = grad, hess = hess %*% 
        t(hess), wts = wts1, nits = nits)
}


