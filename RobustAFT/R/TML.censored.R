TML.censored<-function(formula, delta, data, errors = "Gaussian", initial = "S",  input = NULL, 
              otp = "fixed",  cov=TRUE, cu = NULL, control.S=list(), control.ref=list(), control.tml=list())
{
  if(!(errors %in% c("Gaussian", "logWeibull")))
		stop(gettextf("Errors distribution should be Gaussian or log-Weibull"))
  if(!(initial %in% c("S", "input")))
                stop(gettextf("initial should be S or input"))
  if(!(otp %in% c("fixed", "adaptive")))
		stop(gettextf("'otp' should be fixed or adaptive"))
  	
  call <- match.call()
	if(missing(data))
		data <- environment(formula)
	mf <- match.call(expand.dots = FALSE)
	m <- match(c("formula", "data"), names(mf), 0)
	mf <- mf[c(1, m)]
	mf$drop.unused.levels <- TRUE
	mf[[1]] <- as.name("model.frame")
	mf <- eval(mf, parent.frame())
	mt <- attr(mf, "terms")
	y <- model.response(mf, "any")
	if(length(dim(y)) == 1){
		nm <- rownames(y)
		if(!is.null(nm))
			names(y) <- nm
	}
	X <- if(!is.empty.model(mt))
		model.matrix(mt, mf, contrasts)
	else matrix(, NROW(y), 0)
	n <- length(y)
	
        control.S   <- do.call("TML.censored.control.S", control.S)
	control.ref <- do.call("TML.censored.control.ref", control.ref)
	control.tml <- do.call("TML.censored.control.tml", control.tml)
        N <- control.S$N; q <- control.S$q; sigma0 <- control.S$sigma0
        MAXIT <- control.S$MAXIT; TOL <- control.S$TOL; seed <- control.S$seed
	
	# Gaussian errors
	if(errors == "Gaussian"){

    if (initial == "S") {   
      # Parametric S
      zS  <- SparamG.S(X,y,delta,N,q,sigma0,MAXIT,TOL,ialg=3,seed)
      b.SP <- zS$Bmin
      s.SP <- zS$Smin
	  zR <- RefSG(X, y, delta, b.SP, s.SP, control.ref)
      nit0 <- 0 }

    if (initial == "input") {
      z <- input
      b.SP <- z$theta
      s.SP <- z$sigma
	  zR <- list(Bmin=b.SP, Smin=s.SP)
      nit0 <- 0 }

    B.SP <- zR$Bmin
    S.SP <- zR$Smin 
    nit.ref <- zR$Nit

    # Cut-off values
    
    rs0  <- as.numeric((y - X%*%B.SP)/S.SP)
    ct <- CutoffN(rs0,delta,cu=2.5,zmax=10,gridsize=2000)
    if(otp == "fixed"){
      if(is.null(cu)) cu <- 2.5
      tu <- cu
      tl <- -tu
    }
    else{ 
      tu <- ct$tu
      tl <- -tu
    }
    alpha <- ct$alpha
    
    
    # TML

    # constant for rectangular weight function
    const   <- -tu*dnorm(tu)+tl*dnorm(tl)+pnorm(tu)-pnorm(tl)

    zC    <- TMLG(X,y,delta,B.SP,S.SP,const,tu,control.tml)
    B.WML <- zC$Beta  
    S.WML <- zC$sigma 
    nit.tml <- zC$Nit
    # Covariance matrix of TML

    Beta.tilde  <- B.SP
    sigma.tilde <- S.SP
    Beta.hat    <- B.WML
    sigma.hat   <- S.WML

    rs0         <- as.vector(y-X%*%Beta.tilde)/sigma.tilde
    d.Beta      <- (Beta.hat-Beta.tilde)/sigma.tilde
    d.sigma     <- sigma.hat/sigma.tilde
    Cov         <- NULL
    if (cov) Cov <- Cov.Gauss(d.Beta,d.sigma,sigma.tilde,rs0,delta,X,tu,const)
    
    # Identifying the outliers
    
    Ers0 <- rs0
    for (i in 1:n) {
      if (delta[i]==1) Ers0[i] <- rs0[i]
      if (delta[i]==0) Ers0[i] <- integrate(uf0,lower=rs0[i],upper=10)$val/(1-pnorm(rs0[i])) 
    }
    
    # Create the expected residuals
    Res0 <- (y - X %*% Beta.hat)/sigma.hat
    ERes0 <- Res0
    for (i in 1:n) {
      if (delta[i]==1) ERes0[i] <- Res0[i]
      if (delta[i]==0) ERes0[i] <- integrate(uf0,lower=Res0[i],upper=10)$val/(1-pnorm(Res0[i])) 
    }
  }
  
  # Log-Weibull errors
  if(errors == "logWeibull"){
  
    if (initial == "S") {   
      # Parametric S
      zS  <- SparamW.S(X,y,delta,N,q,sigma0,MAXIT,TOL,ialg=3,seed)
      b.SP <- zS$Bmin
      s.SP <- zS$Smin 
	  zR <- RefSW(X,y,delta,b.SP,s.SP,control.ref)
      nit0 <- 0}

    if (initial == "input") {
      z <- input
      b.SP <- z$tau
      s.SP <- z$v
	  zR <- list(Bmin=b.SP, Smin=s.SP)
      nit0 <- 0 }

    B.SP <- zR$Bmin
    S.SP <- zR$Smin

    # Cut-off values
    rs0  <- as.numeric((y - X%*%B.SP)/S.SP)
    ct <- CutoffW(rs0,delta,cu=1.8554,zmax=10,gridsize=2000)
    if(otp == "fixed"){
      if(is.null(cu)) cu <- 1.855356
      tu <- cu
      tl <- Izero(tu)
    }
    else{
      tl <- ct$tl
      tu <- ct$tu
    }
    alpha <- ct$alpha
    
    
    # TML

    # constant for rectangular weight function
    const <- tl*dlweibul(tl)-tu*dlweibul(tu)+plweibul(tu)-plweibul(tl)
    
    zC    <- TMLW(X,y,delta,B.SP,S.SP,const,tl,tu,control.tml)
    B.WML <- zC$Beta  
    S.WML <- zC$sigma 
    
    # Covariance matrix of TML

    Beta.tilde  <- B.SP
    sigma.tilde <- S.SP
    Beta.hat    <- B.WML
    sigma.hat   <- S.WML

    rs0         <- as.vector(y-X%*%Beta.tilde)/sigma.tilde
    d.Beta      <- (Beta.hat-Beta.tilde)/sigma.tilde
    d.sigma     <- sigma.hat/sigma.tilde
    Cov         <- NULL
    if (cov) Cov <- Cov.LogW(d.Beta,d.sigma,sigma.tilde,rs0,delta,X,tl,tu,const)
    
    # Identifying the outliers
    
    Ers0 <- rs0

    for (i in 1:n) {
      if (delta[i]==1) Ers0[i] <- rs0[i]
      if (delta[i]==0) Ers0[i] <- integrate(uf0w,lower=rs0[i],upper=10)$val/(1-plweibul(rs0[i])) 
    }
    
    # Create the expected residuals
    Res0 <- (y - X %*% Beta.hat)/sigma.hat
    ERes0 <- Res0
    for (i in 1:n) {
      if (delta[i]==1) ERes0[i] <- Res0[i]
      if (delta[i]==0) ERes0[i] <- integrate(uf0w,lower=Res0[i],upper=10)$val/(1-plweibul(Res0[i])) 
    } 
  }
  fitted.values <- X %*% Beta.hat
  wi <- ww(Ers0,tl,tu)
  tn <- length(ERes0) - length((1:n)[wi==0])
  
  res <- list(th0=B.SP, v0=S.SP, nit0=nit0, nit.ref=zR$Nit, th1=B.WML, v1=S.WML, nit.tml=zC$Nit, tl=tl, tu = tu, 
          alpha=alpha, tn=tn, weights=wi, COV=Cov, residuals=ERes0*S.WML, 
          fitted.values=fitted.values, call=call, formula=formula, terms=mt, data=data, 
          errors=errors)
  class(res) <- "TML"
  res
}
