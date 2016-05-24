ssmn.est <- function(y, X, family="sn", method="EM", error =  1e-6, maxit=1000,show.envelope=FALSE)
{

 if(method=="MLE")
 {
    #############################################################################
    ### Maximum Likelihood Estimators of Newton-Raphson procedures
   if(family=="sn")
   {
     start.time  <- Sys.time() #Begin Time
     n          <- length(y)
     if(missing(X)) {X<-matrix(1,n,1)}
     p          <- ncol(X)
     lambda     <- as.numeric(skewness(y))
     beta       <- solve(t(X)%*%X)%*%t(X)%*%y
     sigma2     <- as.numeric(t((y-X%*%beta))%*%(y-X%*%beta)/(n-p))
     t0         <- as.matrix(c(beta,sigma2,lambda),p+2,1)
     LB         <- matrix(c(-Inf*rep(1,p),1e-3,-Inf),p+2,1)
     UB         <- matrix(c(Inf*rep(1,p),Inf,Inf),p+2,1)
     tdir       <- optim(t0,ssmn.logL,method='L-BFGS-B',lower=LB,upper=UB,control=list(maxit=maxit),y=y,X=X,family="sn")
     theta      <- tdir$par
     beta       <- theta[1:p]
     sigma2     <- theta[p+1]
     lambda     <- theta[p+2]
     info       <- ssmn.info(y,X,theta, family)
     logL       <- logL(theta,y,X, family)
     se         <- sqrt(diag(solve(info)))

     pp         <- length(theta)
     AIC        <- (-logL/n)+(pp/n)
     BIC        <- (-logL/n)+((pp/2)*(log(n)/n))

     ttable           <- data.frame(cbind(c(theta),c(se)))
     namesrowAbetas   <- c(); for(i in 1:length(beta)){namesrowAbetas[i] <- paste("beta",i-1,sep="")}
     rownames(ttable) <- c(namesrowAbetas,"sigma2","lambda")
     colnames(ttable) <- c("Estimate","SE")
     end.time         <- Sys.time()
     time.taken       <- end.time - start.time
     obj.out          <- list(ttable = ttable, se=se, criterio=1e-8, time = time.taken, loglik=logL, beta = beta, sigma2 = sigma2, lambda = lambda, iter=tdir$counts[1], n = length(y), aic=AIC, bic=BIC,convergence=tdir$convergence==0)
    }
   if(family=="stn")
   {
     start.time  <- Sys.time() #Begin Time
     n          <- length(y)
     if(missing(X)) {X<-matrix(1,n,1)}
     p          <- ncol(X)
     lambda     <- as.numeric(skewness(y))
     beta       <- solve(t(X)%*%X)%*%t(X)%*%y
     sigma2     <- as.numeric(t((y-X%*%beta))%*%(y-X%*%beta)/(n-p))
     t0         <- as.matrix(c(beta,sigma2,lambda,20),p+3,1)
     LB         <- matrix(c(-Inf*rep(1,p),1e-3,-Inf,1),p+3,1)
     UB         <- matrix(c(Inf*rep(1,p),Inf,Inf,Inf),p+3,1)
     tdir       <- optim(t0,ssmn.logL,method='L-BFGS-B',lower=LB,upper=UB,control=list(maxit=maxit),y=y,X=X,family="stn")
     theta      <- tdir$par
     beta       <- theta[1:p]
     sigma2     <- theta[p+1]
     lambda     <- theta[p+2]
     nu         <- theta[p+3]
     info       <- ssmn.info(y,X,theta, family)
     logL       <- logL(theta,y,X, family)
     aux        <- diag(solve(info))
     se         <- sqrt(diag(solve(info[1:(p+2),1:(p+2)])))
     cse=rep(0,length(theta))
     cse[1:(p+2)]=se
     if (sum(aux<0)==0) {
       se <- sqrt(aux)
       cse=se}

     pp         <- length(theta)
     AIC        <- (-logL/n)+(pp/n)
     BIC        <- (-logL/n)+((pp/2)*(log(n)/n))

     ttable           <- data.frame(cbind(c(theta),c(cse)))
     namesrowAbetas   <- c(); for(i in 1:length(beta)){namesrowAbetas[i] <- paste("beta",i-1,sep="")}
     rownames(ttable) <- c(namesrowAbetas,"sigma2","lambda","nu")
     colnames(ttable) <- c("Estimate","SE")
     end.time         <- Sys.time()
     time.taken       <- end.time - start.time
     obj.out          <- list(ttable = ttable, se=se, criterio=1e-8, time = time.taken, loglik=logL, beta = beta, sigma2 = sigma2, lambda = lambda, nu=nu, iter=tdir$counts[1], n = length(y), aic=AIC, bic=BIC,convergence=tdir$convergence==0)

    }
   if(family=="ssl")
   {
     start.time <- Sys.time() #Begin Time
     n          <- length(y)
     if(missing(X)) {X<-matrix(1,n,1)}
     p          <- ncol(X)
     lambda     <- as.numeric(skewness(y))
     beta       <- solve(t(X)%*%X)%*%t(X)%*%y
     sigma2     <- as.numeric(t((y-X%*%beta))%*%(y-X%*%beta)/(n-p))
     t0         <- as.matrix(c(beta,sigma2,lambda,15),p+3,1)
     LB         <- matrix(c(-Inf*rep(1,p),1e-3,-Inf,1),p+3,1)
     UB         <- matrix(Inf*rep(1,p+3),p+3,1)
     tdir       <- optim(t0,ssmn.logL,method='L-BFGS-B',lower=LB,upper=UB,control=list(maxit=maxit),y=y,X=X,family="ssl")
     theta      <- tdir$par
     beta       <- theta[1:p]
     sigma2     <- theta[p+1]
     lambda     <- theta[p+2]
     nu         <- theta[p+3]
     info       <- ssmn.info(y,X,theta, family)
     logL       <- logL(theta,y,X, family)
     aux        <- diag(solve(info))
     se         <- sqrt(diag(solve(info[1:(p+2),1:(p+2)])))
     cse=rep(0,length(theta))
     cse[1:(p+2)]=se
     if (sum(aux<0)==0) {
       se <- sqrt(aux)
       cse=se}

     pp         <- length(theta)
     AIC        <- (-logL/n)+(pp/n)
     BIC        <- (-logL/n)+((pp/2)*(log(n)/n))

     ttable           <- data.frame(cbind(c(theta),c(cse)))
     namesrowAbetas   <- c(); for(i in 1:length(beta)){namesrowAbetas[i] <- paste("beta",i-1,sep="")}
     rownames(ttable) <- c(namesrowAbetas,"sigma2","lambda","nu")
     colnames(ttable) <- c("Estimate","SE")
     end.time         <- Sys.time()
     time.taken       <- end.time - start.time
     obj.out          <- list(ttable = ttable, se=se, criterio=1e-8, time = time.taken, loglik=logL, beta = beta, sigma2 = sigma2, lambda = lambda, nu = nu, iter=tdir$counts[1], n = length(y), aic=AIC, bic=BIC,convergence=tdir$convergence==0)

    }
   if(family=="scn")
   {
     start.time <- Sys.time() #Begin Time
     n          <- length(y)
     if(missing(X)) {X<-matrix(1,n,1)}
     p          <- ncol(X)
     lambda     <- as.numeric(skewness(y))
     beta       <- solve(t(X)%*%X)%*%t(X)%*%y
     sigma2     <- as.numeric(t((y-X%*%beta))%*%(y-X%*%beta)/(n-p))
     t0         <- as.matrix(c(beta,sigma2,lambda,0.999,0.999),p+4,1)
     LB         <- matrix(c(-Inf*rep(1,p),1e-3,-Inf,1e-3,1e-3),p+4,1)
     UB         <- matrix(c(Inf*rep(1,p),Inf,Inf,1,1),p+4,1)
     tdir       <- optim(t0,ssmn.logL,method='L-BFGS-B',lower=LB,upper=UB,control=list(maxit=maxit),y=y,X=X,family="scn")
     theta      <- tdir$par
     beta       <- theta[1:p]
     sigma2     <- theta[p+1]
     lambda     <- theta[p+2]
     nu         <- theta[p+3]
     gama       <- theta[p+4]
     info       <- ssmn.info(y,X,theta, family)
     logL       <- logL(theta,y,X, family)
     aux        <- diag(solve(info))
     se         <- sqrt(diag(solve(info[1:(p+2),1:(p+2)])))
     cse=rep(0,length(theta))
     cse[1:(p+2)]=se
     if (sum(aux<0)==0) {
       se <- sqrt(aux)
       cse=se}

     pp         <- length(theta)
     AIC        <- (-logL/n)+(pp/n)
     BIC        <- (-logL/n)+((pp/2)*(log(n)/n))

     ttable           <- data.frame(cbind(c(theta),c(cse)))
     namesrowAbetas   <- c(); for(i in 1:length(beta)){namesrowAbetas[i] <- paste("beta",i-1,sep="")}
     rownames(ttable) <- c(namesrowAbetas,"sigma2","lambda","nu","gama")
     colnames(ttable) <- c("Estimate","SE")
     end.time         <- Sys.time()
     time.taken       <- end.time - start.time
     obj.out          <- list(ttable = ttable, se=se, criterio=1e-8, time = time.taken, loglik=logL, beta = beta, sigma2 = sigma2, lambda = lambda, nu=nu, gama=gama, iter=tdir$counts[1], n = length(y), aic=AIC, bic=BIC,convergence=tdir$convergence==0)

    }
   if(family=="sep")
   {
     start.time <- Sys.time() #Begin Time
     n          <- length(y)
     if(missing(X)) {X<-matrix(1,n,1)}
     p          <- ncol(X)
     lambda     <- as.numeric(skewness(y))
     beta       <- solve(t(X)%*%X)%*%t(X)%*%y
     sigma2     <- as.numeric(t((y-X%*%beta))%*%(y-X%*%beta)/(n-p))
     t0         <- as.matrix(c(beta,sigma2,lambda,0.999),p+3,1)
     LB         <- matrix(c(-Inf*rep(1,p),1e-3,-Inf,0.51),p+3,1)
     UB         <- matrix(c(Inf*rep(1,p),Inf,Inf,1),p+3,1)
     tdir       <- optim(t0,ssmn.logL,method='L-BFGS-B',lower=LB,upper=UB,control=list(maxit=maxit),y=y,X=X,family="sep")
     theta      <- tdir$par
     beta       <- theta[1:p]
     sigma2     <- theta[p+1]
     lambda     <- theta[p+2]
     nu         <- theta[p+3]
     info       <- ssmn.info(y,X,theta, family)
     logL       <- logL(theta,y,X, family)
     aux        <- diag(solve(info))
     se         <- sqrt(diag(solve(info[1:(p+2),1:(p+2)])))
     cse=rep(0,length(theta))
     cse[1:(p+2)]=se
     if (sum(aux<0)==0) {
       se <- sqrt(aux)
       cse=se}

     pp         <- length(theta)
     AIC        <- (-logL/n)+(pp/n)
     BIC        <- (-logL/n)+((pp/2)*(log(n)/n))

     ttable           <- data.frame(cbind(c(theta),c(cse)))
     namesrowAbetas   <- c(); for(i in 1:length(beta)){namesrowAbetas[i] <- paste("beta",i-1,sep="")}
     rownames(ttable) <- c(namesrowAbetas,"sigma2","lambda","nu")
     colnames(ttable) <- c("Estimate","SE")
     end.time         <- Sys.time()
     time.taken       <- end.time - start.time
     obj.out          <- list(ttable = ttable, se=se, criterio=1e-8, time = time.taken, loglik=logL, beta = beta, sigma2 = sigma2, lambda = lambda, nu=nu, iter=tdir$counts[1], n = length(y), aic=AIC, bic=BIC,convergence=tdir$convergence==0)

    }
 }

 if(method=="EM")
 {
  if(family=="sn")
  {
   start.time  <- Sys.time() #Begin Time
   # Theta=[beta,sigma2,lambda]
   n           <- length(y)
   if(missing(X)) {X<-matrix(1,n,1)}
   p           <- ncol(X)
   beta        <- solve(t(X)%*%X)%*%t(X)%*%y
   res0        <- y-X%*%beta
   lambda      <- as.numeric(skewness(res0))
   sigma2      <- as.numeric(t((y-X%*%beta))%*%(y-X%*%beta)/(n-p))
   theta0      <- as.matrix(c(beta,sigma2,lambda),p+2,1)
   criterio    <- 1
   cont        <- 0
   while ((criterio > error)&&(cont<maxit))
   {
    cont       <- cont+1
    res        <- y-X%*%beta
    sigma      <- sqrt(sigma2)
    eta        <- lambda*res
    aux        <- eta/sigma
    aux1       <- pmax(aux,-37)
    Wphi       <- dnorm(aux1)/pnorm(aux1)
    t1         <- eta+sigma*Wphi
    t2         <- eta^2+sigma2+sigma*eta*Wphi
    d          <- res^2/sigma2
    beta       <- solve(t(X)%*%(diag(1,n)+lambda^2*diag(1,n))%*%X)%*%t(X)%*%(y-lambda*(t1-lambda*y))
    Qb         <- t(res)%*%res
    sigma2     <- as.numeric((Qb+sum(t2)-2*lambda*t(t1)%*%res+lambda^2*Qb)/(2*n))
    lambda     <- as.numeric(t(t1)%*%res/Qb)
    theta      <- as.matrix(c(beta,sigma2,lambda),p+2,1)
    dif        <- theta-theta0
    criterio   <- (sum(dif^2))^0.5
    theta0     <- theta
   }

    info       <- ssmn.info(y,X,theta, family)
    logL       <- ssmn.logL(theta,y,X, family)
    se         <- sqrt(diag(solve(info)))


    pp         <- length(theta)
    AIC        <- (-logL/n)+(pp/n)
    BIC        <- (-logL/n)+((pp/2)*(log(n)/n))

    ttable           <- data.frame(cbind(c(theta),c(se)))
    namesrowAbetas   <- c(); for(i in 1:length(beta)){namesrowAbetas[i] <- paste("beta",i-1,sep="")}
    rownames(ttable) <- c(namesrowAbetas,"sigma2","lambda")
    colnames(ttable) <- c("Estimate","SE")
    end.time         <- Sys.time()
    time.taken       <- end.time - start.time
    obj.out          <- list(ttable = ttable, se=se, criterio = criterio, time = time.taken, loglik=logL, beta = beta, sigma2 = sigma2, lambda = lambda, iter = cont, n = length(y), aic=AIC, bic=BIC,convergence = (criterio < error))
  }

  #The skew-Student-t normal distribution
  if(family=="stn")
  {
    start.time  <- Sys.time() #Begin Time
    # Theta=[beta,sigma2,lambda,nu]
    n          <- length(y)
    if(missing(X)) {X<-matrix(1,n,1)}
    p          <- ncol(X)
    beta       <- solve(t(X)%*%X)%*%t(X)%*%y
    res0       <- y-X%*%beta
    lambda     <- as.numeric(skewness(res0))
    sigma2     <- as.numeric(t((y-X%*%beta))%*%(y-X%*%beta)/(n-p))
    t0         <- as.matrix(c(beta,sigma2,lambda),p+2,1)
    theta0     <- rbind(t0,20)
    criterio   <- 1
    cont       <- 0
    while ((criterio > 1e-6)&&(cont<maxit))
    {
     cont      <- cont+1
     res       <- y-X%*%beta
     sigma     <- sqrt(sigma2)
     eta       <- lambda*res
     aux       <- eta/sigma
     aux1      <- pmax(aux,-37)
     Wphi      <- dnorm(aux1)/pnorm(aux1)
     t1        <- eta+sigma*Wphi
     t2        <- eta^2+sigma2+sigma*eta*Wphi
     d         <- res^2/sigma2
     nu        <- optimize(ftn,c(1,30),y=y,X=X,theta=theta0)$minimum
     k         <- as.vector((nu+1)/(nu+d))
     beta      <- solve(t(X)%*%(diag(k)+lambda^2*diag(1,n))%*%X)%*%t(X)%*%(diag(k)%*%y-lambda*(t1-lambda*y))
     Qb        <- t(res)%*%res
     Qkb       <- t(res)%*%diag(k)%*%res
     sigma2    <- as.numeric((Qkb+sum(t2)-2*lambda*t(t1)%*%res+lambda^2*Qb)/(2*n))
     lambda    <- as.numeric(t(t1)%*%res/Qb)
     theta     <- as.matrix(c(beta,sigma2,lambda,nu),p+3,1)
     dif       <- theta-theta0
     criterio  <- (sum(dif^2))^0.5
     theta0    <- theta
    }

    info       <- ssmn.info(y,X,theta, family)
    logL       <- ssmn.logL(theta,y,X, family)
    aux        <- diag(solve(info))
    se         <- sqrt(diag(solve(info[1:(p+2),1:(p+2)])))
    cse=rep(0,length(theta))
    cse[1:(p+2)]=se
    if (sum(aux<0)==0) {
      se <- sqrt(aux)
      cse=se}

    pp         <- length(theta)
    AIC        <- (-logL/n)+(pp/n)
    BIC        <- (-logL/n)+((pp/2)*(log(n)/n))

    ttable           <- data.frame(cbind(c(theta),c(cse)))
    namesrowAbetas   <- c(); for(i in 1:length(beta)){namesrowAbetas[i] <- paste("beta",i-1,sep="")}
    rownames(ttable) <- c(namesrowAbetas,"sigma2","lambda","nu")
    colnames(ttable) <- c("Estimate","SE")
    end.time         <- Sys.time()
    time.taken       <- end.time - start.time
    obj.out          <- list(ttable = ttable, se=se, criterio = criterio, time = time.taken, loglik=logL, beta = beta, sigma2 = sigma2, lambda = lambda, nu=nu, iter = cont, n = length(y), aic=AIC, bic=BIC,convergence = (criterio < error))

  }

  #The skew-slash distribution

  if(family=="ssl")
  {
    start.time  <- Sys.time() #Begin Time
    # Theta=[beta,sigma2,lambda,nu]
    n          <- length(y)
    if(missing(X)) {X<-matrix(1,n,1)}
    p          <- ncol(X)
    beta       <- solve(t(X)%*%X)%*%t(X)%*%y
    res0       <- y-X%*%beta
    lambda     <- as.numeric(skewness(res0))
    sigma2     <- as.numeric(t((y-X%*%beta))%*%(y-X%*%beta)/(n-p))
    t0         <- as.matrix(c(beta,sigma2,lambda),p+2,1)
    theta0     <- rbind(t0,10)
    criterio   <- 1
    cont       <- 0
    while ((criterio > 1e-6)&&(cont<maxit))
    {
     cont      <- cont+1
     res       <- y-X%*%beta
     sigma     <- sqrt(sigma2)
     eta       <- lambda*res
     aux       <- eta/sigma
     aux1      <- pmax(aux,-37)
     Wphi      <- dnorm(aux1)/pnorm(aux1)
     t1        <- eta+sigma*Wphi
     t2        <- eta^2+sigma2+sigma*eta*Wphi
     d         <- res^2/sigma2
     nu        <- optimize(fsl,c(0.1,20),y,X,theta0)$minimum
     k         <- as.vector((2*nu+1)/d*(pgamma(1,nu+1.5,scale=2/d)/pgamma(1,nu+0.5,scale=2/d)))
     beta      <- solve(t(X)%*%(diag(k)+lambda^2*diag(1,n))%*%X)%*%t(X)%*%(diag(k)%*%y-lambda*(t1-lambda*y))
     Qb        <- t(res)%*%res
     Qkb       <- t(res)%*%diag(k)%*%res
     sigma2    <- as.numeric((Qkb+sum(t2)-2*lambda*t(t1)%*%res+lambda^2*Qb)/(2*n))
     lambda    <- as.numeric(t(t1)%*%res/Qb)
     theta     <- as.matrix(c(beta,sigma2,lambda,nu),p+3,1)
     dif       <- theta-theta0
     criterio  <- (sum(dif^2))^0.5
     theta0    <- theta
    }

    info       <- ssmn.info(y,X,theta, family)
    logL       <- ssmn.logL(theta,y,X, family)
    aux        <- diag(solve(info))
    se         <- sqrt(diag(solve(info[1:(p+2),1:(p+2)])))
    cse=rep(0,length(theta))
    cse[1:(p+2)]=se
    if (sum(aux<0)==0) {
      se <- sqrt(aux)
      cse=se}

    pp         <- length(theta)
    AIC        <- (-logL/n)+(pp/n)
    BIC        <- (-logL/n)+((pp/2)*(log(n)/n))

    ttable           <- data.frame(cbind(c(theta),c(cse)))
    namesrowAbetas   <- c(); for(i in 1:length(beta)){namesrowAbetas[i] <- paste("beta",i-1,sep="")}
    rownames(ttable) <- c(namesrowAbetas,"sigma2","lambda","nu")
    colnames(ttable) <- c("Estimate","SE")
    end.time         <- Sys.time()
    time.taken       <- end.time - start.time
    obj.out          <- list(ttable = ttable, se=se, criterio = criterio, time = time.taken, loglik=logL, beta = beta, sigma2 = sigma2, lambda = lambda, nu=nu, iter = cont, n = length(y), aic=AIC, bic=BIC,convergence = (criterio < error))

  }

  #The skew-contaminated normal distribution
  if(family=="scn")
  {
    start.time <- Sys.time() #Begin Time
    # Theta=[beta,sigma2,lambda,nu,gama]
    n          <- length(y)
    if(missing(X)) {X<-matrix(1,n,1)}
    p          <- ncol(X)
    beta       <- solve(t(X)%*%X)%*%t(X)%*%y
    res0       <- y-X%*%beta
    lambda     <- as.numeric(skewness(res0))
    sigma2     <- as.numeric(t((y-X%*%beta))%*%(y-X%*%beta)/(n-p))
    t0         <- as.matrix(c(beta,sigma2,lambda),p+2,1)
    theta0     <- rbind(t0,0,1)
    pteta      <- length(theta0)
    criterio   <- 1
    cont       <- 0
    L1         <- rbind(1e-6,1e-6)
    L2         <- rbind(0.9999999,0.9999999)
    while ((criterio > 1e-6)&&(cont<maxit))
    {
     cont      <- cont+1
     res       <- y-X%*%beta
     sigma     <- sqrt(sigma2)
     eta       <- lambda*res
     aux       <- eta/sigma
     aux1      <- pmax(aux,-37)
     Wphi      <- dnorm(aux1)/pnorm(aux1)
     t1        <- eta+sigma*Wphi
     t2        <- eta^2+sigma2+sigma*eta*Wphi
     d         <- res^2/sigma2
     nugama0   <- signif(theta0[c(pteta-1,pteta)],digits=7)
     nugama    <- optim(nugama0,fcn,method='L-BFGS-B',lower=L1,upper=L2,y=y,X=X,theta=theta0)$par
     nu        <- nugama[1]
     gama      <- nugama[2]
     k         <- as.vector((1-nu+nu*gama^(1.5)*exp((1-gama)*d/2))/(1-nu+nu*gama^(0.5)*exp((1-gama)*d/2)))
     beta      <- solve(t(X)%*%(diag(k)+lambda^2*diag(1,n))%*%X)%*%t(X)%*%(diag(k)%*%y-lambda*(t1-lambda*y))
     Qb        <- t(res)%*%res
     Qkb       <- t(res)%*%diag(k)%*%res
     sigma2    <- as.numeric((Qkb+sum(t2)-2*lambda*t(t1)%*%res+lambda^2*Qb)/(2*n))
     lambda    <- as.numeric(t(t1)%*%res/Qb)
     theta     <- as.matrix(c(beta,sigma2,lambda,nu,gama),p+4,1)
     dif       <- theta-theta0
     criterio  <- (sum(dif^2))^0.5
     theta0    <- theta
    }

    info       <- ssmn.info(y,X,theta, family)
    logL       <- ssmn.logL(theta,y,X, family)
    aux        <- diag(solve(info))
    se         <- sqrt(diag(solve(info[1:(p+2),1:(p+2)])))
    cse=rep(0,length(theta))
    cse[1:(p+2)]=se
    if (sum(aux<0)==0) {
      se <- sqrt(aux)
      cse=se}

    pp         <- length(theta)
    AIC        <- (-logL/n)+(pp/n)
    BIC        <- (-logL/n)+((pp/2)*(log(n)/n))

    ttable           <- data.frame(cbind(c(theta),c(cse)))
    namesrowAbetas   <- c(); for(i in 1:length(beta)){namesrowAbetas[i] <- paste("beta",i-1,sep="")}
    rownames(ttable) <- c(namesrowAbetas,"sigma2","lambda","nu","gama")
    colnames(ttable) <- c("Estimate","SE")
    end.time         <- Sys.time()
    time.taken       <- end.time - start.time
    obj.out          <- list(ttable = ttable, se=se, criterio = criterio, time = time.taken, loglik=logL, beta = beta, sigma2 = sigma2, lambda = lambda, nu=nu, gama=gama, iter = cont, n = length(y), aic=AIC, bic=BIC,convergence = (criterio < error))
  }

  #The skew-exponential power distribution
  if(family=="sep")
  {
    start.time <- Sys.time() #Begin Time
    # Theta=[beta,sigma2,lambda,nu]
    n          <- length(y)
    if(missing(X)) {X<-matrix(1,n,1)}
    p          <- ncol(X)
    beta       <- solve(t(X)%*%X)%*%t(X)%*%y
    res0       <- y-X%*%beta
    lambda     <- as.numeric(skewness(res0))
    sigma2     <- as.numeric(t((y-X%*%beta))%*%(y-X%*%beta)/(n-p))
    t0         <- as.matrix(c(beta,sigma2,lambda),p+2,1)
    theta0     <- rbind(t0,1)
    criterio   <- 1
    cont       <- 0
    while ((criterio > 1e-6)&&(cont<maxit))
    {
     cont      <- cont+1
     res       <- y-X%*%beta
     sigma     <- sqrt(sigma2)
     eta       <- lambda*res
     aux       <- eta/sigma
     aux1      <- pmax(aux,-36)
     Wphi      <- dnorm(aux1)/pnorm(aux1)
     t1        <- eta+sigma*Wphi
     t2        <- eta^2+sigma2+sigma*eta*Wphi
     d         <- res^2/sigma2
     nu        <- optimize(fep,c(0.55,1),y,X,theta0)$minimum
     k         <- as.vector(nu*d^(nu-1))
     binv      <- solve(t(X)%*%(diag(k)+lambda^2*diag(1,n))%*%X)
     beta      <- binv%*%t(X)%*%(diag(k)%*%y-lambda*(t1-lambda*y))
     Qb        <- t(res)%*%res
     Qkb       <- t(res)%*%diag(k)%*%res
     sigma2    <- as.numeric((Qkb+sum(t2)-2*lambda*t(t1)%*%res+lambda^2*Qb)/(2*n))
     lambda    <- as.numeric(t(t1)%*%res/Qb)
     theta     <- as.matrix(c(beta,sigma2,lambda,nu),p+3,1)
     dif       <- theta-theta0
     criterio  <- (sum(dif^2))^0.5
     theta0    <- theta
    }

    info       <- ssmn.info(y,X,theta, family)
    logL       <- ssmn.logL(theta,y,X, family)
    aux        <- diag(solve(info))
    se         <- sqrt(diag(solve(info[1:(p+2),1:(p+2)])))
    cse=rep(0,length(theta))
    cse[1:(p+2)]=se
    if (sum(aux<0)==0) {
      se <- sqrt(aux)
      cse=se}

    pp         <- length(theta)
    AIC        <- (-logL/n)+(pp/n)
    BIC        <- (-logL/n)+((pp/2)*(log(n)/n))

    ttable           <- data.frame(cbind(c(theta),c(cse)))
    namesrowAbetas   <- c(); for(i in 1:length(beta)){namesrowAbetas[i] <- paste("beta",i-1,sep="")}
    rownames(ttable) <- c(namesrowAbetas,"sigma2","lambda","nu")
    colnames(ttable) <- c("Estimate","SE")
    end.time         <- Sys.time()
    time.taken       <- end.time - start.time
    obj.out          <- list(ttable = ttable, se=se, criterio = criterio, time = time.taken, loglik=logL, beta = beta, sigma2 = sigma2, lambda = lambda, nu=nu, iter = cont, n = length(y), aic=AIC, bic=BIC,convergence = (criterio < error))
  }
 }

  class(obj.out) <- family
  return(obj.out)
}
