"spCopulaDDP" <- function(y, delta, x=NULL, s, prediction, prior, mcmc, state, status=TRUE, FSA = TRUE, knots,
                                  data=sys.frame(sys.parent()), na.action=na.fail, work.dir=NULL) {
  #########################################################################################
  # call parameters
  #########################################################################################
  m <- mcall <- cl <- match.call()
  
  #########################################################################################
  # data structure and check if FSA is needed
  #########################################################################################
  y <- as.vector(y); n <- length(y);
  delta <- as.vector(delta);
  X <- cbind(rep(1,n), x);
  p <- ncol(X);
  s <- t(s);
  dnn <- .DistMat(s, s);
  if(FSA) {
    ss <- knots$ss; if(is.null(ss)) stop("please specify the knots vector ss if FSA is used");
    ss <- t(as.matrix(ss));
    id <- knots$blockid; if(is.null(id)) stop("please specify the bolock id that each observation is located in if FSA is used ");
    orderindex = order(id);
    if(!(sum(orderindex==(1:n))==n)) stop("please sort the data by ID");
    blocki = c( 0, cumsum(as.vector(table(id))) );
    dnm <- .DistMat(s, ss);
    dmm <- .DistMat(ss, ss);
    #########################################################################################
    # prediction
    #########################################################################################
    s0 <- t(prediction$spred);
    xnew <- prediction$xpred;
    if(is.null(xnew)) stop("please specify xpred");
    if(is.vector(xnew)) {
      npred = length(xnew)
    } else npred = nrow(xnew);
    xpred <- cbind(rep(1,npred), xnew);
    if(!(ncol(xpred)==p)) { stop ("error: ncol(xpred) is not equal to ncol(x)");} 
    ds0n <- .DistMat(s, s0);
    ds0m <- .DistMat(ss, s0);
    ds0block <- matrix(0, n, npred);
    predid <- prediction$predid; 
    if(is.null(predid)) stop("please specify the pred id that each observation is located in if FSA is used ");
    for(i in 1:n){
      for(j in 1:npred){
        ds0block[i,j] = (id[i]==predid[j])+0
      }
    }
  }

  #########################################################################################
  # change working directory (if requested..)
  #########################################################################################
  if(!is.null(work.dir))
  {
    cat("\n Changing working directory to ",work.dir,"\n")
    old.dir <- getwd()  # by default work in current working directory
    setwd(work.dir)
  }
  model.name <- "Spatial copula model for point-referenced time-to-event data"
  
  #########################################################################################
  # prediction
  #########################################################################################
  s0 <- t(prediction$spred);
  xnew <- prediction$xpred;
  if(is.null(xnew)) stop("please specify xpred");
  if(is.vector(xnew)) {
    npred = length(xnew)
  } else npred = nrow(xnew);
  xpred <- cbind(rep(1,npred), xnew);
  if(!(ncol(xpred)==p)) { stop ("error: ncol(xpred) is not equal to ncol(x)");} 
  ds0n <- .DistMat(s, s0);
  
  #########################################################################################
  # initial analysis and mcmc parameters
  #########################################################################################
  nburn <- mcmc$nburn;
  nsave <- mcmc$nsave;
  nskip <- mcmc$nskip;
  ndisplay <- mcmc$ndisplay;
  
  #########################################################################################
  # priors
  #########################################################################################
  fit0 <- survival::survreg(formula = Surv(exp(y), delta) ~ x, dist = "lognormal");
  #fit0=lm(y~x); sfit0=summary(fit0); sig2hat = sfit0$sigma^2; 
  muhat = as.vector(fit0$coefficients);
  sig2hat = fit0$scale^2
  Sighat = as.matrix(fit0$var[(1:p),(1:p)]); Sigscale=30;
  N <- prior$N; if(is.null(N)) N <- 10;
  m0 <- prior$m0; if(is.null(m0)) m0 <- muhat;
  S0 <- prior$S0; if(is.null(S0)) S0 <- Sighat;
  Sig0 <- prior$Sig0; if(is.null(Sig0)) Sig0 <- Sigscale*Sighat; #Sig0 <- diag(rep(1e5,p), nrow=p, ncol=p); #
  k0 <- prior$k0; if(is.null(k0)) k0 <- p+5;
  nua <-prior$nua; nub <- prior$nub;
  if(is.null(nua)) nua=2+1; #nua=2+sig2hat/4; #
  if(is.null(nub)) nub=sig2hat; #nub=sig2hat/4*(nua-1); #
  a0 <-prior$a0; b0 <- prior$b0;
  if(is.null(a0)) a0=2; if(is.null(b0)) b0=2;
  theta0 <- prior$theta0; if(is.null(theta0)) theta0 <- c(1.0, 1.0, 1.0, 1.0);
  spl0 <- prior$spl0; if(is.null(spl0)) spl0 <- round(nburn/2);
  spS0 <- prior$spS0; if(is.null(spS0)) spS0 <- diag(c(0.5, 0.1));
  spadapter <- prior$spadapter; if(is.null(spadapter)) spadapter <- (2.38)^2/2;

  #########################################################################################
  # current state and mcmc specification
  #########################################################################################
  if(status){
    currenty = y;
    mu <- state$mu; if(is.null(mu)) mu <- muhat;
    Sig <- state$Sig; if(is.null(Sig)) Sig <- 25*Sighat;
    beta<- state$beta; if(is.null(beta)) beta = matrix(muhat, p, N);
    sigma2<- state$sigma2; if(is.null(sigma2)) sigma2 = rep(sig2hat/2,N);
    alpha <- state$alpha; if(is.null(alpha)) alpha <- 2;
    K = sample(1:N,n, replace=T);
    V = rbeta(N, 1, alpha); V[N] =1;
    w = V; 
    for (k in 2:N){
      w[k] = max(exp( sum(log(1-V[1:(k-1)]))+log(V[k]) ), 1e-320);
      #w[k] = max( (1 - sum(w[1:(k-1)]))*V[k], 1e-320);
    }
    theta = state$theta; if(is.null(theta)) theta <- c(0.98, 0.2);
  }else{
    K = state$K;
    currenty = state$y;
    V = state$V;
    w = as.vector(state$w);
    beta = state$beta;
    sigma2 = state$sigma2;
    alpha = state$alpha;
    mu = as.vector(state$mu);
    Sig = state$Sig;
    theta = as.vector(state$theta);
    spl0 = 500;
    spS0 = spadapter*state$spSnew;
  }
  
  #########################################################################################
  # calling the c++ code
  #########################################################################################
  if(FSA){
    foo <- .Call("spCopulaDDP_FSA", 
                 nburn_ = nburn, 
                 nsave_ = nsave, 
                 nskip_ = nskip, 
                 ndisplay_ = ndisplay,
                 y_ = currenty, 
                 delta_ = delta, 
                 X_ = as.matrix(t(X)),
                 N_ = N,
                 beta_ = beta, 
                 tau2_ = 1.0/sigma2,
                 K_ = K, 
                 V_ = V,
                 w_ = w,
                 alpha_ = alpha, 
                 mu_ = mu, 
                 Sig_ = Sig,
                 m0_ = m0, 
                 S0_ = S0, 
                 Sig0_ = Sig0, 
                 k0_ = k0,
                 a0_ = a0, 
                 b0_ = b0, 
                 nua_ = nua, 
                 nub_ = nub,
                 xpred_ = as.matrix(xpred), 
                 ds0n_ = ds0n,
                 dnn_ = dnn,
                 theta_ = theta,
                 theta0_ = theta0,
                 spl0_ = spl0,
                 spS0_ = spS0, 
                 spadapter_ = spadapter,
                 dnm_ = dnm, dmm_=dmm, blocki_=blocki,
                 ds0m_= ds0m, ds0block_=ds0block,
                 status_ = status,
                 PACKAGE = "spBayesSurv")
  } else{
    foo <- .Call("spCopulaDDP", 
                 nburn_ = nburn, 
                 nsave_ = nsave, 
                 nskip_ = nskip, 
                 ndisplay_ = ndisplay,
                 y_ = currenty, 
                 delta_ = delta, 
                 X_ = as.matrix(t(X)),
                 N_ = N,
                 beta_ = beta, 
                 tau2_ = 1.0/sigma2,
                 K_ = K, 
                 V_ = V,
                 w_ = w,
                 alpha_ = alpha, 
                 mu_ = mu, 
                 Sig_ = Sig,
                 m0_ = m0, 
                 S0_ = S0, 
                 Sig0_ = Sig0, 
                 k0_ = k0,
                 a0_ = a0, 
                 b0_ = b0, 
                 nua_ = nua, 
                 nub_ = nub,
                 xpred_ = as.matrix(xpred), 
                 ds0n_ = ds0n,
                 dnn_ = dnn,
                 theta_ = theta,
                 theta0_ = theta0,
                 spl0_ = spl0,
                 spS0_ = spS0, 
                 spadapter_ = spadapter,
                 status_ = status,
                 PACKAGE = "spBayesSurv")
  }
  #########################################################################################
  # output
  #########################################################################################
  output <- list(modelname=model.name,
                 beta = foo$beta,
                 sigma2 = foo$sigma2,
                 w = foo$w,
                 alpha = foo$alpha,
                 theta1 = foo$theta1,
                 theta2 = foo$theta2,
                 z = foo$z,
                 y = foo$y,
                 ratey = foo$ratey,
                 ratebeta = foo$ratebeta,
                 ratesigma = foo$ratesigma,
                 rateV = foo$rateV,
                 ratetheta = foo$ratetheta,
                 cpo = foo$cpo,
                 Ypred = foo$Ypred,
                 Zpred = foo$Zpred,
                 V = foo$V,
                 K = foo$K,
                 state=foo$state);
  class(output) <- c("spCopulaDDP")
  output
}