"anovaDDP" <- function(y, delta, x=NULL, prediction, prior, mcmc, state, status=TRUE, 
                       data=sys.frame(sys.parent()), na.action=na.fail, work.dir=NULL) {
    #########################################################################################
    # call parameters
    #########################################################################################
    m <- mcall <- cl <- match.call()
    
    #########################################################################################
    # data structure
    #########################################################################################
    y <- as.vector(y);
    n <- length(y);
    X <- cbind(rep(1,n), x);
    p <- ncol(X);
    
    #########################################################################################
    # change working directory (if requested..)
    #########################################################################################
    if(!is.null(work.dir))
    {
      cat("\n Changing working directory to ",work.dir,"\n")
      old.dir <- getwd()  # by default work in current working directory
      setwd(work.dir)
    }
    model.name <- "ANOVA DDP model for point-referenced time-to-event data"
    
    #########################################################################################
    # prediction
    #########################################################################################
    xnew <- prediction$xpred;
    if(is.null(xnew)) return("please specify xpred")
    if(is.vector(xnew)) {
      npred = length(xnew)
    } else npred = nrow(xnew);
    xpred <- cbind(rep(1,npred), xnew);
    if(!(ncol(xpred)==p)) { return ("error: ncol(xpred) is not equal to ncol(x)");} 

    #########################################################################################
    # initial analysis and priors
    #########################################################################################
    fit0 <- survival::survreg(formula = Surv(exp(y), delta) ~ x, dist = "lognormal");
    #fit0=lm(y~x); sfit0=summary(fit0); sig2hat = sfit0$sigma^2; 
    muhat = as.vector(fit0$coefficients);
    sig2hat = fit0$scale^2
    Sighat = as.matrix(fit0$var[(1:p),(1:p)]); Sigscale=30
    N <- prior$N; if(is.null(N)) N <- 10;
    m0 <- prior$m0; if(is.null(m0)) m0 <- muhat;
    S0 <- prior$S0; if(is.null(S0)) S0 <- Sighat;
    Sig0 <- prior$Sig0; if(is.null(Sig0)) Sig0 <- Sigscale*Sighat; #Sig0 <- diag(rep(1e5,p), nrow=p, ncol=p); 
    k0 <- prior$k0; if(is.null(k0)) k0 <- p+5;
    nua <-prior$nua; nub <- prior$nub;
    if(is.null(nua)) nua=2+1; #nua=2+sig2hat/4; #
    if(is.null(nub)) nub=sig2hat; #nub=sig2hat/4*(nua-1); #
    a0 <-prior$a0; b0 <- prior$b0;
    if(is.null(a0)) a0=2; if(is.null(b0)) b0=2;
    
    #########################################################################################
    # current state and mcmc specification
    #########################################################################################
    nburn <- mcmc$nburn;
    nsave <- mcmc$nsave;
    nskip <- mcmc$nskip;
    ndisplay <- mcmc$ndisplay;
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
    }

    #########################################################################################
    # calling the c++ code
    #########################################################################################
    foo <- .Call("anovaDDP", 
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
                 PACKAGE = "spBayesSurv")
    
    #########################################################################################
    # output
    #########################################################################################
    output <- list(modelname=model.name,
                   beta = foo$beta,
                   sigma2 = foo$sigma2,
                   w = foo$w,
                   alpha = foo$alpha,
                   y = foo$y,
                   cpo = foo$cpo,
                   Ypred = foo$Ypred,
                   V = foo$V,
                   K = foo$K,
                   state=foo$state);
    class(output) <- c("anovaDDP")
    output
  }
