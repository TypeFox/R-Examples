"spCopulaCoxph" <- function(y, delta, x=NULL, s, prediction, prior, mcmc, state, RandomIntervals=F,
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
    X <- t(cbind(x));
    p <- nrow(X);
    s <- t(s);
    dnn <- .DistMat(s, s);
    
    #########################################################################################
    # change working directory (if requested..)
    #########################################################################################
    if(!is.null(work.dir))
    {
      cat("\n Changing working directory to ",work.dir,"\n")
      old.dir <- getwd()  # by default work in current working directory
      setwd(work.dir)
    }
    model.name <- "spatial Copula Cox PH model for time-to-event data"
    
    #########################################################################################
    # prediction
    #########################################################################################
    s0 <- t(prediction$spred);
    xnew <- prediction$xpred;
    if(is.null(xnew)) stop("please specify xpred")
    if(is.vector(xnew)) {
      npred = length(xnew)
    } else npred = nrow(xnew);
    xpred <- cbind(xnew);
    if(!(ncol(xpred)==p)) stop("error: ncol(xpred) is not equal to ncol(x)") 
    ds0n <- .DistMat(s, s0);
    
    #########################################################################################
    # initial analysis and mcmc parameters
    #########################################################################################
    fit0 = survival::survreg( Surv(y, delta) ~ t(X), dist="weibull" );
    thbetaShat0 = fit0$var[c(1,p+2,2:(1+p)),c(1,p+2,2:(1+p))];
    Dtrans = diag(c(-1, -1, rep(-1/fit0$scale,p))); Dtrans[2,-(1:2)]=fit0$coefficients[-1]/fit0$scale;
    thbetaShat = t(Dtrans)%*%thbetaShat0%*%Dtrans;
    ## prior for beta
    Shat = thbetaShat[-(1:2),-(1:2)];
    multShat0 = diag(rep(1e5,p), nrow=p, ncol=p); 
    multShat1 = as.matrix(0.5*Shat);
    
    nburn <- mcmc$nburn;
    nsave <- mcmc$nsave;
    nskip <- mcmc$nskip;
    ndisplay <- mcmc$ndisplay;
    
    #########################################################################################
    # priors
    #########################################################################################
    r0 <- prior$r0; if(is.null(r0)) r0 = 1;
    h0 <- prior$h0; if(is.null(h0)) h0 = as.vector( exp( -fit0$coefficients[1] ) );
    nu0 <- prior$nu0; if(is.null(nu0)) nu0 = 2;
    V0 <- prior$V0; if(is.null(V0)) V0 = nu0*as.vector( exp( -2*fit0$coefficients[1] )*fit0$var[1,1] );
    hl0 <- prior$hl0; if(is.null(hl0)) hl0 <- round(nburn/2);
    hs0 <- prior$hs0; if(is.null(hs0)) hs0 <- 0.5*sqrt( as.vector( exp( -2*fit0$coefficients[1] )*fit0$var[1,1] ) );
    hadapter <- prior$hadapter; if(is.null(hadapter)) hadapter <- (2.38)^2;
    M <- prior$M; if(is.null(M)) M <- 20;
    M1<- M+1;
    d <- prior$d; 
    if(is.null(d)){
      d = as.vector(quantile(y, probs=seq(0,1,length=M1)));
      d = d[-1];
      d[M] = Inf;
    }
    d <- c(0, d);
    if(!(M1==length(d))) stop("error: M is not equal to length(d)");
    mu0 <- prior$mu0; if(is.null(mu0)) mu0 <- as.vector( -fit0$coefficients[-1]/fit0$scale );
    Sig0 <- prior$Sig0; if(is.null(Sig0)) Sig0 <- multShat0;
    l0 <- prior$l0; if(is.null(l0)) l0 <- round(nburn/2);
    S0 <- prior$S0; if(is.null(S0)) S0 <- multShat1;
    adapter <- prior$adapter; if(is.null(adapter)) adapter <- (2.38)^2/p;
    theta0 <- prior$theta0; if(is.null(theta0)) theta0 <- c(1.0, 1.0, 1.0, 1.0);
    spl0 <- prior$spl0; if(is.null(spl0)) spl0 <- round(nburn/2);
    spS0 <- prior$spS0; if(is.null(spS0)) spS0 <- diag(c(0.5,0.1));
    spadapter <- prior$spadapter; if(is.null(spadapter)) spadapter <- (2.38)^2/2;
    
    #########################################################################################
    # current state and mcmc specification
    #########################################################################################
    h = c(0, rep(h0, M));
    beta = as.vector( -fit0$coefficients[-1] );
    theta = state$theta; if(is.null(theta)) theta <- c(0.98, 1);
    
    #########################################################################################
    # calling the c++ code
    #########################################################################################
    if(RandomIntervals){
      foo <- .Call("spCopulaCoxphR", 
                   nburn_ = nburn, 
                   nsave_ = nsave, 
                   nskip_ = nskip, 
                   ndisplay_ = ndisplay,
                   t_ = y,
                   delta_ = delta,
                   X_ = as.matrix(X), 
                   d_ = d,
                   h_ = h,
                   r0_ = r0,
                   h0_ = h0,
                   V0_ = V0,
                   hl0_ = hl0, 
                   hs0_ = hs0, 
                   hadapter_ = hadapter,
                   beta_ = beta, 
                   mu0_ = mu0, 
                   Sig0_ = Sig0,
                   l0_ = l0, 
                   S0_ = S0, 
                   adapter_ = adapter,
                   xpred_ = as.matrix(xpred),
                   ds0n_ = ds0n,
                   dnn_ = dnn,
                   theta_ = theta,
                   theta0_ = theta0,
                   spl0_ = spl0,
                   spS0_ = spS0, 
                   spadapter_ = spadapter,
                   PACKAGE = "spBayesSurv");
      output <- list(modelname=model.name,
                     t = foo$t,
                     z = foo$z,
                     h = foo$h,
                     d = foo$d,
                     beta = foo$beta,
                     hcen = foo$hcen,
                     theta1 = foo$theta1,
                     theta2 = foo$theta2,
                     ratebeta = foo$ratebeta,
                     ratetheta = foo$ratetheta,
                     rateh = foo$rateh,
                     ratehcen = foo$ratehcen,
                     cpo = foo$cpo,
                     Tpred = foo$Tpred,
                     Zpred = foo$Zpred);
    } else{
      foo <- .Call("spCopulaCoxph", 
                   nburn_ = nburn, 
                   nsave_ = nsave, 
                   nskip_ = nskip, 
                   ndisplay_ = ndisplay,
                   t_ = y,
                   delta_ = delta,
                   X_ = as.matrix(X), 
                   d_ = d,
                   h_ = h,
                   r0_ = r0,
                   h0_ = h0,
                   beta_ = beta, 
                   mu0_ = mu0, 
                   Sig0_ = Sig0,
                   l0_ = l0, 
                   S0_ = S0, 
                   adapter_ = adapter,
                   xpred_ = as.matrix(xpred),
                   ds0n_ = ds0n,
                   dnn_ = dnn,
                   theta_ = theta,
                   theta0_ = theta0,
                   spl0_ = spl0,
                   spS0_ = spS0, 
                   spadapter_ = spadapter,
                   PACKAGE = "spBayesSurv");
      output <- list(modelname=model.name,
                     t = foo$t,
                     z = foo$z,
                     h = foo$h,
                     d = foo$d,
                     beta = foo$beta,
                     theta1 = foo$theta1,
                     theta2 = foo$theta2,
                     ratebeta = foo$ratebeta,
                     ratetheta = foo$ratetheta,
                     rateh = foo$rateh,
                     cpo = foo$cpo,
                     Tpred = foo$Tpred,
                     Zpred = foo$Zpred);
    }

    cat("\n\n")
    class(output) <- c("spCopulaCoxph")
    output
  }