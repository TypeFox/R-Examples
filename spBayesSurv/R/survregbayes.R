"survregbayes" <- function (formula, data, na.action, survmodel="PH", 
                            mcmc=list(nburn=3000, nsave=2000, nskip=0, ndisplay=500), 
                            prior=NULL, state=NULL, frailty=NULL, ID=NULL, Proximity=NULL) {
    #########################################################################################
    # call parameters
    #########################################################################################
    Call <- match.call(); # save a copy of the call 
    indx <- match(c("formula", "data", "na.action"),
                  names(Call), nomatch=0) 
    if (indx[1] ==0) stop("A formula argument is required");
    temp <- Call[c(1,indx)]  # only keep the arguments we wanted
    temp[[1]] <- as.name('model.frame')  # change the function called
    
    special <- c("strata", "cluster")
    temp$formula <- if(missing(data)) terms(formula, special)
    else              terms(formula, special, data=data)
    if (is.R()) m <- eval(temp, parent.frame())
    else        m <- eval(temp, sys.parent())

    Terms <- attr(m, 'terms')
    Y <- model.extract(m, "response")
    if (!inherits(Y, "Surv")) stop("Response must be a survival object")
    
    strats <- attr(Terms, "specials")$strata
    cluster<- attr(Terms, "specials")$cluster
    dropx <- NULL
    if (length(cluster)) {
      if (missing(robust)) robust <- TRUE
      tempc <- survival::untangle.specials(Terms, 'cluster', 1:10)
      ord <- attr(Terms, 'order')[tempc$terms]
      if (any(ord>1)) stop ("Cluster can not be used in an interaction")
      cluster <- strata(m[,tempc$vars], shortlabel=TRUE)  #allow multiples
      dropx <- tempc$terms
    }
    if (length(strats)) {
      temp <- survival::untangle.specials(Terms, 'strata', 1)
      dropx <- c(dropx, temp$terms)
      if (length(temp$vars)==1) strata.keep <- m[[temp$vars]]
      else strata.keep <- strata(m[,temp$vars], shortlabel=TRUE)
      strata <- as.numeric(strata.keep)
      nstrata <- max(strata)
    }else{
      nstrata <- 1
      strata <- 0
    }
    
    if (length(dropx)) {
      newTerms <- Terms[-dropx]
      # R (version 2.7.1) adds intercept=T anytime you drop something
      if (is.R()) attr(newTerms, 'intercept') <- attr(Terms, 'intercept')
    } else  newTerms <- Terms
    
    X <- model.matrix(newTerms, m);
    if (is.R()) {
      assign <- lapply(survival::attrassign(X, newTerms)[-1], function(x) x-1)
      xlevels <- .getXlevels(newTerms, m)
      contr.save <- attr(X, 'contrasts')
    }else {
      assign <- lapply(attr(X, 'assign')[-1], function(x) x -1)
      xvars <- as.character(attr(newTerms, 'variables'))
      xvars <- xvars[-attr(newTerms, 'response')]
      if (length(xvars) >0) {
        xlevels <- lapply(m[xvars], levels)
        xlevels <- xlevels[!unlist(lapply(xlevels, is.null))]
        if(length(xlevels) == 0)
          xlevels <- NULL
      } else xlevels <- NULL
      contr.save <- attr(X, 'contrasts')
    }
    
    # drop the intercept after the fact, and also drop strata if necessary
    adrop <- 0  #levels of "assign" to be dropped; 0= intercept
    Xatt <- attributes(X) 
    xdrop <- Xatt$assign %in% adrop  #columns to drop (always the intercept)
    X <- X[, !xdrop, drop=FALSE]
    attr(X, "assign") <- Xatt$assign[!xdrop]
    n <- nrow(X)
    p <- ncol(X)
    #X.scaled <- scale(X, center=rep(0,p), scale=rep(1,p));
    X.scaled <- scale(X);
    X.center = attributes(X.scaled)$`scaled:center`;
    X.scale = attributes(X.scaled)$`scaled:scale`;
    
    #########################################################################################
    # data structure
    #########################################################################################
    t1 = Y[,1]; t2 = Y[,1];
    type <- attr(Y, "type")
    if (type== 'counting') stop ("Invalid survival type")
    exactsurv <- Y[,ncol(Y)] ==1
    if (any(exactsurv)) {
      t1[exactsurv]=Y[exactsurv,1];
      t2[exactsurv]=Y[exactsurv,1];
    }
    if (type=='interval') {
      intsurv <- Y[,3]==3;
      if (any(intsurv)){
        t1[intsurv]=Y[intsurv,1];
        t2[intsurv]=Y[intsurv,2];
      }
    } 
    delta = Y[,ncol(Y)];
    if (!all(is.finite(Y))) {
      stop("Invalid survival times for this distribution")
    } else {
      if (type=='left') delta <- 2- delta;
    }
    
    #########################################################################################
    # check frailty
    #########################################################################################
    if(!is.null(frailty)) {
      if(is.null(ID)) stop("please specify ID");
      orderindex = order(ID); 
      if(!(sum(orderindex==(1:n))==n)) stop("please sort the data by ID");
      blocki = c(0, cumsum(as.vector(table(ID))));
      if(frailty=="CAR") {
        if(is.null(Proximity)) stop("please specify prxoimity matrix");
        W = Proximity;
        D = rowSums(W);
        if (any(D==0)) stop("it seems that some region does not have any neighbers, which is not allowed, pelase check")
      }else{
        W = matrix(0, length(blocki)-1, length(blocki)-1);
        D = rowSums(W);
      }
      
      #########################################################################################
      # initial MLE analysis and mcmc parameters
      #########################################################################################
      ## initial fit
      if(survmodel=="AFT"){
        fit0 <- survival::survreg(formula = Y~X.scaled, dist="loglogistic");
        summary(fit0)
        theta1 = -fit0$coefficients[1];
        theta2 = -log(fit0$scale);
        theta = c(theta1, theta2);
        beta = -fit0$coefficients[-1];
        thbetaShat0 = fit0$var[c(1,p+2,2:(1+p)),c(1,p+2,2:(1+p))];
        Dtrans = diag(c(-1, -1, rep(-1,p)));
        thbetaShat = t(Dtrans)%*%thbetaShat0%*%Dtrans;
        ## prior for theta
        Vhat = thbetaShat[(1:2),(1:2)];
        multVhat0 = 1*Vhat; multVhat1 = 0.5*Vhat;
        ## prior for beta
        Shat = thbetaShat[-(1:2),-(1:2)];
        multShat0 = as.matrix(30*Shat); multShat1 = as.matrix(0.5*Shat);
      }else if(survmodel=="PO"){
        fit0 <- survival::survreg(formula = Y~X.scaled, dist="loglogistic");
        summary(fit0)
        theta1 = -fit0$coefficients[1];
        theta2 = -log(fit0$scale);
        theta = c(theta1, theta2);
        beta = -fit0$coefficients[-1]/fit0$scale;
        thbetaShat0 = fit0$var[c(1,p+2,2:(1+p)),c(1,p+2,2:(1+p))];
        Dtrans = diag(c(-1, -1, rep(-1/fit0$scale,p))); Dtrans[2,-(1:2)]=fit0$coefficients[-1]/fit0$scale;
        thbetaShat = t(Dtrans)%*%thbetaShat0%*%Dtrans;
        ## prior for theta
        Vhat = thbetaShat[(1:2),(1:2)];
        multVhat0 = 1*Vhat; multVhat1 = 0.5*Vhat;
        ## prior for beta
        Shat = thbetaShat[-(1:2),-(1:2)];
        multShat0 = as.matrix(30*Shat); multShat1 = as.matrix(0.5*Shat);
      }else if(survmodel=="PH"){
        fit0 <- survival::survreg(formula = Y~X.scaled, dist="weibull");
        summary(fit0)
        theta1 = -fit0$coefficients[1];
        theta2 = -log(fit0$scale);
        theta = c(theta1, theta2);
        beta = -fit0$coefficients[-1]/fit0$scale;
        thbetaShat0 = fit0$var[c(1,p+2,2:(1+p)),c(1,p+2,2:(1+p))];
        Dtrans = diag(c(-1, -1, rep(-1/fit0$scale,p))); Dtrans[2,-(1:2)]=fit0$coefficients[-1]/fit0$scale;
        thbetaShat = t(Dtrans)%*%thbetaShat0%*%Dtrans;
        ## prior for theta
        Vhat = thbetaShat[(1:2),(1:2)];
        multVhat0 = 1*Vhat; multVhat1 = 0.5*Vhat;
        ## prior for beta
        Shat = thbetaShat[-(1:2),-(1:2)];
        multShat0 = as.matrix(30*Shat); multShat1 = as.matrix(0.5*Shat);
      }else if(survmodel=="AH"){
        fit0 <- survival::survreg(formula = Y~X.scaled, dist="weibull");
        summary(fit0)
        theta1 = -fit0$coefficients[1];
        theta2 = -log(fit0$scale);
        theta = c(theta1, theta2);
        beta = fit0$coefficients[-1]/(fit0$scale-1);
        thbetaShat0 = fit0$var[c(1,p+2,2:(1+p)),c(1,p+2,2:(1+p))];
        Dtrans = diag(c(-1, -1, rep(1/(fit0$scale-1),p))); 
        Dtrans[2,-(1:2)]=-fit0$coefficients[-1]*fit0$scale/(fit0$scale-1)^2;
        thbetaShat = t(Dtrans)%*%thbetaShat0%*%Dtrans;
        ## prior for theta
        Vhat = thbetaShat[(1:2),(1:2)];
        multVhat0 = 1*Vhat; multVhat1 = 0.5*Vhat;
        ## prior for beta
        Shat = thbetaShat[-(1:2),-(1:2)];
        multShat0 = as.matrix(30*Shat); multShat1 = as.matrix(0.5*Shat);
      }
      
      #########################################################################################
      # priors
      # note the priors should be based on scaled data.
      #########################################################################################
      nburn <- mcmc$nburn;
      nsave <- mcmc$nsave;
      nskip <- mcmc$nskip;
      ndisplay <- mcmc$ndisplay;
      maxL <- prior$maxL; if(is.null(maxL)) maxL<-5;
      Ys = rep(0.5, 2^maxL-1);
      a0=prior$a0; if(is.null(a0)) a0=5;
      b0=prior$b0; if(is.null(b0)) b0=1;
      if(a0<=0){
        a0=-1;cpar=state$cpar; 
        if(is.null(cpar)) stop("please specify state$cpar if prior$a0 is missed");
      }
      theta0 <- prior$theta0; if(is.null(theta0)) theta0 <- theta;
      V0 <- prior$V0; if(is.null(V0)) V0 <- multVhat0;
      beta0 <- prior$beta0; if(is.null(beta0)) beta0 <- beta;
      S0 <- prior$S0; if(is.null(S0)) S0 <- multShat0;
      if(any(V0==Inf)){
        V0inv <- diag(c(1e-200,1e-200));
      }else {
        V0inv <- solve(V0);
      }
      S0inv <- solve(S0);
      if(is.null(state$frail)) {
        v <- rep(0, length(blocki)-1);
      } else {
        v <- state$frail; if(length(v)!=(length(blocki)-1)) stop("check the length of frail");
      }
      taua0 = prior$lambdaa0; if(is.null(taua0)) taua0=1;
      taub0 = prior$lambdab0; if(is.null(taub0)) taub0=1;
      
      mcmc = list(nburn=nburn, nsave=nsave, nskip=nskip, ndisplay=ndisplay)
      prior = list(a0=a0, b0=b0, theta0=theta0, V0=V0, beta0=beta0, S0=S0,
                   lambdaa0=taua0, lambdab0=taub0);

      #########################################################################################
      # current state
      #########################################################################################
      cpar=state$cpar; if(is.null(cpar)) cpar=5;
      lambda = state$lambda; if(is.null(lambda)) lambda=2;
      
      #########################################################################################
      # calling the c++ code and # output
      #########################################################################################
      if(survmodel=="AFT"){
        foo <- .Call("frailtyAFTloglogistic", nburn_=nburn, nsave_=nsave, nskip_=nskip, ndisplay_=ndisplay,
                     t1_=t1, t2_=t2, type_=delta, X_=X.scaled, theta_=theta, beta_=beta, cpar_=cpar, Ys_=Ys, maxL_=maxL,
                     a0_=a0, b0_=b0, theta0_=theta0, V0inv_=V0inv, Vhat_=multVhat1, beta0_=beta0, S0inv_=S0inv, 
                     Shat_=multShat1, l0_=min(5000,nsave/2), adapter_=2.38^2, 
                     v_=v, blocki_=blocki, W_=W, lambda_=lambda, a0lambda_=taua0, b0lambda_=taub0, PACKAGE = "spBayesSurv");
      }else if(survmodel=="PO"){
        foo <- .Call("frailtyPOloglogistic", nburn_=nburn, nsave_=nsave, nskip_=nskip, ndisplay_=ndisplay,
                     t1_=t1, t2_=t2, type_=delta, X_=X.scaled, theta_=theta, beta_=beta, cpar_=cpar, Ys_=Ys, maxL_=maxL,
                     a0_=a0, b0_=b0, theta0_=theta0, V0inv_=V0inv, Vhat_=multVhat1, beta0_=beta0, S0inv_=S0inv, 
                     Shat_=multShat1, l0_=min(5000,nsave/2), adapter_=2.38^2, 
                     v_=v, blocki_=blocki, W_=W, lambda_=lambda, a0lambda_=taua0, b0lambda_=taub0, PACKAGE = "spBayesSurv");
      }else if(survmodel=="PH"){
        foo <- .Call("frailtyPHweibull", nburn_=nburn, nsave_=nsave, nskip_=nskip, ndisplay_=ndisplay,
                     t1_=t1, t2_=t2, type_=delta, X_=X.scaled, theta_=theta, beta_=beta, cpar_=cpar, Ys_=Ys, maxL_=maxL,
                     a0_=a0, b0_=b0, theta0_=theta0, V0inv_=V0inv, Vhat_=multVhat1, beta0_=beta0, S0inv_=S0inv, 
                     Shat_=multShat1, l0_=min(5000,nsave/2), adapter_=2.38^2,
                     v_=v, blocki_=blocki, W_=W, lambda_=lambda, a0lambda_=taua0, b0lambda_=taub0, PACKAGE = "spBayesSurv");
        
      }else if(survmodel=="AH"){
        foo <- .Call("frailtyAHweibull", nburn_=nburn, nsave_=nsave, nskip_=nskip, ndisplay_=ndisplay,
                     t1_=t1, t2_=t2, type_=delta, X_=X.scaled, theta_=theta, beta_=beta, cpar_=cpar, Ys_=Ys, maxL_=maxL,
                     a0_=a0, b0_=b0, theta0_=theta0, V0inv_=V0inv, Vhat_=multVhat1, beta0_=beta0, S0inv_=S0inv, 
                     Shat_=multShat1, l0_=min(5000,nsave/2), adapter_=2.38^2, 
                     v_=v, blocki_=blocki, W_=W, lambda_=lambda, a0lambda_=taua0, b0lambda_=taub0, PACKAGE = "spBayesSurv");
      }
      
      #########################################################################################
      # save state
      #########################################################################################
      #### transfer the estimates back to original scales;
      theta.original = foo$theta;
      beta.original = matrix(foo$beta, p, nsave)/matrix(rep(X.scale, nsave), p, nsave);
      
      if(survmodel=="AFT"){
        model.name <- "Accelerated failure time frailty model:";
        theta.original[1,] = foo$theta[1,] - apply(foo$beta, 
                                                   2, function(x) sum(X.center/X.scale*x));
      }else if(survmodel=="PO"){
        model.name <- "Proportional Odds frailty model:";
        theta.original[1,] = foo$theta[1,] - apply(foo$beta*matrix(exp(-foo$theta[2,]),nrow=p, ncol=nsave, byrow=TRUE), 
                                                   2, function(x) sum(X.center/X.scale*x));
      }else if(survmodel=="PH"){
        model.name <- "Proportional hazards frailty model:";
        theta.original[1,] = foo$theta[1,] - apply(foo$beta*matrix(exp(-foo$theta[2,]),nrow=p, ncol=nsave, byrow=TRUE), 
                                                   2, function(x) sum(X.center/X.scale*x));
      }else if(survmodel=="AH"){
        model.name <- "Accelerated hazards frailty model:";
        theta.original[1,] = foo$theta[1,] - apply(foo$beta*matrix(1-exp(-foo$theta[2,]),nrow=p, ncol=nsave, byrow=TRUE), 
                                                   2, function(x) sum(X.center/X.scale*x));
      }
      
      #### coefficients
      coeff1 <- c(apply(beta.original, 1, mean));
      coeff2 <- c(apply(theta.original, 1, mean));
      coeff <- c(coeff1, coeff2);
      names(coeff) = c(colnames(X.scaled),"theta1", "theta2");
      
      #### Save to a list
      output <- list(modelname=model.name,
                     survmodel = survmodel,
                     coefficients=coeff,
                     call=Call,
                     prior=prior,
                     mcmc=mcmc,
                     n=n,
                     p=p,
                     t1=t1,
                     t2=t2,
                     X.scaled=X.scaled,
                     theta = theta.original,
                     beta = beta.original,
                     theta.scaled = foo$theta,
                     beta.scaled = foo$beta,
                     cpar = foo$cpar,
                     maxL = maxL,
                     Ys = foo$Ys,
                     cpo = foo$cpo,
                     pD = foo$pD, 
                     DIC = foo$DIC,
                     ratetheta = foo$ratetheta,
                     ratebeta = foo$ratebeta,
                     rateYs = foo$rateYs,
                     ratec = foo$ratec,
                     v = foo$v,
                     ratev = foo$ratev,
                     lambda = foo$lambda,
                     IDnames = names(table(ID)));
    }else{ 
      ### start non-frality model  
      #########################################################################################
      # initial MLE analysis and mcmc parameters
      #########################################################################################
      ## initial fit
      if(survmodel=="AFT"){
        fit0 <- survival::survreg(formula = Y~X.scaled, dist="loglogistic");
        summary(fit0)
        theta1 = -fit0$coefficients[1];
        theta2 = -log(fit0$scale);
        theta = c(theta1, theta2);
        beta = -fit0$coefficients[-1];
        thbetaShat0 = fit0$var[c(1,p+2,2:(1+p)),c(1,p+2,2:(1+p))];
        Dtrans = diag(c(-1, -1, rep(-1,p)));
        thbetaShat = t(Dtrans)%*%thbetaShat0%*%Dtrans;
        ## prior for theta
        Vhat = thbetaShat[(1:2),(1:2)];
        multVhat0 = 1*Vhat; multVhat1 = 0.5*Vhat;
        ## prior for beta
        Shat = thbetaShat[-(1:2),-(1:2)];
        multShat0 = as.matrix(30*Shat); multShat1 = as.matrix(0.5*Shat);
      }else if(survmodel=="PO"){
        fit0 <- survival::survreg(formula = Y~X.scaled, dist="loglogistic");
        summary(fit0)
        theta1 = -fit0$coefficients[1];
        theta2 = -log(fit0$scale);
        theta = c(theta1, theta2);
        beta = -fit0$coefficients[-1]/fit0$scale;
        thbetaShat0 = fit0$var[c(1,p+2,2:(1+p)),c(1,p+2,2:(1+p))];
        Dtrans = diag(c(-1, -1, rep(-1/fit0$scale,p))); Dtrans[2,-(1:2)]=fit0$coefficients[-1]/fit0$scale;
        thbetaShat = t(Dtrans)%*%thbetaShat0%*%Dtrans;
        ## prior for theta
        Vhat = thbetaShat[(1:2),(1:2)];
        multVhat0 = 1*Vhat; multVhat1 = 0.5*Vhat;
        ## prior for beta
        Shat = thbetaShat[-(1:2),-(1:2)];
        multShat0 = as.matrix(30*Shat); multShat1 = as.matrix(0.5*Shat);
      }else if(survmodel=="PH"){
        fit0 <- survival::survreg(formula = Y~X.scaled, dist="weibull");
        summary(fit0)
        theta1 = -fit0$coefficients[1];
        theta2 = -log(fit0$scale);
        theta = c(theta1, theta2);
        beta = -fit0$coefficients[-1]/fit0$scale;
        thbetaShat0 = fit0$var[c(1,p+2,2:(1+p)),c(1,p+2,2:(1+p))];
        Dtrans = diag(c(-1, -1, rep(-1/fit0$scale,p))); Dtrans[2,-(1:2)]=fit0$coefficients[-1]/fit0$scale;
        thbetaShat = t(Dtrans)%*%thbetaShat0%*%Dtrans;
        ## prior for theta
        Vhat = thbetaShat[(1:2),(1:2)];
        multVhat0 = 1*Vhat; multVhat1 = 0.5*Vhat;
        ## prior for beta
        Shat = thbetaShat[-(1:2),-(1:2)];
        multShat0 = as.matrix(30*Shat); multShat1 = as.matrix(0.5*Shat);
      }else if(survmodel=="AH"){
        fit0 <- survival::survreg(formula = Y~X.scaled, dist="weibull");
        summary(fit0)
        theta1 = -fit0$coefficients[1];
        theta2 = -log(fit0$scale);
        theta = c(theta1, theta2);
        beta = fit0$coefficients[-1]/(fit0$scale-1);
        thbetaShat0 = fit0$var[c(1,p+2,2:(1+p)),c(1,p+2,2:(1+p))];
        Dtrans = diag(c(-1, -1, rep(1/(fit0$scale-1),p))); 
        Dtrans[2,-(1:2)]=-fit0$coefficients[-1]*fit0$scale/(fit0$scale-1)^2;
        thbetaShat = t(Dtrans)%*%thbetaShat0%*%Dtrans;
        ## prior for theta
        Vhat = thbetaShat[(1:2),(1:2)];
        multVhat0 = 1*Vhat; multVhat1 = 0.5*Vhat;
        ## prior for beta
        Shat = thbetaShat[-(1:2),-(1:2)];
        multShat0 = as.matrix(30*Shat); multShat1 = as.matrix(0.5*Shat);
      }
      
      #########################################################################################
      # priors
      # note the priors should be based on scaled data.
      #########################################################################################
      nburn <- mcmc$nburn;
      nsave <- mcmc$nsave;
      nskip <- mcmc$nskip;
      ndisplay <- mcmc$ndisplay;
      maxL <- prior$maxL; if(is.null(maxL)) maxL<-5;
      Ys = rep(0.5, 2^maxL-1);
      a0=prior$a0; if(is.null(a0)) a0=5;
      b0=prior$b0; if(is.null(b0)) b0=1;
      if(a0<=0){
        a0=-1;cpar=state$cpar; 
        if(is.null(cpar)) stop("please specify state$cpar if prior$a0 is missed");
      }
      theta0 <- prior$theta0; if(is.null(theta0)) theta0 <- theta;
      V0 <- prior$V0; if(is.null(V0)) V0 <- multVhat0;
      beta0 <- prior$beta0; if(is.null(beta0)) beta0 <- beta;
      S0 <- prior$S0; if(is.null(S0)) S0 <- multShat0;
      if(any(V0==Inf)){
        V0inv <- diag(c(1e-200,1e-200));
      }else {
        V0inv <- solve(V0);
      }
      S0inv <- solve(S0);
      
      mcmc = list(nburn=nburn, nsave=nsave, nskip=nskip, ndisplay=ndisplay)
      prior = list(a0=a0, b0=b0, theta0=theta0, V0=V0, beta0=beta0, S0=S0);
      
      #########################################################################################
      # current state
      #########################################################################################
      cpar=state$cpar; if(is.null(cpar)) cpar=5;
      
      #########################################################################################
      # calling the c++ code and # output
      #########################################################################################
      if(survmodel=="AFT"){
        foo <- .Call("AFTloglogistic", nburn_=nburn, nsave_=nsave, nskip_=nskip, ndisplay_=ndisplay,
                     t1_=t1, t2_=t2, type_=delta, X_=X.scaled, theta_=theta, beta_=beta, cpar_=cpar, Ys_=Ys, maxL_=maxL,
                     a0_=a0, b0_=b0, theta0_=theta0, V0inv_=V0inv, Vhat_=multVhat1, beta0_=beta0, S0inv_=S0inv, 
                     Shat_=multShat1, l0_=min(5000,nsave/2), adapter_=2.38^2, PACKAGE = "spBayesSurv");
      }else if(survmodel=="PO"){
        foo <- .Call("POloglogistic", nburn_=nburn, nsave_=nsave, nskip_=nskip, ndisplay_=ndisplay,
                     t1_=t1, t2_=t2, type_=delta, X_=X.scaled, theta_=theta, beta_=beta, cpar_=cpar, Ys_=Ys, maxL_=maxL,
                     a0_=a0, b0_=b0, theta0_=theta0, V0inv_=V0inv, Vhat_=multVhat1, beta0_=beta0, S0inv_=S0inv, 
                     Shat_=multShat1, l0_=min(5000,nsave/2), adapter_=2.38^2, PACKAGE = "spBayesSurv");
      }else if(survmodel=="PH"){
        foo <- .Call("PHweibull", nburn_=nburn, nsave_=nsave, nskip_=nskip, ndisplay_=ndisplay,
                     t1_=t1, t2_=t2, type_=delta, X_=X.scaled, theta_=theta, beta_=beta, cpar_=cpar, Ys_=Ys, maxL_=maxL,
                     a0_=a0, b0_=b0, theta0_=theta0, V0inv_=V0inv, Vhat_=multVhat1, beta0_=beta0, S0inv_=S0inv, 
                     Shat_=multShat1, l0_=min(5000,nsave/2), adapter_=2.38^2, PACKAGE = "spBayesSurv");
        
      }else if(survmodel=="AH"){
        foo <- .Call("AHweibull", nburn_=nburn, nsave_=nsave, nskip_=nskip, ndisplay_=ndisplay,
                     t1_=t1, t2_=t2, type_=delta, X_=X.scaled, theta_=theta, beta_=beta, cpar_=cpar, Ys_=Ys, maxL_=maxL,
                     a0_=a0, b0_=b0, theta0_=theta0, V0inv_=V0inv, Vhat_=multVhat1, beta0_=beta0, S0inv_=S0inv, 
                     Shat_=multShat1, l0_=min(5000,nsave/2), adapter_=2.38^2, PACKAGE = "spBayesSurv");
      }
      
      #########################################################################################
      # save state
      #########################################################################################
      #### transfer the estimates back to original scales;
      theta.original = foo$theta;
      beta.original = matrix(foo$beta, p, nsave)/matrix(rep(X.scale, nsave), p, nsave);
      
      if(survmodel=="AFT"){
        model.name <- "Accelerated failure time model:";
        theta.original[1,] = foo$theta[1,] - apply(foo$beta, 
                                                   2, function(x) sum(X.center/X.scale*x));
      }else if(survmodel=="PO"){
        model.name <- "Proportional Odds model:";
        theta.original[1,] = foo$theta[1,] - apply(foo$beta*matrix(exp(-foo$theta[2,]),nrow=p, ncol=nsave, byrow=TRUE), 
                                                   2, function(x) sum(X.center/X.scale*x));
      }else if(survmodel=="PH"){
        model.name <- "Proportional hazards model:";
        theta.original[1,] = foo$theta[1,] - apply(foo$beta*matrix(exp(-foo$theta[2,]),nrow=p, ncol=nsave, byrow=TRUE), 
                                                   2, function(x) sum(X.center/X.scale*x));
      }else if(survmodel=="AH"){
        model.name <- "Accelerated hazards model:";
        theta.original[1,] = foo$theta[1,] - apply(foo$beta*matrix(1-exp(-foo$theta[2,]),nrow=p, ncol=nsave, byrow=TRUE), 
                                                   2, function(x) sum(X.center/X.scale*x));
      }
      
      #### coefficients
      coeff1 <- c(apply(beta.original, 1, mean));
      coeff2 <- c(apply(theta.original, 1, mean));
      coeff <- c(coeff1, coeff2);
      names(coeff) = c(colnames(X.scaled),"theta1", "theta2");
      
      #### Save to a list
      output <- list(modelname=model.name,
                     survmodel = survmodel,
                     coefficients=coeff,
                     call=Call,
                     prior=prior,
                     mcmc=mcmc,
                     n=n,
                     p=p,
                     t1=t1,
                     t2=t2,
                     X.scaled=X.scaled,
                     theta = theta.original,
                     beta = beta.original,
                     theta.scaled = foo$theta,
                     beta.scaled = foo$beta,
                     cpar = foo$cpar,
                     maxL = maxL,
                     Ys = foo$Ys,
                     cpo = foo$cpo,
                     pD = foo$pD, 
                     DIC = foo$DIC,
                     ratetheta = foo$ratetheta,
                     ratebeta = foo$ratebeta,
                     rateYs = foo$rateYs,
                     ratec = foo$ratec);
    }
    class(output) <- c("survregbayes")
    output
  }

#### print, summary, plot
"print.survregbayes" <- function (x, digits = max(3, getOption("digits") - 3), ...) 
{
  cat(x$modelname,"\nCall:\n", sep = "")
  print(x$call)
  
  cat("\nPosterior Inference of Parameters:\n")
  print.default(format(x$coefficients, digits = digits), print.gap = 2, 
                quote = FALSE)
  
  cat("\nLPML:", sum(log(x$cpo)))
  cat("\nDIC:", x$DIC)
  cat("\nn =",x$n)
  invisible(x)
}

"plot.survregbayes" <- function (x, xpred, tgrid=NULL, CI=0.95, PLOT=FALSE, ...) {
  if(is(x,"survregbayes")){
    if(is.null(xpred)) {
      stop("please specify xpred")
    }else{
      if(is.vector(xpred)) xpred=matrix(xpred, nrow=1);
      if(ncol(xpred)!=x$p) stop("please make sure the number of columns matches!");
    };
    if(is.null(tgrid)) tgrid = seq(0, max(c(x$t1), na.rm=T), length.out=200);
    X.center = attributes(x$X.scaled)$`scaled:center`;
    X.scale = attributes(x$X.scaled)$`scaled:scale`;
    xpred = cbind(xpred);
    nxpred = nrow(xpred);
    for(i in 1:nxpred) xpred[i,] = (xpred[i,]-X.center)/X.scale;
    if(x$survmodel=="AFT"){
      estimates <- .Call("AFTloglogistic_plots", tgrid, xpred, x$theta.scaled, x$beta.scaled, x$Ys, x$maxL, CI, PACKAGE = "spBayesSurv");
    }else if(x$survmodel=="PO"){
      estimates <- .Call("POloglogistic_plots", tgrid, xpred, x$theta.scaled, x$beta.scaled, x$Ys, x$maxL, CI, PACKAGE = "spBayesSurv");
    }else if(x$survmodel=="PH"){
      estimates <- .Call("PHweibull_plots", tgrid, xpred, x$theta.scaled, x$beta.scaled, x$Ys, x$maxL, CI, PACKAGE = "spBayesSurv");
    }else if(x$survmodel=="AH"){
      estimates <- .Call("AHweibull_plots", tgrid, xpred, x$theta.scaled, x$beta.scaled, x$Ys, x$maxL, CI, PACKAGE = "spBayesSurv");
    }
    if(PLOT){
      par(mfrow = c(ceiling(nxpred/2),2));
      for(i in 1:nxpred){
        par(cex=1.5,mar=c(4.1,4.1,1,1),cex.lab=1.4,cex.axis=1.1)
        plot(tgrid, estimates$Shat[,i], "l", lwd=3, xlab="time", ylab="survival", main=paste(i));
        polygon(x=c(rev(tgrid),tgrid),
                y=c(rev(estimates$Shatlow[,i]),estimates$Shatup[,i]),
                border=NA,col="lightgray");
        lines(tgrid, estimates$Shat[,i], lty=3, lwd=3, col=1);
      }
    }
  }
  estimates$tgrid=tgrid;
  invisible(estimates)
}

"summary.survregbayes" <- function(object, CI.level=0.95, ...) {
  ans <- c(object[c("call", "modelname")])
  
  ### CPO
  ans$cpo <- object$cpo
  
  ### Median information
  mat <- as.matrix(object$beta)
  coef.p <- object$coefficients[(1:object$p)];
  coef.m <- apply(mat, 1, median)    
  coef.sd <- apply(mat, 1, sd)
  limm <- apply(mat, 1, function(x) as.vector(coda::HPDinterval(coda::as.mcmc(x), prob=CI.level)))
  coef.l <- limm[1,]
  coef.u <- limm[2,]
  
  coef.table <- cbind(coef.p, coef.m, coef.sd, coef.l , coef.u)
  dimnames(coef.table) <- list(names(coef.p), c("Mean", "Median", "Std. Dev.", 
                                                  paste(CI.level*100, "%HPD-Low", sep=""),
                                                  paste(CI.level*100, "%HPD-Upp", sep="")))
  ans$coeff <- coef.table
  
  ### Baseline Information
  mat <- as.matrix(object$theta)
  coef.p <- object$coefficients[-(1:object$p)];
  coef.m <- apply(mat, 1, median)    
  coef.sd <- apply(mat, 1, sd)
  limm <- apply(mat, 1, function(x) as.vector(coda::HPDinterval(coda::as.mcmc(x), prob=CI.level)))
  coef.l <- limm[1,]
  coef.u <- limm[2,]
  
  coef.table <- cbind(coef.p, coef.m, coef.sd, coef.l , coef.u)
  dimnames(coef.table) <- list(names(coef.p), c("Mean", "Median", "Std. Dev.", 
                                                paste(CI.level*100, "%HPD-Low", sep=""),
                                                paste(CI.level*100, "%HPD-Upp", sep="")))
  ans$basepar <- coef.table
  
  ### Precision parameter
  if(object$prior$a0<=0){
    ans$prec <- NULL
  }else{
    mat <- object$cpar
    coef.p <- mean(mat)    
    coef.m <- median(mat)    
    coef.sd <- sd(mat)
    limm <- as.vector(coda::HPDinterval(coda::mcmc(mat), prob=CI.level))
    coef.l <- limm[1]
    coef.u <- limm[2]
    
    coef.table <- cbind(coef.p, coef.m, coef.sd, coef.l , coef.u)
    dimnames(coef.table) <- list(names(coef.p), c("Mean", "Median", "Std. Dev.", 
                                                  paste(CI.level*100, "%HPD-Low", sep=""),
                                                  paste(CI.level*100, "%HPD-Upp", sep="")))
    ans$prec <- coef.table
  }
  
  ### frailty variance parameter
  ans$n <- object$n
  ans$p <- object$p
  ans$LPML <- sum(log(object$cpo))
  ans$DIC <- object$DIC
  
  ### acceptance rates
  ans$ratetheta = object$ratetheta;
  ans$ratebeta = object$ratebeta;
  ans$rateYs = object$rateYs;
  ans$ratec = object$ratec;
  
  class(ans) <- "summary.survregbayes"
  return(ans)
}


"print.summary.survregbayes"<-function (x, digits = max(3, getOption("digits") - 3), ...) 
{
  cat(x$modelname,"\nCall:\n", sep = "")
  print(x$call)
  
  cat("\nPosterior Inference of Regression Parameters\n")
  cat("(Adaptive M-H acceptance rate: ", x$ratebeta, "):\n", sep="")
  print.default(format(x$coeff, digits = digits), print.gap = 2, 
                quote = FALSE)
  
  cat("\nPosterior Inference of Baseline Parameters\n")
  cat("(Adaptive M-H acceptance rate: ", x$ratetheta, "):\n", sep="")
  print.default(format(x$basepar, digits = digits), print.gap = 2, 
                quote = FALSE)
  cat("(Adaptive M-H acceptance rate for conditional probabilities: ", x$rateYs, ")\n", sep="")
  
  if (!is.null(x$prec)) {
    cat("\nPosterior Inference of Precision Parameter\n")
    cat("(Adaptive M-H acceptance rate: ", x$ratec, "):\n", sep="")
    print.default(format(x$prec, digits = digits), print.gap = 2, 
                  quote = FALSE)
  }
  
  cat("\nLog pseudo marginal likelihood: LPML", x$LPML, sep="=")
  cat("\nDeviance Information Criterion: DIC", x$DIC, sep="=")
  cat("\nNumber of subjects:", x$n, sep="=")     
  invisible(x)
}

