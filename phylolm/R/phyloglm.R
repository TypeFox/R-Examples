phyloglm <- function(formula, data=list(), phy, method=c("logistic_MPLE","logistic_IG10","poisson_GEE"),
                     btol = 10, log.alpha.bound = 4, start.beta=NULL, start.alpha=NULL, 
                     boot = 0, full.matrix = TRUE)
{
  ### initialize	
  if (!inherits(phy, "phylo")) stop("object \"phy\" is not of class \"phylo\".")	
  if (is.null(phy$edge.length)) stop("the tree has no branch lengths.")
  if (is.null(phy$tip.label)) stop("the tree has no tip labels.")
  method = match.arg(method)
  phy = reorder(phy,"pruningwise")
  ## save original branch lengths for likelihood calculation
  original.edge.length = phy$edge.length	
  n <- length(phy$tip.label)
  N <- dim(phy$edge)[1]
  ROOT <- n + 1L
  anc <- phy$edge[, 1]
  des <- phy$edge[, 2]
  mf = model.frame(formula=formula,data=data)
  if (nrow(mf)!=length(phy$tip.label))
    stop("the number of rows in the data does not match the number of tips in the tree.")
  if (is.null(rownames(mf)))
    warning("the data has no names, order assumed to be the same as tip labels in the tree.\n")
  else {
    ordr = match(phy$tip.label,rownames(mf))
    if (sum(is.na(ordr))>0)
      stop("the row names in the data do not match the tip labels in the tree.\n")
    mf = mf[ordr,,drop=F]
  }
  X = model.matrix(attr(mf, "terms"), data=mf); 
  y = model.response(mf);
  dk = ncol(X) # number of predictors, including intercept, = dimension of beta
  
  dis = pruningwise.distFromRoot(phy) # distance from root to each tip
  
  ## check condition and initialize for Logistic models
  if (method %in% c("logistic_MPLE","logistic_IG10")) {    
    if ( sum(!(y %in% c(0,1))) )
      stop("The model by Ives and Garland requires a binary response (dependent variable).")
    if (var(y)==0) stop("the response (dependent variable) is always 0 or always 1.")  
    
    btouch = 0
    proposedBetaSD = 0.05
    
    ### preparing for generalized tree-structure      
    D = max(dis[1:n]) - dis[1:n]
    D = D - mean(D)
    externalEdge = (des <= n)
    phy$edge.length[externalEdge] <- phy$edge.length[externalEdge] + D[des[externalEdge]]
    # phy is now ultrametric
    times <- pruningwise.branching.times(phy)
    names(times) <- (n+1):(n+phy$Nnode)
    Tmax <- max(times)
    #plot(phy); add.scale.bar(); cat("Tmax=",Tmax,"\ntimes are:",times,"\n")
    intern = which(phy$edge[,2] > n)
    lok = rep(-1,N)
    lok[intern] = des[intern]-n
  } 
  
  ## check condition and initialize for Poisson models
  if (method == "poisson_GEE") {
    if ( (!isTRUE(all(y == floor(y)))) )
      stop("The Poisson regression requires an integer response (dependent variable).")
    if ( sum(y<0) )  
      stop("The Poisson regression requires a positive response (dependent variable).")
  }
  
  ## transform branch lengths fot poisson_GEE model
  transf.branch.lengths_poisson_GEE <- function(beta) {
    if (dk > 1) g = X%*%beta else g = rep(1,n)*beta
    mu = as.vector(exp(g))
    root.edge = 0 
    diag = sqrt(mu/dis[1:n])
    edge.length = phy$edge.length
    return(list(edge.length,root.edge,diag))
  }
    
  ## transform branch lengths fot Logistic models
  transf.branch.lengths <- function(B,lL) {
    if (dk > 1) g = X%*%B else g = rep(1,n)*B
    mu = as.vector(1/(1+exp(-g)))
    p = mean(mu)					
    alpha = 1/exp(lL)
    edge.length = numeric(N)
    distFromRoot <-  exp(-2*alpha*times)
    tmp=.C("transbranchlengths_IvesGarland2010", as.integer(N),as.integer(des),
      as.integer(anc-n), as.integer(lok), as.double(distFromRoot),
      as.integer(externalEdge),as.double(mu), as.double(p), as.double(alpha),
      as.double(D), el=as.double(1:N), di=as.double(1:n))
    edge.length=tmp$el
    diag=tmp$di		
    root.edge = min(distFromRoot)	
    if (any(is.nan(edge.length)))
      stop("edge.length[i] is NaN. Please reduce btol and/or log.alpha.bound.")
    return(list(edge.length,root.edge,diag))	
  }

  ## use three-point structure to do computation
  three.point.compute <- function(trans,y,X) {
    ole=4 + 2*dk + dk*dk
    tmp=.C("threepoint", as.integer(N),as.integer(n),as.integer(phy$Nnode),
      as.integer(1), as.integer(dk), as.integer(ROOT), as.double(trans[[2]]), as.double(trans[[1]]),
      as.integer(des), as.integer(anc), as.double(as.vector(y)), as.double(as.vector(X)),
      result=double(ole))$result  # result=as.double(1:ole)
    # tmp has, in this order:
    # logdetV, 1'V^{-1}1, y'V^{-1}1, y'V^{-1}y, X'V^{-1}1, X'V^{-1}X, X'V^{-1}y
    return(list(vec11=tmp[2], y1=tmp[3], yy=tmp[4],
                X1=tmp[5:(4+dk)], XX=matrix(tmp[(5+dk):(ole-dk)], dk,dk),
                Xy=tmp[(ole-dk+1):ole],logd=tmp[1]))
  }
  
  ## overview of what follows:
  ## B: coefficients
  ## lL = log(1/alpha), alpha: phylogenetic signal
  ## plogregfunct = main optimizing function for IG10 method. It calls:
  ##        plogreglLfunct to optimize alpha given beta (GEE approx)
  ##        plogregBfunct  to optimize the beta coefficients given alpha
  ## BSE = SE for beta parameters, covBSE=variance-covariance matrix for estimated beta's.
  ## npllh = function to calculate negative penalized log likelihood
  
  plogregfunct <- function(startB,startlL) {
    convergeflag = 0
    clL = startlL # current -log(alpha)
    cB  = startB  # current beta coefficients
    diflL = 100   # difference between old and current value
    difB  = 100
    counter = 0   # number of iterations of optimizing beta then alpha
    ttozero = 10^6	
    optss<-list(reltol=.Machine$double.eps^0.5, maxit=100000, parscale=1) # control prm for subplex optimization
		
    while (((diflL>10^-6)|(difB>10^-6)|(ttozero>10^-1))&(counter<20)) {
      counter = counter+1 
      oldlL = clL
      oldB = cB		
      olddiflL = diflL
      olddifB = difB
      #cat("counter:",counter,"starting -log(alpha) and beta:",clL,cB,"\n")
      #if (counter==10) cat("  arg... counter reached 10...\n")

      ## optimize alpha conditional of beta:
      opt <- optim(par = clL, fn = function(par){plogreglLfunct(cB,par)}, method = "L-BFGS-B")
      #opt<-subplex(par=clL, fn = function(par){plogreglLfunct(cB,par)}, control=optss)
      clL = as.numeric(opt$par)
      #cat("  optimized, new -log(alpha):",clL,"\n")
      diflL = (clL-oldlL)^2
      if (counter>=10)
        clL = (clL+oldlL)/2

      ## optimize beta conditional on alpha:
      opt <- optim(par = cB, fn = function(par){plogregBfunct(par,clL)}, method = "L-BFGS-B", control = list(factr=1e12))
      #opt<-subplex(par=cB, fn = function(par){plogregBfunct(par,clL)}, control=optss)
      cB = as.vector(opt$par)
      #cat("  optimized, new beta:",cB,"\n")
      ttozero = as.numeric(opt$value)
      if (ttozero > 10^-2) {
        Btemp = rnorm(dk,startB,proposedBetaSD*pmax(abs(startB),rep(0.1,dk)))
        opt <- optim(par = Btemp, fn = function(par){plogregBfunct(par,clL)}, method = "L-BFGS-B", control = list(factr=1e12))
        #opt<-subplex(par= Btemp, fn = function(par){plogregBfunct(par,clL)}, control=optss)
        Btemp = as.vector(opt$par)
        newttozero = as.numeric(opt$value)
        if (newttozero < ttozero) {
          cB = Btemp
          ttozero = newttozero 
        }
        #cat("    further optimized, new beta:",cB,"\n")
      } 
      difB = sum((cB-oldB)*(cB-oldB))
      if (counter>=10)
        cB = (cB+oldB)/2
    }
    if (counter >= 19) 
      if ((max(abs(c(oldlL - clL,oldB - cB)))>0.1)|(ttozero > 0.5)) convergeflag = 1
    return(list(B=cB,lL=clL,convergeflag = convergeflag))	
  }

  plogregBfunct <- function(B,lL) {
    ## returns L2 norm of penalized score. We want this to be 0
    if (dk > 1) g = X%*%B else g = rep(1,n)*B
    if (any(abs(g) >= btol)) {
      btouch <<- 1
      return(1e6)
    }
    mu = as.vector(1/(1+exp(-g)))
				
    temp  = transf.branch.lengths(B,lL) # edge lengths, root edge, and
    dia  = temp[[3]]   # diagonal terms for tips
    comp = three.point.compute(temp[1:2],(y-mu)/dia,mu*(1-mu)*X/dia)
		
    #compn = three.point.compute.old(temp[1:2],(y-mu)/dia,mu*(1-mu)*X/dia)
    #cat("    (beta) new vs. old 3point: logd, vec11, y1, yy, Xy:\n")
    #print(cbind(comp$logd,  compn$logd))
    #print(cbind(comp$vec11, compn$vec11))
    #print(cbind(comp$y1,    compn$y1))
    #print(cbind(comp$yy,    compn$yy))
    #print(rbind(comp$Xy,    compn$Xy))
    #print(rbind(as.vector(comp$XX),as.vector(compn$XX)))
		
    logdetC = comp$logd + 2*sum(log(dia)) - sum(log(mu*(1-mu)))	
    if (logdetC < -100*log(10)) return(1e6)
		
    Z = comp$Xy
		
    ## Firth correction, i.e. penalty
    if (dk == 1)
      FirthC = (1-2*mu)/2
    else {
      Dx = 0.1
      infoM = comp$XX
      invInfoM = solve(infoM)
      FirthC = rep(NA,dk)
      for (i in 1:dk) {
        ## increase
        dB = B
        dB[i] = dB[i]+Dx
        g = X%*%dB
        if (any(abs(g) >= btol)) return(1e6)
        mu = as.vector(1/(1+exp(-g)))
        ttemp = transf.branch.lengths(dB,lL)
        tdiag = ttemp[[3]]
        tcomp = three.point.compute(ttemp[1:2],(y-mu)/tdiag,mu*(1-mu)*X/tdiag)
        dinfoMp = tcomp$XX
				
        ## decrease
        dB = B
        dB[i] = dB[i]-Dx
        g = X%*%dB
        if (any(abs(g) >= btol)) return(1e6)
        mu = as.vector(1/(1+exp(-g)))
        ttemp = transf.branch.lengths(dB,lL)
        tdiag = ttemp[[3]]
        tcomp = three.point.compute(ttemp[1:2],(y-mu)/tdiag,mu*(1-mu)*X/tdiag)
        dinfoMm = tcomp$XX

        DinfoM = (dinfoMp - dinfoMm)/Dx/2
        FirthC[i] = sum(diag(invInfoM%*%DinfoM))/2
      }
    }
    tozero = Z + FirthC
    return(sum(tozero^2))
  }
  
  plogreglLfunct <- function(B,lL) {
    ## returns sum-of-square-type negative log-likelihood, which we want small:
    ## log|V|/2 + 1/2 (y-mu)' V^{-1} (y-mu).
    g = X%*%B # actually works in case dk=1.
    mu = as.vector(1/(1+exp(-g)))
	
    if (abs(lL - log(Tmax)) >= log.alpha.bound) return(1e10)
    temp = transf.branch.lengths(B,lL)
    dia  = temp[[3]]
    comp = three.point.compute(temp[1:2],(y-mu)/dia ,mu*(1-mu)*X/dia)
    
    LL = (comp$logd + 2*sum(log(dia)) + comp$yy)/2
    if (!is.finite(LL)) LL = 1e10
    return(LL)
  }
	
  plogregBSEfunct <- function(B,lL) { # standard errors for beta coefficients
    g = X%*%B
    mu = as.vector(1/(1+exp(-g)))
    temp = transf.branch.lengths(B,lL)
    dia  = temp[[3]]
    comp = three.point.compute(temp[1:2],(y-mu)/dia, mu*(1-mu)*X/dia )
    infoM = comp$XX # inverse of information matrix
    covBSE = solve(infoM)
    BSE = sqrt(diag(covBSE))
    return(list(BSE = BSE, covBSE = covBSE, info = infoM))
  }

  ## function to calculate the penalized log-likelihood
  npllh <- function(par) {
    if (abs(par[dk+1] - log(Tmax)) >= log.alpha.bound) return(1e10)
    g = X %*% par[1:dk] # beta = first dk components of 'par'ameters
    if (any(abs(g) >= btol)) {
      btouch <<- 1
      return(1e10)
    }
    mu = as.vector(1/(1+exp(-g)))
    temp = transf.branch.lengths(par[1:dk],par[dk+1]) # -log(alpha) = last value in par
    dia = temp[[3]]
    comp = three.point.compute(temp[1:2],numeric(n), mu*(1-mu)*X/dia) # using y=zeros here
    infoM = comp$XX
    llk <- .C("logistreglikelihood", as.integer(N),as.integer(n),as.integer(phy$Nnode),
              as.integer(ROOT), as.double(original.edge.length), as.integer(des), as.integer(anc),
              as.integer(as.vector(y)), as.double(as.vector(mu)),as.integer(dk), as.double(exp(-par[dk+1])),
              loglik=double(1))$loglik
    if (dk==1) pllik = llk + log(abs(infoM))/2
    else       pllik = llk + log(det(infoM))/2
    -pllik
  }
  llh <- function(mu, alpha) { # log-likelihood only (no penalty)
    .C("logistreglikelihood", as.integer(N),as.integer(n),as.integer(phy$Nnode),
       as.integer(ROOT), as.double(original.edge.length), as.integer(des), as.integer(anc),
       as.integer(as.vector(y)), as.double(as.vector(mu)),as.integer(dk), as.double(alpha),
       loglik=double(1))$loglik
  }
	
  ## GEE method for poisson_GEE regression
  ## Iterate: beta_{t+1} = beta_t + I^{-1}[(AX)'R^{-1}(Y-mu)]
  ## I = (AX)'R^{-1}(AX)
  ## R is the correlation matrix
  
  iterate_beta <- function(beta) {
    difbeta = 1
    maxint = 10000
    count = 0
    curbeta = beta
   
    while ((difbeta > 1e-10 )&&(count < maxint)) {
      mu = as.vector(exp(X%*%curbeta))      
      temp = transf.branch.lengths_poisson_GEE(curbeta)      
      dia = temp[[3]]
      #print(dia)
      if (sum(which(mu==0))>0) break
      comp = three.point.compute(temp[1:2],(y-mu)/dia,mu*X/dia)
      invI = solve(comp$XX)
      newbeta = curbeta + invI%*%comp$Xy      
      count = count + 1
      difbeta = sum(abs(newbeta - curbeta))
      curbeta = newbeta              
    }
   
    mu = as.vector(exp(X%*%curbeta))
    r = (y-mu)/sqrt(mu)    
    phi = sum(r^2)/(n-dk)
    covBSE = phi*invI
    BSE = sqrt(diag(covBSE))    
    if (difbeta > 1e-10) convergeflag = 1 else convergeflag = 0       
    return(list(beta = as.vector(curbeta),BSE = BSE,covBSE = covBSE,phi = phi,convergeflag = convergeflag))  
  }
	
  ## Starting values
  if (is.null(start.beta)) {
    if (method %in% c("logistic_MPLE","logistic_IG10")) {    
      fit = glm(y~X-1,family=binomial)
      ## logistf(y~X-1) for regular logistic regression with 
      ## Firth correction. But doesn't work with intercept only.
      startB = fit$coefficients
      if (any(abs(X%*%startB) >= btol)) {
        warning("The estimated coefficients in the absence of phylogenetic signal lead\n  to some linear predictors beyond 'btol'. Increase btol?\n  Starting from beta=0 other than intercept.")
        startB = numeric(dk)
        iint = match( "(Intercept)", colnames(X))
        if (!is.na(iint))
          startB[iint] = log(sum(y==1)/sum(y==0))
        if (any(abs(X%*%startB) >= btol))
          startB[iint] = 0 # all beta's are zero -> linear predictors are all 0.
      }
    }
    if (method == "poisson_GEE") {
      fit = glm(y~X-1,family=poisson)
      start.beta = fit$coefficients
    }
  } else {
      if (length(start.beta)!=dk)
        stop(paste("start.beta shoudl be of length",dk))
      if (method %in% c("logistic_MPLE","logistic_IG10")) {
        startB = as.vector(start.beta)
      if (any(abs(X%*%startB) >= btol))
        stop("With these starting beta values, some linear predictors are beyond 'btol'.\n  Increase btol or choose new starting values for beta.")
      }  
    }
  
  if (method %in% c("logistic_MPLE","logistic_IG10")) {
    if (is.null(start.alpha))
      startlL = log(Tmax) # i.e. alpha = 1/Tmax
    else {
      if (length(start.alpha)!=1) stop("start.alpha should be a single positive value")
      if (start.alpha<=0) stop("start.alpha should be a positive value")
      startlL = -log(start.alpha)
      if (abs(startlL - log(Tmax)) >= log.alpha.bound) {
        tmp = 'start.alpha is outside the bounds, which are\n  exp(+/-log.alpha.bound)/Tmax: '
        tmp = paste(tmp,signif(exp(-log.alpha.bound)/Tmax,3),",",
                    signif(exp(log.alpha.bound)/Tmax,3),' (Tmax=',Tmax,').',
                    '\n  Change start.alpha or increase log.alpha.bound.',sep="")
        stop(tmp)
      }
    }  
  }
  
  
  ### Estimation
  if (method %in% c("logistic_MPLE","logistic_IG10")) {
    if (method == "logistic_IG10"){
      plogreg = plogregfunct(startB,startlL)
      lL = plogreg$lL
      B = plogreg$B
      convergeflag = plogreg$convergeflag
    }
    if (method == "logistic_MPLE") {
      opt <- optim(par=c(startB,startlL), fn=npllh, method="L-BFGS-B", control=list(factr=1e12))
      #optss<-list(reltol=.Machine$double.eps^0.5, maxit=100000, parscale=10)
      #opt<-subplex(par=cB, fn = function(par){robfunct(par)}, control=optss)  	
      B  = opt$par[1:dk]
      lL = opt$par[dk+1]
      convergeflag = opt$convergence
    }
    
    if ((lL - log(Tmax) + 0.02) > log.alpha.bound) {
      warn = paste("the estimate of 'alpha' (",1/exp(lL),
                   ") reached the lower bound (",1/Tmax/exp(log.alpha.bound),
                   ").\n This may reflect a flat likelihood at low alpha values near 0,\n",
                   " meaning that the phylogenetic correlation is estimated to be maximal\n",
                   " under the model in Ives and Garland (2010).", sep="")
      warning(warn)	
    }
    if ((lL - log(Tmax) - 0.02) < - log.alpha.bound) {
      warn = paste("the estimate of 'alpha' (",1/exp(lL),
                   ") reached the upper bound (",exp(log.alpha.bound)/Tmax,
                   ").\n This may simply reflect a flat likelihood at large alpha values,\n",
                   " meaning that the phylogenetic correlation is estimated to be negligible.",sep="")
      warning(warn)	
    }
    if (btouch == 1) 
      warning("the boundary of the linear predictor has been reached during the optimization procedure.
You can increase this bound by increasing 'btol'.")    
    
    plogregBSE = plogregBSEfunct(B,lL)
    
    results <- list(coefficients = B,
                    alpha = 1/exp(lL),
                    sd = plogregBSE$BSE,
                    vcov = plogregBSE$covBSE,
                    convergence = convergeflag
    )
  }
  
  if (method == "poisson_GEE") {
    res = iterate_beta(as.vector(start.beta))        
    results <- list(coefficients = res$beta,
                    scale = res$phi,
                    sd = res$BSE,
                    vcov = res$covBSE,
                    convergence = res$convergeflag
    )    
  }
  
  if (results$converge) warning("phyloglm failed to converge.\n")
  
  names(results$coefficients) = colnames(X)
  colnames(results$vcov) = colnames(X)
  rownames(results$vcov) = colnames(X)
  results$linear.predictors = as.vector(X %*% results$coefficients)
  names(results$linear.predictors) = names(y)
  
  if (method %in% c("logistic_MPLE","logistic_IG10")) {
    if (max(abs(results$linear.predictors)) + 0.01 > btol)
      warning("the linear predictor reaches its bound for one (or more) tip.")  
    results$fitted.values = as.vector(1/(1+exp(-results$linear.predictors)))
    results$mean.tip.height = Tmax
    results$logLik    = llh(results$fitted.values, results$alpha)
    results$penlogLik = results$logLik + log(det(as.matrix(plogregBSE$info)))/2
    results$aic       = -2*results$logLik + 2*(dk+1)
  }  
  
  if (method == "poisson_GEE") {
    results$fitted.values = as.vector(exp(-results$linear.predictors))    
    results$logLik    = NA
    results$penlogLik = NA
    results$aic       = NA
  }    
  
  names(results$fitted.values ) = names(y)
  results$residuals = y - results$fitted.values
  results$y = y
  results$n = n
  results$d = dk  
  results$formula = formula
  results$call = match.call()
  results$method = method
  results$X = X
  results$boot = boot
  
  if ((boot>0)&&(method %in% c("logistic_MPLE","logistic_IG10"))) {
    
    # Turn off warnings
    options(warn=-1)
    
    # simulate all bootstrap data sets
    bootobject <- rbinTrait(n = boot, phy = phy, beta = results$coefficients,
                            alpha = results$alpha, X = X, model = "LogReg")
    responsecolumn <- attr(attr(mf, "terms"), "response")
    
    # analyze these bootstrapped data
    ncoeff = length(results$coefficients)
    bootmatrix <- matrix(NA, boot, ncoeff + 1)
    colnames(bootmatrix) <- c(names(results$coefficients), "alpha")
    
    for (i in 1:boot){
      y = bootobject[,i]
      if (method == "logistic_IG10") {
        bootfit <- try(plogregfunct(startB,startlL), silent=TRUE)
        if (!inherits(bootfit, 'try-error')){
          bootmatrix[i, 1:ncoeff] <- bootfit$B
          bootmatrix[i, ncoeff + 1] <- 1/exp(bootfit$lL)
        }
      }
      
      if (method == "logistic_MPLE") {
        bootfit <- try(optim(par=c(startB,startlL), fn=npllh, method="L-BFGS-B", control=list(factr=1e12)), 
                       silent=TRUE)
        if (!inherits(bootfit, 'try-error')){
          bootmatrix[i, 1:ncoeff] <- bootfit$par[1:dk]
          bootmatrix[i, ncoeff + 1] <- 1/exp(bootfit$par[dk+1])
        }
      }
    }
    
    # summarize bootstrap estimates
    ind.na <- which(is.na(bootmatrix[,1]))
    # indices of replicates that failed: phyloglm had an error
    if (length(ind.na)>0) {
      bootmatrix <- bootmatrix[-ind.na,]
      numOnes <- range(apply(bootobject[,ind.na],2,sum))
    }
    bootmean <- apply(bootmatrix, 2, mean)
    bootsd <- apply(bootmatrix, 2, sd)
    bootconfint95 <- apply(bootmatrix, 2, quantile, probs = c(.025, .975))
    
    bootmeanAlog <- mean(log(bootmatrix[, ncoeff + 1]))
    bootsdAlog <- sd(log(bootmatrix[, ncoeff + 1]))	
    
    results$bootmean = bootmean
    results$bootsd = bootsd
    results$bootconfint95 = bootconfint95
    results$bootmeanAlog = bootmeanAlog
    results$bootsdAlog = bootsdAlog
    results$bootnumFailed = length(ind.na)
    if (full.matrix) results$bootstrap = bootmatrix
    
    ### Turn on warnings
    options(warn=0)
  }
  
  class(results) = "phyloglm"
  results
}


################################################
print.phyloglm <- function(x, digits = max(3, getOption("digits") - 3), ...){
  cat("Call:\n")
  print(x$call)
  if (x$method %in% c("logistic_MPLE","logistic_IG10")) {
    aiclogLik = c(x$aic,x$logLik,x$penlogLik)
    names(aiclogLik) = c("AIC","logLik","Pen.logLik")
    print(aiclogLik, digits = digits)  
  }  
  cat("\nParameter estimate(s) from ")
  if (x$method=="logistic_IG10") cat("GEE approximation:\n")
  if (x$method=="logistic_MPLE") cat("MPLE:\n")
  if (x$method=="poisson_GEE") cat("poisson_GEE:\n")
  if (x$method %in% c("logistic_MPLE","logistic_IG10")) cat("alpha:",x$alpha,"\n")
  cat("\nCoefficients:\n")
  print(x$coefficients)
}
################################################
summary.phyloglm <- function(object, ...) {
  se <- object$sd
  zval <- coef(object) / se
  if (object$boot == 0) 
    TAB <- cbind(Estimate = coef(object), StdErr = se, z.value = zval,
               p.value = 2*pnorm(-abs(zval)))
  else 
    TAB <- cbind(Estimate = coef(object), StdErr = se, z.value = zval,
                 p.value = 2*pnorm(-abs(zval)), 
#                  bootMean = object$bootmean[1:object$d], bootStdErr = object$bootsd[1:object$d],
                 lowerbootCI = object$bootconfint95[1,1:object$d], 
                 upperbootCI = object$bootconfint95[2,1:object$d])
  if (object$method %in% c("logistic_MPLE","logistic_IG10")) {
    res <- list(call=object$call,
                coefficients=TAB,
                residuals = object$residuals,
                alpha=object$alpha,
                aic=object$aic,   	
                logLik=object$logLik,
                penlogLik=object$penlogLik,
                d = object$d,
                method=object$method,
                mean.tip.height=object$mean.tip.height,
                bootNrep = ifelse(object$boot>0, object$boot - object$bootnumFailed, 0)
                )
    if (res$bootNrep>0) {
      res$bootmean = object$bootmean
      res$bootsd = object$bootsd
      res$bootconfint95 = object$bootconfint95
      res$bootmeanAlog <- object$bootmeanAlog
    }
  }
  
  if (object$method == "poisson_GEE") {
    res <- list(call=object$call,
                coefficients=TAB,
                residuals = object$residuals,
                scale=object$scale,                
                method=object$method)  
  }
  
  class(res) = "summary.phyloglm"
  res
}
################################################
print.summary.phyloglm <- function(x, digits = max(3, getOption("digits") - 3), ...){
  cat("\nCall:\n")
  print(x$call)
  if (x$method %in% c("logistic_MPLE","logistic_IG10")) {
    aiclogLik = c(x$aic,x$logLik,x$penlogLik)
    names(aiclogLik) = c("AIC","logLik","Pen.logLik")
    print(aiclogLik, digits = digits)  
  }  
  cat("\nMethod:",x$method)
  if (x$method %in% c("logistic_MPLE","logistic_IG10")) cat("\nMean tip height:",x$mean.tip.height)
  cat("\nParameter estimate(s):\n")
  if (x$method %in% c("logistic_MPLE","logistic_IG10")) {
    cat("alpha:",x$alpha,"\n")
    if (x$bootNrep > 0) {
      cat("      bootstrap mean: ",exp(x$bootmeanAlog)," (on log scale, then back transformed)","\n",sep="")
      cat("      so possible ",ifelse(x$bootmeanAlog>log(x$alpha),"upward","downward")," bias.","\n", sep="")
      cat("      bootstrap 95% CI: (",x$bootconfint95[1,x$d+1],",",x$bootconfint95[2,x$d+1],")\n", sep="")
    }
  }
  if (x$method == "poisson_GEE") cat("scale:",x$scale,"\n")
  cat("\nCoefficients:\n")
  printCoefmat(x$coefficients, P.values=FALSE, has.Pvalue=FALSE)
  if (x$method %in% c("logistic_MPLE","logistic_IG10"))
    cat("\nNote: Wald-type p-values for coefficients, conditional on alpha=",
      x$alpha,"\n",sep="")
  if (x$method == "poisson_GEE") 
    cat("\nNote: Wald-type p-values for coefficients, conditional on scale=",
      x$scale,"\n",sep="")
  if (x$bootNrep > 0)
    cat("      Parametric bootstrap results based on",x$bootNrep,"fitted replicates\n")
  cat("\n")
}
################################################
residuals.phyloglm <-function(object,type=c("response"), ...){
  type <- match.arg(type)
  r <- object$residuals 
  r	 
}
################################################
plot.phyloglm <-function(x, ...){
  plot(fitted(x), x$residuals, xlab = "Fitted value", ylab = "Residuals", ...)
  abline(h=0, lty = 2)
}
################################################
vcov.phyloglm <- function(object, ...){
  vcov = object$vcov
  vcov
}
################################################
logLik.phyloglm <- function(object, ...){
  res = list(logLik = object$logLik, df = object$d+1)
  class(res) = "logLik.phylolm"
  res
}
print.logLik.phyloglm <- function (x, ...) {
  cat("'log Lik.' ",x$logLik," (df=",x$df,")\n", sep = "")
}
AIC.logLik.phyloglm <- function(object, k=2, ...) {
  return(k*object$df - 2*object$logLik)
}
AIC.phyloglm <- function(object, k=2, ...) {
  return(AIC(logLik(object),k))
}
################################################
