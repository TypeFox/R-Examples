#Function: ghap.blmm
#License: GPLv3 or later
#Modification date: 13 Apr 2016
#Written by: Yuri Tani Utsunomiya
#Contact: ytutsunomiya@gmail.com
#Descriptiont: Uses Gibbs sampling to calculate the posterior mean of mixed model parameters

ghap.blmm<-function(
  fixed,             #Formula for the fixed part of the model
  random,            #Column name with the random part of the model
  weights=NULL,      #Weights for records
  ordinal=FALSE,     #Type of response variable (FALSE = continuous, TRUE = ordinal)
  env.eff=FALSE,     #Should permanent environmental effects be included
  data,              #Data frame containing model data
  K,                 #Covariance matrix of random effects
  vu=4,              #Prior degrees of freedom for correlated random variance
  vp=4,              #Prior degrees of freedom for permanent environmental effect variance
  ve=4,              #Prior degrees of freedom for residual variance
  R2=0.5,            #Starting guess for variance explained by random effects
  nchain=1,          #Number of independent Markov chains
  nsim=500,          #Number of MCMC simulations
  burnin=0,          #Initial number of MCMC simulations
  thin=1,            #Thinning interval
  ncores=1,          #Number of cores to split independent chains
  verbose=TRUE
){
  
  #Log message
  if(verbose==TRUE){
    if(ordinal == TRUE){
      cat("\nAssuming a threshold/liability model for ordered categorical data.\n")
    }else{
      cat("\nAssuming continuous data.\n")
    }
    cat("Assembling design matrices... ")
  }
  
  #Build design matrices
  if (ordinal == TRUE) {
    mf <- model.frame(fixed, data = data)
    y <- model.response(mf)
    y <- factor(y, ordered = TRUE)
    y.levels <- levels(y)
    y.nlevels <- length(y.levels)
    y.obs <- as.numeric(y)
    y <- rep(0,length(y))
    if (length(y) <= y.nlevels) {
      stop("The number of records cannot be smaller than or equal to the number of categories.\n")
    }
    X <- model.matrix(mf,mf)
    random.fm <- as.formula(paste("~ 0 +",random))
    Z <- model.matrix(random.fm, data = data)
    colnames(Z) <- gsub(pattern = random, replacement = "", x = colnames(Z))
    Z <- Z[,colnames(K)]
    if(is.null(weights) == FALSE){
      stop("Residual weights can only be defined for continuous responses.")
    }else{
      weights <- rep(1,times=length(y))
    }
  }
  if (ordinal == FALSE) {
    mf <- model.frame(fixed, data = data)
    y <- model.response(mf)
    y.levels <- NULL
    y.nlevels <- NULL
    y.obs <- NULL
    X <- model.matrix(mf,mf)
    random.fm <- as.formula(paste("~ 0 +",random))
    Z <- model.matrix(random.fm, data = data)
    colnames(Z) <- gsub(pattern = random, replacement = "", x = colnames(Z))
    Z <- Z[,colnames(K)]
    if(is.null(weights) == FALSE){
      if(length(weights) != length(y)){
        stop("The vector of residual weights must have the same length as the number of records.")
      }
    }else{
      weights <- rep(1,times=length(y))
    }
    w <- sqrt(1/weights)
    y <- w*y
    X <- w*X
    Z <- w*Z
  }
  
  #Check data
  if(identical(colnames(K),colnames(Z)) == FALSE){
    stop("Names declared in random effects and covariance matrix do not match!")
  }
  
  #Decorrelate random effects
  svdK <- svd(K)
  L <- svdK$u %*% diag(sqrt(svdK$d))
  ZL <- Z%*%L
  
  #Log message
  if(verbose==TRUE){
    cat("Done.\n")
    cat(length(y),"records will be fitted to an intercept,",ncol(X)-1,"fixed effects and",ncol(Z),"random effects.\n")
    if(env.eff==TRUE){
      cat("Additionally,",ncol(Z),"permanent environmental effects will be included.\n")
    }
    if(nchain > 1){
      cat("\nGibbs sampling will be based on",nchain,"parallel Markov chains of:\n")
    }else{
      cat("\nGibbs sampling will be based on a single Markov chain of:\n")
    }
    cat(nsim,"Monte Carlo simulations with a thinning interval of",thin,"samples.\n")
    cat(burnin,"samples will be collected and discarded before the main simulation.\n")
  }
  
  #Gibbs sampler
  gibbs.FUN <- function(mychain,y,y.obs,y.levels,y.nlevels,Xmat,Z,L,ZL,K,svdK,weights,w,vu,vp,ve,R2,nsim,burnin,thin,verbose,ordinal,env.eff){
    
    #Starting values and priors
    if (ordinal == TRUE) {
      #Starting values
      iter.b <- rep(0,times=ncol(Xmat))
      iter.u <- rep(0,times=ncol(Z))
      iter.vare <- 1
      iter.varu <- (R2*iter.vare)/(1-R2)
      iter.thr <- c(-Inf,qnorm(cumsum(table(y.obs)/length(y.obs))))
      iter.dev <- 0
      sum.y <- rep(0,times=length(y))
      sum.thr <- rep(0,times=y.nlevels+1)
      eta <- as.vector(Xmat%*%iter.b)
      #Prior scale parameters
      Su <- iter.varu*(vu-2)/vu
      vuSu <- vu*Su
      posterior.vu <- vu+ncol(Z)
      if(env.eff==TRUE){
        iter.varu <- 0.5*iter.varu
        iter.varp <- iter.varu
        iter.p <- iter.u
        Sp <- iter.varp*(vp-2)/vp
        vpSp <- vp*Sp
        posterior.vp <- vp+ncol(Z)
      }
    }else{
      #Starting values
      vary <- var(y)
      iter.b <- lm.fit(x=Xmat,y=y)
      iter.b <- as.vector(iter.b$coefficients)
      iter.u <- rep(0,times=ncol(ZL))
      iter.vare <- (1-R2)*vary
      iter.varu <- R2*vary
      iter.dev <- 0
      eta <- as.vector(Xmat%*%iter.b)
      #Prior scale parameters
      Su <- iter.varu*(vu-2)/vu
      Se <- iter.vare*(ve-2)/ve
      vuSu <- vu*Su
      veSe <- ve*Se
      posterior.vu <- vu+ncol(Z)
      posterior.ve <- ve+length(y)
      if(env.eff==TRUE){
        iter.varu <- 0.5*iter.varu
        iter.varp <- iter.varu
        iter.p <- iter.u
        Sp <- iter.varp*(vp-2)/vp
        vpSp <- vp*Sp
        posterior.vp <- vp+ncol(Z)
      }
    }
    
    #Initialize iteration parameters
    i.thin <- 0
    sum.sim <- 0
    sum.b <- rep(0,times=ncol(Xmat))
    sum.u <- rep(0,times=ncol(Z))
    sum.vare <- 0
    sum.varu <- 0
    sum.dev <- 0
    if(env.eff==TRUE){
      sum.p <- sum.u
      sum.varp <- sum.varu
    }
    
    #Log message
    if(verbose==TRUE){
      cat("\nExecuting main algorithm...\n")
    }
    
    #Gibbs sampler
    for(i in 1:(nsim+burnin)){
      
      #Increment thinning interval index
      if(i > burnin){
        i.thin <- i.thin + 1
      }
      
      #Sampling for latent variable
      if (ordinal == TRUE) {
        rmin <- pnorm(iter.thr[y.obs],mean = eta,sd = 1)
        rmax <- pnorm(iter.thr[y.obs+1],mean = eta,sd = 1)
        p <- runif(n=length(y.obs), min=rmin, max=rmax)
        z <- qnorm(p, mean = eta, sd=1)
        z[z == -Inf] <- NA; z[is.na(z)] <- min(z, na.rm=TRUE)
        z[z == Inf] <- NA; z[is.na(z)] <- max(z, na.rm=TRUE)
        y <- z
        for(l in 2:max(y.obs)){
          rmin <- max(y[y.obs == l-1],iter.thr[l-1])
          rmax <- min(y[y.obs == l],iter.thr[l+1])
          iter.thr[l] <- runif(n=1,min = rmin,max = rmax)
        }
      }
      
      #Sampling for fixed effects
      for(j in 1:ncol(Xmat)){
        eta <- eta - Xmat[,j]*iter.b[j]
        xpw <- sum(Xmat[,j]*(y-eta))
        xpx <- sum(Xmat[,j]^2)
        a <- xpw/xpx
        v2 <- iter.vare/xpx
        iter.b[j] <- rnorm(n=1,mean=a,sd=sqrt(v2))
        if(i > burnin & i.thin == thin){
          sum.b[j] <- sum.b[j] + iter.b[j]
        }
        eta <- eta + Xmat[,j]*iter.b[j]
      }
      
      #Sampling for correlated random effects
      for(j in 1:ncol(Z)){
        eta <- eta - ZL[,j]*iter.u[j]
        xpw <- sum(ZL[,j]*(y-eta))
        xpx <- sum(ZL[,j]^2) + (iter.vare/iter.varu)
        a <- xpw/xpx
        v2 <- iter.vare/xpx
        iter.u[j] <- rnorm(n=1,mean=a,sd=sqrt(v2))
        if(i > burnin & i.thin == thin){
          sum.u[j] <- sum.u[j] + iter.u[j]
        }
        eta <- eta + ZL[,j]*iter.u[j]
      }
      
      #Sampling for independent random effects
      if(env.eff==TRUE){
        for(j in 1:ncol(Z)){
          eta <- eta - Z[,j]*iter.p[j]
          xpw <- sum(Z[,j]*(y-eta))
          xpx <- sum(Z[,j]^2) + (iter.vare/iter.varp)
          a <- xpw/xpx
          v2 <- iter.vare/xpx
          iter.p[j] <- rnorm(n=1,mean=a,sd=sqrt(v2))
          if(i > burnin & i.thin == thin){
            sum.p[j] <- sum.p[j] + iter.p[j]
          }
          eta <- eta + Z[,j]*iter.p[j]
        }
      }
      
      #Sampling for variance
      iter.varu <- (sum(iter.u^2) + vuSu)/rchisq(n=1,df=posterior.vu)
      if(env.eff==TRUE){
        iter.varp <- (sum(iter.p^2) + vpSp)/rchisq(n=1,df=posterior.vp)
      }
      res <- (1/w)*(y-eta)
      if (ordinal == FALSE) {
        iter.vare <- (sum(res^2) + veSe)/rchisq(n=1,df=posterior.ve)
      }
      if(i > burnin & i.thin == thin){
        sum.varu <- sum.varu + iter.varu
        if(env.eff==TRUE){
          sum.varp <- sum.varp + iter.varp
        }
        sum.vare <- sum.vare + iter.vare
        if(ordinal == TRUE){
          sum.y <- sum.y + y
          sum.thr <- sum.thr + iter.thr
        }
      }
      
      #dev
      if(i > burnin & i.thin == thin){
        iter.dev <- (-length(y)/2)*log(2*pi*iter.vare)
        iter.dev <- iter.dev + (sum(res^2)/(2*iter.vare))
        iter.dev <- -2*iter.dev
        sum.dev <- sum.dev + iter.dev
      }
      
      #Update
      if(i > burnin & i.thin == thin){
        sum.sim <- sum.sim + 1
        i.thin <- 0
      }
      
      #Log message
      if(verbose==TRUE){
        if(i <= burnin){
          cat("Burn-in: ", i,"\r")
        }else{
          cat("Main simulation: ", i-burnin,"\r")
        }
      }
      
    }
    
    #Output
    results <- NULL
    if(ordinal == TRUE){
      results$sum.sim <- sum.sim
      results$liability <- sum.y/sum.sim; names(results$liability) <- data[,random]
      results$thresholds <- as.vector(sum.thr/sum.sim)
      results$b <- sum.b/sum.sim; names(results$b) <- colnames(Xmat)
      results$u <- as.vector(L%*%(sum.u/sum.sim)); names(results$u) <- colnames(K)
      if(env.eff==TRUE){
        results$p <- sum.p/sum.sim; names(results$p) <- colnames(Z)
      }
      results$varu <- sum.varu/sum.sim
      if(env.eff==TRUE){
        results$varp <- sum.varp/sum.sim
      }
      results$vare <- sum.vare/sum.sim
      if(env.eff==TRUE){
        results$h2 <- results$varu/(results$varu+results$vare+results$varp)
        results$H2 <- (results$varu+results$varp)/(results$varu+results$vare+results$varp)
      }else{
        results$h2 <- results$varu/(results$varu+results$vare)
      }
      results$dev <- sum.dev/sum.sim
    }else{
      results$sum.sim <- sum.sim
      results$b <- sum.b/sum.sim; names(results$b) <- colnames(Xmat)
      results$u <- as.vector(L%*%(sum.u/sum.sim)); names(results$u) <- colnames(K)
      if(env.eff==TRUE){
        results$p <- sum.p/sum.sim; names(results$p) <- colnames(Z)
      }
      results$varu <- sum.varu/sum.sim
      if(env.eff==TRUE){
        results$varp <- sum.varp/sum.sim
      }
      results$vare <- sum.vare/sum.sim
      if(env.eff==TRUE){
        results$h2 <- results$varu/(results$varu+results$vare+results$varp)
        results$H2 <- (results$varu+results$varp)/(results$varu+results$vare+results$varp)
      }else{
        results$h2 <- results$varu/(results$varu+results$vare)
      }
      results$dev <- sum.dev/sum.sim
    }
    
    return(results)
    
  }
  
  if(nchain > 1){
    
    tmp <- mclapply(gibbs.FUN,X=1:nchain,mc.cores=ncores,y=y,y.obs=y.obs,y.levels=y.levels,y.nlevels=y.nlevels,Xmat=X,
                    Z=Z,L=L,ZL=ZL,K=K,svdK=svdK,weights=weights,w=w,vu=vu,vp=vp,ve=ve,R2=R2,
                    nsim=nsim,burnin=burnin,thin=thin,verbose=FALSE,ordinal=ordinal,env.eff=env.eff)
    
    if(verbose==TRUE){
      cat("\nAlgorithm finished.\n")
      cat("Assembling results... ")
    }
    
    if(ordinal == TRUE){
      liability <- tmp[[1]]$liability
      thresholds <- tmp[[1]]$thresholds
      b <- tmp[[1]]$b
      u <- tmp[[1]]$u
      varu <- tmp[[1]]$varu
      vare <- tmp[[1]]$vare
      h2 <- tmp[[1]]$h2
      sum.sim <- tmp[[1]]$sum.sim
      if(env.eff==TRUE){
        p <- tmp[[1]]$p
        varp <- tmp[[1]]$varp
        H2 <- tmp[[1]]$H2
      }
      dev <- tmp[[1]]$dev
      for(i in 2:nchain){
        liability <- rbind(liability,tmp[[i]]$liability)
        thresholds <- rbind(thresholds,tmp[[i]]$thresholds)
        b <- rbind(b,tmp[[i]]$b)
        u <- rbind(u,tmp[[i]]$u)
        varu <- rbind(varu,tmp[[i]]$varu)
        vare <- rbind(vare,tmp[[i]]$vare)
        h2 <- rbind(h2,tmp[[i]]$h2)
        sum.sim <- rbind(sum.sim,tmp[[i]]$sum.sim)
        if(env.eff==TRUE){
          p <- rbind(p,tmp[[i]]$p)
          varp <- rbind(varp,tmp[[i]]$varp)
          H2 <- rbind(H2,tmp[[i]]$H2)
        }
        dev <- rbind(dev,tmp[[i]]$dev)
      }
      results <- NULL
      results$nchain <- nchain
      results$nsim <- nsim
      results$thin <- thin
      results$eff.nsim <- apply(sum.sim,MARGIN=2,sum)
      results$liability <- apply(liability,MARGIN=2,mean)
      results$thresholds <- apply(thresholds,MARGIN=2,mean)
      results$b <- apply(b,MARGIN=2,mean)
      results$u <- apply(u,MARGIN=2,mean)
      if(env.eff==TRUE){
        results$p <- apply(p,MARGIN=2,mean)
      }
      results$varu <- apply(varu,MARGIN=2,mean)
      if(env.eff==TRUE){
        results$varp <- apply(varp,MARING=2,mean)
      }
      results$vare <- apply(vare,MARGIN=2,mean)
      results$h2 <- apply(h2,MARGIN=2,mean)
      if(env.eff==TRUE){
        results$H2 <- apply(H2,MARGIN=2,mean)
      }
      results$k <- as.vector(svdK$v%*%diag(1/svdK$d)%*%t(svdK$u)%*%results$u); names(results$k) <- colnames(K)
      results$y <- y.obs; names(results$y) <- data[,random]
      results$weights <- weights; names(results$weights) <- data[,random]
      if(env.eff==TRUE){
        results$residuals <- as.vector(results$liability - X%*%results$b - Z%*%results$u - Z%*%results$p); names(results$residuals) <- data[,random]
      }else{
        results$residuals <- as.vector(results$liability - X%*%results$b - Z%*%results$u); names(results$residuals) <- data[,random]
      }
      results$dev <- apply(dev,MARGIN=2,mean)
      results$interchain$b.sd <- apply(b,MARGIN=2,sd)
      results$interchain$u.sd <- apply(u,MARGIN=2,sd)
      if(env.eff==TRUE){
        results$interchain$p.sd <- apply(p,MARGIN=2,sd)
      }
      results$interchain$varu.sd <- apply(varu,MARGIN=2,sd)
      if(env.eff==TRUE){
        results$interchain$varp.sd <- apply(varp,MARGIN=2,sd)
      }
      results$interchain$vare.sd <- apply(vare,MARGIN=2,sd)
      results$interchain$h2.sd <- apply(h2,MARGIN=2,sd)
      if(env.eff==TRUE){
        results$interchain$H2.sd <- apply(H2,MARGIN=2,sd)
      }
      results$interchain$dev.sd <- apply(dev,MARGIN=2,sd)
    }else{
      b <- tmp[[1]]$b
      u <- tmp[[1]]$u
      varu <- tmp[[1]]$varu
      vare <- tmp[[1]]$vare
      h2 <- tmp[[1]]$h2
      sum.sim <- tmp[[1]]$sum.sim
      if(env.eff==TRUE){
        p <- tmp[[1]]$p
        varp <- tmp[[1]]$varp
        H2 <- tmp[[1]]$H2
      }
      dev <- tmp[[1]]$dev
      for(i in 2:nchain){
        b <- rbind(b,tmp[[i]]$b)
        u <- rbind(u,tmp[[i]]$u)
        varu <- rbind(varu,tmp[[i]]$varu)
        vare <- rbind(vare,tmp[[i]]$vare)
        h2 <- rbind(h2,tmp[[i]]$h2)
        sum.sim <- rbind(sum.sim,tmp[[i]]$sum.sim)
        if(env.eff==TRUE){
          p <- rbind(p,tmp[[i]]$p)
          varp <- rbind(varp,tmp[[i]]$varp)
          H2 <- rbind(H2,tmp[[i]]$H2)
        }
        dev <- rbind(dev,tmp[[i]]$dev)
      }
      results <- NULL
      results$nchain <- nchain
      results$nsim <- nsim
      results$thin <- thin
      results$eff.nsim <- apply(sum.sim,MARGIN=2,sum)
      results$b <- apply(b,MARGIN=2,mean)
      results$u <- apply(u,MARGIN=2,mean)
      if(env.eff==TRUE){
        results$p <- apply(p,MARGIN=2,mean)
      }
      results$varu <- apply(varu,MARGIN=2,mean)
      if(env.eff==TRUE){
        results$varp <- apply(varp,MARGIN=2,mean)
      }
      results$vare <- apply(vare,MARGIN=2,mean)
      results$h2 <- apply(h2,MARGIN=2,mean)
      if(env.eff==TRUE){
        results$H2 <- apply(H2,MARGIN=2,mean)
      }
      results$k <- as.vector(svdK$v%*%diag(1/svdK$d)%*%t(svdK$u)%*%results$u); names(results$k) <- colnames(K)
      results$y <- (1/w)*y; names(results$y) <- data[,random]
      results$weights <- weights; names(results$weights) <- data[,random]
      if(env.eff==TRUE){
        results$residuals <- (1/w)*as.vector(results$y - X%*%results$b - Z%*%results$u - Z%*%results$p); names(results$residuals) <- data[,random]
      }else{
        results$residuals <- (1/w)*as.vector(results$y - X%*%results$b - Z%*%results$u); names(results$residuals) <- data[,random]
      }
      results$dev <- apply(dev,MARGIN=2,mean)
      results$interchain$b.sd <- apply(b,MARGIN=2,sd)
      results$interchain$u.sd <- apply(u,MARGIN=2,sd)
      if(env.eff==TRUE){
        results$interchain$p.sd <- apply(p,MARGIN=2,sd)
      }
      
      results$interchain$varu.sd <- apply(varu,MARGIN=2,sd)
      if(env.eff==TRUE){
        results$interchain$varp.sd <- apply(varp,MARGIN=2,sd)
      }
      results$interchain$vare.sd <- apply(vare,MARGIN=2,sd)
      results$interchain$h2.sd <- apply(h2,MARGIN=2,sd)
      if(env.eff==TRUE){
        results$interchain$H2.sd <- apply(H2,MARGIN=2,sd)
      }
      results$interchain$dev.sd <- apply(dev,MARGIN=2,sd)
    }
  }else{
    
    tmp <- mclapply(FUN=gibbs.FUN,X=nchain,mc.cores = ncores,
                    y=y,y.obs=y.obs,y.levels=y.levels,y.nlevels=y.nlevels,Xmat=X,
                    Z=Z,L=L,ZL=ZL,K=K,svdK=svdK,weights=weights,w=w,vu=vu,vp=vp,ve=ve,R2=R2,
                    nsim=nsim,burnin=burnin,thin=thin,verbose=verbose,ordinal=ordinal,env.eff=TRUE)
    
    if(verbose==TRUE){
      cat("\nAlgorithm finished.\n")
      cat("Assembling results... ")
    }
    
    if(ordinal == TRUE){
      results <- NULL
      results$nchain <- nchain
      results$nsim <- nsim
      results$thin <- thin
      results$eff.nsim <- tmp[[1]]$sum.sim
      results$liability <- tmp[[1]]$liability
      results$thresholds <- tmp[[1]]$thresholds
      results$b <- tmp[[1]]$b
      results$u <- tmp[[1]]$u
      if(env.eff==TRUE){
        results$p <- tmp[[1]]$p
      }
      results$varu <- tmp[[1]]$varu
      if(env.eff==TRUE){
        results$varp <- tmp[[1]]$varp
      }
      results$vare <- tmp[[1]]$vare
      results$h2 <- tmp[[1]]$h2
      if(env.eff==TRUE){
        results$H2 <- tmp[[1]]$H2
      }
      results$k <- as.vector(svdK$v%*%diag(1/svdK$d)%*%t(svdK$u)%*%results$u); names(results$k) <- colnames(K)
      results$y <- y.obs; names(results$y) <- data[,random]
      results$weights <- weights; names(results$weights) <- data[,random]
      if(env.eff==TRUE){
        results$residuals <- as.vector(results$liability - X%*%results$b - Z%*%results$u - Z%*%results$p); names(results$residuals) <- data[,random]
      }else{
        results$residuals <- as.vector(results$liability - X%*%results$b - Z%*%results$u); names(results$residuals) <- data[,random]
      }
      results$dev <- tmp[[1]]$dev
    }else{
      results <- NULL
      results$nchain <- nchain
      results$nsim <- nsim
      results$thin <- thin
      results$eff.nsim <- tmp[[1]]$sum.sim
      results$b <- tmp[[1]]$b
      results$u <- tmp[[1]]$u
      if(env.eff==TRUE){
        results$p <- tmp[[1]]$p
      }
      results$varu <- tmp[[1]]$varu
      if(env.eff==TRUE){
        results$varp <- tmp[[1]]$varp
      }
      results$vare <- tmp[[1]]$vare
      results$h2 <- tmp[[1]]$h2
      if(env.eff==TRUE){
        results$H2 <- tmp[[1]]$H2
      }
      results$k <- as.vector(svdK$v%*%diag(1/svdK$d)%*%t(svdK$u)%*%results$u); names(results$k) <- colnames(K)
      results$y <- (1/w)*y; names(results$y) <- data[,random]
      results$weights <- weights; names(results$weights) <- data[,random]
      if(env.eff==TRUE){
        results$residuals <- (1/w)*as.vector(results$y - X%*%results$b - Z%*%results$u - Z%*%results$p); names(results$residuals) <- data[,random]
      }else{
        results$residuals <- (1/w)*as.vector(results$y - X%*%results$b - Z%*%results$u); names(results$residuals) <- data[,random]
      }
      results$dev <- tmp[[1]]$dev
    }
  }
  
  #Posterior deviance
  results$pdev <- ((-length(results$y)/2)*log(2*pi*results$vare)) + ((sum(results$residuals^2)/(2*results$vare)))
  results$pdev <- -2*results$pdev
  
  if(verbose==TRUE){
    cat("Done.\n")
    cat("Results are based on an effective number of simulations of ",results$eff.nsim,".\n\n",sep="")
  }
  class(results) <- "GHap.blmm"
  return(results)
  
}