#####################
## Generalized CWM ##
#####################

cwm2 <- function(formulaY, data, Y, Xnorm, Xmult, Xpois, Xbin,colXm,colXn,colXb,colXp,Xbtrials, n, m, Xmod, k, modelXnorm, familyY, 
                    method, initialization, start.z, iter.max, threshold, loglikplot, seed,maxR,eps) 
  # colXm number of factors in data
  # colXn number of Xnorm variables in data
  # m vector with number of levels of factors
{
  t_df <- vY <- NULL       # only for the t distribution  
  nuY <- NULL       # only for the Gamma distribution   
  familyYname <- familyY$family  
  if(!is.null(seed)) set.seed(seed) 
  if(initialization=="manual" |initialization=="mclust") maxR <- 1
  
  
  
  
  # Xnorm Parameters definition ---------------------------------------------------------
  if(colXn>0){
  }  
  # Parameters definition ---------------------------------------------------------
  prior   <- numeric(k) # weights
  
  
  if(!is.null(data)){ 
    if(familyYname=="student.t"){ 
    vY    <- array(1,c(n,k),dimnames=list(1:n,paste("comp.",1:k,sep="")))
    t_df     <- rep(30,k)   # only for the t distribution
    }
    if (familyYname=="binomial" && !is.factor(Y)) mY <- rowSums(Y)
    else mY <- rep(1,nrow(Y))
  }
  result <- list()
  # Do Estimation ----------------------------------------------------------------
  for (r in 1:maxR){
    if (maxR>1) cat(paste("\nrandom init.",r))
    z <- .postInit(initialization=initialization,Y=Y,Xnorm=Xnorm,Xpois=Xpois,Xbin=Xbin,n=n,k=k,start.z=start.z)
    
    # EM algorithm --------------------------------------------------
    # Preliminary definition of convergence criterions
    check     <- 0
    iteration <- 1
    loglik    <- NULL
    aloglik   <- c(0,0)
    a         <- c(0,0)
    
    while(check<1){
      
      # M-step --------------------------------------------------
      # Weights #
      prior <- colMeans(z)
      nj    <- colSums(z) #for multinomial variables
      
      # Posterior probabilities
      norm  <- .PX_norm(colXn=colXn,Xnorm=Xnorm,modelXnorm=modelXnorm,z=z,k=k,n=n,eps=eps)
      bin   <- .PX_bin(k=k,X=Xbin,weights=z,Xbtrials=Xbtrials,n=n)
      pois  <- .PX_Poisson(k=k,X=Xpois,weights=z,n=n) 
      multi <- .PX_multi(colXm=colXm,Xmod=Xmod,z=z,nj=nj,m=m,k=k,n=n)
      
      # Y|x pag 45 of McCullag & Nelder (1989) -----------------------------------------------------
      if (!is.null(data))
        l <- do.call(paste0(".familyY.",familyYname),
                   list(familyY=familyY,k=k,Y=Y,mY=mY,n=n,data=data,z=z,method=method,vY=vY,t_df=t_df,formulaY=formulaY))
      else l <- list(PY=1)
      # Global - Observed-data log-likelihood ## 
      dens              <- matrix(rep(prior,n),n,k,byrow=TRUE)*norm$PX*multi$PX*pois$PX*bin$PX*l$PY
      loglik[iteration] <- sum(log(rowSums(dens)))
        
      # Aitkane's Acceleration-Based Stopping Criterion # #
      
      if(iteration>2 & k > 1){
        if (!is.na(loglik[iteration-1]-loglik[iteration-2]) && !is.na(loglik[iteration])){
          if(abs(loglik[iteration-1]-loglik[iteration-2])>0){
            a[iteration-1]      <- (loglik[iteration]-loglik[iteration-1])/(loglik[iteration-1]-loglik[iteration-2])
            aloglik[iteration]  <- loglik[iteration-1]+(1/(1-a[iteration-1])*(loglik[iteration]-loglik[iteration-1]))
            if(abs(aloglik[iteration]-loglik[iteration])<threshold) 
              check <- 1
          }
          else
            check <- 1    
        }
      }
      if(iteration==iter.max | k==1) check <- 1
      
      cat("*")
      iteration <- iteration + 1
      
      # E-Step ---------------------------------------------------------------------
      dens[dens==0] <- 1E-320
      z      <- dens/matrix(rep(rowSums(dens) ,k),ncol=k)                     # (n x k)
      
      if(!is.null(data)) if(familyYname=="student.t"){
        t_df <- l$t_df
        normzvY <- l$zvY/matrix(rep(colSums(l$zvY),n),n,k,byrow=TRUE) 
        for(j in 1:k)
          vY[,j] <- (t_df[j]+1)/(t_df[j]+(Y-fitted(l$lmodelY[[j]]))^2/l$sig[j])
      }
    }
    
    finalloglik <- loglik[iteration-1] 
    
    # EM-algorithm is finished ---------------------------------------------- 
    
    # Check the number of parameters -------------------------------
    ####################################################################### Sistemare
    
    df_model <- length(coef(l$lmodelY[[1]]))*k
    if(!is.null(data)){
      if(familyYname=="gaussian" | familyYname=="Gamma")  df_model <- df_model + k
      if(familyYname=="student.t")  df_model <- df_model + 2*k
    }
    if (!is.null(modelXnorm)) df_model <- df_model + ncovpar(modelname=strtrim(paste0(modelXnorm,"II"), 3), p=colXn, G=k)    
    
    df <- df_model + (k-1) +colXb*k+colXp*k+ colXn*k + k*(sum(m)-colXm)
    if (df >= n) warning(paste0(df," parameters with only ",n," observations."),call. =FALSE)
    
    
    # Classification Matrix --------------------------------------------------
    cluster <- apply(z,1,which.max)  
    
    # Information criteria --------------------------------------------------
    IC <- list()
    IC$AIC   <- 2*finalloglik - df*2
    IC$BIC  <- 2*finalloglik - df*log(n)
    IC$AIC3  <- 2*finalloglik - df*3  
    IC$AICc  <- IC$AIC - (2*df*(df+1))/(n-df-1)
    IC$AICu  <- IC$AICc - n*log(n/(n-df-1))
    IC$CAIC  <- 2*finalloglik - df*(1+log(n))
    IC$AWE   <- 2*finalloglik - 2*df*(3/2+log(n))  
    z.const    <- (z<10^(-322))*10^(-322)+(z>10^(-322))*z   # vincolo per evitare i NaN nel calcolo di tau*log(tau)
    hard.z     <- (matrix(rep(apply(z,1,max),k),n,k,byrow=F)==z)*1
    ECM     <- sum(hard.z*log(z.const))
    IC$ICL <- IC$BIC+ECM
## GLM parameters 
    GLModel <- NULL
    if(!is.null(data)){
      GLModel <- list()
      #coef <- SEreg(l$lmodelY)    
      for (i in seq_len(k)) {
        parameters  = list(model=l$lmodelY[[i]])
        if (familyYname=="gaussian" |familyYname=="student.t") parameters <- c(parameters, list(sigma = sqrt(l$VarY[1,i])))
        if (familyYname=="student.t") parameters <- c(parameters, list(t_df = t_df[i]))
        if (familyYname=="Gamma") parameters <- c(parameters, list(nuY =l$nuY[i]))
        GLModel[[i]] <- parameters
      }
      names(GLModel) <- paste0("comp.",1:k)
    }
## Concomitant
    concomitant <- NULL
    if (!is.null(Xnorm)){
      concomitant <- list()
       normal <- list(
        d    = norm$PX,
        mu       = norm$mu,
        Sigma    = norm$Sigma,
        model    = modelXnorm
      )
      concomitant <- c(concomitant, normal= normal)
    }
    if (!is.null(Xmult)){
      multinomial <- list(
        d         = multi$PX,
        probs     = multi$alpha 
      )
      concomitant <- c(concomitant, multinomial= multinomial)
    }
    if (!is.null(Xpois)){
        poisson <- list(
        d         = pois$PX,
        lambda    = pois$lambda
      )
      concomitant <- c(concomitant, poisson= poisson)
    }
    if (!is.null(Xbin)){
      binomial <- list(
        d    = bin$PX,
        p    = bin$p
      )
      concomitant <- c(concomitant, binomial= binomial)
    }
    result[[r]] <- list(
        posterior = z,
        iter      = iteration,
        k         = k,
        size      = table(cluster),
        cluster   = cluster,
        logLik    = finalloglik,
        df        = df,      
        prior     = prior,
        IC        = IC,
        converged = TRUE,  #Logical, TRUE if EM algorithm converged.
        GLModel     = GLModel,
        concomitant = concomitant
    )
  }
  result <- result[[which.max(sapply(1:maxR,function(r) result[[r]]$IC$BIC))]]
  class(result) <- "cwm"
  return(result)
}


