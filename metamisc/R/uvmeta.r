#Multivariate meta-analyse: http://blogs.sas.com/content/iml/2012/10/31/compute-the-log-determinant-of-a-matrix/ (directly take log-determinant)

#TODO: allow data transformations
uvmeta <- function(r, vars, model="random", method="MOM", labels, na.action,
                   pars, verbose=FALSE, ...) 
  UseMethod("uvmeta")

uvmeta.default <- function(r, vars, model="random", method="MOM", labels, na.action, 
                           pars, verbose=FALSE, ...)
{
  calcProfile <- function (mleObj, pars) {
    levels = c(pars$level, 0.50)
    quantiles = c(((1-pars$level)/2),0.50,((1-(1-pars$level)/2)))
    pci = array(NA,dim=c(length(coef(mleObj)),3))
    colnames(pci) = paste(quantiles*100,"%",sep=" ")
    for (i in 1:length(levels)) {
      pcint <- confint(mleObj,level=levels[i], quietly=T)
      
      if (length(coef(mleObj))>1) {
        cols.select <- which(colnames(pcint) %in% colnames(pci))
        pci[,colnames(pcint)[cols.select]] <- pcint[,cols.select]
      } else {
        cols.select <- which(names(pcint) %in% colnames(pci))
        pci[,names(pcint)[cols.select]] <- pcint[cols.select]
      }
            
    }
    return(pci)
  }
  
  pars.default <- list(level = 0.95,
                       hp.mu.mean = 0, 
                       hp.mu.var = 1000,
                       n.chains=4, #JAGS (# chains)
                       n.adapt=5000, #JAGS
                       n.init=5000,  #JAGS burn-in
                       n.iter=10000) #JAGS
  
  if (length(r)!=length(vars)) {
    stop("The vectors 'r' and 'vars' have different lengths!")
  }
  
  ds <- as.data.frame(cbind(as.vector(r),as.vector(vars)))
  colnames(ds) <- c("theta","v")
  
  if (!missing(labels)) {
    if (length(labels) != length(r))
      stop("The vectors 'labels' and 'r' have different lengths!")
    rownames(ds) = labels
  } 
  
  if (missing(na.action)) 
    na.action <- "na.fail"
  if (length(na.action)) 
    ds <- do.call(na.action, list(ds))
  
  if (!missing(pars)) {
    for (i in 1:length(pars)) {
      element <- ls(pars)[i]
      pars.default[[element]] <- pars[[i]]
    }
  }
  
  quantiles <- c((1-pars.default$level)/2, 0.50, (1-((1-pars.default$level)/2)))
  est <- NA 
  
  #############################################################################
  # Start analyses
  #############################################################################
  numstudies = dim(ds)[1]
  dfr = numstudies-1
  
  if(numstudies < 3) {
    warning("There are very few primary studies!")
  }
  
  
  if (method == "MOM") { 
    results = as.data.frame(array(NA,dim=c(4, 5)))
    colnames(results) = c("Estimate","Var",paste(quantiles*100,"%",sep=""))
    rownames(results) = c("mu","tausq","Q","Isq")
    
    # FIXED EFFECTS MODEL
    w = 1/ds$v
    weighted_Tbar = sum(ds$theta*w)/sum(w)
    var_T = 1/sum(w)
    se_T = sqrt(var_T)
    Q = sum(w*(ds$theta-weighted_Tbar)**2)
    results["Q",] = c(Q,NA,rep(NA,length(quantiles)))
  
    # RANDOM EFFECTS MODEL
    # Between-study variance
    if (model=="random" & Q > dfr) {
      re_C =  sum(w) - sum(w**2)/sum(w)
      between_study_var = (Q - dfr)/re_C
    } else {
      between_study_var = 0
    }

    re_v = vars + between_study_var # Within-study plus between-study variance
    re_w = 1/re_v # Updated weights
    re_weighted_Tbar =  sum(ds$theta*re_w)/sum(re_w) # Combined effect
    re_var_T  = 1/sum(re_w)     # Variance of the combined effect   
    re_se_T = sqrt(re_var_T)     # Standard error of combined effect
    
    if (model=="random") {
      results["mu",] = c(re_weighted_Tbar,re_var_T,re_weighted_Tbar+qnorm(quantiles)*sqrt(re_var_T))
      results["tausq",] = c(between_study_var,NA,rep(NA,3))
      
      # Calculate I2 and its confidence limits
      Isq <- (results["Q",]-dfr)/results["Q",]
      Isq[which(Isq>1)] <- 1
      Isq[which(Isq<0)] <- 0
      results["Isq",] = Isq
    } else if (model=="fixed") {
      results["mu",] = c(weighted_Tbar,var_T,weighted_Tbar+qnorm(quantiles)*sqrt(var_T))
      results["tausq",] = c(0,0,rep(NA,3))
    }
    pred.int <- results["mu","Estimate"] + qt(quantiles,df=(numstudies-2))*sqrt(results["tausq","Estimate"]+results["mu","Var"])
    names(pred.int) <- paste(quantiles*100,"%",sep="")
    
    est <- list(results=results, model=model,df=dfr,numstudies=numstudies, pred.int=pred.int)
    
  } else if (method=="ml" | method=="pl") {
    results = as.data.frame(array(NA,dim=c(4, length(quantiles)+2)))
    colnames(results) = c("Estimate","Var",paste(quantiles*100,"%",sep=""))
    rownames(results) = c("mu","tausq","Q","Isq")
    
    #mle.loglik <- function( theta, tausq, ds) {
    #  loglik <- -0.5*sum(log(2*pi*(ds$v+tausq)))-0.5*sum(((ds$theta-theta)**2)/(ds$v+tausq))
    #  return (-loglik) #return negative log-likelihood
    #}
    mle.loglik.random <- function(theta, tausq, ds) { #random effects
          loglik <- sum(dnorm(x=ds$theta,mean=theta, sd=sqrt(tausq+ds$v),log=T))
          return (-loglik)
    }
    mle.loglik.fixed <- function(theta, ds) { #fixed effects
      loglik <- sum(dnorm(x=ds$theta,mean=theta, sd=sqrt(ds$v),log=T))
      return (-loglik)
    }
    
    # first apply fixed-effects analysis
    mle.fixed <- mle2(minuslogl=mle.loglik.fixed,start=list(theta=0),data=list(ds=ds))
    Q <- sum((ds$theta-coef(mle.fixed)["theta"])**2/ds$v) #use theta of the fixed-effects analysis
    results["Q",] = c(Q,NA,rep(NA,length(quantiles)))
    
    if (model=="random") {
      mle.random <- mle2(minuslogl=mle.loglik.random,start=list(theta=0, tausq=0),data=list(ds=ds),method="L-BFGS-B",lower=list(theta=-Inf,tausq=0))
      
      if (method=="pl") {
        profile = calcProfile(mle.random, pars.default)  #Use profile likelihood confidence intervals
        results["mu",] = c(coef(mle.random)["theta"],diag(vcov(mle.random))["theta"],profile[1,])
        results["tausq",] = c(coef(mle.random)["tausq"],diag(vcov(mle.random))["tausq"],profile[2,])
      } else {
        results["mu",] = c(coef(mle.random)["theta"],diag(vcov(mle.random))["theta"],coef(mle.random)["theta"]+qnorm(quantiles)*sqrt(diag(vcov(mle.random))["theta"]))
        results["tausq",] = c(coef(mle.random)["tausq"],diag(vcov(mle.random))["tausq"],rep(NA,length(quantiles)))
      }
      
      # Calculate I2 
      Isq <- (results["Q",]-dfr)/results["Q",]
      Isq[which(Isq>1)] <- 1
      Isq[which(Isq<0)] <- 0
      results["Isq",] = Isq
      loglik = -attr(mle.random,"min")
    } else {
      profile = calcProfile(mle.fixed,pars.default)  #Use profile likelihood confidence intervals
      results["mu",] = c(coef(mle.fixed)["theta"],diag(vcov(mle.fixed))["theta"],profile)
      results["tausq",] = c(0,0,rep(0,length(quantiles)))
      loglik = -attr(mle.fixed,"min")
    }
    pred.int <- results["mu","Estimate"] + qt(quantiles,df=(numstudies-2))*sqrt(results["tausq","Estimate"]+results["mu","Var"])
    names(pred.int) <- paste(quantiles*100,"%",sep="")    
    est <- list(results=results, model=model,df=dfr,numstudies=numstudies, pred.int=pred.int, loglik=loglik)
  } else if (method == "reml") {
    
  }
  else if (method == "bayes") { 
    #require("coda")
    #require("rjags")
    
    quiet = !verbose
    
    #Start with fixed effects model to calculate Q and I square statistic
    modelfile <-  system.file(package="metamisc", "model", "uvmeta_fixef.bug")
   
    jags <- jags.model(modelfile,
                       data = list('r' = ds$theta,
                                   'vars' = ds$v,
                                   'hp.mu.mean' = pars.default$hp.mu.mean,
                                   'hp.mu.prec' = 1/pars.default$hp.mu.var,
                                   'k' = numstudies), #prior precision matrix
                       n.chains = pars.default$n.chains,
                       n.adapt = pars.default$n.adapt,
                       quiet = quiet)
    update(jags, pars.default$n.init) #initialize burn-in
    samples <- coda.samples(jags, c('mu','Q','Isq','theta.new'),n.iter=pars.default$n.iter)
    results <- summary(samples,quantiles=quantiles) 
    
    results.overview = as.data.frame(array(NA,dim=c(dim(results[[1]])[1], length(quantiles)+2)))
    colnames(results.overview) = c("Estimate","Var",paste(quantiles*100,"%",sep=""))
    rownames(results.overview) = rownames(results[[2]])
    results.overview[,1] = (results[[1]])[,"Mean"]
    results.overview[,2] = (results[[1]])[,"SD"]**2
    for (i in 1:length(quantiles)) {
      results.overview[,(i+2)] = (results[[2]])[,i]
    }
    results.overview["tausq",] = 0
    results.overview = results.overview[c("mu","tausq","Q","Isq"),]
    
    if (model=="random") {
      modelfile <- system.file(package="metamisc", "model", "uvmeta_ranef.bug")
      jags <- jags.model(modelfile,
                         data = list('r' = ds$theta,
                                     'vars' = ds$v,
                                     'k' = numstudies,
                                     'hp.mu.mean' = pars.default$hp.mu.mean,
                                     'hp.mu.prec' = 1/pars.default$hp.mu.var), #prior precision matrix
                         n.chains = pars.default$n.chains,
                         n.adapt = pars.default$n.adapt,
                         quiet = quiet)
      update(jags, pars.default$n.init) #initialize
      samples <- coda.samples(jags, c('mu','tausq','theta.new'),n.iter=pars.default$n.iter)
      
      results <- summary(samples,quantiles=quantiles) 
      
      #Update 'mu' and 'tausq'
      results.overview[c("mu","tausq"),1] = (results[[1]])[c("mu","tausq"),"Mean"]
      results.overview[c("mu","tausq"),2] = (results[[1]])[c("mu","tausq"),"SD"]**2
      
    }
    
    # Calculate prediction interval
    pred.int=(results[[2]])["theta.new",]
    
    # Calculate deviance
    m.deviance <- dic.samples(jags, n.iter=pars.default$n.iter) # Deviance Information Criterion
    pD <- sum(m.deviance$deviance) # deviance information criterion
    popt <- pD + sum(m.deviance$penalty) #penalized expected deviance
    
    
    est <- list(results=results.overview, model=model,df=dfr,numstudies=numstudies,pred.int=pred.int, pD=pD, popt=popt)
  } else {
    stop("Invalid meta-analysis method!")
  }
  attr(est$results,"level") = pars.default$level
  est$data <- ds
  est$na.action <- na.action
  est$method <- method
  est$call <- match.call()
  class(est) <- "uvmeta"
  return(est)
}

plot.uvmeta <- function(x, ...) {
  
  level <- attr(x$results,"level")
  quantiles <- c((1-level)/2, (1-((1-level)/2)))

  ci <- x$data[,"theta"]+t(qnorm(quantiles)*matrix(rep(sqrt(x$data[,"v"]),length(quantiles)),nrow=(length(quantiles)), ncol=dim(x$data)[1],byrow=T))
  
  xlim <- c(min(ci),max(ci))
  ylim <- c(2,(x$numstudies+5))

  par(mfrow=c(1,1), mar=( c(5, 12, 4, 4) + 0.1))
  
  loc = c((x$numstudies+4):3)
    
  lcol = "#EBEBEB"
  plot(-500,-500,xlim=xlim,ylim=ylim,xlab="", ylab="",yaxt="n",...)
  axis(2,at=c((x$numstudies+4):5,3),labels=c(rownames(x$data),"Pooled estimate"),las=1)
  
  abline(v=0.00,lty=2,col=lcol)
  for (i in 1:x$numstudies) {
    yloc = loc[i]
    points(x$data[i,"theta"],yloc,pch=15)
    lines(c(min(ci[i,]),max(ci[i,])),c(yloc,yloc))
    for (j in 1:dim(ci)[2])
      lines(c(ci[i,j],ci[i,j]),c((yloc-0.1),(yloc+0.1)),pch=3)
  }
  
  ci.bounds <- x$results["mu",paste(quantiles*100,"%", sep="")]
  lines(rep(x$results["mu","Estimate"],2), c(2.9,3.1),pch=3)
  lines(c(min(ci.bounds),max(ci.bounds)),c(3,3))
  points(x$results["mu","Estimate"],3,pch=23,bg="white")
  
  box()
}



print.uvmeta <- function(x, ...)
{
  out <- (x$results)
  text.model <- if (x$model=="fixed") "Fixed" else "Random"
  text.method <- if(x$method=="bayes") "credibility" else "confidence"
  cat(paste(text.model,"effects estimates with corresponding", text.method, "intervals:\n\n"))
	print(out)
  if (x$model=="random") {
    cat(paste("\n\nPrediction interval for mu:\n\n"))
    print(x$pred.int)
  }
  if(x$method=="ml" | x$method=="pl") { #display MLE
    cat(paste("\nLog-likelihood: ", round(x$loglik,2),"\n"))
  } else if (x$method=="bayes") {
    cat(paste("\nDeviance information criterion (DIC): ", round(x$pD,3),"\nPenalized expected deviance: ", round(x$popt,3),"\n"))
  }
  
	out
}


summary.uvmeta <- function(object, ...)
{
    cat("Call:\n")
    print(object$call)
    if (object$model=="fixed")  cat(paste("\nFixed effects summary:\t",round(object$results["mu","Estimate"],5))," (SE: ",round(sqrt(object$results["mu","Var"]),5), ")",sep="")
    if (object$model=="random") {
        cat(paste("\nRandom effects summary:\t",round(object$results["mu","Estimate"],5))," (SE: ",round(sqrt(object$results["mu","Var"]),5), ")",sep="")
        cat(paste("\n\nTau squared: \t\t",round(object$results["tausq","Estimate"],5),sep=""))
    }
    Q_p = 1-pchisq(object$results["Q","Estimate"],df=object$df)
    cat(paste("\nCochran's Q statistic: \t",round(object$results["Q","Estimate"],5)," (p-value: ",round(Q_p,5),")",sep=""))
    cat(paste("\nI-square index: \t", round(object$results["Isq","Estimate"]*100,3)," %\n",sep=""))
}



