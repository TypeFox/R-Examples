#'@title  Monte Carlo simulation for one or two RR variables
#' 
#' @description Simulate and analyse bivariate data including either one RR variable (either correlation, logistic, or linear regression model) or two RR variables (only correlations). Useful for power analysis, parametric bootstraps or for testing the effects of noncompliance on the stability of estimates.
#' @param numRep number of replications
#' @param n sample size
#' @param pi true proportion of carriers of sensitive attribute (for 2 RR variables: \code{vector})
#' @param model either one or two RR model (as \code{vector}), see \code{\link{RRuni}}
#' @param p randomization probability (for 2 RR variables: a \code{list}). See \code{\link{RRuni}}  for details.
#' @param cor true Pearson-correlation used for data generation (for \code{\link{RRcor}}). Can also be used to generate data with two dichotomous RR variables.
#' @param b.log true regression coefficient in logistic regression (for \code{\link{RRlog}})
#   @param b.lin true regression coefficient in linear regression (for \code{\link{RRlin}})
#   @param sig.lin error variance for linear regression model, only used if \code{b.lin} not 0
#' @param complyRates vector with two values giving the proportions of participants who adhere to the instructions in the subset with or without the sensitive attribute, respectively (for 2 RR variables: a \code{list})
#' @param sysBias probability of responding 'yes' (coded as 1 in the RR variable) in case of non-compliance for carriers and non-carriers, respectively. See \code{\link{RRgen}}
#' @param method vector specifying which RR methods to be used in each replication. For a single RR variable, the methods \code{\link{RRuni}}, \code{\link{RRcor}},\code{\link{RRlog}}, and \code{\link{RRlin}} are available. For 2 RR variables, only \code{\link{RRcor}} is available.
#' @param alpha significance threshold for testing the logistic regression parameter \code{beta}
#' @param groupRatio proportion of participants in group 1. Only for two-group models (e.g., \code{"SLD"}) (for 2 RR variables: \code{vector})
#' @param MLest concerns \code{\link{RRuni}}: whether to use \code{optim} to get ML instead of moment estimates (only relevant if pi is outside of [0,1])
#' @param getPower whether to compute power for \code{method="RRcor"} (performs an additional bootstrap assuming independence)
#' @param nCPU integer: how many processors to use? (use 'max' for automatic detection on Windows)
#' @return A list containing
#'  \item{parEsts}{matrix containing the estimated parameters}
#'  \item{results}{matrix  with mean parameters, standard errors, and number of samples to which the respective method could not be fitted}
#'  \item{power}{vector with the estimated power of the selected randomized response procedures}
#' @details  For a single RR variable:
#' 
#' The parameter \code{b.log} is the slope-coefficient for the true, latent values in a logistic regression model that is used for data generation. 
#' 
#' The argument \code{cor} is used for data generation for linear models. The directly measured covariate is sampled from a normal distribution with shifted means, depending on the true state on the sensitive attribute (i.e., the true, underlying values on the RR variable). For dichotomous RR variables, this corresponds to the assumption of an ordinary t-test, where the dependent variable is normally distributed within groups with equal variance. The difference in means is chosen in a way, to obtain the point-biserial correlation defined by \code{cor}.
#' 
#' For two RR variables:
#' 
#' \code{cor} has to be used. In case of two dichotomous RR variables, the true group membership of individuals is sampled from a 2x2 cross table. Within this table, probabilities are chosen in a way, to obtain the point-tetrachoric correlation defined by \code{cor}
#' 
#' Note, that for the FR model with multiple response categories (e.g., from 0 to 4), the specified \code{cor} is not the exact target of the sampling procedure. It assumes a normal distribution for each true state, with constant differences between the groups (i.e., it assumes an interval scaled variable).
#' @examples # Not run: Simulate data according to the Warner model
#' # mcsim <-  RRsimu(numRep=100, n=300, pi=.4, 
#' #                  model="Warner", p=.2, cor=.3)
#' # print(mcsim)
#' @import doParallel
#' @import parallel
#' @import foreach
#' @export
RRsimu <- function(numRep, n, pi, model, p, cor=0, b.log=0, 
                   complyRates=c(1,1), sysBias=c(0,0), 
                   method=c("RRuni","RRcor","RRlog","RRlin"), 
                   alpha=0.05, groupRatio=0.5, MLest=FALSE, 
                   getPower=TRUE, nCPU=1){
  try(stopCluster(cl.tmp), silent = T)
  
  modelNames <- c("Warner","UQTknown","UQTunknown","Mangat","Kuk",
                  "FR","Crosswise","direct","CDM","CDMsym","SLD","custom")
  model <- pmatch(model, modelNames, duplicates.ok=T )
  model <- modelNames[model]
  
  checkzero <- c(cor!=0, b.log!=0)
  if(sum(checkzero) >1){
    if (cor!= 0){
      if (abs(cor) > 1)
        warning("True correlation cor must be in the interval [-1, 1].")
      b.log <- 0
    }
    warning(paste("Only the argument 'cor' is used for data generation (linear model)! Currently used: cor =", cor,"; b.log =", b.log))
  }
  
  nn <- c("pi.true","pi.RRuni","pi2.true","pi2.RRuni","piSE.RRuni",
          "cor.true","cor.RRcor",
          "beta.true","betaSE.true","beta.RRlog","betaSE.RRlog",
          "lincoef.true","lincoefSE.true","lincoef.RRlin","lincoefSE.RRlin",
          "beta.deltaG2.RRlog",
          
          "par2.true","par2.RRuni","par2SE.RRuni")
  nvar <- length(nn)
  
  # check: one or two RR variables
  model <- model[model!="direct"]
  if (length(model)<1 | length(model) > 2)
    stop("Definition of either one or two RR variables required (e.g., model='Warner').")
  twoRR <- length(model)==2
  if (twoRR &&  b.log != 0)
    warning("For two RR variables, only the argument cor is used to simulate data.")
  if (cor != 0){
    if (!twoRR){
      # if (class(complyRates)!= "list") complyRates <- list(complyRates)
      # if only one RR variable included: transform correlation to mean difference
      diffMean <- sqrt(cor^2/(pi[1]*(1-pi[1])*(1-cor^2))) * sign(cor)
      if (model %in% c("custom", "FR") & length(pi)>2){
        ppi <- pi
        if (min(pi)==0){
          ppi <- pi + .001
        }
        diffMean <-  sign(cor)*mean(sqrt(cor^2/(ppi*(1-ppi)*(1-cor^2)))) * sqrt(pi%*% pi)
      }
    }else{
      # two RR variables: think about data generation !
      if (missing(complyRates)) complyRates <- list(c(1,1),c(1,1))
      if (missing(groupRatio)) groupRatio <- c(.5,.5)
      if(length(pi)!=2 ||min(pi)<=0 ||max(pi)>=1) stop("Please provide a vector pi with two values in (0,1)")
    }
  }
  
  ### CDM and CDMsym not available in RRcor
  if (any(model %in% c("CDM","CDMsym")) && any( c("RRlin","RRcor") %in% method)){
    method <- method[!(method == "RRcor")]
    warning("Models 'CDM' and 'CDMsym' are not available in RRcor because both define a third group of cheaters (see vignette('RRreg') ).")
  }
  
  ############################################### ONE RR   ###############################################
  getParEsts.oneRR <- function(replic){
    parEsts <- matrix(NA,nrow=replic,ncol=nvar)
    colnames(parEsts) <- nn
    parEsts <- as.data.frame(parEsts)
    cnt <- 1 ; max.tries <- 3  #; trys <- 0
    #     okcor <- T; oklog <- T; oklin <- T; 
    while(cnt <= replic){
      
      # generate normal data dependent on true value on RR variable
      # k : loop through true states
      if (cor != 0){
        # generate continuous and discrete variable (point-biserial correlation)
        dat <- RRgen(n=n,pi.true=pi,model=model,p=p,complyRates=complyRates,
                     sysBias=sysBias, groupRatio=groupRatio)
        n.true <- table(dat$true)
        dat$cov <- rep(0,n)
        for (k in 0:(length(n.true)-1) ){
          dat[dat$true==k,]$cov <-  rnorm(n.true[k+1],k*diffMean,1)
        }
        # generate from logistic regression model
      }else if (b.log != 0){
        cov <- rnorm(n)
        modelpred <- log(pi/(1-pi)) + b.log*cov
        probsens <- 1/(1+exp(-modelpred))
        true <- rbinom(n, 1, probsens)
        dat <- RRgen(n=n,model=model,p=p,complyRates=complyRates,
                     sysBias=sysBias, groupRatio=groupRatio, trueState=true)
        dat$cov <- cov
      }else{  ## use RRlin with b.lin=0 anyways
        #         true <- rbinom(n, 1, pi)
        #         dat <- RRgen(model=model, p=p, complyRates=complyRates,
        #                         sysBias=sysBias, groupRatio=groupRatio, trueState=true)
        #         dat$cov <- b.lin * true + rnorm(n, sd=sig.lin)
        dat <- RRgen(n=n,pi.true=pi,model=model,p=p,complyRates=complyRates,
                     sysBias=sysBias, groupRatio=groupRatio)
        dat$cov <- rnorm(n)
      }
      
      if (is.null(dat$group)) 
        group <- rep(1,nrow(dat))
      else 
        group <- dat$group
      
      fit.success <- c(RRcor=FALSE, RRlog=FALSE, RRlin=FALSE)
      for(trys in 1:max.tries){
        suppressWarnings(try(rm(log1), silent=TRUE))
        suppressWarnings(try(rm(linmod), silent=TRUE))
        
        # analyse with RRcor
        if (!fit.success[1] && "RRcor" %in% method){
          cor1 <- RRcor(x=dat[,c("response","cov")],
                        models=c(model,"d"),p.list= list(p,1),
                        group=group )  
          if (!is.na(cor1$r["cov","response"])){
            parEsts[cnt,"cor.true"] <- cor(dat$true,dat$cov)
            parEsts[cnt,"cor.RRcor"] <- cor1$r["cov","response"]
            fit.success[1] <- TRUE
          }
        }else{
          fit.success[1] <- TRUE
        }
        
        # analyse with RRuni
        if(trys == 1){
          uni1 <- RRuni(response=dat$response,model=model,p=p,group=group,MLest=MLest)  
          if ("RRuni" %in% method){
            parEsts[cnt,"pi.true"] <- table(dat$true)["1"]/n
            parEsts[cnt,"pi.RRuni"] <- uni1$pi[ifelse(model %in% c("custom","FR"),2,1)]
            parEsts[cnt,"piSE.RRuni"] <- uni1$piSE[ifelse(model %in% c("custom","FR"),2,1)]
            
            if (model %in% c("CDMsym","CDM")){
              parEsts[cnt,"pi.true"] <- table(dat[dat$comply==1,"true"])["1"]/n
              parEsts[cnt,"par2.true"] = 1- colMeans(dat)["comply"]
              parEsts[cnt,"par2.RRuni"] = uni1$gamma
              #         parEsts[cnt,"par2SE.RRuni"] = uni1$gammaSE
            }else if (model== "SLD"){
              parEsts[cnt,"par2.true"] = colMeans(dat[dat$true==1,])["comply"]
              parEsts[cnt,"par2.RRuni"] = uni1$t
              #         parEsts[cnt,"par2SE.RRuni"] = uni1$tSE
            } 
          }
        }
        
        # analyse with RRlog
        if (!fit.success[2] && "RRlog" %in% method){
          #         oklog <- F
          parEsts[cnt,"beta.true"] <- NA
          # adjust dependent variable for CDM: having sensitive attribute and commiting in RR design!
          if (model %in% c("CDMsym","CDM")){
            dat$true <- dat$true & dat$comply
          }
          if(trys == 1)
            try({
              glm1 <- glm(true~cov,dat,family=binomial(link = "logit"))
              parEsts[cnt,"beta.true"] <- glm1$coefficients["cov"]
              parEsts[cnt,"betaSE.true"] <- summary(glm1)$coef[2,2]    
            }, silent=T)
          
          try({log1 <- RRlog(response~cov,data=dat,model=model,p=p,group=group, 
                          LR.test=T, fit.n=1)}, silent=T)
          if (exists("log1", envir= environment()) && 
              !is.null(log1$coefficients["cov"]) && !is.na(log1$coefficients["cov"])){
            if(abs(log1$coefficients["cov"])>3){
              # print(log1)
              try({log1 <- RRlog(response~cov,data=dat,model=model,p=p,group=group, 
                                  LR.test=T, fit.n=3)}, silent=T)
              # print(log1)
            }
            parEsts[cnt,"beta.RRlog"] <- log1$coefficients["cov"]
            suppressWarnings(try(parEsts[cnt,"betaSE.RRlog"] <- sqrt(log1$vcov["cov","cov"]), silent=TRUE))
            try(parEsts[cnt,"beta.deltaG2.RRlog"] <- -2*log1$deltaLogLik["cov"], silent=TRUE)
            fit.success[2] <- TRUE
          }
        }else{
          fit.success[2] <- TRUE
        } 
        
        
        # analyse with RRlin
        if (!fit.success[3] && "RRlin" %in% method){
          #         oklin <- F
          parEsts[cnt,"lincoef.true"] <- NA
          # adjust dependent variable for CDM: having sensitive attribute and commiting in RR design!
          if (model %in% c("CDMsym","CDM")){
            dat$true <- dat$true & dat$comply
          }
          lm1 <- lm(cov~true,dat)
          parEsts[cnt,"lincoef.true"] <- coef(lm1)["true"]
          parEsts[cnt,"lincoefSE.true"] <- sqrt(vcov(lm1)["true","true"])
          
          try(linmod <- RRlin(cov~response,data=dat,models=model,p.list=list(p),
                               group=group, fit.n=1), silent=T)
          if (exists("linmod", envir= environment()) && 
              !is.null(linmod$beta["response"]) && !is.na(linmod$beta["response"])){
            parEsts[cnt,"lincoef.RRlin"] <- linmod$beta["response"]
            parEsts[cnt,"lincoefSE.RRlin"] <- sqrt(vcov(linmod)["response","response"])
            parEsts[cnt,"lincoef.deltaG2.RRlin"] <- -2*logLik(linmod)["response"]
            fit.success[3] <- TRUE
          }
        }else{
          fit.success[3] <- TRUE
        }
        if(all(fit.success))
          break
      }
      cnt <- cnt+1
    }
    parEsts
  }
  
  ############################################### TWO RR   ###############################################
  getParEsts.twoRR <- function(replic){
    parEsts <- matrix(NA,nrow=replic,ncol=nvar)
    colnames(parEsts) <- nn
    parEsts <- as.data.frame(parEsts)
    for (i in 1:(replic)){
      # generate data for two discrete, dichotomous RR variables
      # for a fourfold table, with marginal proportions pi[1] and pi[2]
      # a b
      # c d
      cov <- cor * sqrt(pi[1]*(1-pi[1])*pi[2]*(1-pi[2]))
      a <- pi[1]*pi[2] + cov
      b <- pi[1]*(1-pi[2]) - cov
      c <- (1-pi[1])*pi[2] - cov
      d <- (1-pi[1])*(1-pi[2]) + cov
      if(any(c(a,b,c,d)<0))
        stop("The given correlation coefficient cor=",cor," (phi coefficient) is not observable with n=",n, " and pi=c(",pi[1],", ",pi[2],").")
      dat <- expand.grid(0:1,0:1)[sample(1:4, size=n, replace=TRUE, prob=c(a,b,c,d)),]
      colnames(dat) <- c("true1", "true2")
      RR1 <- RRgen(model =model[1], pi.true=.3, p=p[[1]], complyRates = complyRates[[1]],
                   sysBias = sysBias, groupRatio = groupRatio[1], trueState = dat$true1)
      colnames(RR1) <- paste0(colnames(RR1),1)
      RR2 <- RRgen(model =model[2], pi.true=.3, p=p[[2]], complyRates = complyRates[[2]],
                   sysBias = sysBias, groupRatio = groupRatio[2], trueState = dat$true2)
      colnames(RR2) <- paste0(colnames(RR2),2)
      dat <- cbind(RR1, RR2)
      group <- dat[,grep("^group",colnames(dat))]
      
      # analyse with RRcor and RRlog
      if ("RRcor" %in% method){
        cor1 <- RRcor(x=dat[,c("response1","response2")],
                      models=model,p.list= p,
                      group=group)
        parEsts[i,"cor.true"] <- cor(dat$true1,dat$true2)
        parEsts[i,"cor.RRcor"] <- cor1$r["response1","response2"]
      }
    }
    parEsts
  }
  
  ###############################################
  
  # multi core processing
  if (nCPU!=1){
    #     require(doParallel, quietly=T)
    if (nCPU=="max"){
      try(nCPU <-  as.numeric(Sys.getenv('NUMBER_OF_PROCESSORS')))
      if (nCPU=="max") nCPU <- 2    
    }
    cl.tmp =  makeCluster(nCPU) 
    registerDoParallel(cl.tmp, cores=nCPU)
    #   try(rm(parEsts), silent=TRUE)
    if (!twoRR){
      parEsts <- foreach(k=1:nCPU , .combine= rbind,.packages='RRreg') %dopar% { 
        getParEsts.oneRR(ceiling(numRep/nCPU)) }
    }else{
      parEsts <- foreach(k=1:nCPU , .combine= rbind,.packages='RRreg') %dopar% { 
        getParEsts.twoRR(ceiling(numRep/nCPU)) }
    }
    stopCluster(cl.tmp)
  }else{
    if(!twoRR)
      parEsts <- getParEsts.oneRR(numRep)
    else
      parEsts <- getParEsts.twoRR(numRep)
  }
  
  # significance testing of beta coefficients
  #parEsts[,"beta.perctSignif.RRlog"] <- parEsts[,"beta.prob.RRlog"]  <alpha
  parEsts <- as.matrix(parEsts)
  if ("RRlog" %in% method){
    cutoff <- max(parEsts[,"beta.true"], na.rm=T)
    parEsts[parEsts[,"beta.RRlog"]>max(10, 5*cutoff),"beta.RRlog"] <- NA
  }
  
  # summarize simulation results
  NAs <-  colSums(is.na(parEsts))
  results <- data.frame(Mean=colMeans(parEsts,na.rm =T),
                        #                         SE=bs.se(parEsts,)
                        SD=apply(parEsts,2,sd,na.rm =T)
  )
  # use online estimates between quantiles(.01, .99)
  for(i in 1:nrow(results)){
    if(!is.na(results[i,1])){
      sel <- parEsts[,i] >= quantile(parEsts[,i], probs = .01, na.rm=TRUE) &
        parEsts[,i] <= quantile(parEsts[,i], probs = .99, na.rm=TRUE)
      results[i,1] <- mean(parEsts[sel,i],na.rm =T)
      results[i,2] <- sd(parEsts[sel,i],na.rm =T)
    }
  }
  ####################################### COMPUTE POWER ###########################################
  power <- rep(-1,4) 
  names(power) <-  c("RRuni(z-test)", "RRcor(par_bootstrap)","RRlog(LR-test)","RRlin(Wald-test)")
  if(!twoRR && "RRuni" %in% method){
    z_val <- parEsts[,"pi.RRuni"] / parEsts[,"piSE.RRuni"]
    power[1] <- sum( pnorm(z_val,lower.tail=F) < alpha, na.rm =T)/ (nrow(parEsts)-NAs["pi.RRuni"])
  }
  if ("RRcor" %in% method && getPower){
    simH0 <- RRsimu(numRep=numRep, n=n, pi=pi, model = model, cor=0,
                    p=p, MLest=TRUE, complyRates =complyRates,
                    sysBias = sysBias, groupRatio=groupRatio, 
                    method="RRcor",getPower=F, nCPU=nCPU)
    crit <- quantile(ecdf(abs(simH0$parEsts[,"cor.RRcor"])), 1-alpha)
    power[2] <- sum(abs(parEsts[,"cor.RRcor"]) > crit, na.rm=T)/
      (nrow(parEsts)-NAs["cor.RRcor"])
  }
  if ("RRlog" %in% method){
    prob <- pchisq(parEsts[,"beta.deltaG2.RRlog"], 1, lower.tail=F)
    power[3] <- sum(prob < alpha, na.rm =T)/ (nrow(parEsts)-NAs["beta.deltaG2.RRlog"])
  }
  if (!twoRR &&"RRlin" %in% method){
    z_val <- parEsts[,"lincoef.RRlin"] / parEsts[,"lincoefSE.RRlin"]
    power[4] <- sum( pnorm( abs(z_val), lower.tail=F) < alpha/2, na.rm =T) / (nrow(parEsts)-NAs["lincoef.RRlin"])
  }
  power <- power[power != -1]
  
  
  results <- cbind(results,NAs)
  results <- results[!is.na(results[,1]),]
  parEsts <- parEsts[,NAs!=numRep]
  # filtering for
  # sim$parEsts <- sim$parEsts[sim$parEsts[,"beta.deltaG2.RRlog"]<200,]
  
  sim <- list(parEsts = parEsts, results = results, power = power,
              model=model,p=p,n=n,numRep=numRep, alpha=alpha,
              complyRates=complyRates, sysBias=sysBias, groupRatio=groupRatio,
              method=method)
  class(sim) <- "RRsimu"
  sim
}


# Print results of RRsimu

#' @export
print.RRsimu <- function(x,...){
  if (length(x$model) ==1)
    cat( paste0(x$model," ; n= ",x$n,
                "; randomization probability: ", 
                gsub(", ",", ",toString(round(x$p,3))),"\n" ))
  if (length(x$model) ==2){
    cat( paste0("\n",x$model[1],
                " ; n= ",x$n,
                "; randomization probability: ",
                gsub(", ",", ",toString(round(x$p[[1]],3)))))
    cat( paste0("\n",x$model[2],
                " ; n= ",x$n,
                "; randomization probability: ",
                gsub(", ",", ",toString(round(x$p[[2]],3))),"\n" ))
  }  
  cat(paste0("\nBootstrapped mean and SD of parameters (",x$numRep," replications):\n"))
  print( round(x$results[,1:2],5))
  cat(paste0("\nPower on alpha=",x$alpha," level:\n"))
  print( round(x$power,5))
  sel <- x$results[,3] > x$numRep* .05
  if(any(sel)){
    warning("Fitting routines did fail in \n  ",
            paste(x$results[,3]*100/x$numRep, collapse=","), 
                  "% of bootstrap samples for:\n      ", 
            paste(rownames(x$results)[sel], collapse="; "), ", respectively.", 
            "\n  This might be due to extreme prevalence rates (e.g., pi=.01).")
  }
  
  # if("RRlog" %in% x$method){
  #     warning("Fitting routines did fail in more than 5% of bootstrap samples for:\n      ", 
  #             paste(rownames(x$results)[sel], collapse="; "), 
  #             "\n  This might be due to extreme prevalence rates (e.g., pi=.01).")
  # }
}

# Plot Results of Monte Carlo Simulation
#
# Generic method to plot histograms of the parameter estimates of a Monte Carlo simulation for RR models. Mean parameter estimates are shown in red.
# @param x an \code{\link{RRsimu}} object
# @param ... ignored
# 
#' @export
plot.RRsimu <- function(x,...){
  parEsts <- x$parEsts
  results <- x$results
  nvar <- ncol(parEsts)
  nn <- colnames(parEsts)
  
  # how many plots?
  root <- sqrt(nvar)
  mfrow <- c( floor(root),1+floor(nvar/floor(root))) 
  plot.new()
  par(mfrow=mfrow)
  
  for (i in 1:nvar){
    if(!all(is.na(parEsts[,i]))){
      hist(x= parEsts[,i],freq=FALSE,nclass=max(25, round(x$numRep/25)),
           main=colnames(parEsts)[i],col="grey",
           xlab=paste0("mean=",round(results[i,1],4),",sd=",round(results[i,2],4)))
      abline(v=results[i,1],col="red")
      if (nn[i] %in% c("pi.RRuni","pi.RRlog")){
        dd <- function(x) dnorm(x,mean=results[i,1],sd=results["piSE.RRuni",1])
        curve(dd,col="blue",n=500,add=T)
      }
      else if (nn[i] %in% c("par2.RRuni")){
        dd <- function(x) dnorm(x,mean=results[i,1],sd=results["par2SE.RRuni",1])
        curve(dd,col="blue",n=500,add=T)
      }
      else if(nn[i]=="beta.RRlog"){
        dd <- function(x) dnorm(x,mean=results[i,1],sd=results["betaSE.RRlog",1])
        curve(dd,col="blue",n=500,add=T)
      }
      # asymptotic distribution of beta in RRlog. should be chisq for diffMean=0
      else if (nn[i]=="beta.deltaG2.RRlog"){
        wr <- function(x) dchisq(x,1)
        curve(wr, col = "blue", add = TRUE,n=500)
      }
      #     else if (nn[i]=="beta.prob.RRlog"){
      #       curve(dunif, col = "blue", add = TRUE,n=500)
      #     }
      
    }
  }
  par(mfrow=c(1,1))
}