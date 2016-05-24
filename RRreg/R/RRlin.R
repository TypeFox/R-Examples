# Randomized response as predictor in linear regression

#' Linear randomized response regression
#' 
#' Linear regression for a continuous criterion, using randomized-response (RR) variables as predictors.
#' @param formula a continuous criterion is predicted by one or more categorical RR variables defined by \code{models}. If the number of predictors exceeds the number defined by the vector \code{models}, the remaining predictors are treated as non-randomized variables (e.g., direct questions). Interactions including any of the RR variables cannot be included.
#' @param data an optional data frame, list or environment, containing the variables in the model. 
#' @param models character vector specifying RR model(s) in order of appearance in formula. Available models: \code{"Warner"}, \code{"UQTknown"}, \code{"UQTunknown"}, \code{"Mangat"}, \code{"Kuk"}, \code{"FR"}, \code{"Crosswise"}, \code{"CDM"}, \code{"CDMsym"}, \code{"SLD"}, \code{"custom"} (custom: a randomization matrix must be specified in the corresponding element of \code{p.list}, where the entry \code{p[i,j]} defines the probability of responding i (i-th row) given a true state of j (j-th column)).
#' @param p.list list of randomization probabilities for RR models in the same order as specified in \code{models}. Note, that the randomization probabilities p must be provided in a \code{\link{list}}, e.g., \code{list(p=c(.2, .3))}. See \code{\link{RRuni}} for details.
# LR.test if true, regression coefficients are tested by a likelihood ratio test by stepwise exclusion of each predictor
#' @param group vector or matrix specifying group membership by the indices 1 and 2. Only for multigroup RR models, e.g., \code{UQTunknown}, \code{CDM} or \code{SLD}
#' @param Kukrep defines the number of repetitions in Kuk's card playing method
#' @param bs.n Number of samples used for the non-parametric bootstrap
#' @param nCPU Number of cores used for the bootstrap
#' @param maxit maximum number of iterations in optimization routine
#' @param fit.n number of fitting runs with random starting values 
#' @param pibeta approximate ratio of probabilities pi to regression weights beta (to adjust scaling). Can be used for speeding-up and fine-tuning ML estimation (i.e., choosing a smaller value for larger beta values).
# @param ... ignored
#' @author Daniel W. Heck
#' @seealso \code{vignette('RRreg')} or \url{https://dl.dropboxusercontent.com/u/21456540/RRreg/index.html} for a detailed description of the RR models and the appropriate definition of \code{p} 
#' @return Returns an object \code{RRlin} which can be analysed by the generic method \code{\link{summary}}
#' @references van den Hout, A., & Kooiman, P. (2006). Estimating the linear regression model with categorical covariates subject to randomized response. \emph{Computational Statistics & Data Analysis, 50}, 3311-3323. 
#' @examples
#' # generate two RR predictors
#' dat <- RRgen(n=500, pi=.4, model="Warner", p=.3)
#' dat2 <- RRgen(n=500, pi=c(.4,.6), model="FR", p=c(.1,.15))
#' dat$FR <- dat2$response
#' dat$trueFR <- dat2$true
#' 
#' # generate a third predictor and continuous dependent variables
#' dat$nonRR <- rnorm(500, 5, 1)
#' dat$depvar <- 2*dat$true - 3*dat2$true + 
#'                        .5*dat$nonRR +rnorm(500, 1, 7) 
#' 
#' # use RRlin and compare to regression on non-RR variables
#' linreg <- RRlin(depvar~response+FR+nonRR, data=dat,
#'                 models=c("Warner","FR"),
#'                 p.list=list(.3, c(.1,.15)), fit.n=1)
#' summary(linreg)
#' summary(lm(depvar~true +trueFR+nonRR, data=dat))
#' @rdname RRlin
#' @import doParallel
#' @import parallel
#' @import foreach
#' @export
RRlin <- function(formula, data, models, p.list, group=NULL, 
                  Kukrep=1, bs.n=0, nCPU=1, maxit=1000, 
                  fit.n = 3, pibeta=0.05) {
  #   UseMethod("RRlin")
  
  # formula interface: from formula to design matrix
  # @export
  # RRlin.formula <- function(formula, group, data, models,  ...){
  mf <- model.frame(formula=formula, data=data,na.action=na.omit)
  x <- model.matrix(attr(mf, "terms"), data=mf)
  
  y <- as.matrix(model.response(mf, "numeric"))
  dimnames(y) <- list(NULL, all.vars(formula)[1])
  
  M <- length(models)  
  intercept <-  ifelse(attr(attr(mf,"terms"),"intercept")==1, T,F)
  if (intercept){
    x <- x[,-1,drop=F]
  }
  w <- x[,c(1:M), drop = F]
  # check, whether DQ/nonRR variables are included in regression
  # U = number of nonRR variables
  if(ncol(x) > M){
    u <- x[, (M+1):ncol(x), drop = F]
    U <- ncol(u)
  }else{
    u <- NULL
    U <- 0
  }
  
  ### check for interactions of RR variables
#   if(length(all.vars(formula))-1 != ncol(x)-ifelse(intercept, 1, 0))
#     warning("No interactions including any RR-variables are allowed (the interaction term will not be corrected for the additional randomness due to the RR procedure). Interactions including only directly observed variables are allowed.")
  
  if (missing(group)){
    group <- matrix(1, length(y), M, dimnames = list(NULL,paste0("g",1:M)))
  }
  if (!missing(data) && class(substitute(group)) == "name") {
    try({data1 <- as.data.frame(data)
         gg <-  as.matrix(eval(substitute(group),data1, parent.frame()),dimnames=list(NULL,group))
         group <- gg
    },silent=T)
    if (class(substitute(group)) == "name"){
      stop(paste0("'group' (",substitute(group),") could not be found in 'data' (",substitute(data),")."))
    }
  }
  if(class(group) != "matrix"){
    group <- matrix(group,ncol=1)
  }
  
  
  # send to default function:
  #   ywu <- list(y,w,u, intercept)
  #   est <- RRlin.default(y=y, w=w, u=u, group=group, data=data, models=models, intercept=intercept, ...)
  #   est <- RRlin.default(formula=ywu, group=group, data=data, models=models,  ...)
  
  
  
  # @export
  # RRlin.default <- function(y, w, u, group, data=NULL, models, p.list, intercept,bs.n=0, nCPU=1, ...){
  # RRlin.default <- function(formula, group, data, models, p.list, bs.n=0, nCPU=1, ...){
  #   
  #   y <- formula[[1]]; w <- formula[[2]] ; u <- formula[[3]]
  #   intercept <- formula[[4]]
  
  # INPUT
  modelNames <- c("Warner","UQTknown","UQTunknown","Mangat","Kuk","FR","Crosswise","CDM","CDMsym","SLD","custom")
  models <- pmatch(models, modelNames, duplicates.ok=T )
  models <- modelNames[models]
  #   M <- length(models)
  RRcheck.lin(y, w, u, models, p.list, group)
  
  ##################
  # get estimates of second parameter for 2-group models
  par2 <- rep(NA,M)
  names(par2) <- models
  for (m in 1: M){
    if (models[m] %in% c("UQTunknown","SLD")){
      sel <- group[,m] != 0
      ytemp <- w[sel,m]
      grouptemp <- group[sel,m]
      uni <- RRuni(response=ytemp,model=models[m],p=p.list[[m]],group=grouptemp,MLest=T)
      if (models[m]=="SLD"){
        par2[m] <- uni$t
      }else if (models[m]=="UQTunknown"){
        par2[m] <- uni$piUQ
      }
      if (par2[m] < 1e-5) par2[m] <- 1e-5
      if (par2[m] > 1-1e-5) par2[m] <- 1-1e-5
    }
  }
  
  ###################
  # construct missclassification matrix/array PW
  # for multiple group models: list/array of PW matrices needed
  #   M <- length(models)
  N <- length(y)   # number observations
  gcombs <- unique(group)  # carefull - can be matrix (>1 RR variables) or vector (1 RR var)
  G <- nrow(gcombs)  # number of group-combinations
  responses <- apply(w,1,paste,collapse=":")   # all responses of persons
  PWarray <- getPWarray(M, G, models, p.list, group, par2, Kukrep=Kukrep)
  
  patterns <- rownames(PWarray[,,1])  #   possible response patterns
  J <- length(patterns)    # number of possible response patterns
  K <- ncol(PWarray[,,1])  # number of possible true states
  
  # group index for correct comparison of PWarray[,,group] and PWpp below
  gidx <- numeric(N)
  for (g in 1:G){
    # select all participants in group combination gcomb[g,]
    sel <- apply(group,1, function(x) all.equal(x,gcombs[g,])) == T
    gidx[sel] <- g
  }
  
  # missclass-probabilities per person (corresponding rows from PWarray, that match response pattern)
  # slow loop
  #   PWpp <- matrix(0,N,K)
  #   for (i in 1:N){
  #     cnt <- cnt + ifelse(PWpp[i,] == PWarray[ match(responses[i],patterns),,gidx[i]],1,0)
  #   }
  # the same vectorized
  PWpp <- t(apply(cbind(responses,gidx),1, 
                  function(x) PWarray[ match(x[1],patterns),,as.numeric(x[2])]))
  # image(PWpp)
  
  #   # design matrix X (intercept in first column)
  if ( !missing(u)){
    x <- cbind(w,u)
  }else{
    x <- w
  }
  if (intercept){
    nam <- colnames(x)
    x <- cbind(1,x)
    colnames(x) <- c("(Intercept)",nam)
  }
  
  
  ####################
  # STARTING VALUES ##
  # standard linear regression
  startmod <- lm(y~x-1)
  beta <- coef(startmod)
  names(beta) <- substr(names(beta),2,1000)
  if (intercept){
    RRvariables <- names(beta)[1:M+1]
  }else{
    RRvariables <- names(beta)[1:M]
  }
  rdf <- startmod$df.residual
  r <- startmod$residuals
  sigma <- sqrt( sum(r^2)/rdf)
  
  pi.est <- matrix(nrow=K, ncol=G,dimnames=list(colnames(PWarray[,,1]),NULL))
  n.star <- matrix(nrow=J, ncol=G,dimnames=list(patterns,NULL))
  # already above: gidx
  #   groupidx <- apply(group,1,paste,collapse="")
  #   groupidx <- as.numeric(factor(groupidx,labels=1:G))
  for (g in 1:G){
    sel <- gidx==g
    # observed prevalence without randomization
    for (j in 1:J){
      n.star[j,g] <-  sum(patterns[j]==responses[sel])
    }
    # correction for RR procedure
    if ( "SLD" %in% models || "UQTunknown" %in% models ||
           ("Kuk" %in% models && Kukrep > 1) ){
      # random start values
      pi.est[,g] <- runif(K)
    }else{
      pi.est[,g] <- solve(PWarray[,,g]) %*% n.star[,g] /sum(sel)
      pi.est[,g] <- pi.est[,g]
    }
  }
  
  
  ############### REPEATED FITTING ##############
  est <- list() ; est.new <- list()
  for (ff in 1:fit.n){
    ## adjust for negative estimates of pi (or larger than 1): set to 0 or 1
    if (ff == 1){
      pi.estR <- rowMeans(pi.est)
      beta1 <- beta
      sigma1 <- sigma
    }else{
      pi.estR <- rowMeans(pi.est)* runif(length(pi.est[,1]),.5,1.5)
      beta1 <- beta*runif(length(beta), .3, 3)
      sigma1 <- sigma*runif(1,.3,3)
    }
    
    pi.estR[pi.estR<0] <- 0
    pi.estR[pi.estR>1] <- 1
    pi.estR <- pi.estR/sum(pi.estR)
    
    phi <- c(beta1,sigma=sigma1,pi=pi.estR)
    phi <- phi[-length(phi)]  # last pi is a fixed parameter: sum(pi)=1 !
    nbeta <- length(beta)
    npi <- length(pi.estR)
    
    # EM algorithm
    
    # NEWTON RAPHSON like opimization
    try(est.new <- optim(par=phi,fn=RRlin.ll, method="L-BFGS-B",
                         lower=c(rep(-Inf,nbeta),0,  rep(0,npi-1)),
                         upper=c(rep( Inf,nbeta),Inf,rep(1,npi-1)),
                         control=list(fnscale=-1, maxit=maxit,
                                      parscale=c(rep(1,nbeta),1,rep(pibeta/max(abs(beta)),npi-1))),
                         hessian=T,  
                         y=y, u=u, gidx=gidx, U=U, intercept=intercept,
                         PWpp=PWpp, nbeta=nbeta, n.star=n.star, PWarray=PWarray, M=M), silent=T)
    if (is.null(est$value)){
      est <- est.new
      start <- phi
    }
    if(!is.null(est.new$value) && est.new$value >est$value){
      est <- est.new
      start <- phi
    }
  }
  
  #   if(est$convergence !=0){
  #     warning("The L-BFGS-B fitting algorithm did not converge. If model$convergence==1, try to fit the model with a larger number of maximum iterations, e.g., maxit=500.")
  #   }
  
  # calculate last value of pi: fixed by sum
  phi <- est$par
  pi <- c(phi[(nbeta+2):length(phi)] , 0)
  pi[npi] <- 1-sum(pi)
  names(pi) <- colnames(PWarray[,,1])
  
  res <- list(models=models,p.list=p.list, PWarray=PWarray, par2=par2, N=N,
              beta=phi[1:nbeta], sigma=phi[nbeta+1],pi=pi, 
              hessian=est$hessian, convergence=est$convergence, 
              message = est$message, counts=est$counts,
              RRvariables=RRvariables, call=match.call(),
              start=start, intercept=intercept, logLik=est$value, bs.n=bs.n) 
  res$npar <- length(phi)
  res$vcov <- est$hessian
  try({
    res$vcov <- solve(-est$hessian)
    dimnames(res$hessian) <- list(names(phi) , names(phi) )
    dimnames(res$vcov) <- list(names(phi) , names(phi) )
  }, silent=T)
  
  res$df <- nrow(unique(cbind(y,x))) -length(phi)
  
  ######################################################################################
  # LR test for each parameter
  #   if (LR.test){
  #     deltaLogLik <- rep(NA,nbeta)
  #     
  #     if (intercept){
  #       # fit model without intercept
  #       try({model2 <- RRlin(formula=y~cbind(w,u)-1, group=group, Kukrep=Kukrep,
  #                            models=models, p.list=p.list)
  #            deltaLogLik[1] <- model2$logLik - est$value})
  #       
  #       # LR test for RR models
  #       if (M==1){
  #         warning('LR test not functional for a single RR variable.')
  #         ## find logLik without RR variables
  # #         model2 <- lm(y ~ u)
  # #         print(model2)
  # #         sss <- summary(model2)$sigma
  # #         print(sss)
  # #         ll0 <- sum(dnorm(y,mean=predict(model2), sd=sss, log=T))
  # #         print(ll0)
  # #         deltaLogLik[2] <- ll0 - est$value
  #       }else{
  #         for (i in 1:M){  
  #           # fit without the RR variable number i
  #           model2 <- RRlin(y~ cbind(w[,-i,drop=F],u), group=group[,-i,drop=F], 
  #                           models=models[-i], p.list=p.list[-i], Kukrep=Kukrep)
  #           deltaLogLik[1+i] <- model2$logLik - est$value
  #         }
  #       }
  #       if (U>0){
  #         for (i in 1:U){  
  #           # fit without the RR variable number i
  #           model2 <- RRlin(y~ cbind(w,u[,-i,drop=F]), group=group, 
  #                           models=models, p.list=p.list, Kukrep=Kukrep)
  #           deltaLogLik[1+M+i] <- model2$logLik - est$value
  #         }
  #       }
  #     }else{
  #       warning('LR test only available for models including an intercept')
  #     }  
  #     
  #     res$prob <- 1-pchisq( -2*deltaLogLik,1)
  #     res$deltaLogLik <- deltaLogLik
  #     names(res$prob) <- names(phi)[1:nbeta]
  #     names(res$deltaLogLik) <- names(phi)[1:nbeta]
  #   }
  
  # nonparametric bootstrap
  if(bs.n>0){
    # use function for parallel computation
    npboot <- function(R){
      # initialize matrices/vectors for results
      bs.beta <- matrix(nrow=R, ncol=nbeta)
      bs.pi <- matrix(nrow=R, ncol=length(pi))
      bs.sigma <- rep(NA, R) 
      bs.logLik <- bs.sigma
      cnt <- 1
      while( cnt <= R){
        # sampling: identical group sizes in bootstrap sample!
        sel <- c()
        for (g in 1:G){
          sel <- c(sel, sample((1:N)[gidx==g], table(gidx)[g], T) ) 
        }
        list(y,sel,models,y,group,w,p.list,u)
        
        # fit model to selected nonparametric BS sample
        try({model.bs <-  RRlin(formula=y[sel,,drop=F]~cbind(w[sel,,drop=F],u[sel,,drop=F]),
                                group=group[sel,,drop=F], Kukrep=Kukrep,
                                models=models, p.list=p.list, fit.n=2)
             bs.beta[cnt,] <- model.bs$beta
             bs.sigma[cnt] <- model.bs$sigma
             bs.pi[cnt,] <- model.bs$pi
             bs.logLik[cnt] <- model.bs$logLik
        })
        if (!is.na(bs.beta[cnt,1])) 
          cnt <- cnt+1
      }
      # return results as list
      bs.res <- list(beta=bs.beta, pi=bs.pi, sigma=bs.sigma, logLik=bs.logLik)
    }
    if (nCPU==1){
      bs.res <- npboot(R=bs.n)
      res$bs.beta <- bs.res$beta; res$bs.sigma <- bs.res$sigma; 
      res$bs.pi <- bs.res$pi; res$bs.logLik <- bs.res$logLik; 
    }else{
      #       require(doParallel, quietly=T)
      if (nCPU=="max"){
        try(nCPU <-  as.numeric(Sys.getenv('NUMBER_OF_PROCESSORS')))
        if (nCPU=="max") nCPU <- 2
      }
      cl.tmp =  makeCluster(nCPU) 
      registerDoParallel(cl.tmp, cores=nCPU)
      bs.res <- foreach(k=1:nCPU , .packages='RRreg') %dopar% { 
        npboot(R=ceiling(bs.n/nCPU)) }
      
      #       res$bs.beta <- bs.res$beta; res$bs.sigma <- bs.res$sigma; 
      #       res$bs.pi <- bs.res$pi; res$bs.logLik <- bs.res$logLik; 
      for (i in 1:nCPU){
        res$bs.beta <- rbind(res$bs.beta, bs.res[[i]]$beta); 
        res$bs.sigma <- c(res$bs.sigma, bs.res[[i]]$sigma); 
        res$bs.pi <- rbind(res$bs.pi, bs.res[[i]]$pi); 
        res$bs.logLik <- c(res$bs.logLik, bs.res[[i]]$logLik);
      }
    }   
  }
  
  # index models which mix DQ and RR procedure
  for (m in 1:M){
    models[m] <- ifelse ( 0 %in% gcombs[,m],
                          paste0(models[m],"+DQ"),
                          models[m])
  }
  res$models <- models  
  res$call <- match.call()
  res$formula <- formula
  class(res) <- "RRlin"
  res
}


#' @export
print.RRlin <- function(x, ...)
{
  cat("Call:\n")
  print(x$call)
  cat("\nRegression coefficients (beta):\n")
  print(round(x$beta,5))
  cat(paste0("\nResidual standard error (sigma): ",round(x$sigma,5),"\n"))
  cat("\nPrevalence estimates (pi):\n")
  print(round(x$pi,5))
}


#' @export
summary.RRlin <- function(object, ...)
{
  se <- sqrt(diag(object$vcov))
  wald_chi <- (object$beta/se[1:length(object$beta)])^2
  bcoef <- cbind(Estimate = object$beta,
                 StdErr = se[1:length(object$beta)],
                 "Wald test"=wald_chi,
                 "Pr(>Chi2,df=1)"= 1-pchisq(wald_chi, df=1) )
  #   if(object$LR.test){
  #     bcoef <- cbind(bcoef,
  #                    "deltaG2"= -2*object$deltaLogLik,
  #                    "Pr(>deltaG2)" = object$prob)
  #   }
  
  pcoef <- cbind(Estimate = object$pi[-length(object$pi)],
                 StdErr = se[(length(object$beta)+2):length(se)]
  )
  #    rownames(pcoef) <- paste("pi =",rownames(pcoef))
  fitInfo <- c(N=object$N, logLik= object$logLik, df=object$df)
  p.char <- unlist( as.character(lapply(object$p.list, round, 3)))
  
  res <- list(call=object$call, 
              RRvariables = object$RRvariables, p.char = p.char,models=object$models,
              bcoef=bcoef, pcoef=pcoef, sigma=object$sigma, 
              intercept=object$intercept, N=object$N, bs.n=object$bs.n
  )
  if (object$bs.n>0){
    bs.tab <- matrix(c(colMeans(object$bs.beta,na.rm=T), 
                       mean(object$bs.sigma,na.rm=T), 
                       colMeans(object$bs.pi,na.rm=T), 
                       apply(object$bs.beta,2,sd,na.rm=T),
                       sd(object$bs.sigma,na.rm=T),
                       apply(object$bs.pi,2,sd,na.rm=T)
    ), ncol=2)
    colnames(bs.tab) <- c("Mean","SE")
    rownames(bs.tab) <- c(names(object$beta),"sigma",paste("pi =",names(object$pi))) #,"logLik")
    res$bs.tab <- bs.tab
  }
  class(res) <- "summary.RRlin"
  res
}


#' @export
print.summary.RRlin <- function(x, ...){
  cat("Call:\n")
  print(x$call)
  cat("\nRandomized response variables:\n")
  TAB <- cbind(Variable = x$RRvariables,
               Model = x$models,
               p = x$p.char)
  rownames(TAB)=1:length(x$models)
  print(TAB, quote = FALSE)
  cat("\nCoefficients (beta):\n")
  printCoefmat( round(x$bcoef,5))
  cat(paste0("\nResidual standard error (sigma): ",round(x$sigma,3), "; N=",x$N,"\n"))
#   cat("\nPrevalence estimates for combinations of RR responses:\n")
#   printCoefmat(round(x$pcoef,5))
  if(x$bs.n>0){
    cat(paste0("\n\nResults of nonparametric bootstrap (",x$bs.n," samples):\n"))
    printCoefmat(round(x$bs.tab,4))
  }
}

#' @export
logLik.RRlin <- function(object, ...){
  return(object$logLik)
}

#' @export
vcov.RRlin <- function(object, ...){
  return(object$vcov)
}

### bootstrap methods
# bs.se <- function(ests,true){
#   if (class(ests) == 'matrix'){
#     sqrt(apply((t(ests)-true)^2,1,sum,na.rm=TRUE) /nrow(ests))
#   }else{
#     sqrt(sum((ests-true)^2,na.rm=TRUE)/(length(ests)))
#   }        
# }