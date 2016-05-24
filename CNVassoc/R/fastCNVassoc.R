fastCNVassoc <- function(probs, formula, data, model = "additive",
    family = "binomial", nclass = 3, colskip = 5, tol = 1e-06, max.iter = 30, verbose = FALSE, multicores=0)
{
    cl <- match.call()
    if (missing(data)) 
        data <- environment(formula)
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "na.action"), names(mf), 0L) #offset
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    mt <- attr(mf, "terms")
    model.type <- charmatch(model, c("additive","multiplicative"))
    if (is.na(model.type))
      stop(" argument 'model' must be either 'multiplicative' or 'additive'")
    if (model.type!=1)
      stop(" only 'additive' model is implemented")
    if (!family %in% c("weibull","binomial"))
      stop(" only 'biomial' and 'weibull' models are implemented by now")
    nIndiv <- NROW(mf)
    mm <- model.matrix(mt, mf, contrasts)
    y <- model.response(mf)
    nVar <- NCOL(mf) - 1
    if (nVar != 0) {
        Xcov <- mm[, -1, drop = FALSE]
        nCov <- NCOL(Xcov)
    }
    else {
        Xcov <- NULL
        nCov <- 0
    }
    N<-as.integer(nIndiv)
    
    # RESPONSE
    if (family!="weibull")
      Y<-y
    else {
      Y<-y[,1]
      CENS<-y[,2]
    }
    
    # PARAM INI
    if (nVar > 0){
      if (family=="binomial")
        fit0<-glm(Y~Xcov,family="binomial")
      if (family=="weibull")
        fit0<-survreg(Surv(Y,CENS)~Xcov)
    }else{
      if (family=="binomial")
        fit0<-glm(Y~1,family="binomial")
      if (family=="weibull")
        fit0<-survreg(Surv(Y,CENS)~1)
    }
    if (family!="weibull")
      PARAM<-c(0,coef(fit0))
    else
      PARAM<-c(0,unlist(translate(fit0)))
    PARAM[1:2]<-PARAM[2:1] 
    
    # OTHER INPUTS
    TOL <- tol
    MAX_ITER <- as.integer(max.iter)
    VERBOSE <- as.integer(0)
    P <- as.integer(nclass)
    Q <- as.integer(nCov)
    Y <- as.double(Y)
    if (nCov>0)
      C <- t(Xcov)
    else
      C <- NULL

    # PROBS
    if (is.character(probs)){
      cat("Scanning probs data...\n")
      read.st <- system.time(PROBS <- scan(probs,what="character",sep="\n",quiet=TRUE))[3]
      cat("Done! Took ",read.st,"seconds\n\n")
      nCNVs<-length(PROBS)
    } else {
      if (is.matrix(probs) || is.data.frame(probs))
        PROBS<-as.matrix(probs)  
      else
        stop(" 'probs' must be a character (file) or a matrix or a data.frame")        
      nCNVs<-nrow(PROBS)
    }
    
    # RUN
    if (multicores==0){  # no use of parallel package
      ans<-lapply(1:nCNVs, function(i){
        if (verbose)
          cat("Iteration ",i,"\n")      
        if (is.character(probs)){
          PROBS.i<-strsplit(PROBS[i],split=" ")[[1]]
        }else
          PROBS.i<-PROBS[i,]
        if (colskip>0)
          PROBS.i<-PROBS.i[-(1:colskip)]
        WW<-matrix(as.double(PROBS.i),nrow=nclass)
        if (nCov > 0){
          if (family=="binomial")
            fit <- try(.Call("NRlogisticcov", N, PARAM, P, Y, WW, Q, C, TOL, MAX_ITER, VERBOSE))
          if (family=="weibull")
            fit <- try(.Call("NRweibullcov", N, PARAM, P, Y, CENS, WW, Q, C, TOL, MAX_ITER, VERBOSE))
        } else{ 
          if (family=="binomial")
            fit <- try(.Call("NRlogistic", N, PARAM, P, Y, WW, TOL, MAX_ITER, VERBOSE))
          if (family=="weibull")
            fit <- try(.Call("NRweibull", N, PARAM, P, Y, CENS, WW, TOL, MAX_ITER, VERBOSE))          
        }
        if (inherits(fit,"try-error")) 
          return(rep(NA,3))
        logor<-fit[[1]][2]
        se<-try(sqrt(diag(solve(-fit[[4]])))[2],silent=TRUE)
        if (inherits(se,"try-error")) 
          se<-NA
        iter <- fit[[5]]
        return(c(logor,se,iter))
      })
      ans<-matrix(unlist(ans),ncol=3,byrow=TRUE)
    } else {            # use of parallel package
      requireNamespace("parallel", quietly=TRUE)
      ans<-parallel::mclapply(1:nCNVs, function(i){
        if (verbose)
          cat("Iteration",i,"\n")      
        if (is.character(probs)){
          PROBS.i<-strsplit(PROBS[i],split=" ")[[1]]
        } else
          PROBS.i<-PROBS[i,]
        if (colskip>0)
          PROBS.i<-PROBS.i[-(1:colskip)]
        WW<-matrix(as.double(PROBS.i),nrow=nclass)
        if (nCov > 0){
          if (family=="binomial")
            fit <- try(.Call("NRlogisticcov", N, PARAM, P, Y, WW, Q, C, TOL, MAX_ITER, VERBOSE))
          if (family=="weibull")
            fit <- try(.Call("NRweibullcov", N, PARAM, P, Y, CENS, WW, Q, C, TOL, MAX_ITER, VERBOSE))
        }else{ 
          if (family=="binomial")
            fit <- try(.Call("NRlogistic", N, PARAM, P, Y, WW, TOL, MAX_ITER, VERBOSE))
          if (family=="weibull")
            fit <- try(.Call("NRweibull", N, PARAM, P, Y, CENS, WW, TOL, MAX_ITER, VERBOSE))          
        }
        if (inherits(fit,"try-error")) 
          return(rep(NA,3))
        logor<-fit[[1]][2]
        se<-try(sqrt(diag(solve(-fit[[4]])))[2],silent=TRUE)
        if (inherits(se,"try-error")) 
          se<-NA        
        iter <- fit[[5]]
        return(c(logor,se,iter))
      },mc.cores=multicores)
      ans<-matrix(unlist(ans),ncol=3,byrow=TRUE)
    }

   zscore <- ans[,1]/ans[,2]
   pval <- 2*(1-pnorm(abs(zscore)))
   if (nrow(ans)==1)
     ans <- rbind(c(ans[1,1:2], zscore, pval, ans[1,3]))  
   else
     ans <- cbind(ans[,1:2], zscore, pval, ans[,3])
   colnames(ans)<-c("beta","se","zscore","pvalue","iter")
   ans<-data.frame("variant"=1:nrow(ans),ans)
   ans

}


