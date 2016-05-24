######################################################
### Authors: Simone Padoan.
### Emails: simone.padoan@unibocconi.it,
### moreno.bevilacqua@uv.cl
### Institutions: Department of Decision Sciences,
### University Bocconi of Milan and
### Departamento de Estad?stica
### Universidad de Valparaiso
### File name: Statistics.r
### Description:
### This file contains a set of procedures
### for the computation of composite likelihood-based
### statistics and tests.
### Last change: 28/03/2013.
######################################################

### Procedures are in alphabetical order.

### Statistical hypothesis testing for nested models
HypoTest <- function(object1, object2, ..., statistic)
  {
      # check the type of test:
      Istest <- function(statistic)
      {
          Istest <- switch(statistic,
                           Rao=1,
                           Wald=2,
                           WilksCB=3,
                           WilksPSS=4,
                           WilksRJ=5,
                           WilksS=6)
          return(Istest)
      }
      # check if correlation models are nested:
      Isnested <- function(model1, model2)
      {
          result <- FALSE
          if(identical(model1, model2))
              result <- TRUE
          if(identical(model1, "stable"))
              if(identical(model2, "exponential"))
                  result <- TRUE
          if(identical(model1, "matern"))
              if(identical(model2, "exponential"))
                  result <- TRUE
          if(identical(model1, "matern_cauchy"))
              if(identical(model2, "exp_cauchy"))
                  result <- TRUE
          if(identical(model1, "matern_exp"))
              if(identical(model2, "exp_exp"))
                  result <- TRUE
          return(result)
      }
      # compute gradient, sensitivity and variability under the null:
      UnderNull <- function(model1, model2)
      {
          fixed <- model1$fixed
          if(!is.null(fixed)){
              start <- model2$fixed[!names(model2$fixed)%in%names(fixed)]
              namesparam <- names(start)
              fixed <- as.list(model1$fixed)
              start <- as.list(c(start,model2$param))}
          else{
              fixed <- NULL
              start <- as.list(c(model2$param,model2$fixed))
              namesparam <- names(model2$fixed)}
          # derive the initial parameters:
          initparam <- InitParam(model2$coordx,model2$coordy,model2$coordt,model2$corrmodel,
                                 model2$data,model2$distance,"Fitting",fixed,model2$grid,
                                 model2$likelihood,model2$margins,model2$srange[2],model2$trange[2],
                                 model2$model,NULL,NULL,NULL,FALSE,model2$numrep,start,NULL,NULL,
                                 model2$threshold,model2$type,model2$type,TRUE,model2$vartype,
                                 NULL,model2$winconst,model2$winstp)
          # set useful quantities for the computation of the Godambe matrix
          numparam <- initparam$numparam
          dimmat <- numparam^2
          dmat <- numparam*(numparam+1)/2
          score <- double(numparam)
          sensmat <- double(dmat)# H - sensitivity matrix
          varimat <- double(dmat)# J - variability matrix
          eps <- (.Machine$double.eps)^(1/3)
          spacetime <- length(model2$coordt)>1
          param <- c(initparam$param,initparam$fixed)
          paramcorr <- param[initparam$namescorr]
          nuisance <- param[initparam$namesnuis]
          # compute the gradient and the components of the Godambe matrix:
          GD=.C('GodambeMat',as.double(model2$coordx),as.double(model2$coordy),
             as.integer(initparam$corrmodel),as.double(model2$data),as.integer(initparam$distance),
             as.double(eps),as.integer(initparam$flagcorr),as.integer(initparam$flagnuis),
             as.integer(model2$grid),as.integer(initparam$likelihood),
             as.integer(initparam$model),as.integer(numparam),
             as.integer(initparam$numparamcorr),as.double(paramcorr),as.double(nuisance),
             score=score,sensmat=sensmat,as.integer(spacetime),as.double(model2$threshold),
             as.integer(initparam$type),varimat=varimat,as.integer(initparam$vartype),
             as.double(initparam$winconst),as.double(initparam$winstp),
             PACKAGE='CompRandFld',DUP=TRUE,NAOK=TRUE)
             
             sensmat<-GD$sensmat; varimat<-GD$varimat;score<-GD$score
          # define the sensitivity and variability matrices:
          H <- matrix(rep(0,dimmat),ncol=numparam)# H - sensitivity matrix
          J <- matrix(rep(0,dimmat),ncol=numparam)# J - variability matrix
          namesmat <- c(initparam$namesnuis[as.logical(initparam$flagnuis)],
                        initparam$namescorr[as.logical(initparam$flagcorr)])
          names(score) <- namesmat
          score <- score[initparam$namesparam]
          dimnames(H) <- list(namesmat,namesmat)
          dimnames(J) <- list(namesmat,namesmat)
          if(numparam>1){
              H[lower.tri(H,diag=TRUE)] <- sensmat
              H <- t(H)
              H[lower.tri(H,diag=TRUE)] <- sensmat
              J[lower.tri(J,diag=TRUE)] <- varimat
              J <- t(J)
              J[lower.tri(J,diag=TRUE)] <- varimat
              H <- H[initparam$namesparam,initparam$namesparam]
              J <- J[initparam$namesparam,initparam$namesparam]}
          else{
              H[1,1] <- sensmat
              J[1,1] <- varimat}
          # Delete global variables:
          .C('DeleteGlobalVar', PACKAGE='CompRandFld', DUP=TRUE, NAOK=TRUE)
          return(list(score=score, sensmat=H, varimat=J))
      }
      # compute the statistic:
      StatiTest <- function(df, model1, model2, statistic)
      {
          namesparam <- names(model1$param)[!names(model1$param)%in%names(model2$param)]
          # compute the Wald-type statistic:
          if(statistic=='Wald'){
              theta <- model1$param[namesparam]-model2$fixed[namesparam]
              varcov <- model1$varcov[namesparam,namesparam]# Restricted variance-covariance matrix
              W <- t(theta)%*%solve(varcov)%*%theta
              nu <- df}
          if(statistic=="WilksCB"){
              W <- 2*(model1$logCompLik-model2$logCompLik)
              theta <- model1$param[namesparam]-model2$fixed[namesparam]
              G <- solve(model1$varcov)[namesparam,namesparam]
              H <- model1$sensmat[namesparam,namesparam]
              W <- W*((t(theta)%*%G%*%theta)/(t(theta)%*%H%*%theta))
              nu <- df}
          if(any(statistic==c("WilksPSS","WilksRJ","WilksS","Rao"))){
              null <- UnderNull(model1,model2)
              H <- null$sensmat# H - sensitivity matrix
              Hi <- solve(H)# inverse of H
              J <- null$varimat# J - variability matrix
              Ji <- solve(J)# inverse of J
              G <- H%*%Ji%*%H# compute the Godambe matrix
              varcov <- solve(G)# variance-covariance matrix
              if(statistic=="WilksPSS"){
                  W <- 2*(model1$logCompLik-model2$logCompLik)
                  varcov <- varcov[namesparam,namesparam]
                  G <- try(solve(varcov),silent=TRUE)
                  score <- null$score[namesparam]# select the score related to the null
                  Hi <- Hi[namesparam,namesparam]# select the inversed of sens related to the null
                  Rao <- t(score)%*%Hi%*%G%*%Hi%*%score# compute the Rao statistic
                  W <- W*Rao/(t(score)%*%Hi%*%score)
                  nu <- df}
              if(statistic=="WilksRJ"){
                  W <- 2*(model1$logCompLik-model2$logCompLik)
                  lambda <- eigen(solve(Hi[namesparam, namesparam])%*%varcov[namesparam, namesparam])$values
                  W <- W/mean(lambda)
                  nu <- df}
              if(statistic=="WilksS"){
                  W <- 2*(model1$logCompLik-model2$logCompLik)
                  lambda <- eigen(solve(Hi[namesparam, namesparam])%*%varcov[namesparam, namesparam])$values
                  slambda <- sum(lambda)
                  nu <- slambda^2/sum(lambda^2)
                  W <- nu*W/slambda}
              if(statistic=="Rao"){
                  varcov <- varcov[namesparam,namesparam]
                  G <- try(solve(varcov),silent=TRUE)
                  score <- null$score[namesparam]# select the score related to the null
                  Hi <- Hi[namesparam,namesparam]# select the inversed of sens related to the null
                  W <- t(score)%*%Hi%*%G%*%Hi%*%score# compute the Rao statistic
                  nu <- df}}
          return(list(W=W,nu=nu))
      }
      ### START THE MAIN BODY OF THE PROCEDURE
      # check fitted models:
      if(any(missing(object1), missing(object2)))
         stop("Models one and two must be specified\n")
      # check statistics:
      if(missing(statistic))
         stop('Insert the type of statistic use in the hypothesis test\n')
      if(is.null(Istest(statistic)))
          stop("The name of test does not match with one those available\n")
      # check if there are multipl fitted models:
      objects <- as.list(substitute(list(...)))[-1]
      objects <- sapply(objects,function(x) deparse(x))
      if(!length(objects)) objects <- NULL
      # build a sequence of fitted models:
      models <- c(deparse(substitute(object1)),
                  deparse(substitute(object2)),
                  objects)
      nummod <- length(models)# number of models
      numparam <- NULL# parameters size
      numstat <- length(statistic)# number of statistics
      lmodels <- vector("list", nummod)# define a list of models
      W <- double(nummod - 1)# statistics
      pvalue <- double(nummod - 1)# p-values
      df <- double(nummod - 1)# degrees of freedoms (of the tests)
      nu <- double(nummod - 1)# adjusted df
      # Compute hypothesis tests:
      for(i in 1:nummod)
      {   # consider the ith model:
          model <- get(models[i], envir=parent.frame())
          if(!inherits(model, "FitComposite"))
              stop("use HypoTest only with 'FitComposite' objects\n")
          numparam <- c(numparam, length(model$param))
          lmodels[[i]] <- model
          if(!is.matrix(lmodels[[i]]$varcov))
              stop("one of the fitted models does not have a valid variance-covariance matrix\n")
          if(i>1){
            j <- i-1
            if((!all(names(lmodels[[j]]) %in% names(lmodels[[i]]))) &&
               (!identical(lmodels[[j]]$model,lmodels[[i]]$model)) &&
               (!Isnested(lmodels[[j]]$corrmodel,lmodels[[i]]$corrmodel)))
              stop("models are not nested\n")
            # Define the degrees of freedom:
            df[j] <- length(lmodels[[j]]$param)-length(lmodels[[i]]$param)
            if(df[j] <= 0) stop("model are not nested\n")
            stat <- StatiTest(df[j],lmodels[[j]],lmodels[[i]],statistic)
            nu[j] <- stat$nu
            W[j] <- stat$W
            pvalue[j] <- pchisq(W[j], df = nu[j], lower.tail = FALSE)
          }
      }
      ### END THE MAIN BODY OF THE PROCEDURE
      #print a table with the hypothesis testing:
      table <- data.frame(numparam, c(NA, df), c(NA, nu),c(NA, W), c(NA, pvalue))
      dimnames(table) <- list(models, c("Num.Par", "Diff.Par", "Df","Chisq", "Pr(>chisq)"))
      structure(table, heading = c("Statistical Hypothesis Test Table\n"),
                class = c("data.frame"))
  }
