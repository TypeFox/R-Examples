################################################################################

## slightly modified from glmnet's cvcompute function

.relax.cvcompute <- function (mat, weights, foldid, nlams) 
{
    wisum = tapply(weights, foldid, sum)
    nfolds = max(foldid)
    outmat = matrix(NA, nfolds, ncol(mat))
    good = matrix(0, nfolds, ncol(mat))
    mat[is.infinite(mat)] = NA
    for (i in seq(nfolds)) {
        mati = mat[foldid == i, , drop = FALSE]
        wi = weights[foldid == i]
        outmat[i, ] = apply(mati, 2, weighted.mean, w = wi, na.rm = TRUE)
        good[i, seq(nlams[i])] = 1
    }
    N = apply(good, 2, sum)
    list(cvraw = outmat, weights = wisum, N = N)
}

## slightly modified from glmnet's cv.elnet function

.relax.cv.elnet <- function(outlist,
                           lambda,
                           x,
                           y,
                           weights,
                           offset,
                           foldid,
                           type.measure,
                           grouped) {
  typenames=c(deviance="Mean-Squared Error",mse="Mean-Squared Error",mae="Mean Absolute Error")
  if(type.measure=="default")type.measure="mse"
  if(!match(type.measure,c("mse","mae","deviance"),FALSE)){
    warning("Only 'mse', 'deviance' or 'mae'  available for Gaussian models; 'mse' used")
    type.measure="mse"
  }
     if(!is.null(offset))y=y-drop(offset)
     predmat=matrix(NA,length(y),length(lambda))
    nfolds=max(foldid)
    nlams=double(nfolds)
    for(i in seq(nfolds)){
      which=foldid==i
      fitobj=outlist[[i]]

      if(inherits(fitobj, "relaxnet.intercept.only")) {

        predmat[which, ] <- fitobj
        ## just fill it in with intercept value stored in fitobj

        nlami <- length(lambda)
        
      } else {
        
        ## the beta mat, intercepts (a0) and lambda should already
        ## have been subset so that only the lower lambda values
        ## are processed (no need to repeat the values that were
        ## already done in the main model)

        fitobj$offset=FALSE
        preds <- predict(fitobj,
                         x[which, rownames(fitobj$beta),
                           drop=FALSE])
        nlami=length(fitobj$lambda)
        predmat[which,seq(nlami)]=preds
      }
      nlams[i]=nlami
    }

  N=length(y) - apply(is.na(predmat),2,sum)
  cvraw=switch(type.measure,
    "mse"=(y-predmat)^2,
    "deviance"=(y-predmat)^2,
    "mae"=abs(y-predmat)
    )
   if( (length(y)/nfolds <3)&&grouped){
    warning("Option grouped=FALSE enforced in cv.relaxnet, since < 3 observations per fold",call.=FALSE)
    grouped=FALSE
  }
 if(grouped){
   cvob=.relax.cvcompute(cvraw,weights,foldid,nlams)
  cvraw=cvob$cvraw;weights=cvob$weights;N=cvob$N
 }

  cvm=apply(cvraw,2,weighted.mean,w=weights,na.rm=TRUE)
  cvsd=sqrt(apply(scale(cvraw,cvm,FALSE)^2,2,weighted.mean,w=weights,na.rm=TRUE)/(N-1))
  list(cvm=cvm,cvsd=cvsd,name=typenames[type.measure])
}

################################################################################

## slightly modified from glmnet's cv.lognet function

.relax.cv.lognet <- function(outlist,
                            lambda,
                            x,
                            y,
                            weights,
                            offset,
                            foldid,
                            type.measure,
                            grouped){
  typenames=c(mse="Mean-Squared Error",mae="Mean Absolute Error",deviance="Binomial Deviance",auc="AUC",class="Misclassification Error")
  if(type.measure=="default")type.measure="deviance"
  if(!match(type.measure,c("mse","mae","deviance","auc","class"),FALSE)){
    warning("Only 'deviance', 'class', 'auc', 'mse' or 'mae'  available for binomial models; 'deviance' used")
    type.measure="deviance"
  }

###These are hard coded in the Fortran, so we do that here too
  prob_min=1e-5
  prob_max=1-prob_min
  ###Turn y into a matrix
  nc = dim(y)
  if (is.null(nc)) {
    y = as.factor(y)
    ntab = table(y)
    nc = as.integer(length(ntab))
    y = diag(nc)[as.numeric(y), ]
  }
  N=nrow(y)
  nfolds=max(foldid)
  if( (N/nfolds <10)&&type.measure=="auc"){
    warning("Too few (< 10) observations per fold for type.measure='auc' in cv.lognet; changed to type.measure='deviance'. Alternatively, use smaller value for nfolds",call.=FALSE)
    type.measure="deviance"
  }
  if( (N/nfolds <3)&&grouped){
    warning("Option grouped=FALSE enforced in cv.relaxnet, since < 3 observations per fold",call.=FALSE)
    grouped=FALSE
  }


  if(!is.null(offset)){
    is.offset=TRUE
    offset=drop(offset)
  }else is.offset=FALSE
  predmat=matrix(NA,nrow(y),length(lambda))
  nlams=double(nfolds)
  for(i in seq(nfolds)){
    which=foldid==i
    fitobj=outlist[[i]]

    if(inherits(fitobj, "relaxnet.intercept.only")) {

      if(is.offset)
        stop("Internal Error:\n",
             "haven't dealt with offsets for binomial intercept models yet")
      
      predmat[which, ] <- 1 / (1 + exp(-fitobj))
      ## just fill it in with intercept value stored in fitobj
      ## check this again to make sure I'm doing this right
      
      nlami <- length(lambda)
        
    } else {

      ## the beta mat, intercepts (a0) and lambda should already
      ## have been subset so that only the lower lambda values
      ## are processed (no need to repeat the values that were
      ## already done in the main model)


      if(is.offset) off_sub=offset[which]
      preds <- predict(fitobj,
                       x[which, rownames(fitobj$beta), drop=FALSE],
                       offset=off_sub, type="response")
      nlami=length(fitobj$lambda)
      predmat[which,seq(nlami)]=preds
    }
    
    nlams[i]=nlami
  }
  ##If auc we behave differently
  if(type.measure=="auc") {
    cvraw=matrix(NA,nfolds,length(lambda))
    good=matrix(0,nfolds,length(lambda))
    for(i in seq(nfolds)){
      good[i,seq(nlams[i])]=1
      which=foldid==i
      for(j in seq(nlams[i])){
        cvraw[i,j]=auc.mat(y[which,],predmat[which,j],weights[which])
      }
    }
    N=apply(good,2,sum)
    weights=tapply(weights,foldid,sum)

  } else {

    ##extract weights and normalize to sum to 1
    ywt=apply(y,1,sum)
    y=y/ywt
    weights=weights*ywt

    N=nrow(y) - apply(is.na(predmat),2,sum)
    cvraw=switch(type.measure,
      "mse"=(y[,1]-(1-predmat))^2 +(y[,2]-predmat)^2,
      "mae"=abs(y[,1]-(1-predmat)) +abs(y[,2]-predmat),
      "deviance"= {
        predmat=pmin(pmax(predmat,prob_min),prob_max)
        lp=y[,1]*log(1-predmat)+y[,2]*log(predmat)
        ly=log(y)
        ly[y==0]=0
        ly=drop((y*ly)%*%c(1,1))
        2*(ly-lp)},
      "class"=y[,1]*(predmat>.5) +y[,2]*(predmat<=.5)
      )
    if(grouped){
      cvob=.relax.cvcompute(cvraw,weights,foldid,nlams)
      cvraw=cvob$cvraw;weights=cvob$weights;N=cvob$N
    }
  }
  cvm=apply(cvraw,2,weighted.mean,w=weights,na.rm=TRUE)
  cvsd=sqrt(apply(scale(cvraw,cvm,FALSE)^2,2,weighted.mean,w=weights,na.rm=TRUE)/(N-1))
  list(cvm=cvm,cvsd=cvsd,name=typenames[type.measure])
}



################################################################################

## cross-validation (on lambda only) for relaxnet models
## adapted from the cv.glmnet function from package glmnet

cv.relaxnet <- function(x, y, family = c("gaussian", "binomial"),
                        nlambda = 100,
                        alpha = 1,

                        relax = TRUE,
                        relax.nlambda = 100,
                        relax.max.vars = min(nrow(x), ncol(x)) * 0.8,

                        lambda = NULL,
                        relax.lambda.index = NULL,
                        relax.lambda.list = NULL,

                        ##type.measure=c("mse","deviance","class","auc","mae"),
                        ## just set it in code for now
                        
                        nfolds = 10, ## set min for this to 3
                        foldid,
                        
                        multicore = FALSE,
                        mc.cores,
                        mc.seed = 123,
                              
                        ...) {

  start.time <- Sys.time()

  ## can do a lot of the argument checking in relaxnet function
  
  family = match.arg(family)

  if(family == "gaussian") type.measure <- "mse"
  if(family == "binomial") type.measure <- "deviance"

  if(relax && ncol(x) == 1) {

    warning("x has only one column, setting relax to FALSE")
    relax <- FALSE
  }

  if(!relax) { ## just run cv.glmnet

    cv.glmnet.result <- cv.glmnet(x = x, y = y, family = family,
                                  nlambda = nlambda,
                                  alpha = alpha,
                                  lambda = lambda,
                                  nfolds = nfolds,
                                  foldid = foldid,
                                  type.measure = type.measure,
                                  ...)


    relaxnet.fit <- list(main.glmnet.fit = cv.glmnet.result$glmnet.fit,
                         relax = relax,
                         relax.glmnet.fits = NA,
                         relax.num.vars = NA,
                         relax.lambda.index = NA,
                         total.time = NA,
                         main.fit.time = NA,
                         relax.fit.times = NA)

    class(relaxnet.fit) <- "relaxnet"
    
    min.cvm <- min(cv.glmnet.result$cvm)

    end.time = Sys.time()
    
    obj <- list(call = match.call(),
                relax = relax,
                lambda = cv.glmnet.result$lambda,
                cvm = cv.glmnet.result$cvm,
                cvsd = cv.glmnet.result$cvsd,
                cvup = cv.glmnet.result$cvup,
                cvlo = cv.glmnet.result$cvlo,
                nzero = cv.glmnet.result$nzero,
                name = cv.glmnet.result$cvname,
                relaxnet.fit = relaxnet.fit,
                relax.cvstuff.list = NA,
                relax.lambda.list.trunc = NA,
                which.model.min = "main",
                overall.lambda.min = cv.glmnet.result$lambda.min,
                min.cvm = min.cvm,

                ## maybe take these next three out
                main.lambda.min = cv.glmnet.result$lambda.min,
                main.lambda.1se = cv.glmnet.result$lambda.1se,
                main.min.cvm = min.cvm,

                total.time = as.double(difftime(end.time, start.time,
                                              units = "secs")),
                full.data.fit.time = NA,
                cv.fit.times = NA)


    class(obj) <- "cv.relaxnet"

    return(obj)
  }

  lambda.gap <- 10e-10 ## necessary since the lambda values are sometimes
  ## perturbed slightly, probably somewhere in the fortran code
  
  ## if(missing(type.measure))type.measure="default"
  ## else type.measure=match.arg(type.measure)
  
  ## just do this for now
  ## offset <- NULL
  ## weights <- rep(1, nrow(x))

  if(!is.null(lambda) && length(lambda) < 2)
    stop("Need more than one value of lambda for cv.relaxnet")

  if(!is.null(relax.lambda.index) && is.null(lambda))
    stop("if you are specifying a relax.lambda.index",
         "you must also specify a lambda")

  N=nrow(x)

  if(missing(foldid)) {

    ## check nfolds here
    
    foldid <- sample(rep(seq(nfolds), length=N))

  } else {

    ## check foldid here
    
    nfolds <- length(table(foldid))
  }

  if(nfolds < 3)
    stop("nfolds must be bigger than 3; nfolds=10 recommended")

  ## multicore stuff -- adapted from multiPIM package

  if(multicore) {

    tryCatch(library(parallel), error = function(e) {

      stop("multicore is TRUE, but unable to load package parallel.\n",
           "Error message was:\n",
           e$message)
    })

    if(missing(mc.cores))
      stop("if multicore = TRUE, mc.cores must.be.specified")

    if(mode(mc.cores) != "numeric" || length(mc.cores) != 1
       || mc.cores < 1 || mc.cores %% 1 != 0)
      stop("mc.cores should be a single integer giving the number of\n",
           "cores/CPUs to be used for multicore processing")

    ## set mc.cores to nfolds if it's greater

    if(mc.cores > nfolds) mc.cores <- nfolds

    ## check mc.seed

    if(mode(mc.seed) != "numeric"
       || length(mc.seed) != 1
       || (mc.seed %% 1) != 0)
      stop("mc.seed must be a single integer")

    ## set RNG kind to l'ecuyer, saving current RNG kinds

    prev.RNG.kinds <- RNGkind("L'Ecuyer-CMRG")

    ## set the seed

    set.seed(mc.seed)
  }

  ## add this once you add the weights argument
  ##if(missing(weights))weights=rep(1.0,N)else weights=as.double(weights)
  
  ##Fit the model once to get dimensions etc of output

  y=drop(y) # we dont like matrix responses unless we need them

  relaxnet.fit <- relaxnet(x, y, family = family,
                              nlambda = nlambda,
                              alpha = alpha,

                              relax = relax,
                              relax.nlambda = relax.nlambda,
                              relax.max.vars = relax.max.vars,

                              lambda = lambda,
                              relax.lambda.index = relax.lambda.index,
                              relax.lambda.list = relax.lambda.list,

                              ...)
  
  is.offset=relaxnet.fit$main.glmnet.fit$offset


  lambda <- relaxnet.fit$main.glmnet.fit$lambda
  relax.lambda.index <- relaxnet.fit$relax.lambda.index

  relax.lambda.list <- lapply(relaxnet.fit$relax.glmnet.fits,
                              function(fit) fit$lambda)

  ## if(inherits(glmnet.object,"multnet")){
  ##   nz=predict(glmnet.object,type="nonzero")
  ##   nz=sapply(nz,function(x)sapply(x,length))
  ##   nz=ceiling(apply(nz,1,median))
  ## }
  ## else
  nz <- sapply(predict(relaxnet.fit$main.glmnet.fit, type="nonzero"),
               length)

  outlist <- vector("list", length = nfolds)

  cv.fit.times = rep(as.double(NA), nfolds)
  
  ##Now fit the nfold models and store them

    
  run.one.fold <- function(which) { ## which is logical index
                                        # identifying the current fold

    if(is.matrix(y)){
      y.sub <- y[!which,]
    } else {
      y.sub <- y[!which]
    }

    ## subset the offset here if doing that

    time1 <- Sys.time()

    fold.relaxnet.result <- relaxnet(x[!which, , drop=FALSE], y.sub,
                                     family = family,
                                     alpha = alpha,
                                     
                                     relax = relax, 
                             
                                     lambda = lambda,
                                     relax.lambda.index = relax.lambda.index,
                                     relax.lambda.list = relax.lambda.list,
                                     ...)
    time2 <- Sys.time()

    return(list(fold.relaxnet.result = fold.relaxnet.result,
                fold.time = as.double(difftime(time2, time1, units = "secs"))))
  }

  which.list = lapply(1:nfolds, function(i) foldid == i)

  if(multicore) {

    ## needed for reproducibility

    mc.reset.stream()

    ## start parallel jobs

    fold.results <- mclapply(which.list, run.one.fold,
                             mc.preschedule = TRUE,
                             mc.set.seed = TRUE,
                             mc.cores = mc.cores)

  } else {

    fold.results <- lapply(which.list, run.one.fold)

  }

  outlist <- lapply(fold.results, function(x) x$fold.relaxnet.result)

  cv.fit.times <- sapply(fold.results, function(x) x$fold.time)
  
  ###What to do depends on the type.measure and the model fit

  fun <- paste("cv", class(relaxnet.fit$main.glmnet.fit)[[1]], sep=".")

  ## Do main fits first
  
  current.outlist <- lapply(outlist, function(fit) fit$main.glmnet.fit)

  ## just do this for now
  offset <- NULL
  weights <- rep(1, nrow(x))
  grouped <- TRUE

  cvstuff <- do.call(fun, list(current.outlist,
                               lambda,
                               x,
                               y,
                               weights,
                               offset,
                               foldid,
                               type.measure,
                               grouped))

  cvm=cvstuff$cvm
  cvsd=cvstuff$cvsd
  cvname=cvstuff$name

  relax.num.models <- length(relax.lambda.index)
  
  relax.cvstuff.list <- vector("list", length = relax.num.models)

  fun <- paste(".relax.", fun, sep = "")
  
  relax.lambda.list.trunc <- vector("list", length = relax.num.models)

  for(i in 1:relax.num.models) {

    current.outlist <- lapply(outlist, function(fit) fit$relax.glmnet.fits[[i]])

    ## before processing, remove duplicate lambda values which were already
    ## done in the main model -- remove from the fit objects and the lambda.list

    relax.lambda.start.val <- lambda[relax.lambda.index[i]]

    relax.lambda.list.trunc[[i]] <-
      relax.lambda.list[[i]][relax.lambda.list[[i]] < relax.lambda.start.val]
    
    ## remove from relax.fit objects (if it's not an intercept only model)

    for(j in 1:nfolds) {

      if(!inherits(current.outlist[[j]],
                   "relaxnet.intercept.only")) {

        index <-
          current.outlist[[j]]$lambda < (relax.lambda.start.val - lambda.gap)
        ## figure out a better way to take care of the lambda.gap

        current.outlist[[j]]$beta <- current.outlist[[j]]$beta[, index,
                                                               drop = FALSE]

        ## also need to subset lambda and intercepts (a0)

        current.outlist[[j]]$lambda <- current.outlist[[j]]$lambda[index]
        current.outlist[[j]]$a0 <- current.outlist[[j]]$a0[index]
      }
    }

    relax.cvstuff.list[[i]] <-
      do.call(fun,
              list(current.outlist,
                   relax.lambda.list.trunc[[i]],
                   x,
                   y,
                   weights,
                   offset,
                   foldid,
                   type.measure,
                   grouped
                   ))[c("cvm", "cvsd")]
  }

  end.time <- Sys.time()

  lamin <- if(type.measure=="auc") {
             getmin(lambda,-cvm,cvsd)
           } else {
             getmin(lambda,cvm,cvsd)
           }

  relax.mins <- sapply(relax.cvstuff.list,
                       function(x) min(x$cvm, na.rm = TRUE))

  min.relax <- min(relax.mins)
  min.main <- min(cvm, na.rm = TRUE)

  if(min.main <= min.relax) {

    which.model.min <- "main"
    overall.lambda.min <- lamin$lambda.min
    min.cvm <- min.main
    
  } else {
  
    which.model.min <- which.min(relax.mins)
    lambda.min.index <-
      which.min(relax.cvstuff.list[[which.model.min]]$cvm)
    overall.lambda.min <-
      relax.lambda.list.trunc[[which.model.min]][lambda.min.index]
    min.cvm <- min.relax
  }

  if(multicore) {

    ## reset RNG kinds (seed will probably be reset)

    RNGkind(prev.RNG.kinds[1], prev.RNG.kinds[2])
  }
  
  obj <- list(call = match.call(),
              relax = relax,
              lambda=lambda,
              cvm=cvm,
              cvsd=cvsd,
              cvup=cvm+cvsd,
              cvlo=cvm-cvsd,
              nzero=nz,
              name=cvname,
              relaxnet.fit = relaxnet.fit,
              relax.cvstuff.list = relax.cvstuff.list,
              relax.lambda.list.trunc = relax.lambda.list.trunc,
              which.model.min = which.model.min,
              ## either "main" or a number corresponding to a relax model,
              overall.lambda.min = overall.lambda.min,
              min.cvm = min.cvm,

              ## maybe take these next three out
              main.lambda.min = lamin$lambda.min,
              main.lambda.1se = lamin$lambda.1se,
              main.min.cvm = min.main,

              total.time = as.double(difftime(end.time, start.time,
                                              units = "secs")),
              full.data.fit.time = relaxnet.fit$total.time,
              cv.fit.times = cv.fit.times)


  class(obj) <- "cv.relaxnet"

  obj
}
