## meant to be called from within widenet


runCVrelaxnetOnAlphaVal <- function(alphaVal, x, y, family, foldid,
                                    multicore, mc.cores, mc.seed,
                                    ...) {

  cv.relaxnet(x, y, family,
              alpha = alphaVal,
              foldid = foldid,
              multicore = multicore,
              mc.cores = mc.cores,
              mc.seed = mc.seed,
              ...)
}


runRelaxnetOnAlphaVal <- function(alphaVal, x, y, family, ...) {

  relaxnet(x, y, family, alpha = alphaVal, ...)
}


runRelaxnetOnLambdaStuff <- function(lambda.stuff, x, y, family, ...) {

  relaxnet(x, y, family,
           alpha = lambda.stuff$alpha,
           lambda = lambda.stuff$lambda,
           relax.lambda.index = lambda.stuff$relax.lambda.index,
           relax.lambda.list = lambda.stuff$relax.lambda.list,
           ...)
}


## widenet function, main function from widenet package

widenet <- function(x, y, family = c("gaussian", "binomial"),
                    order = 1:3,
                    alpha = 1,

                    nfolds = 10,
                    foldid,

                    screen.method = c("none", "cor", "ttest"),
                    screen.num.vars = 50,

                    ## maybe add this later
                    ## screen.glmnet.alpha = 1, ## ignored if num.vars missing
                    ## screen.glmnet.nlambda = 100, ## ignored if num.vars missing

                    multicore = FALSE,
                    mc.cores,
                    mc.seed = 123,
                    ...) {

  start.time <- Sys.time()

  family <- match.arg(family)

  screen.method <- match.arg(screen.method)

  if(any(alpha <= 0) || any(alpha > 1))
    stop("all alpha values must be > 0 and <= 1")

  if(!all(order %in% 1:3))
    stop("order must have values in c(1,2,3)")

  order <- names(table(order))
  
  ## check for binary columns if any order is > 1

  colsBinary <- NULL
  
  if(any(order > 1)) {

    colsBinary <- check.binary(x)

    if(any(colsBinary == 1))
      stop("x has at least one column which contains a single repeated value")
  }

  ## deal with nfolds and foldid

  if(missing(foldid)) {

    ## check nfolds here
    
    foldid <- sample(rep(seq(nfolds), length = nrow(x)))

  } else {

    ## check foldid here
    
    nfolds <- length(table(foldid))
  }


  ## multicore stuff -- copied from cv.relaxnet function

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


  ## see if we need to screen, if not just run cv.relaxnet for each order
  ## if yes, need to screen separately within cross-validation folds
  
  if(screen.method == "none") {

    if(any(order > 1)) {

      numBinary <- sum(colsBinary == 2)
      numNotBinary <- sum(colsBinary == 3)

      if( (numBinary + numNotBinary) != length(colsBinary) )
        stop("Internal Error #1")
    }

    relaxnet.fit.list <- NULL

    ## results.list <- vector("list", length = length(order))
    ## names(results.list) <- paste("order", order, sep = "")


    if(1 %in% order) {

      order1fits <- lapply(alpha, runCVrelaxnetOnAlphaVal,
                           x, y, family,
                           foldid = foldid,
                           multicore = multicore,
                           mc.cores = mc.cores,
                           mc.seed = mc.seed,
                           ...)

      names(order1fits) <- paste("alpha", alpha, sep = "")
    
      relaxnet.fit.list <- c(relaxnet.fit.list, order1 = list(order1fits))
    }

    if(any(c(2, 3) %in% order)) {

      x2 <- expand.to.order.2(x, colsBinary, numBinary, numNotBinary)

      ## remove any columns resulting from interactions between
      ## dummmy vars from same factor
      
      x2.colsBinary <- check.binary(x2)

      x2.full <- cbind(x, x2[, x2.colsBinary != 1, drop = FALSE])

      if(2 %in% order) {

        order2fits <- lapply(alpha, runCVrelaxnetOnAlphaVal,
                             x2.full, y, family,
                             foldid = foldid,
                             multicore = multicore,
                             mc.cores = mc.cores,
                             mc.seed = mc.seed,
                             ...)

        names(order2fits) <- paste("alpha", alpha, sep = "")

        relaxnet.fit.list <- c(relaxnet.fit.list, order2 = list(order2fits))
      }
    }

    if(3 %in% order) {

      x3 <- expand.to.order.3(x, x2, colsBinary, numBinary, numNotBinary)

      x3.colsBinary <- check.binary(x3)
      
      x3.full <- cbind(x2.full, x3[, x3.colsBinary != 1, drop = FALSE])

################################################################################
##      rm(x3)  ## check if I'll need it later  
################################################################################

      order3fits <- lapply(alpha, runCVrelaxnetOnAlphaVal,
                           x3.full, y, family,
                           foldid = foldid,
                           multicore = multicore,
                           mc.cores = mc.cores,
                           mc.seed = mc.seed,
                           ...)

      names(order3fits) <- paste("alpha", alpha, sep = "")
    
      relaxnet.fit.list <- c(relaxnet.fit.list, order3 = list(order3fits))
    }

    min.cvm.mat <- sapply(relaxnet.fit.list,
                          function(x) sapply(x, function(y) y$min.cvm))

    dim(min.cvm.mat) <- c(length(alpha), length(order))

    dimnames(min.cvm.mat) <- list(alpha = alpha, order = order)

    win.indices <- which(min.cvm.mat == min(min.cvm.mat), arr.ind = TRUE)

    alpha.min <- alpha[win.indices[1]]

    order.min <- order[win.indices[2]]

    end.time <- Sys.time()
    
    obj <- list(call = match.call(),
                order = order,
                alpha = alpha,
                screen.method = screen.method,
                screened.in.index = NA,
                colsBinary = colsBinary,
                cv.relaxnet.results = relaxnet.fit.list,
                min.cvm.mat = min.cvm.mat,
                which.order.min = order.min,
                which.alpha.min = alpha.min,
                total.time = as.double(difftime(end.time, start.time, units = "secs")))

    class(obj) <- "widenet"

    return(obj)


  }

################################################################################

  ## now do screening part

  if(screen.num.vars < 2)
      stop("screen.num.vars must be >= 2")

  fun <- paste("screen.", screen.method, sep = "")

  screened.in.index <- do.call(fun, list(x = x,
                                         y = y,
                                         family = family,
                                         num.vars = screen.num.vars))

  ## if any order is greater than 1, check for binary columns

  if(any(order > 1)) {

    colsBinary.screened <- colsBinary[screened.in.index]

    numBinary.screened <- sum(colsBinary.screened == 2)
    numNotBinary.screened <- sum(colsBinary.screened == 3)

    if( (numBinary.screened + numNotBinary.screened)
        != length(colsBinary.screened) )
      stop("Internal Error #2")
  }

  relaxnet.fit.list <- NULL
  
  ## results.list <- vector("list", length = length(order))
  ## names(results.list) <- paste("order", order, sep = "")


  x.screened <- x[, screened.in.index, drop = FALSE]

  if(1 %in% order) {

    order1fits <- lapply(alpha, runRelaxnetOnAlphaVal,
                         x.screened, y, family, ...)

    names(order1fits) <- paste("alpha", alpha, sep = "")

    lambda.stuff.order1 <-
      mapply(function(fit, alpha.val) {

        return(list(alpha = alpha.val,
                    lambda = fit$main.glmnet.fit$lambda,
                    relax.lambda.index = fit$relax.lambda.index,
                    relax.lambda.list = lapply(fit$relax.glmnet.fits,
                                               function(fit) fit$lambda)))},
        order1fits, alpha, SIMPLIFY = FALSE)

    relaxnet.fit.list <- c(relaxnet.fit.list, order1 = list(order1fits))
  }

  if(any(c(2, 3) %in% order)) {

    x2.screened <- expand.to.order.2(x.screened, colsBinary.screened,
                            numBinary.screened, numNotBinary.screened)

    x2.colsBinary.screened <- check.binary(x2.screened)

    x2.screened.full <- cbind(x.screened,
                              x2.screened[, x2.colsBinary.screened != 1,
                                          drop = FALSE])

    if(2 %in% order) {

      order2fits <- lapply(alpha, runRelaxnetOnAlphaVal,
                           x2.screened.full, y, family, ...)

      names(order2fits) <- paste("alpha", alpha, sep = "")

      lambda.stuff.order2 <-
        mapply(function(fit, alpha.val) {

          return(list(alpha = alpha.val,
                      lambda = fit$main.glmnet.fit$lambda,
                      relax.lambda.index = fit$relax.lambda.index,
                      relax.lambda.list = lapply(fit$relax.glmnet.fits,
                        function(fit) fit$lambda)))},
               order2fits, alpha, SIMPLIFY = FALSE)

      relaxnet.fit.list <- c(relaxnet.fit.list, order2 = list(order2fits))
    }
  }

  if(2 %in% order && !(3 %in% order))
    rm(x2.screened, x2.screened.full)

  if(3 %in% order) {

    x3.screened <- expand.to.order.3(x.screened, x2.screened,
                                     colsBinary.screened,
                                     numBinary.screened, numNotBinary.screened)

    x3.colsBinary.screened <- check.binary(x3.screened)

    x3.screened.full <- cbind(x2.screened.full,
                              x3.screened[, x3.colsBinary.screened != 1,
                                          drop = FALSE])

    rm(x3.screened, x2.screened, x2.screened.full)


    order3fits <- lapply(alpha, runRelaxnetOnAlphaVal,
                         x3.screened.full, y, family, ...)

    rm(x3.screened.full)
    
    names(order3fits) <- paste("alpha", alpha, sep = "")

    lambda.stuff.order3 <-
      mapply(function(fit, alpha.val) {

        return(list(alpha = alpha.val,
                    lambda = fit$main.glmnet.fit$lambda,
                    relax.lambda.index = fit$relax.lambda.index,
                    relax.lambda.list = lapply(fit$relax.glmnet.fits,
                      function(fit) fit$lambda)))},
             order3fits, alpha, SIMPLIFY = FALSE)
    
    relaxnet.fit.list <- c(relaxnet.fit.list, order3 = list(order3fits))
  }

  run.one.fold <- function(which) { ## which is logical index
    ## identifying the current fold

    ## subset y
  
    if(is.matrix(y)){
      y <- y[!which,]
    } else {
      y <- y[!which]
    }

    ## subset x
  
    x <- x[!which, , drop = FALSE]
  
    ## do the screening
  
    fun <- paste("screen.", screen.method, sep = "")
    
    screened.in.index <- do.call(fun, list(x = x,
                                           y = y,
                                           family = family,
                                           num.vars = screen.num.vars))

    ## if any order is greater than 1, check for binary columns

    if(any(order > 1)) {

      colsBinary.screened <- colsBinary[screened.in.index]

      numBinary.screened <- sum(colsBinary.screened == 2)
      numNotBinary.screened <- sum(colsBinary.screened == 3)

      if( (numBinary.screened + numNotBinary.screened)
         != length(colsBinary.screened) )
        stop("Internal Error #2")
    }
  
    relaxnet.fit.list <- NULL
  
    ## results.list <- vector("list", length = length(order))
    ## names(results.list) <- paste("order", order, sep = "")


    x.screened <- x[, screened.in.index, drop = FALSE]

    if(1 %in% order) {

      order1fits <- lapply(lambda.stuff.order1, runRelaxnetOnLambdaStuff,
                           x.screened, y, family, ...)

      names(order1fits) <- paste("alpha", alpha, sep = "")
    
      relaxnet.fit.list <- c(relaxnet.fit.list, order1 = list(order1fits))
    }

    if(any(c(2, 3) %in% order)) {

      x2.screened <- expand.to.order.2(x.screened, colsBinary.screened,
                                       numBinary.screened, numNotBinary.screened)

      x2.colsBinary.screened <- check.binary(x2.screened)

      x2.screened.full <- cbind(x.screened,
                                x2.screened[, x2.colsBinary.screened != 1,
                                            drop = FALSE])

      if(2 %in% order) {

        order2fits <- lapply(lambda.stuff.order2, runRelaxnetOnLambdaStuff,
                             x2.screened.full, y, family, ...)

        names(order2fits) <- paste("alpha", alpha, sep = "")
    
        relaxnet.fit.list <- c(relaxnet.fit.list, order2 = list(order2fits))
      }
    }

    if(2 %in% order && !(3 %in% order))
      rm(x2.screened, x2.screened.full)

    if(3 %in% order) {

      x3.screened <- expand.to.order.3(x.screened, x2.screened,
                                       colsBinary.screened,
                                       numBinary.screened, numNotBinary.screened)

      x3.colsBinary.screened <- check.binary(x3.screened)
    
      x3.screened.full <- cbind(x2.screened.full,
                                x3.screened[, x3.colsBinary.screened != 1,
                                            drop = FALSE])

      rm(x3.screened, x2.screened, x2.screened.full)

      order3fits <- lapply(lambda.stuff.order3, runRelaxnetOnLambdaStuff,
                           x3.screened.full, y, family, ...)

      names(order3fits) <- paste("alpha", alpha, sep = "")
    
      relaxnet.fit.list <- c(relaxnet.fit.list, order3 = list(order3fits))
    }

    return(list(relaxnet.fit.list = relaxnet.fit.list,
                screened.in.index = screened.in.index))
  }

  ## adapted from cv.relaxnet function
  
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


  ## What to do depends on the type.measure and the model fit

  fun1 <- paste(".widenet.cv", class(relaxnet.fit.list[[1]][[1]]$main.glmnet.fit)[[1]],
               sep=".")

  fun2 <- paste(".relax", fun1, sep = "")

  lambda.list <- NULL
  
  if(1 %in% order) lambda.list <- c(lambda.list,
                                    list(lambda.stuff.order1))

  if(2 %in% order) lambda.list <- c(lambda.list,
                                    list(lambda.stuff.order2))

  if(3 %in% order) lambda.list <- c(lambda.list,
                                    list(lambda.stuff.order3))


  lambda.gap <- 10e-10 ## necessary since the lambda values are sometimes
  ## perturbed slightly,
  

  ## set type.measure

  if(family == "binomial") {

    type.measure <- "deviance"

  } else if(family == "gaussian") {

    type.measure <- "mse"

  }

  cv.relaxnet.results <- vector("list", length = length(order))
  names(cv.relaxnet.results) <- paste("order", order, sep = "")

  for(i in 1:length(order)) {
    cv.relaxnet.results[[i]] <- vector("list", length = length(alpha))
    names(cv.relaxnet.results[[i]]) <- paste("alpha", alpha, sep = "")
  }

  for(ii in 1:length(order)) {

    for(jj in 1:length(alpha)) {

      current.outlist <- lapply(fold.results, function(fold.result) {
        
        fold.result$relaxnet.fit.list[[ii]][[jj]]$main.glmnet.fit})

      screened.in.indices <- lapply(fold.results, function(fold.result) {

        fold.result$screened.in.index})
          
      ## just do this for now
      offset <- NULL
      weights <- rep(1, nrow(x))
      grouped <- TRUE

      cvstuff <- do.call(fun1, list(current.outlist,
                                    lambda.list[[ii]][[jj]]$lambda,
                                    x,
                                    y,
                                    weights,
                                    offset,
                                    foldid,
                                    type.measure,
                                    grouped,
                                    order = order[ii],
                                    colsBinary = colsBinary,
                                    screened.in.indices = screened.in.indices))

      cvm=cvstuff$cvm
      cvsd=cvstuff$cvsd
      cvname=cvstuff$name

      relax.lambda.index <- lambda.list[[ii]][[jj]]$relax.lambda.index

      relax.num.models <- length(relax.lambda.index)

      relax.cvstuff.list <- vector("list", length = relax.num.models)

      relax.lambda.list.trunc <- vector("list", length = relax.num.models)

      for(i in 1:relax.num.models) {

        current.outlist <-
          lapply(fold.results, function(fold.result) {

            fold.result$relaxnet.fit.list[[ii]][[jj]]$relax.glmnet.fits[[i]]})

        ## before processing, remove duplicate lambda values which were already
        ## done in the main model -- remove from the fit objects and the lambda.list

        relax.lambda.start.val <-
          lambda.list[[ii]][[jj]]$lambda[relax.lambda.index[i]]

        relax.lambda.list.trunc[[i]] <-
          lambda.list[[ii]][[jj]]$relax.lambda.list[[i]][
                    lambda.list[[ii]][[jj]]$relax.lambda.list[[i]]
                                            < relax.lambda.start.val]
    
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
          do.call(fun2,
                  list(current.outlist,
                       relax.lambda.list.trunc[[i]],
                       x,
                       y,
                       weights,
                       offset,
                       foldid,
                       type.measure,
                       grouped,
                       order = order[ii],
                       colsBinary = colsBinary,
                       screened.in.indices = screened.in.indices
                       ))[c("cvm", "cvsd")]

      }

      lambda <- lambda.list[[ii]][[jj]]$lambda
      
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

      nz <- sapply(predict(relaxnet.fit.list[[ii]][[jj]]$main.glmnet.fit,
                           type="nonzero"),
                   length)
      
      cv.relaxnet.results[[ii]][[jj]] <-
        list(call = "this object generated by widenet function",
             relax = relaxnet.fit.list[[ii]][[jj]]$relax,
             lambda=lambda,
             cvm=cvm,
             cvsd=cvsd,
             cvup=cvm+cvsd,
             cvlo=cvm-cvsd,
             nzero=nz,
             name=cvname,
             relaxnet.fit = relaxnet.fit.list[[ii]][[jj]],
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
             
             total.time = NA,
             full.data.fit.time = NA,
             cv.fit.times = NA)

      class(cv.relaxnet.results[[ii]][[jj]]) <- "cv.relaxnet"

    }
  }
   
  if(multicore) {

    ## reset RNG kinds (seed will probably be reset)

    RNGkind(prev.RNG.kinds[1], prev.RNG.kinds[2])
  }

  min.cvm.mat <- sapply(cv.relaxnet.results,
                        function(x) sapply(x, function(y) y$min.cvm))

  dim(min.cvm.mat) <- c(length(alpha), length(order))

  dimnames(min.cvm.mat) <- list(alpha = alpha, order = order)

  win.indices <- which(min.cvm.mat == min(min.cvm.mat), arr.ind = TRUE)

  alpha.min <- alpha[win.indices[1]]

  order.min <- order[win.indices[2]]

  end.time <- Sys.time()

  obj <- list(call = match.call(),
              order = order,
              alpha = alpha,
              screen.method = screen.method,
              screened.in.index = screened.in.index,
              colsBinary = colsBinary,
              cv.relaxnet.results = cv.relaxnet.results,
              min.cvm.mat = min.cvm.mat,
              which.order.min = order.min,
              which.alpha.min = alpha.min,
              total.time = as.double(difftime(end.time, start.time, units = "secs")))

  class(obj) <- "widenet"

  return(obj)
}
