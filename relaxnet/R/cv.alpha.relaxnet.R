################################################################################

## cross-validation (on both lambda and alpha) for relaxnet models
## adapted from the cv.glmnet function from package glmnet

cv.alpha.relaxnet <- function(x, y, family = c("gaussian", "binomial"),
                              nlambda = 100,
                              alpha = c(.1, .3, .5, .7, .9),

                              relax = TRUE,
                              relax.nlambda = 100,
                              relax.max.vars = min(nrow(x), ncol(x)) * 0.8,

                              lambda = NULL,
                              relax.lambda.index = NULL,
                              relax.lambda.list = NULL,

                              ##type.measure=c("mse","deviance","class","auc","mae"),
                              ## just set it in code for now
                        
                              nfolds = 10,
                              foldid,
                              
                              multicore = FALSE,
                              mc.cores,
                              mc.seed = 123,
                              
                              ...) {

  start.time <- Sys.time()

  family = match.arg(family)

  ## do something with this if including type.measure as an arg
  
  ## if(family == "gaussian") type.measure <- "mse"
  ## if(family == "binomial") type.measure <- "deviance"

  N=nrow(x)

  if(missing(foldid)) {

    ## check nfolds here

    foldid <- sample(rep(seq(nfolds), length=N))

  }## else {

    ## check foldid

##  }

  
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

    ## set mc.cores to length(alpha) if it's greater

    if(mc.cores > length(alpha)) mc.cores <- length(alpha)

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
  
  if(relax && ncol(x) == 1) {

    warning("x has only one column, setting relax to FALSE")
    relax <- FALSE
  }

  run.one.alpha.val <- function(alpha.val) {

    cv.relaxnet(x, y, family,
                nlambda,
                alpha.val,
                relax,
                relax.nlambda,
                relax.max.vars,
                lambda,
                relax.lambda.index,
                relax.lambda.list,
                foldid = foldid,
                ...)
  }
  
  if(multicore) {

    ## needed for reproducibility

    mc.reset.stream()

    ## start parallel jobs

    cv.relaxnet.results <- mclapply(alpha, run.one.alpha.val,
                                    mc.preschedule = TRUE,
                                    mc.set.seed = TRUE,
                                    mc.cores = mc.cores)
  } else {

    cv.relaxnet.results <- lapply(alpha, run.one.alpha.val)
  }

  alpha.min.index <- which.min(sapply(cv.relaxnet.results, function(cv.result) cv.result$min.cvm))

  if(multicore) {

    ## reset RNG kinds (seed will probably be reset)

    RNGkind(prev.RNG.kinds[1], prev.RNG.kinds[2])
  }
  
  end.time <- Sys.time()
  
  obj <- list(call = match.call(),
              relax = relax,
              alpha = alpha,
              cv.relaxnet.results = cv.relaxnet.results,
              which.alpha.min = alpha[alpha.min.index],
              total.time = as.double(difftime(end.time, start.time,
                                              units = "secs")))

  class(obj) <- "cv.alpha.relaxnet"

  obj
}
