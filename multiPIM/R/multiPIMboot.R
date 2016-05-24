################################################################################
################################################################################
################################ multiPIMboot.R ################################
################################################################################
################################################################################

## Author: Stephan Ritter

## Function definition for the multiPIMboot function.
## multiPIM is run once on the original data, then bootstrapped.
## multicore option allows for parallel processing over the bootstrap replicates

multiPIMboot <- function(Y, A, W = NULL,
                         times = 5000,
                         id = 1:nrow(Y),
                         multicore = FALSE,
                         mc.num.jobs,
                         mc.seed = 123,
                         estimator = c("TMLE", "DR-IPCW", "IPCW", "G-COMP"),
                         g.method = "main.terms.logistic", g.sl.cands = NULL,
                         g.num.folds = NULL, g.num.splits = NULL,
                         Q.method = "sl", Q.sl.cands = "default",
                         Q.num.folds = 5, Q.num.splits = 1,
                         Q.type = NULL,
                         adjust.for.other.As = TRUE,
                         truncate = 0.05,
                         return.final.models = TRUE,
                         na.action,
                         verbose = FALSE,
                         extra.cands = NULL,
                         standardize = TRUE,
                         ...) {

  if( !( (mode(times) == "numeric") && (length(times) == 1) &&
        !is.na(times) && ((times %% 1) == 0) && (times >= 2)) )
    stop("times must be a single integer greater than or equal to 2")

  if(length(id) != nrow(Y))
    stop("id must have length equal to nrow(Y)")

  if(!(identical(multicore, FALSE) || identical(multicore, TRUE)))
    stop("argument multicore must be either TRUE or FALSE")
  
  if(!(identical(verbose, FALSE) || identical(verbose, TRUE)))
    stop("argument verbose must be either TRUE or FALSE")

  if(multicore) {

    if(!requireNamespace("parallel", quietly = TRUE)) {

      stop("multicore is TRUE, but unable to load package parallel")

    }

    ## check mc.num.jobs

    if(missing(mc.num.jobs))
      stop("if multicore = TRUE, mc.num.jobs must.be.specified")

    if(mode(mc.num.jobs) != "numeric" || length(mc.num.jobs) != 1
       || mc.num.jobs < 1 || mc.num.jobs %% 1 != 0)
      stop("mc.num.jobs should be a single integer giving the number of\n",
           "cores/CPUs to be used")

    ## set mc.num.jobs to times if it's greater

    if(mc.num.jobs > times) mc.num.jobs <- times
    
    ## find out number of bootstrap samples per job
    ## (the final job may have fewer samples than this)

    samples.per.job <- (times + mc.num.jobs - 1) %/% mc.num.jobs

    ## find out how many jobs really need to be run (this is important for
    ## preventing errors when mc.num.jobs is not much less than times)

    mc.num.jobs <- times %/% samples.per.job

    if(times %% samples.per.job != 0) mc.num.jobs <- mc.num.jobs + 1

    ## if final job would have zero samples, subtract 1 from mc.num.jobs

    if( (mc.num.jobs - 1) * samples.per.job == times)
      mc.num.jobs = mc.num.jobs - 1

    ## check mc.seed

    if(mode(mc.seed) != "numeric"
       || length(mc.seed) != 1
       || (mc.seed %% 1) != 0)
      stop("mc.seed must be a single integer")

    ## set RNG kind to l'ecuyer, saving current RNG kinds
    
    prev.RNG.kinds <- RNGkind("L'Ecuyer-CMRG")

    ## set the seed

    set.seed(mc.seed)
    
  } else { ## multicore is false

    ## do everything in one job
    
    mc.num.jobs <- 1
    samples.per.job <- times
  }
    
  if(verbose) cat("Starting Main Run\n")

  main.run <- multiPIM(Y, A, W,
                       estimator,
                       g.method, g.sl.cands,
                       g.num.folds, g.num.splits,
                       Q.method, Q.sl.cands,
                       Q.num.folds, Q.num.splits,
                       Q.type,
                       adjust.for.other.As,
                       truncate,
                       return.final.models,
                       na.action,
                       check.input = TRUE,
                       verbose = FALSE,
                       extra.cands,
                       standardize)

  ## Prepare for getting bootstrap samples

  unique.ids <- unique(id)

  num.ids <- length(unique.ids)

  id.index.list <- split(1:length(id), id)

  bootstr.distr <- array(0, dim = c(times, dim(main.run$param.estimates)))

  bootstrapped.W <- NULL
  
  ## set Q.type so that there will be no checking of Y
  ## to determine which type to use

  Q.type <- main.run$Q.type

  run.one.job <- function(job.num) {

    ## find out how many samples this job need to run

    num.samples <- ifelse(job.num != mc.num.jobs, samples.per.job,
                          (times - (job.num - 1) * samples.per.job))

    ## instantiate an array to store results for this job

    job.result <- array(0, dim = c(num.samples, dim(main.run$param.estimates)))

    ## loop over samples

    for(i in 1:num.samples) {

      unique.id.index.sample <- sample(1:num.ids, num.ids, replace = T)

      boot.index.vec <- unlist(id.index.list[unique.id.index.sample],
                               use.names = FALSE)

      if(!is.null(W)) bootstrapped.W <- W[boot.index.vec, , drop = FALSE]

      if(verbose) {

        if(multicore) {
          cat("Starting bootstrap run number", i, "of", num.samples, "\n",
              "for job number", job.num, "\n")
        } else {
          cat("Starting bootstrap run number", i, "of", num.samples, "\n")
        }
      }

      job.result[i, , ] <- multiPIM(Y[boot.index.vec, , drop = FALSE],
                                    A[boot.index.vec, , drop = FALSE],
                                    bootstrapped.W,
                                    estimator,
                                    g.method, g.sl.cands,
                                    g.num.folds, g.num.splits,
                                    Q.method, Q.sl.cands,
                                    Q.num.folds, Q.num.splits,
                                    Q.type,
                                    adjust.for.other.As,
                                    truncate,
                                    return.final.models = FALSE,
                                    na.action,
                                    check.input = FALSE,
                                    verbose = FALSE,
                                    extra.cands,
                                    standardize)$param.estimates
    }

    job.result
  }

  if(multicore) {

    ## this is needed for reproducibility
    
    parallel::mc.reset.stream()

    ## start parallel jobs, then collect results

    jobs <- lapply(1:mc.num.jobs, function(x) parallel::mcparallel(run.one.job(x),
                                                       name = x))

    results.list <- parallel::mccollect(jobs)

    for(job.num in 1:mc.num.jobs) {

      begin.index <- samples.per.job * (job.num - 1) + 1
      end.index <- ifelse(job.num != mc.num.jobs, samples.per.job * job.num,
                          times)
      bootstr.distr[begin.index:end.index, , ] <- results.list[[job.num]]
    }
  } else {

    bootstr.distr <- run.one.job(1)

  }

  if(multicore) {

    ## reset RNG kinds (seed will probably be reset)

    RNGkind(prev.RNG.kinds[1], prev.RNG.kinds[2])
  }
    
  dimnames(bootstr.distr) <- c(sample.number = list(1:times),
                               dimnames(main.run$param.estimates))

  main.run$call <- match.call()
  main.run$boot.param.array <- bootstr.distr

  return(main.run)

} ## end multiPIMboot function
