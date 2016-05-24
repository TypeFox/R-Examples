# This file contains gof methods for ergm, btergm, mtergm, Siena, and other 
# network models. It also contains the compute.goflist function, which accepts 
# simulations and target networks and applies gof statistics to them.


setMethod("getformula", signature = className("btergm", "btergm"), 
    definition = function(x) x@formula)

setMethod("getformula", signature = className("mtergm", "btergm"), 
    definition = function(x) x@formula)

setMethod("getformula", signature = className("ergm", "ergm"), 
    definition = function(x) x$formula)


# function which reduces a statistic x nsim matrix and computes summary stats
# input: two matrices of a certain type of statistics (simulated and observed)
# goal: get rid of empty rows at the end, e.g., where dsp(37) or so is usually
# not observed; return value: object containing the summary statistics
reduce.matrix <- function(sim, obs) {
  
  numsim <- ncol(as.matrix(sim))
  numobs <- ncol(as.matrix(obs))
  xlabels <- rownames(obs)
  
  # if geodist statistic: put aside the last 'Inf' row
  if (is.null(rownames(sim)) || rownames(sim)[nrow(sim)] == "Inf") {
    geo <- TRUE
    inf.sim <- sim[nrow(sim), ]  # put aside this row for now and reuse later
    sim <- sim[-nrow(sim), ]
    if (class(obs) == "matrix") {
      inf.obs <- obs[nrow(obs), ]
      obs <- matrix(obs[-nrow(obs), ], ncol = numobs)
    } else {
      inf.obs <- obs[length(obs)]
      obs <- matrix(obs[-length(obs)], ncol = numobs)
    }
  } else {
    geo <- FALSE
  }
  
  # find first empty row for simulation matrix
  sim.rs <- rowSums(sim)
  sim.remove <- length(sim.rs)  # at which row index can removal start?
  for (i in (length(sim.rs) - 1):1) {
    if (sim.rs[i] == 0 && sim.remove == (i + 1)) {
      sim.remove <- i  # remember which is the first empty row
    }
  }
  
  if (class(obs) != "matrix") {  # one network is compared
    obs.remove <- length(obs)
    for (i in (length(obs) - 1):1) {
      if (obs[i] == 0 && obs.remove == (i + 1)) {
        obs.remove <- i
      }
    }
  } else {  # several networks are compared
    obs.rs <- rowSums(obs)
    obs.remove <- length(obs.rs)
    for (i in (length(obs.rs) - 1):1) {
      if (obs.rs[i] == 0 && obs.remove == (i + 1)) {
        obs.remove <- i
      }
    }
  }
  rem <- max(c(obs.remove, sim.remove), na.rm = TRUE)  # which one is longer?
  
  # remove unnecessary rows
  if (class(obs) != "matrix") {  # get rid of empty observations or rows of obs
    obs <- matrix(obs[-(rem:length(obs))], ncol = numobs)
  } else {
    obs <- matrix(obs[-(rem:nrow(obs)), ], ncol = numobs)
  }
  sim <- matrix(sim[-(rem:nrow(sim)), ], ncol = numsim)  # same for sim stats
  
  if (nrow(obs) < rem) {
    for (i in (nrow(obs) + 1):rem) {
      obs <- rbind(obs, rep(0, ncol(obs)))
    }
  }
  if (nrow(sim) < rem) {
    for (i in (nrow(sim) + 1):rem) {
      sim <- rbind(sim, rep(0, ncol(sim)))
    }
  }
  
  # for geodist, add Inf row again
  if (geo == TRUE) {
    sim <- rbind(sim, inf.sim)
    obs <- rbind(obs, inf.obs)
    rownames(sim) <- c(1:(nrow(sim) - 1), "Inf")
  }
  
  # create final object which will contain the raw simulations + the comparison
  reducedobject <- list()
  reducedobject$sim <- as.data.frame(sim)
  
  rownames(sim) <- NULL
  rownames(obs) <- NULL
  
  # compute means, p values, etc. and put them in a data frame
  x <- matrix()
  if (ncol(obs) == 1 || ncol(sim) == 1) {  # compute all the summary statistics
    x.obs <- obs
    x.mean <- apply(sim, 1, mean)
    x.min <- apply(sim, 1, min)
    x.max <- apply(sim, 1, max)
    x.median <- apply(sim, 1, median)
    zval <- (x.mean - x.obs) / sd(x.mean)
    pval <- 2 * pnorm(abs(zval), lower.tail = FALSE)
    x.pval <- pval
    x <- data.frame(x.obs, x.mean, x.median, x.min, x.max, x.pval)
    colnames(x) <- c("obs", "sim: mean", "median", "min", "max", "Pr(>z)")
  } else {  # for several target networks, compute distribution
    x.obs.mean <- apply(obs, 1, mean)
    x.obs.min <- apply(obs, 1, min)
    x.obs.max <- apply(obs, 1, max)
    x.obs.median <- apply(obs, 1, median)
    x.mean <- apply(sim, 1, mean)
    x.min <- apply(sim, 1, min)
    x.max <- apply(sim, 1, max)
    x.median <- apply(sim, 1, median)
    x.pval <- numeric()
    for (i in 1:nrow(sim)) {
      tryCatch(
        expr = {
          x.pval[i] <- t.test(obs[i, ], sim[i, ])$p.value  # compare group means
        }, 
        error = function(e) {
          x.pval[i] <- 1  # if both are 0, a t.test cannot be computed...
        }, 
        finally = {}
      )
      
      if (is.nan(x.pval[i])) {  # geodist contains "Inf"
        x.pval[i] <- 1
      }
    }
    x <- data.frame(x.obs.mean, x.obs.median, x.obs.min, x.obs.max, x.mean, 
        x.median, x.min, x.max, x.pval)
    colnames(x) <- c("obs: mean", "median", "min", "max", 
        "sim: mean", "median", "min", "max", "Pr(>z)")
  }
  if (geo == TRUE) {
    rownames(x) <- c(xlabels[c(1:(nrow(x) - 1))], "Inf")
  } else {
    rownames(x) <- xlabels[1:nrow(x)]
  }
  rownames(reducedobject$sim) <- rownames(x)
  
  reducedobject$comparison <- as.data.frame(x)

  return(reducedobject)
}


# evaluate statistics and create gof object
compute.goflist <- function(simulations, target, statistics, parallel = "no", 
    ncpus = 1, cl = NULL, verbose = TRUE) {
  
  # prepare parallel processing
  stop.parallel <- FALSE
  if (parallel[1] == "snow" && is.null(cl) && ncpus > 1) {
    cl <- makeCluster(ncpus)
    stop.parallel <- TRUE
  }
  
  # go through statistics functions and compute results
  goflist <- list()
  for (z in 1:length(statistics)) {
    args <- length(formals(statistics[[z]]))  # determine number of arguments
    if (args == 1) {  # either sim or obs
      tryCatch(
        expr = {
          label <- suppressMessages(attributes(statistics[[z]](
              simulations[[1]]))$label)
          if (verbose == TRUE) {
            message(paste("Processing statistic:", label))
          }
          
          if (parallel[1] == "no") {
            simulated <- suppressMessages(sapply(simulations, statistics[[z]]))
            observed <- suppressMessages(sapply(target, statistics[[z]]))
          } else if (parallel[1] == "multicore") {
            test <- suppressMessages(statistics[[z]](simulations[[1]]))
            if (class(test) == "numeric" && length(test) == 1) {
              simulated <- suppressMessages(unlist(mclapply(simulations, 
                  statistics[[z]], mc.cores = ncpus)))
              observed <- suppressMessages(unlist(mclapply(target, 
                  statistics[[z]], mc.cores = ncpus)))
            } else {  # no mcsapply available because different length vectors
              simulated <- suppressMessages(mclapply(simulations, 
                  statistics[[z]], mc.cores = ncpus))
              observed <- suppressMessages(mclapply(target, statistics[[z]], 
                  mc.cores = ncpus))
              max.length.sim <- max(sapply(simulated, length), na.rm = TRUE)
              max.length.obs <- max(sapply(observed, length), na.rm = TRUE)
              max.length <- max(max.length.sim, max.length.obs, na.rm = TRUE)
              simulated <- sapply(simulated, function(x) {
                c(x, rep(0, max.length - length(x)))
              })
              observed <- sapply(observed, function(x) {
                c(x, rep(0, max.length - length(x)))
              })
            }
          } else {
            clusterEvalQ(cl, library("ergm"))
            clusterEvalQ(cl, library("xergm.common"))
            simulated <- suppressMessages(parSapply(cl = cl, simulations, 
                statistics[[z]]))
            observed <- suppressMessages(parSapply(cl = cl, target, 
                statistics[[z]]))
          }
          
          # if simulations have different dimensions, convert list to matrix
          if (class(simulated) == "list") {
            lengths <- sapply(simulated, length)
            l <- max(lengths)
            index <- which(lengths == l)[1]
            rn <- names(simulated[[index]])
            simulated <- sapply(simulated, function(x) {
                c(x, rep(0, l - length(x)))
            })
            rownames(simulated) <- rn
          }
          if (class(observed) == "list") {
            lengths <- sapply(observed, length)
            l <- max(lengths)
            index <- which(lengths == l)[1]
            rn <- names(observed[[index]])
            observed <- sapply(observed, function(x) {
                c(x, rep(0, l - length(x)))
            })
            rownames(observed) <- rn
          }
          
          gofobject <- list()
          gofobject$label <- label
          if (class(simulated) == "matrix") {  # boxplot-type GOF
            reduced <- reduce.matrix(simulated, observed)
            gofobject$type <- "boxplot"
            gofobject$stats <- reduced$comparison
            gofobject$raw <- Matrix(as.matrix(reduced$sim))
            class(gofobject) <- "boxplot"
          } else if (class(simulated) == "numeric") {  # density-type GOF
            gofobject$type <- "univariate"
            gofobject$obs <- observed
            gofobject$sim <- simulated
            class(gofobject) <- "univariate"
          }
          goflist[[length(goflist) + 1]] <- gofobject
          names(goflist)[length(goflist)] <- label
        }, 
          error = function(e) {
          if (verbose == TRUE) {
            cat(paste("  Skipping statistic for the following reason:", e))
          }
        }, 
        finally = {}
      )
    } else if (args == 2) {  # sim and obs; ROCPR-type GOF
      tryCatch(
        expr = {
          label <- "Tie prediction"
          if (verbose == TRUE) {
            message(paste("Processing statistic:", label))
          }
          gofobject <- statistics[[z]](simulations, target)
          goflist[[length(goflist) + 1]] <- gofobject
          names(goflist)[length(goflist)] <- label
        }, 
        error = function(e) {
          cat(paste("  Skipping this statistic for the following reason:", 
              e))
        }, 
        finally = {}
      )
    }
  }
  class(goflist) <- "gof"
  if (stop.parallel == TRUE) {
    stopCluster(cl)
  }
  return(goflist)
}


# gof method for ergm and btergm objects
gof.btergm <- function(object, target = NULL, formula = getformula(object), 
    nsim = 100, MCMC.interval = 1000, MCMC.burnin = 10000, parallel = c("no", 
    "multicore", "snow"), ncpus = 1, cl = NULL, statistics = c(dsp, esp, deg, 
    ideg, geodesic, rocpr, walktrap.modularity), verbose = TRUE, ...) {
  
  if (nsim < 2) {
    stop("The 'nsim' argument must be greater than 1.")
  }
  
  if (is.function(statistics)) {
    statistics <- c(statistics)
  }
  
  # prepare parallel processing; translate options into statnet arguments
  if (is.null(ncpus) || ncpus == 0) {
    ncpus <- 1
  }
  if (!parallel[1] %in% c("no", "multicore", "snow")) {
    parallel <- "no"
    warning("'parallel' argument not recognized. Using 'no' instead.")
  }
  #statnet.stop.parallel <- FALSE
  #if (parallel[1] == "no") {
  #  statnet.parallel.type <- NULL
  #  statnet.parallel <- 0
  #} else if (parallel[1] == "multicore") {
  #  statnet.parallel.type <- "PSOCK"
  #  statnet.parallel <- makeForkCluster(nnodes = ncpus)
  #  statnet.stop.parallel <- TRUE
  #} else if (parallel[1] == "snow") {
  #  statnet.parallel.type <- "PSOCK"
  #  if (!is.null(cl)) {
  #    statnet.parallel <- cl
  #  } else {
  #    statnet.parallel <- ncpus
  #  }
  #}
  if (verbose == TRUE) {
    if (parallel[1] == "no") {
      parallel.msg <- "on a single computing core."
    } else if (parallel[1] == "multicore") {
      parallel.msg <- paste("using multicore forking on", ncpus, "cores.")
    } else if (parallel[1] == "snow") {
      parallel.msg <- paste("using parallel processing on", ncpus, "cores.")
    }
    message("\nStarting GOF assessment ", parallel.msg, "...")
  }
  
  # call tergmprepare and integrate results as a child environment in the chain
  if (class(object)[1] == "btergm") {
    env <- tergmprepare(formula = formula, offset = object@offset, 
        verbose = verbose)
    parent.env(env) <- environment()
    offset <- object@offset
  } else {
    env <- tergmprepare(formula = formula, offset = FALSE, verbose = FALSE)
    parent.env(env) <- environment()
    offset <- FALSE
  }
  
  # check and rearrange target network(s)
  if (is.null(target)) {
    if (verbose == TRUE) {
      message(paste("\nNo 'target' network(s) provided. Using networks on the",
          "left-hand side of the model formula as observed networks.\n"))
    }
    target <- env$networks
  } else if (class(target) == "network" || class(target) == "matrix") {
    target <- list(target)
    if (verbose == TRUE) {
      message("\nOne observed ('target') network was provided.\n")
    }
  } else if (class(target) == "list") {
    if (verbose == TRUE) {
      message(paste("\n", length(target), "observed ('target') networks were",
          "provided.\n"))
    }
  } else {
    stop("'target' must be a network, matrix, or list of matrices or networks.")
  }
  
  # extract coefficients from object
  if (class(object)[1] == "btergm" && offset == TRUE) {
    coefs <- c(coef(object), -Inf)  # -Inf for offset matrix
  } else {
    coefs <- coef(object)
  }
  
  # adjust formula at each step, and simulate networks
  sim <- list()
  degen <- list()
  for (index in 1:env$time.steps) {
    i <- index  # index 'i' is used in formula construction in 'env'!
    # simulations for statnet-style and rocpr GOF
    if (verbose == TRUE) {
      if ("btergm" %in% class(object) || "mtergm" %in% class(object)) {
        f.i <- gsub("\\[\\[i\\]\\]", paste0("[[", index, "]]"), 
            paste(deparse(env$form), collapse = ""))
        f.i <- gsub("\\s+", " ", f.i)
        if ("btergm" %in% class(object)) {
          f.i <- gsub("^networks", env$lhs.original, f.i)
        }
      } else if ("ergm" %in% class(object)) {
        f.i <- paste(deparse(formula), collapse = "")
        f.i <- gsub("\\s+", " ", f.i)
      } else {
        stop(paste("Unknown object type:", class(object)))
      }
      message(paste("Simulating", nsim, 
          "networks from the following formula:\n", f.i, "\n"))
    }
    # parallel processing in simulate.formula is SLOW!
    #sim[[index]] <- simulate.formula(env$form, nsim = nsim, coef = coefs, 
    #    constraints = ~ ., control = control.simulate.formula(MCMC.interval = 
    #    MCMC.interval, MCMC.burnin = MCMC.burnin, parallel = statnet.parallel, 
    #    parallel.type = statnet.parallel.type))
    sim[[index]] <- simulate.formula(env$form, nsim = nsim, coef = coefs, 
        constraints = ~ ., control = control.simulate.formula(MCMC.interval = 
        MCMC.interval, MCMC.burnin = MCMC.burnin))
  }
  
  # check basis network(s)
  if (verbose == TRUE) {
    if (env$time.steps == 1) {
      message("One network from which simulations are drawn was provided.\n")
    } else {
      message(paste(env$time.steps, "networks from which simulations are",
          "drawn were provided.\n"))
    }
  }
  
  # unpack nested lists of simulations and target statistics
  simulations <- list()
  for (i in 1:length(sim)) {
    for (j in 1:length(sim[[i]])) {
      simulations[[length(simulations) + 1]] <- sim[[i]][[j]]
    }
  }
  rm(sim)
  
  # if NA in target networks, put them in the base network, too, and vice-versa
  if (length(env$networks) == length(target)) {
    for (i in 1:env$time.steps) {
      env$networks[[i]] <- as.matrix(env$networks[[i]])
      target[[i]] <- as.matrix(target[[i]])
      if (nrow(target[[i]]) != nrow(env$networks[[i]])) {
        stop(paste0("Dimensions of observed network and target do not match ", 
            "at t = ", i, ": observed network has ", nrow(env$networks[[i]]), 
            " rows while target has ", nrow(target[[i]]), " rows."))
      }
      if (ncol(target[[i]]) != ncol(env$networks[[i]])) {
        stop(paste0("Dimensions of observed network and target do not match ", 
            "at t = ", i, ": observed network has ", ncol(env$networks[[i]]), 
            " columns while target has ", ncol(target[[i]]), " columns."))
      }
      env$networks[[i]][is.na(as.matrix(target[[i]]))] <- NA
      env$networks[[i]] <- network::network(env$networks[[i]], 
          directed = env$directed, bipartite = env$bipartite)
      target[[i]][is.na(as.matrix(env$networks[[i]]))] <- NA
      target[[i]] <- network::network(target[[i]], directed = env$directed, 
          bipartite = env$bipartite)
    }
  }
  
  # data preparation
  sptypes <- c("dgCMatrix", "dgTMatrix", "dsCMatrix", "dsTMatrix", "dgeMatrix")
  
  directed <- logical()
  twomode <- logical()
  for (i in 1:length(target)) {
    if (class(target[[i]]) == "network") {
      directed[i] <- is.directed(target[[i]])
      twomode[i] <- is.bipartite(target[[i]])
      target[[i]] <- Matrix(as.matrix(target[[i]]))
    } else if (class(target[[i]]) == "matrix") {
      directed[i] <- is.mat.directed(target[[i]])
      twomode[i] <- !is.mat.onemode(target[[i]])
      target[[i]] <- Matrix(target[[i]])
    } else if (class(target[[i]]) %in% sptypes) {
      # OK
      directed[i] <- is.mat.directed(target[[i]])
      twomode[i] <- !is.mat.onemode(target[[i]])
    }
  }
  simulations <- lapply(simulations, function(x) Matrix(as.matrix(x)))
  goflist <- compute.goflist(simulations = simulations, target = target, 
      statistics = statistics, parallel = parallel, ncpus = ncpus, cl = cl, 
      verbose = verbose)
  #if (statnet.stop.parallel == TRUE) {
  #  stopCluster(statnet.parallel)
  #}
  return(goflist)
}

# generic methods for goodness-of-fit assessment
setMethod("gof", signature = className("btergm", "btergm"), 
    definition = gof.btergm)

setMethod("gof", signature = className("ergm", "ergm"), 
    definition = gof.btergm)

setMethod("gof", signature = className("mtergm", "btergm"), 
    definition = gof.btergm)


# gof method for sienaFit objects
gof.sienaFit <- function(object, period = NULL, parallel = c("no", 
    "multicore", "snow"), ncpus = 1, cl = NULL, structzero = 10, 
    statistics = c(esp, deg, ideg, geodesic, rocpr, walktrap.modularity), 
    groupName = object$f$groupNames[[1]], varName = NULL, 
    outofsample = FALSE, sienaData = NULL, sienaEffects = NULL, 
    nsim = NULL, verbose = TRUE, ...) {
  
  # check correct specification of all data
  if (outofsample == FALSE) {
    if (is.null(varName)) {
      varName <- object$f$depNames[[1]]
    }
    if (is.null(object$sims) || any(sapply(object$sims, is.null))) {
      stop(paste("The object does not contain any simulations. Please", 
          "re-estimate the model again with argument 'returnDeps = TRUE'."))
    }
    if (verbose == TRUE) {
      if (!is.null(sienaData)) {
        message(paste("'sienaData' argument is ignored because", 
            "outofsample = FALSE."))
      }
      if (!is.null(sienaEffects)) {
        message(paste("'sienaEffects' argument is ignored because", 
            "outofsample = FALSE."))
      }
      if (!is.null(nsim)) {
        message(paste("'nsim' argument is ignored because", 
            "outofsample = FALSE."))
      }
    }
  } else {  # check if all data are there for out-of-sample prediction
    if (is.null(sienaData) || class(sienaData) != "siena") {
      stop(paste("For out-of-sample prediction, the 'sienaData'", 
          "argument should provide a 'siena' object (as created by the", 
          "'sienaDependent' function), and this object should contain two", 
          "networks: the base network to simulate from and the target network", 
          "to be predicted."))
    }
    if (sienaData$observations > 2) {
      stop(paste("For out-of-sample prediction, the 'sienaData' argument", 
          "should contain only two networks: the base network to simulate", 
          "from and the target network to be predicted."))
    }
    if (is.null(sienaEffects) || class(sienaEffects) != "sienaEffects") {
      stop(paste("For out-of-sample GOF assessment, a 'sienaEffects' object", 
          "must be provided. This object must contain the same effects that", 
          "were used in the original estimation, and it must be based on", 
          "the 'sienaData' object."))
    }
    if (is.null(varName)) {
      varName <- names(sienaData$depvars)[1]
    }
    if (is.null(nsim)) {
      if (verbose == TRUE) {
        message("'nsim' not given. Using a default of 100 simulations.")
      }
      nsim <- 100
    }
    if (verbose == TRUE) {
      if (!is.null(period)) {
        message("The 'period' argument is ignored because outofsample = TRUE.")
      }
    }
  }
  if (verbose == TRUE) {
    message(paste("Variable name:", varName))
    if (outofsample == FALSE) {
      message(paste("Group name:", groupName))
    }
  }
  
  # interpret period argument
  if (outofsample == FALSE) {
    if (is.null(period)) {
      base <- 1:(attr(object$f[[1]]$depvars[[1]], "netdims")[3] - 1)
    } else {
      if (period < 2) {
        stop("The 'period' argument should have a minimum value of 2.")
      } else if (period > object$observations) {
        stop(paste("The 'period' argument should not be larger than the number", 
            "of observed time periods."))
      }
      base <- period - 1
    }
  } else {
    base <- 1
  }
  
  # out of sample: extract, inject and fix coefficients and simulate networks
  if (outofsample == TRUE) {
    coefs <- object$theta[object$BasicRateFunction == FALSE]
    numtheta.new <- length(sienaEffects$initialValue[sienaEffects$include])
    sienaEffects$initialValue[sienaEffects$include][2:numtheta.new] <- coefs
    sim_model <- RSiena::sienaAlgorithmCreate(projname = "sim_model", 
        cond = FALSE, useStdInits = FALSE, nsub = 0 , n3 = nsim, simOnly = TRUE)
    sim_ans <- siena07(sim_model, data = sienaData, effects = sienaEffects, 
        batch = TRUE, silent = TRUE, returnDeps = TRUE)
    object <- sim_ans  # replace object by new simulations
    if (is.null(groupName)) {
      groupName <- names(sim_ans$f)[1]
      if (verbose == TRUE) {
        message(paste("Group name:", groupName))
      }
    }
  }
  
  # save the target object in a list and remove/handle structural zeros
  #dv <- eval(parse(text = varName))
  target <- list()
  for (i in base) {
    if (outofsample == TRUE) {
      dv.temp <- sienaData$depvars[[varName]][, , i + 1]
    } else {
      #dv.temp <- dv[, , i + 1]
      dv.temp <- object$f[[groupName]]$depvars[[varName]][, , i + 1]
    }
    rownames(dv.temp) <- 1:nrow(dv.temp)
    colnames(dv.temp) <- 1:ncol(dv.temp)
    dv.temp <- suppressMessages(handleMissings(dv.temp, na = structzero, 
        method = "remove", logical = FALSE))
    target[[length(target) + 1]] <- dv.temp
  }
  
  # define function for extracting all simulations and making them conformable
  simSiena <- function(i, myobject, mybase, mygroup, mydvname, target) {
    l <- list(length(target))
    count <- 0
    for (j in mybase) {
      count <- count + 1
      s <- RSiena::sparseMatrixExtraction(i = i, obsData = myobject$f, 
          sims = myobject$sims, period = j, groupName = mygroup, 
          varName = mydvname)
      rownames(s) <- 1:nrow(s)
      colnames(s) <- 1:ncol(s)
      l[[count]] <- Matrix(adjust(as.matrix(s), target[[count]]))
    }
    return(l)
  }
  
  # extract the simulations, possibly in parallel
  if (is.null(ncpus) || ncpus == 0) {
    ncpus <- 1
  }
  if (!parallel[1] %in% c("no", "multicore", "snow")) {
    parallel <- "no"
    warning("'parallel' argument not recognized. Using 'no' instead.")
  }
  if (parallel[1] == "snow") {
    if (is.null(cl)) {
      cl <- makeCluster(ncpus)
    }
    if (verbose == TRUE) {
      message(paste("Using snow parallelization with", ncpus, "cores."))
    }
    clusterEvalQ(cl, library("xergm.common"))
    clusterEvalQ(cl, library("Matrix"))
    sim <- parLapply(cl, 1:length(object$sims), simSiena, 
        myobject = object, mybase = base, mygroup = groupName, 
        mydvname = varName, target = target)
  } else if (parallel[1] == "multicore") {
    if (verbose == TRUE) {
      message(paste("Using multicore parallelization with", ncpus, "cores."))
    }
    sim <- mclapply(1:length(object$sims), simSiena, myobject = object, 
        mybase = base, mygroup = groupName, mydvname = varName, 
        target = target, mc.cores = ncpus)
  } else {
    if (verbose == TRUE) {
      message(paste0("Parallelization is switched off. Extracting ", 
          "simulations sequentially."))
    }
    sim <- lapply(1:length(object$sims), simSiena, myobject = object, 
        mybase = base, mygroup = groupName, mydvname = varName, target = target)
  }
  
  # unpack nested lists of simulations
  simulations <- list(length(base) * length(object$sims))
  count <- 1
  for (t in 1:length(base)) {
    for (i in 1:length(object$sims)) {
      simulations[[count]] <- sim[[i]][[t]]
      count <- count + 1
    }
  }
  rm(sim)
  
  # apply auxiliary functions and return list of comparisons
  goflist <- compute.goflist(simulations = simulations, target = target, 
      statistics = statistics, parallel = parallel, ncpus = ncpus, cl = cl, 
      verbose = verbose)
  return(goflist)
}

setMethod("gof", signature = className("sienaFit", "RSiena"), 
    definition = gof.sienaFit)


# gof method for dyadic-independence models with custom data and coefficients
gof.network <- function(object, covariates, coef, target = NULL, 
    nsim = 100, mcmc = FALSE, MCMC.interval = 1000, MCMC.burnin = 10000, 
    parallel = c("no", "multicore", "snow"), ncpus = 1, cl = NULL, 
    statistics = c(dsp, esp, deg, ideg, geodesic, rocpr, walktrap.modularity), 
    verbose = TRUE, ...) {
  
  if (nsim < 2) {
    stop("The 'nsim' argument must be greater than 1.")
  }
  
  # check dependent network
  nw <- object
  if (class(nw) == "network") {
    directed <- network::is.directed(nw)
    bipartite <- network::is.bipartite(nw)
  } else if (class(nw) == "matrix") {
    directed <- !isSymmetric(nw)
    if (nrow(nw) == ncol(nw)) {
      bipartite <- FALSE
    } else {
      bipartite <- TRUE
    }
    nw <- network::network(nw, bipartite = bipartite, directed = directed)
  } else {
    stop("'object' must be a network object or a matrix.")
  }
  time.steps <- 1
  num.vertices <- nrow(as.matrix(nw))
  
  # check and rearrange target network(s)
  if (is.null(target)) {
    if (verbose == TRUE) {
      message(paste("\nNo 'target' network(s) provided. Using networks on the",
          "left-hand side of the model formula as observed networks.\n"))
    }
    target <- nw
  } else if (class(target) == "network" || class(target) == "matrix") {
    # do nothing
    if (verbose == TRUE) {
      message("\nOne observed ('target') network was provided.\n")
    }
  } else if (class(target) == "list") {
    if (verbose == TRUE) {
      message(paste("\n", length(target), "observed ('target') networks were",
          "provided. Using the first network\n"))
    }
    target <- target[[1]]
    if (class(target) != "matrix" && class(target) != network) {
      stop("First target network was not a matrix or network object.")
    }
  } else {
    stop("'target' must be a network, matrix, or list of matrices or networks.")
  }
  
  # check predictors and assemble formula
  if (class(covariates) != "list") {
    stop("Covariates must be provided as a list of matrices.")
  }
  numcov <- length(covariates)
  if (numcov + 1 != length(coef)) {
    stop(paste("The 'coef' vector ought to have a coefficient for edges", 
        "plus the same number of coefficients as there are covariates.", 
        "Right now, there are", length(coef), "coefficients and", numcov, 
        "covariates."))
  }
  rhs <- "edges"
  for (i in 1:numcov) {
    if (!class(covariates[[i]]) %in% c("network", "matrix")) {
      stop(paste("Covariate", i, "is not a matrix or network object."))
    }
    if (nrow(as.matrix(covariates[[i]])) != nrow(as.matrix(nw))) {
      stop(paste("Number of row nodes of covariate", i, "is not compatible."))
    }
    if (ncol(as.matrix(covariates[[i]])) != ncol(as.matrix(nw))) {
      stop(paste("Number of column nodes of covariate", i, 
          "is not compatible."))
    }
    rhs <- paste(rhs, "+ edgecov(covariates[[", i, "]])")
  }
  form <- as.formula(paste("nw ~", rhs))
  
  # simulations for statnet-style and rocpr GOF
  message(paste("Simulating", nsim, 
      "networks from the following formula:\n", 
      gsub("\\s+", " ", paste(deparse(form), collapse = "")), "\n"))
  if (mcmc == TRUE) {  # TODO: IMPLEMENT PARALLEL PROCESSING HERE
    simulations <- simulate.formula(form, nsim = nsim, coef = coef, 
        constraints = ~ ., 
        control = control.simulate.formula(MCMC.interval = MCMC.interval, 
        MCMC.burnin = MCMC.burnin))
  } else {
    dat <- sapply(covariates, function(x) c(as.matrix(x)))
    dat <- cbind(rep(1, nrow(dat)), dat)
    prob <- plogis(coef %*% t(dat))
    simval <- t(sapply(prob, function(x) rbinom(nsim, 1, x)))
    simulations <- apply(simval, 2, function(x) Matrix(x, 
        nrow = num.vertices, byrow = FALSE))
  }
  
  # if NA in target networks, put them in the base network, too, and vice-versa
  nw <- as.matrix(nw)
  nw[is.na(as.matrix(target))] <- NA
  nw <- network::network(nw, directed = directed, bipartite = bipartite)
  target <- as.matrix(target)
  target[is.na(as.matrix(nw))] <- NA
  target <- list(Matrix(target))
  
  # data preparation
  sptypes <- c("dgCMatrix", "dgTMatrix", "dsCMatrix", "dsTMatrix", "dgeMatrix")
  
  directed <- logical()
  twomode <- logical()
  for (i in 1:length(target)) {
    if (class(target[[i]]) == "network") {
      directed[i] <- is.directed(target[[i]])
      twomode[i] <- is.bipartite(target[[i]])
      target[[i]] <- Matrix(as.matrix(target[[i]]))
    } else if (class(target[[i]]) == "matrix") {
      directed[i] <- is.mat.directed(target[[i]])
      twomode[i] <- !is.mat.onemode(target[[i]])
      target[[i]] <- Matrix(target[[i]])
    } else if (class(target[[i]]) %in% sptypes) {
      # OK
      directed[i] <- is.mat.directed(target[[i]])
      twomode[i] <- !is.mat.onemode(target[[i]])
    }
  }
  goflist <- compute.goflist(simulations = simulations, target = target, 
      statistics = statistics, parallel = parallel, ncpus = ncpus, cl = cl, 
      verbose = verbose)
  return(goflist)
}

setMethod("gof", signature = className("network", "network"), 
    definition = gof.network)

setMethod("gof", signature = className("matrix", "base"), 
    definition = gof.network)
