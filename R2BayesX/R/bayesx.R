bayesx <- function(formula, data, weights = NULL, subset = NULL, 
  offset = NULL, na.action = NULL, contrasts = NULL, 
  control = bayesx.control(...), model = TRUE,
  chains = NULL, cores = NULL, ...)
{
  ## multiple core processing
  if(!is.null(cores)) {
    setseed <- round(runif(cores) * .Machine$integer.max)
    control$dir.rm <- if(is.null(control$outfile)) TRUE else control$dir.rm
    bayesx_parallel <- function(j) {
      control$setseed <- setseed[j]
      bayesx(formula, data, weights, subset, offset, na.action,
        contrasts, control, model, chains = chains, cores = NULL)
    }
    rval <- parallel::mclapply(1:cores, bayesx_parallel, mc.cores = cores)
    if(inherits(rval[[1]][[1]], "bayesx")) {
      k <- length(rval[[1]])
      rval2 <- list(); mn <- NULL
      for(i in 1:cores) {
        mn <- c(mn, paste("Core", i,
          if(k > 1) paste("Chain", 1:k, sep = "_") else NULL,
          sep = "_"))
        rval2 <- c(rval2, rval[[i]])
      }
      names(rval2) <- mn
      class(rval2) <- "bayesx"
      return(rval2)
    } else {
      names(rval) <- paste("Core", 1:cores, sep = "_")
      class(rval) <- "bayesx"
      return(rval)
    }
  }

  ## multiple chain processing
  if(!is.null(chains)) {
    if(!is.null(control$setseed))
      set.seed(control$setseed)
    setseed <- round(runif(chains) * .Machine$integer.max)
    outfile <- if(nout <- is.null(control$outfile)) {
      file.path(tempfile(), paste("Chain", 1:chains, control$model.name, sep = "_"))
    } else {
      if(length(control$outfile) < 2)
        file.path(path.expand(control$outfile), paste("Chain", 1:chains, control$model.name, sep = "_"))
      else
        path.expand(control$outfile)
    }
    if(length(unique(outfile)) != chains)
      stop(paste("there must be", chains, "direcories supplied within outfile for parallel computing!"))
    bayesx_chain <- function(j) {
      control$setseed <- setseed[j]
      control$outfile <- outfile[j]
      control$dir.rm <- if(nout) TRUE else control$dir.rm
      bayesx(formula, data, weights, subset, offset, na.action,
        contrasts, control, model, chains = NULL, cores = NULL)
    }
    rval <- list()
    for(j in 1:chains) {
      rval[[j]] <- bayesx_chain(j)
    }
    names(rval) <- paste("Chain", 1:chains, sep = "_")
    class(rval) <- "bayesx"
    return(rval)
  }

  ## special handling of control arguments
  args <- list(...)
  if(is.function(args$family))
    args$family <- args$family()$family
  if(!is.null(args$family) || !is.null(args$method)) {
    if((!is.null(args$family) && control$family != args$family) || 
       (!is.null(args$method) && control$method != args$method))
      control <- do.call("bayesx.control", args)
  }    

  res <- list()
  res$formula <- formula
  weights <- deparse(substitute(weights), backtick = TRUE, width.cutoff = 500L)
  offset <- deparse(substitute(offset), backtick = TRUE, width.cutoff = 500L)
  subset <- deparse(substitute(subset), backtick = TRUE, width.cutoff = 500L)

  ## setup files for bayesx
  res$bayesx.setup <- parse.bayesx.input(formula, data,
    weights, subset, offset, na.action, contrasts, control)

  ## write prg file
  res$bayesx.prg <- write.bayesx.input(res$bayesx.setup)

  ## model.frame
  if(!model) {
    res$bayesx.setup$data <- NULL
    res$bayesx.setup <- rmbhmf(res$bayesx.setup)
  }

  ## now estimate with BayesX
  res$bayesx.run <- run.bayesx(file.path(res$bayesx.prg$file.dir, 
    prg.name = res$bayesx.prg$prg.name), verbose = res$bayesx.setup$verbose)

  if(is.null(res$bayesx.setup$hlevel))
    tm <- res$bayesx.prg$model.name
  else
    tm <- s4hm(res$bayesx.prg$file.dir, control$model.name)

  ## read bayesx output files
  if(control$read) {
    if(!grepl("ERROR:", res$bayesx.run$log[length(res$bayesx.run$log)]) || 
      res$bayesx.run$log[1L] == 0) {
      if(length(tm) > 1L) {
        res.h <- read.bayesx.output(res$bayesx.prg$file.dir, tm)
        for(k in 1:length(res.h)) {
          res.h[[k]]$model.fit$method = "HMCMC"
          res.h[[k]]$model.fit$model.name <- tm[k] 
        }
        res.h <- res.h[length(res.h):1L]
        res.h[[1L]]$call <- match.call()
        res.h[[1L]][names(res)] <- res[names(res)]
        class(res.h) <- "bayesx"
      } else {
        res.h <- NULL
        res <- c(res, read.bayesx.model.output(res$bayesx.prg$file.dir, tm))
        res$call <- match.call()
        res$terms <- res$bayesx.setup$term.labels
        res$model.fit$family <- control$family
      }
    } else warning("an error occured during runtime of BayesX!")
  } else {
    cat("All output files are stored in:", res$bayesx.prg$file.dir, "\n")
    return(invisible(res$bayesx.prg$file.dir))
  }

  ## remove output folder
  if((is.null(res$bayesx.setup$outfile)) && res$bayesx.setup$dir.rm) {
    wd <- getwd()
    try(Sys.chmod(res$bayesx.prg$file.dir, mode = "7777"), silent = TRUE)
    setwd(res$bayesx.prg$file.dir)
    files <- list.files()
    for(k in 1L:length(files)) {
      ok <- try(suppressWarnings(file.remove(files[k])), silent = TRUE)
      if(inherits(ok, "try-error"))
        try(suppressWarnings(system(paste("rm --force", files[k]))), silent = TRUE)
    }
    ok <- try(suppressWarnings(file.remove(res$bayesx.prg$file.dir)), silent = TRUE)
    if(inherits(ok, "try-error"))
      try(suppressWarnings(system(paste("rm --force", res$bayesx.prg$file.dir))), silent = TRUE)
    setwd(wd)
  }

  if(!is.null(res.h)) {
    res.h <- res.h[sort(names(res.h))]
    class(res.h) <- "bayesx"
    return(res.h)
  } else {
    class(res) <- "bayesx"
    return(res)
  }
}


rmbhmf <- function(x) {
  if(!is.null(x$h.random)) {
    for(k in 1L:length(x$h.random))
      x$h.random[[k]]$data <- NULL
    x$h.random[[k]] <- rmbhmf(x$h.random[[k]])
  }
  return(x)
}
