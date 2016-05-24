#  R package rjags file R/jags.R
#  Copyright (C) 2006-2013 Martyn Plummer
#
#  This program is free software; you can redistribute it and/or
#  modify it under the terms of the GNU General Public License version
#  2 as published by the Free Software Foundation.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/
#

.quiet.messages <- function(quiet)
{
    .Call("quietMessages", quiet, PACKAGE="rjags")
}

print.jags <- function(x, ...)
{
  cat("JAGS model:\n\n")

  model <- x$model()
  for (i in 1:length(model)) {
    cat(model[i],"\n",sep="")
  }

  data <- x$data()
  full <- !sapply(lapply(data, is.na), any)
  if (any(full)) {
    cat("Fully observed variables:\n", names(data)[full], "\n")
  }
  part <- !full & !sapply(lapply(data, is.na), all)
  if (any(part)) {
    cat("Partially observed variables:\n", names(data)[part], "\n")
  }
}

jags.model <- function(file, data=NULL, inits,
                       n.chains = 1, n.adapt=1000, quiet=FALSE)
{
    if (missing(file)) {
        stop("Model file name missing")
    }
    if (is.character(file)) {
      modfile <- file
      ## Check file exists and can be opened in text mode
      con <- try(file(modfile, "rt"))
      if (inherits(con, "try-error")) {
        stop(paste("Cannot open model file \"", modfile, "\"", sep=""))
      }
      close(con)
      ## Need this for print method and for recompile function
      model.code <- readLines(file, warn=FALSE)
    }
    else if (inherits(file, "connection")) {
        modfile <- tempfile()
        ## JAGS library requires a physical file, so we need to copy
        ## the contents of the connection to a temporary file
        model.code <- readLines(file, warn=FALSE)
        writeLines(model.code, modfile)
    }
    else {
        stop("'file' must be a character string or connection")
    }

    if (quiet) {
        .quiet.messages(TRUE)
        on.exit(.quiet.messages(FALSE), add=TRUE)
    }

    p <- .Call("make_console", PACKAGE="rjags")
    .Call("check_model", p, modfile, PACKAGE="rjags")
    if (!is.character(file)) {
        unlink(modfile) #Remove temporary copy
    }

    varnames <- .Call("get_variable_names", p, PACKAGE="rjags")
    if (missing(data) || is.null(data)) {
        data <- list()
    }
    else if (is.environment(data)) {
        ##Get a list of numeric objects from the supplied environment
        data <- mget(varnames, envir=data, mode="numeric",
                     ifnotfound=list(NULL))
        ##Strip null entries
        data <- data[!sapply(data, is.null)]
    }
    else if (is.list(data)) {
        v <- names(data)
        if (is.null(v) && length(v) != 0) {
            stop("data must be a named list")
        }
        if (any(nchar(v)==0)) {
            stop("unnamed variables in data list")
        }
        if (any(duplicated(v))) {
            stop("Duplicated names in data list: ",
                 paste(v[duplicated(v)], collapse=" "))
        }
        relevant.variables <- v %in% varnames
        data <- data[relevant.variables]
        unused.variables <- setdiff(v, varnames)
        for (i in seq(along=unused.variables)) {
            warning("Unused variable \"", unused.variables[i], "\" in data")
        }
        ### Check for data frames
        df <- which(as.logical(sapply(data, is.data.frame)))
        for (i in seq(along=df)) {
            if (all(sapply(data[[df[i]]], is.numeric))) {
                #Turn numeric data frames into matrices
                data[[df[i]]] <- as.matrix(data[[df[i]]])
            }
            else {
                stop("Data frame with non-numeric elements provided as data: ",
                     names(data)[df[i]])
            }
        }
    }
    else {
        stop("data must be a list or environment")
    }

    .Call("compile", p, data, as.integer(n.chains), TRUE, PACKAGE="rjags")

### Setting initial values

    if (!missing(inits) && !is.null(inits))  {

        checkParameters <- function(inits) {
            if(!is.list(inits)) {
                return("inits parameter must be a list")
            }

            inames <- names(inits)
            if (is.null(inames) || any(nchar(inames) == 0)) {
                return("No variable names supplied for the initial values")
            }

            dupinames <- duplicated(inames)
            if (any(dupinames)) {
                return(paste("Duplicated initial values for variable(s): ",
                             paste0(unique(inames[dupinames]), collapse = ", ")
                             )
                       )
            }

            if (any(inames==".RNG.name")) {
                rngname <- inits[[".RNG.name"]]
                if (!is.character(rngname) || length(rngname) != 1) {
                    return("Incorrect .RNG.name value")
                }
                inits[[".RNG.name"]] <- NULL
            }

            ## Strip null initial values, but give a warning
            null.inits <- sapply(inits, is.null)
            if (any(null.inits)) {
                warning(paste("NULL initial value supplied for variable(s) ",
                        paste(inames[null.inits], collapse=", "), sep=""))
                inits <- inits[!null.inits]
            }

            num_vals <- sapply(inits, is.numeric)
            if (any(!num_vals)) {
                return(paste("Non-numeric initial values supplied for variable(s) ",
                       paste(inames[!num_vals], collapse=", "), sep=""))
            }

            return ("ok")
        }

        setParameters <- function(inits, chain) {
            if (!is.null(inits[[".RNG.name"]])) {
                .Call("set_rng_name", p, inits[[".RNG.name"]],
                      as.integer(chain), PACKAGE="rjags")
                inits[[".RNG.name"]] <- NULL
            }
            .Call("set_parameters", p, inits, as.integer(chain),
                  PACKAGE="rjags")
        }

        init.values <- vector("list", n.chains)

        if (is.function(inits)) {
            if (any(names(formals(inits)) == "chain")) {
                for (i in 1:n.chains) {
                    init.values[[i]] <- inits(chain=i)
                }
            }
            else {
                for (i in 1:n.chains) {
                    init.values[[i]] <- inits()
                }
            }
        }
        else if (is.list(inits)) {

            if ( !is.null(names(inits)) ) {
                ## Replicate initial values for all chains
                for (i in 1:n.chains) {
                    init.values[[i]] <- inits
                }
            }
            else {
                if (length(inits) != n.chains) {
                    stop("Length mismatch between inits and n.chains")
                }
                init.values <- inits
            }
        }

        for (i in 1:n.chains) {
            msg  <- checkParameters(init.values[[i]])
            if (!identical(msg, "ok")) {
                stop("Invalid parameters for chain ", i, ":\n", msg);
            }
            setParameters(init.values[[i]], i)
            unused.inits <- setdiff(names(init.values[[i]]), varnames)
            unused.inits <- setdiff(unused.inits,
                                    c(".RNG.seed", ".RNG.state", ".RNG.name"))
            for (j in seq(along=unused.inits)) {
                warning("Unused initial value for \"", unused.inits[j],
                        "\" in chain ", i)
            }
        }
    }

    .Call("initialize", p, PACKAGE="rjags")

    model.state <- .Call("get_state", p, PACKAGE="rjags")
    model.data <- .Call("get_data", p, PACKAGE="rjags")
    model <- list("ptr" = function() {p},
                  "data" = function() {model.data},
                  "model" = function() {model.code},
                  "state" = function(internal=FALSE)
                  {
                      if(!internal) {
                          for(i in 1:n.chains) {
                              model.state[[i]][[".RNG.state"]] <- NULL
                              model.state[[i]][[".RNG.name"]] <- NULL
                          }
                      }
                      return(model.state)
                  },
                  "nchain" = function()
                  {
                      .Call("get_nchain", p, PACKAGE="rjags")
                  },
                  "iter" = function()
                  {
                      .Call("get_iter", p, PACKAGE="rjags")
                  },
                  "sync" = function() {

                      model.state <<- .Call("get_state", p, PACKAGE="rjags")
                  },
                  "recompile" = function() {
                      ## Clear the console
                      .Call("clear_console", p, PACKAGE="rjags")
                      p <<- .Call("make_console", PACKAGE="rjags")
                      ## Write the model to a temporary file so we can re-read it
                      mf <- tempfile()
                      writeLines(model.code, mf)
                      .Call("check_model", p, mf, PACKAGE="rjags")
                      unlink(mf)
                      ## Re-compile
                      .Call("compile", p, data, n.chains, FALSE, PACKAGE="rjags")
                      ## Re-initialize
                      if (!is.null(model.state)) {
                          if (length(model.state) != n.chains) {
                              stop("Incorrect number of chains in saved state")
                          }
                          for (i in 1:n.chains) {
                              statei <- model.state[[i]]
                              rng <- statei[[".RNG.name"]]
                              if (!is.null(rng)) {
                                  .Call("set_rng_name", p, rng, i, PACKAGE="rjags")
                                  statei[[".RNG.name"]] <- NULL
                              }
                              .Call("set_parameters", p, statei, i, PACKAGE="rjags")
                          }
                          .Call("initialize", p, PACKAGE="rjags")
                          ## Redo adaptation
                          adapting <- .Call("is_adapting", p, PACKAGE="rjags")
                          if(n.adapt > 0 && adapting) {
                              cat("Adapting\n")
                              .Call("update", p, n.adapt, PACKAGE="rjags")
                              if (!.Call("check_adaptation", p, PACKAGE="rjags")) {
                                  warning("Adaptation incomplete");
                              }
                          }
                          model.state <<- .Call("get_state", p, PACKAGE="rjags")
                      }
                      invisible(NULL)
                  })
    class(model) <- "jags"

    if (n.adapt > 0) {
        pb <- if(quiet) NULL else getOption("jags.pb")
        ok <- adapt(model, n.adapt, end.adaptation=FALSE, progress.bar=pb)
        if (ok) {
            .Call("adapt_off", p, PACKAGE="rjags")
        }
        else {
            warning("Adaptation incomplete")
        }
    }
    return(model)
}

parse.varname <- function(varname) {

  ## Try to parse string of form "a" or "a[n,p:q,r]" where "a" is a
  ## variable name and n,p,q,r are integers

  v <- try(parse(text=varname, n=1), silent=TRUE)
  if (!is.expression(v) || length(v) != 1)
    return(NULL)

  v <- v[[1]]
  if (is.name(v)) {
    ##Full node array requested
    return(list(name=deparse(v)))
  }
  else if (is.call(v) && identical(deparse(v[[1]]), "[") && length(v) > 2) {
    ##Subset requested
    ndim <- length(v) - 2
    lower <- upper <- numeric(ndim)
    if (any(nchar(sapply(v, deparse)) == 0)) {
      ##We have to catch empty indices here or they will cause trouble
      ##below
      return(NULL)
    }
    for (i in 1:ndim) {
      index <- v[[i+2]]
      if (is.numeric(index)) {
        ##Single index
        lower[i] <- upper[i] <- index
      }
      else if (is.call(index) && length(index) == 3 &&
               identical(deparse(index[[1]]), ":") &&
               is.numeric(index[[2]]) && is.numeric(index[[3]]))
        {
          ##Index range
          lower[i] <- index[[2]]
          upper[i] <- index[[3]]
        }
      else return(NULL)
    }
    if (any(upper < lower))
      return (NULL)
    return(list(name = deparse(v[[2]]), lower=lower, upper=upper))
  }
  return(NULL)
}

parse.varnames <- function(varnames)
{
  names <- character(length(varnames))
  lower <- upper <- vector("list", length(varnames))
  for (i in seq(along=varnames)) {
    y <- parse.varname(varnames[i])
    if (is.null(y)) {
      stop(paste("Invalid variable subset", varnames[i]))
    }
    names[i] <- y$name
    if (!is.null(y$lower)) {
      lower[[i]] <- y$lower
    }
    if (!is.null(y$upper)) {
      upper[[i]] <- y$upper
    }
  }
  return(list(names=names, lower=lower, upper=upper))
}


jags.samples <-
  function(model, variable.names, n.iter, thin=1, type="trace", ...)
{
    if (class(model) != "jags")
      stop("Invalid JAGS model")

    if (!is.character(variable.names) || length(variable.names) == 0)
      stop("variable.names must be a character vector")

    if (!is.numeric(n.iter) || length(n.iter) != 1 || n.iter <= 0)
      stop("n.iter must be a positive integer")
    if (!is.character(type))
      stop("type must be a character vector")

    pn <- parse.varnames(variable.names)
    status <- .Call("set_monitors", model$ptr(), pn$names, pn$lower, pn$upper,
                    as.integer(thin), type, PACKAGE="rjags")
    if (!any(status)) stop("No valid monitors set")
    update.jags(model, n.iter, ...)
    ans <- .Call("get_monitored_values", model$ptr(), type, PACKAGE="rjags")
    for (i in seq(along=ans)) {
        class(ans[[i]]) <- "mcarray"
        attr(ans[[i]], "varname") <- names(ans)[i]
    }
    for (i in seq(along=variable.names)) {
        if (status[i]) {
            .Call("clear_monitor", model$ptr(), pn$names[i], pn$lower[[i]],
                  pn$upper[[i]], type, PACKAGE="rjags")
        }
    }
    return(ans)
}

list.samplers <- function(object)
{
    if (!inherits(object, "jags")) {
        stop("not a jags model object")
    }
    .Call("get_samplers", object$ptr(), PACKAGE="rjags")
}

list.factories <- function(type)
{
    type = match.arg(type, c("sampler","monitor","rng"))
    as.data.frame(.Call("get_factories", type, PACKAGE="rjags"))
}

set.factory <- function(name, type, state)
{
    if (!is.character(name) || length(name) != 1)
        stop("invalid name")
    if (!is.character(type) || length(type) != 1)
        stop("invalid name")
    if (length(state) != 1)
        stop("invalid state")

    type <- match.arg(type, c("sampler","rng","monitor"))
    .Call("set_factory_active", name, type, as.logical(state), PACKAGE="rjags")
}

coda.names <- function(basename, dim)
{
    ## Utility function used to get the names of the individual elements
    ## of a node array

    if (prod(dim) == 1)
      return(basename)

    ##Default lower and upper limits
    ndim <- length(dim)
    lower <- rep(1, ndim)
    upper <- dim

    ##If the node name is a subset, we try to parse it to get the
    ##names of its elements. For example, if basename is "A[2:3]"
    ##we want to return names "A[2]", "A[3]" not "A[2:3][1]", "A[2:3][2]".
    pn <- parse.varname(basename)
    if (!is.null(pn) && !is.null(pn$lower) && !is.null(pn$upper)) {
        if (length(pn$lower) == length(pn$upper)) {
            dim2 <- pn$upper - pn$lower + 1
            if (isTRUE(all.equal(dim[dim!=1], dim2[dim2!=1],
                                 check.attributes=FALSE))) {
                basename <- pn$name
                lower <- pn$lower
                upper <- pn$upper
                ndim <- length(dim2)
            }
        }
    }

    indices <- as.character(lower[1]:upper[1])
    if (ndim > 1) {
        for (i in 2:ndim) {
            indices <- outer(indices, lower[i]:upper[i], FUN=paste, sep=",")
        }
    }
    paste(basename,"[",as.vector(indices),"]",sep="")
}

nchain <- function(model)
{
    if (!inherits(model, "jags"))
      stop("Invalid JAGS model object in nchain")

    .Call("get_nchain", model$ptr(), PACKAGE="rjags")
}

coda.samples <- function(model, variable.names=NULL, n.iter, thin=1,
                         na.rm = TRUE, ...)
{
    start <- model$iter() + thin
    out <- jags.samples(model, variable.names, n.iter, thin, type="trace", ...)

    ans <- vector("list", nchain(model))
    for (ch in 1:nchain(model)) {
        ans.ch <- vector("list", length(out))

        vnames.ch <- NULL
        for (i in seq(along=out)) {

            varname <- names(out)[[i]]
            d <- dim(out[[i]])
            if (length(d) < 3) {
                stop("Invalid dimensions for sampled output")
            }
            vardim <- d[1:(length(d)-2)]
            nvar <- prod(vardim)
            niter <- d[length(d) - 1]
            nchain <- d[length(d)]

            values <- as.vector(out[[i]])
            var.i <- matrix(NA, nrow=niter, ncol=nvar)
            for (j in 1:nvar) {
                var.i[,j] <- values[j + (0:(niter-1))*nvar + (ch-1)*niter*nvar]
            }
            vnames.ch <- c(vnames.ch, coda.names(varname, vardim))
            ans.ch[[i]] <- var.i
        }

        ans.ch <- do.call("cbind", ans.ch)
        colnames(ans.ch) <- vnames.ch
        ans[[ch]] <- mcmc(ans.ch, start=start, thin=thin)
    }

    if (isTRUE(na.rm)) {
        ## Drop variables that are missing for all iterations in at least
        ## one chain
        all.missing <- sapply(ans, function(x) {apply(is.na(x), 2, any)})
        drop.vars <- if (is.matrix(all.missing)) {
            apply(all.missing, 1, any)
        }
        else {
            any(all.missing)
        }
        ans <- lapply(ans, function(x) return(x[, !drop.vars, drop=FALSE]))
    }

    mcmc.list(ans)
}

load.module <- function(name, path, quiet=FALSE)
{
    if (name %in% list.modules()) {
        ## This is a stop-gap measure as JAGS 2.1.0 does allow you
        ## to load the same module twice. This should be fixed in
        ## later versions.
        return(invisible()) #Module already loaded
    }

    if (missing(path)) {
        path = getOption("jags.moddir")
        if (is.null(path)) {
            stop("option jags.moddir is not set")
        }
    }
    if (!is.character(path) || length(path) != 1)
        stop("invalid path")
    if (!is.character(name) || length(name) != 1)
        stop("invalid name")

    file <- file.path(path, paste(name, .Platform$dynlib.ext, sep=""))
    if (!file.exists(file)) {
        stop("File not found: ", file)
    }
    if (!isDLLLoaded(file)) {
        ## We must avoid calling dyn.load twice on the same DLL This
        ## may result in the DLL being unloaded and then reloaded,
        ## which will invalidate pointers to the distributions,
        ## functions and factories in the module.
        dyn.load(file)
    }
    ok <- .Call("load_module", name, PACKAGE="rjags")
    if (!ok) {
        stop(paste("module", name, "not found\n"))
    }
    else if (!quiet) {
        message("module ", name, " loaded")
    }
    invisible()
}

unload.module <- function(name, quiet=FALSE)
{
    if (!is.character(name) || length(name) != 1)
        stop("invalid name")

    ok <- .Call("unload_module", name, PACKAGE="rjags")
    if (!ok) {
        warning(paste("module", name, "not loaded"))
    }
    else if (!quiet) {
        cat("Module", name, "unloaded\n", sep=" ")
    }
    invisible()
}

list.modules <- function()
{
    .Call("get_modules", PACKAGE="rjags");
}

isDLLLoaded <- function(file)
{
    dll.list <- getLoadedDLLs()
    for (i in seq(along=dll.list)) {
        if (dll.list[[i]]["path"][1] == file)
            return(TRUE)
    }
    return(FALSE)
}

parallel.seeds <- function(factory, nchain)
{
    .Call("parallel_seeds", factory, nchain, PACKAGE="rjags")
}
