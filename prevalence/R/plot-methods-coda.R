## to do: gelman.plot does not par mfrow

setMethod("densplot", "prev",
  function(x, exclude_fixed = TRUE, ...) {
    ## check inputs
    checkInput(exclude_fixed, "exclude_fixed", class = "logical")

    ## guess which function generated 'x'
    multi <- is.null(x@par$SE)

    ## calculate number of plots
    if (multi) {
      if (exclude_fixed) {
        is_fixed <- sapply(x@mcmc, function(x) var(unlist(x)) == 0)
        is_fixed <- head(is_fixed, -1)
        n <- length(x@mcmc) - sum(is_fixed) - 1
        N <- which(!is_fixed)
      } else {
        n <- length(x@mcmc) - 1
        N <- seq(n)
      }

    } else {
      if (exclude_fixed) {
        is_fixed <- sapply(x@mcmc, function(x) var(unlist(x)) == 0)
        N <- which(!is_fixed)
        n <- length(N)

      } else {
        n <- length(x@mcmc)
        N <- seq(n)
      }
    }

    ## define 'ask'
    ask_old <- par("ask")
    ask_new <- prod(par("mfrow")) < n
    devAskNewPage(ask_new)
    on.exit(devAskNewPage(ask_old))

    ## density plots
    for (i in N)
      densplot(x@mcmc[[i]],
               main = paste("Density of", names(x@mcmc)[i]),
               ask = FALSE, ...)
  }
)

setMethod("traceplot", "prev",
  function(x, exclude_fixed = TRUE, ...) {
    ## check inputs
    checkInput(exclude_fixed, "exclude_fixed", class = "logical")

    ## guess which function generated 'x'
    multi <- is.null(x@par$SE)

    ## calculate number of plots
    if (multi) {
      if (exclude_fixed) {
        is_fixed <- sapply(x@mcmc, function(x) var(unlist(x)) == 0)
        is_fixed <- head(is_fixed, -1)
        n <- length(x@mcmc) - sum(is_fixed) - 1
        N <- which(!is_fixed)
      } else {
        n <- length(x@mcmc) - 1
        N <- seq(n)
      }

    } else {
      if (exclude_fixed) {
        is_fixed <- sapply(x@mcmc, function(x) var(unlist(x)) == 0)
        N <- which(!is_fixed)
        n <- length(N)

      } else {
        n <- length(x@mcmc)
        N <- seq(n)
      }
    }

    ## define 'ask'
    ask_old <- par("ask")
    ask_new <- prod(par("mfrow")) < n
    devAskNewPage(ask_new)
    on.exit(devAskNewPage(ask_old))

    ## trace plots
    for (i in N)
      traceplot(x@mcmc[[i]],
                main = paste("Trace of", names(x@mcmc)[i]),
                ask = FALSE, ...)
  }
)

setMethod("gelman.plot", "prev",
  function(x, ...) {
    ## guess which function generated 'x'
    multi <- is.null(x@par$SE)

    ## calculate number of plots
    if (multi) {
      is_fixed <- sapply(x@mcmc, function(x) var(unlist(x)) == 0)
      is_fixed <- head(is_fixed, -1)
      n <- length(x@mcmc) - sum(is_fixed) - 1
      N <- which(!is_fixed)

    } else {
      is_fixed <- sapply(x@mcmc, function(x) var(unlist(x)) == 0)
      N <- which(!is_fixed)
      n <- length(N)
    }

    ## define 'ask'
    ask_old <- par("ask")
    ask_new <- prod(par("mfrow")) < n
    devAskNewPage(ask_new)
    on.exit(devAskNewPage(ask_old))

    ## gelman plots
    for (i in N)
      gelman.plot(x@mcmc[[i]],
                  main = paste("BGR plot of", names(x@mcmc)[i]),
                  ask = TRUE, auto.layout = FALSE, ...)
  }
)

setMethod("autocorr.plot", "prev",
  function(x, exclude_fixed = TRUE, chain = 1, ...) {
    ## check inputs
    checkInput(exclude_fixed, "exclude_fixed", class = "logical")
    checkInput(chain, "chain", class = "integer", min = 1)

    ## check number of chains
    if (chain > length(x@mcmc$TP))
      stop(paste("'x' only has", length(x@mcmc$TP), "chains"))

    ## guess which function generated 'x'
    multi <- is.null(x@par$SE)

    ## calculate number of plots
    if (multi) {
      if (exclude_fixed) {
        is_fixed <- sapply(x@mcmc, function(x) var(unlist(x)) == 0)
        is_fixed <- head(is_fixed, -1)
        n <- length(x@mcmc) - sum(is_fixed) - 1
        N <- which(!is_fixed)
      } else {
        n <- length(x@mcmc) - 1
        N <- seq(n)
      }

    } else {
      if (exclude_fixed) {
        is_fixed <- sapply(x@mcmc, function(x) var(unlist(x)) == 0)
        N <- which(!is_fixed)
        n <- length(N)

      } else {
        n <- length(x@mcmc)
        N <- seq(n)
      }
    }

    ## define 'ask'
    ask_old <- par("ask")
    ask_new <- prod(par("mfrow")) < n
    devAskNewPage(ask_new)
    on.exit(devAskNewPage(ask_old))

    ## autocorrelation plots
    for (i in N)
      autocorr.plot(x@mcmc[[i]][[chain]],
                    main = paste("Autocorrelation of", names(x@mcmc)[i]),
                    ask = TRUE, auto.layout = FALSE, ...)
  }
)