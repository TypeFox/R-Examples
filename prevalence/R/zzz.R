##= Define S4 classes =====================================================
setOldClass("prevModel")  # virtual S3 class
setOldClass("mcmc.list")  # virtual S3 class
setClass("prev",
  representation(
    par = "list",
    model = "prevModel",
    mcmc = "list",
    diagnostics = "list"))

##= Define S4 methods =====================================================
setMethod("show", "prev",
  function(object)
    print(object)
)

setMethod("print", "prev",
  function(x, conf.level = 0.95, dig = 3, ...) {
    ## guess which function generated 'x'
    multi <- is.null(x@par$SE)

    ## get summary statistics from 'summary()'
    stats <- summary(x, conf.level)
    if (!is.list(stats)) stats <- list(stats)

    summary_row <- x@par$nchains + 1
    out <- t(sapply(stats, function(x) x[summary_row, c(1:4, 6:7)]))

    if (multi) {
      h <- log(length(x@par$x), 2)
      method <-
        ifelse(x@par$prior[[1]][[1]] == "TP", "covariance", "conditional")

      rownames(out) <-
        switch(method,
               conditional = c(" TP",
                               paste0(rep(c("SE", "SP"), times = h),
                                      rep(seq(h), each = 2))),
               covariance = get_nodes(h))
    }

    ## get BGR statistic
    BGR <- x@diagnostics$BGR

    ## if multinomial, get bayesP
    if (multi) bayesP <- x@diagnostics$bayesP

    ## print 'out' dataframe
    print(round(out, dig), ...)

    ## print diagnostic information
    # if only one node: mpsrf == NULL
    if (is.null(BGR$mpsrf)) {
      cat("\nBGR statistic = ", round(BGR[[1]][, 1], 4),
          " (upper CL = ", round(BGR[[1]][, 2], 4), ")\n", sep = "")
      cat("BGR values substantially above 1 indicate lack of convergence\n")

    # if multiple nodes
    } else {
      cat("\nMultivariate BGR statistic = ",
          round(BGR$mpsrf, 4), "\n", sep = "")
      cat("BGR values substantially above 1 indicate lack of convergence\n")
    }

    if (multi && method == "conditional") {
      cat("Bayes-P statistic =", round(bayesP, 2), "\n")
      cat("Bayes-P values substantially different from 0.5",
          "indicate lack of convergence\n")
    }
  }
)

setMethod("summary", "prev",
  function(object, conf.level = 0.95) {
    ## derive lower and upper confidence level
    if (sum(object@par$x) == 0) {
      p <- c(0, conf.level)

    } else if (ifelse(length(object@par$x) == 1,
                 object@par$x == object@par$n,
                 sum(object@par$x) == length(object@par$x))) {
      p <- c(1 - conf.level, 1)

    } else {
      p <- c((1 - conf.level) / 2,
              1 - (1 - conf.level) / 2)
    }
    ciLabel <- paste0(100 * p, "%")

    ## guess which function generated 'object'
    multi <- is.null(object@par$SE)

    if (multi) {
      nodes <- names(object@mcmc)[-length(names(object@mcmc))]

    } else {
      nodes <- names(object@mcmc)
    }

    stat_list <- vector("list", length(nodes))
    names(stat_list) <- nodes

    for (node in seq_along(nodes)) {
      ## define 'stats' matrix
      n <- object@par$nchains
      stats <- matrix(ncol = 8, nrow = n + 1)
      colnames(stats) <- list("mean", "median", "mode", "sd", "var",
                              ciLabel[1], ciLabel[2], "samples")
      dimnames(stats)[[1]] <- c(paste(rep("chain", n), seq(n)), "all chains")

      ## extract mcmc samples for this node
      mcmc <- object@mcmc[[node]]

      ## calculate summary statistics per chain
      for (i in seq(object@par$nchains)) {
        stats[i, 1] <- mean(mcmc[[i]], na.rm = TRUE)
        stats[i, 2] <- median(mcmc[[i]], na.rm = TRUE)

        if (var(mcmc[[i]]) > 0) {
          d <- density(mcmc[[i]], na.rm = TRUE)
          stats[i, 3] <- d$x[which.max(d$y)]
        } else {
          stats[i, 3] <- mcmc[[i]][1]
        }

        stats[i, 4] <- sd(mcmc[[i]], na.rm = TRUE)
        stats[i, 5] <- var(mcmc[[i]], na.rm = TRUE)
        stats[i, 6] <- quantile(mcmc[[i]], probs = p[1], na.rm = TRUE)
        stats[i, 7] <- quantile(mcmc[[i]], probs = p[2], na.rm = TRUE)
        stats[i, 8] <- length(mcmc[[i]])
      }

      ## calculate overall summary statistics
      y <- unlist(mcmc)
      i <- i + 1

      stats[i, 1] <- mean(y, na.rm = TRUE)
      stats[i, 2] <- median(y, na.rm = TRUE)

      if (var(y) > 0) {
        d <- density(y, na.rm = TRUE)
        stats[i, 3] <- d$x[which.max(d$y)]
      } else {
        stats[i, 3] <- y[1]
      }

      stats[i, 4] <- sd(y, na.rm = TRUE)
      stats[i, 5] <- var(y, na.rm = TRUE)
      stats[i, 6] <- quantile(y, probs = p[1], na.rm = TRUE)
      stats[i, 7] <- quantile(y, probs = p[2], na.rm = TRUE)
      stats[i, 8] <- length(y)

      stat_list[[node]] <- stats
    }

    ## return resulting 'stat' list
    return(stat_list)
  }
)

setMethod("plot", "prev",
  function(x, y = NULL, ...) {
    ## define 'y' if missing
    if (missing(y)) y <- "TP"

    ## define 'ask'
    ask_old <- par("ask")
    ask_new <- ifelse(prod(par("mfrow")) == 4, FALSE, TRUE)
    devAskNewPage(ask_new)
    on.exit(devAskNewPage(ask_old))

    ## guess which function generated 'x'
    multi <- is.null(x@par$SE)
    if (multi) {
      h <- log2(length(x@par$x))

      if (length(x@mcmc) == 1 + length(get_nodes(h))) {
        choices <- get_nodes(h)
      } else {
      choices <-
        c("TP", paste0(rep(c("SE", "SP"), each = h), seq(h)))
      }

      y <- match.arg(y, choices)
      mcmc <- x@mcmc[[y]]

    } else {
      choices <- c("TP", "SE", "SP")
      y <- match.arg(y, choices)
      mcmc <- x@mcmc[[y]]
    }

    ## 4 plots
    of_y <- ifelse(ask_new, paste("of", y), "")
    cex.main <- ifelse(ask_new, 1.2, 1)
    line <- ifelse(ask_new, 1.5, 1)

    densplot(mcmc, main = "", ylab = "density", ask = FALSE)
    title(main = paste("Density", of_y),
          cex.main = cex.main, line = line)
    traceplot(mcmc, main = "", ylab = "prevalence", ask = FALSE)
    title(main = paste("Trace", of_y),
          cex.main = cex.main, line = line)
    gelman.plot(mcmc, ask = TRUE, auto.layout = FALSE)
    title(main = expression(symbol("\250")),
          col.main = "white", cex.main = 5)
    title(main = paste("BGR plot", of_y),
          cex.main = cex.main, line = line)
    autocorr.plot(mcmc[[1]], ask = TRUE, auto.layout = FALSE)
    title(main = expression(symbol("\250")),
          col.main = "white", cex.main = 5)
    title(main = paste("Autocorrelation", of_y),
          cex.main = cex.main, line = line)

    if (!ask_new)
      title(y, outer = TRUE, line = -1.5)
  }
)

setMethod("as.matrix", "prev",
  function(x, iters = FALSE, chains = FALSE) {
    ## convert MCMC to matrix
    mx <- sapply(x@mcmc, unlist)

    ## add iteration numbers
    if (iters)
      mx <- cbind(ITER = rep(seq(x@par$update), x@par$nchains), mx)

    ## add chain numbers
    if (chains)
      mx <- cbind(CHAIN = rep(seq(x@par$nchains), each = x@par$update), mx)

    ## return matrix
    return(mx)
  }
)
