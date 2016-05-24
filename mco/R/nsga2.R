##
## nsga2.R - Interface to nsga2.c
##
## Authors:
##  Heike Trautmann  <trautmann@statistik.uni-dortmund.de>
##  Detlef Steuer    <detlef.steuer@hsu-hamburg.de>
##  Olaf Mersmann    <olafm@statistik.uni-dortmund.de>
##

nsga2 <- function(fn, idim, odim, ...,
                  constraints=NULL, cdim=0,
                  lower.bounds=rep(-Inf, idim),
                  upper.bounds=rep(Inf, idim),
                  popsize=100, generations=100,
                  cprob=0.7, cdist=5,
                  mprob=0.2, mdist=10,
                  vectorized=FALSE) {
  ff <- if (vectorized)  {
    function(x) fn(x, ...)
  } else {
    function(x) apply(x, 1, fn, ...)
  }

  cf <- if (vectorized) {
    function(x) constraints(x, ...)
  } else {
    function(x) apply(x, 1, constraints, ...)
  }

  ## Make sure popsize is a multiple of 4
  if (popsize %% 4 != 0)
    stop("Population size must be a multiple of 4")

  ## Check bounding box
  if (any(is.infinite(lower.bounds)) || any(is.infinite(upper.bounds))) {
    warning("While it is possible to optimize an unconstrained problem, it is not recommended. Please consider adding finite upper and lower bounds.")
    sane_lower <- -.Machine$double.xmax / 4
    sane_upper <- .Machine$double.xmax / 4
    lower.bounds <- ifelse(lower.bounds == -Inf, sane_lower,
                           ifelse(lower.bounds == Inf, sane_upper,
                                  lower.bounds))

    upper.bounds <- ifelse(upper.bounds == -Inf, sane_lower,
                           ifelse(upper.bounds == Inf, sane_upper,
                                  upper.bounds))
  }
  ## Lag generations:
  ## C source expects each element of generations to be the number of
  ## generations to go forward before saving the next result set. This is
  ## unintuitive to specify, so we compute the lagged differences here.
  if (length(generations) > 1)
    generations <- c(generations[1], diff(generations))

  if (!all(generations > 0))
    stop("Cannot go back in time! Your generations argument must be sorted!")

  ## Set cdim = 0 if no cfn was given:
  if (is.null(constraints)) cdim <- 0

  res <- .Call(do_nsga2,
               ff, cf, sys.frame(),
               as.integer(odim),
               as.integer(cdim),
               as.integer(idim),
               lower.bounds, upper.bounds,
               as.integer(popsize), as.integer(generations),
               cprob, as.integer(cdist),
               mprob, as.integer(mdist))
  if (1 == length(res)) {
    res <- res[[1]]
    names(res) <- c("par", "value", "pareto.optimal")
    class(res) <- c("nsga2", "mco")
  } else {
    for (i in 1:length(res)) {
      names(res[[i]]) <- c("par", "value", "pareto.optimal")
      class(res[[i]]) <- c("nsga2", "mco")
    }
    class(res) <- "nsga2.collection"
  }
  return (res)
}

plot.nsga2 <- function(x, ...) {
  v <- x$value
  o <- x$pareto.optimal
  d <- ncol(v)
  col <- ifelse(o, "red", "blue")
  pch <- ifelse(o, 4, 19)
  if (d <= 2) {
    plot(v, col=col, pch=pch, ...)
    ov <- v[o,]
    ov <- ov[order(ov[,1]),]
    lines (ov, col="red", type="s")
  } else if (d == 3) {
    if (require("scatterplot3d")) {
      scatterplot3d::scatterplot3d(v, color=ifelse(o, "red", "blue"))
    } else {
      pairs(v, col=col, pch=pch, ...)
    }
  } else {
    pairs(v, col=col, pch=pch, ...)
  }
}

plot.nsga2.collection <- function(x, ...) {
  oask <- devAskNewPage(TRUE)
  on.exit(devAskNewPage(oask))
  sapply(x, plot)
  return;
}

paretoSet.nsga2 <- function(x, ...)
  x$par[x$pareto.optimal,]

paretoFront.nsga2 <- function(x, ...)
  x$value[x$pareto.optimal,,drop=FALSE]
