
as.rvsummary <- function (x, ...) {
  UseMethod("as.rvsummary")
}

is.rvsummary <- function (x) {
  inherits(x, "rvsummary")
}

print.rvsummary <- function (x, digits=3, ...) # METHOD
{
  s <- summary(x)
  for (i in which(sapply(s, is.numeric))) {
    s[[i]] <- round(s[[i]], digits=digits)
  }
  print(s)
}

as.rvsummary.default <- function (x, ...)  # NOEXPORT
{
  as.rvsummary(as.rv(x), ...)
}

as.rvsummary.rv <- function (x, quantiles=(0:200/200), ...)  # NOEXPORT
{
  y <- if (is.logical(x)) {
    as.rvsummary.rvlogical(x, ...)
  } else if (is.integer(x)) {
    as.rvsummary.rvinteger(x, quantiles=quantiles, ...)
  } else {
    as.rvsummary.rvnumeric(x, quantiles=quantiles, ...)
  }
  return(y)
}





as.rvsummary.rvsummary <- function (x, ...)  # NOEXPORT
{
  return(x)
}

as.rvsummary.rvnumeric <- function (x, quantiles=(0:200/200), ...) # NOEXPORT
{
  ms <- .rvmeansd(x, names.=c("mean", "sd", "NAS", "n.sims"))
  for (name in names(ms)) {
    rvattr(x, name) <- ms[[name]]
  }
  for (i in seq_along(x)) {
    a <- attributes(x[[i]])
    Q <- quantile(x[[i]], probs=quantiles, na.rm=TRUE)
    attributes(Q) <- a
    x[[i]] <- Q
  } 
  structure(x, class=c("rvsummary_numeric", "rvsummary"), quantiles=quantiles)
}

as.rvsummary.rvinteger <- function (x, quantiles=(0:200/200), ...) # NOEXPORT
{
  ms <- .rvmeansd(x, names.=c("mean", "sd", "NAS", "n.sims"))
  for (name in names(ms)) {
    rvattr(x, name) <- ms[[name]]
  }
  for (i in seq_along(x)) {
    a <- attributes(x[[i]])
    Q <- quantile(x[[i]], probs=quantiles, na.rm=TRUE)
    attributes(Q) <- a
    x[[i]] <- Q
  } 
  structure(x, class=c("rvsummary_integer", "rvsummary"), quantiles=quantiles)
}



as.rvsummary.rvlogical <- function (x, ...) # NOEXPORT
{
  ms <- .rvmeansd(x, names.=c("mean", "sd", "NAS", "n.sims"))
  for (name in names(ms)) {
    rvattr(x, name) <- ms[[name]]
  }
  for (i in seq_along(x)) {
    a <- attributes(x[[i]])
    x[[i]] <- ms[["mean"]][i]
    attributes(x[[i]]) <- a
  } 
  structure(x, class=c("rvsummary_logical", "rvsummary"))
}

as.rvsummary.rvfactor <- function (x, ...) # NOEXPORT
{
  levels <- levels(x)
  llev <- length(levels)
  num.levels <- seq_len(llev)
  #
  S <- sims(x)
  a <- apply(S, 2, function (x) table(c(x, num.levels))) # ensure that all levels are tried
  if (is.null(dim(a))) {
    dim(a) <- c(ncol(S), llev)
  }
  a <- (a-1) # And now subtract the extra counts from the matrix that was obtained.
  ns <- rvnsims(x)
  if (any(naS <- is.na(S))) {
    NAS <- (colMeans(naS)*100)
  } else {
    NAS <- rep.int(0, length(x))
  }
  nax <- if (is.null(dim(x))) NULL else names(x)
  M <- a
  rownames(M) <- levels
  remaining <- (ns-colSums(M))
  if (any(remaining>0)) {
    stop("Impossible: levels won't sum up to 0")
  } 
  P <- t(M/ns)  # compute proportions in each category and transpose
  for (i in seq_along(x)) {
    a <- attributes(x[[i]])
    x[[i]] <- P[i,]
    attributes(x[[i]]) <- a
    ###names(x[[i]]) <- levels
  } 
  rvattr(x, "n.sims") <- as.list(ns)
  rvattr(x, "NAS") <- as.list(NAS)
  structure(x, class=c("rvsummary_rvfactor", "rvsummary"))
}

as.rvsummary.data.frame <- function (x, quantiles=rvpar("summary.quantiles.numeric"), ...)
{
  name <- names(x)
  rnames <- rownames(x)
  q.columns <- (regexpr("^([0-9.]+)%", name)>0)
  q.names <- name[q.columns]
  d.quantiles <- (as.numeric(gsub("^([0-9.]+)%", "\\1", name[q.columns]))/100)
  lx <- nrow(x)
  ms <- list()
  ms$n.sims <- if ("sims" %in% name) { x[["sims"]] } else { rep(Inf, lx) }
  ms$NAS    <- if ("NA%" %in% name) { x[["NA%"]] } else { rep(0L, lx) }
  ms$mean   <- if ("mean" %in% name) { x[["mean"]] } else { rep(NA, lx) }
  ms$sd     <- if ("sd" %in% name) { x[["sd"]] } else { rep(NA, lx) }
  if (length(d.quantiles)==0) {
    if (length(quantiles)>0) {
      x <- lapply(seq_along(ms$mean), function (i) qnorm(quantiles, mean=ms$mean[i], sd=ms$sd[i]))
    } else {
      x <- as.list(rep.int(NA, nrow(x)))
    }
    d.quantiles <- quantiles
  } else {
    x <- as.matrix(x[q.columns])
    x <- split(x, row(x))
  }
  for (name in names(ms)) {
    rvattr(x, name) <- ms[[name]]
  }
  structure(x, class=c("rvsummary_numeric", "rvsummary"), quantiles=d.quantiles, names=rnames)
}

as.double.rvsummary <- function (x, ...)
{
  if (is.null(attr(x, "quantiles"))) {
    stop("Cannot coerce to double.")
  }
  return(x)
}

print.rvsummary_rvfactor <- function (x, all.levels=FALSE, ...) # METHOD
{
  print(summary(x, all.levels=all.levels, ...))
}

as.data.frame.rvsummary <- function (x, ...) {
  S <- summary(x, ...)
  rownames(S) <- S[["name"]]
  S[["name"]] <- NULL
  return(S)
}

summary.rvsummary <- function (object, ...)
{
  # supposed to be called AFTER other methods (e.g. summary.rvsummary_rvnumeric)
  # assumes that the 'summary' slot exists in the rvsummary object;
  # this function adds:
  #  1. name, NA%, n.sims
  #  2. Rhat, n.eff
  #  3. dimnames columns
  #
  x <- object
  xdim <- dim(x)
  xdimnames <- dimnames(x)
  Summary <- attr(object, "summary")
   if (length(.names <- names(x))>0) {
    Summary <- cbind(name=.names, Summary)
  }
  Col <- NULL
  n.sims. <- rvnsims(x)
  NAS <- unlist(rvattr(x, "NAS"))
  if (all(NAS==0)) {
    Col <- data.frame(sims=n.sims.)
  } else {
    Col <- data.frame("NA%"=NAS, sims=n.sims.)
  }
  if (!all(is.na(Rhats <- rvRhat(x)))) {
    Col <- cbind(Col, Rhat=Rhats)
  }
  if (!all(is.na(n.effs <- rvneff(x)))) {
    Col <- cbind(Col, n.eff=n.effs)
  }
  Summary <- cbind(Summary, Col)
  if (!is.null(unlist(xdimnames))) {
    # 'is.null(unlist(dimnames))' since we might have e.g. list(NULL, NULL) 
    sud <- rvpar("summary.dimnames")
    if (is.null(sud) || isTRUE(sud)) {
      .f <- function (i) {
        X <- .dimind(dim.=xdim, MARGIN=i)
        na <- xdimnames[[i]]
        if (!is.null(na)) { na <- na[X] }
        return(na)
      }
      da <- lapply(seq_along(xdim), .f)
      names(da) <- names(xdimnames)
      if (is.null(names(da))) {
        names(da) <- if (length(xdim)==2) c("row", "col") else paste("d", seq_along(da), sep="")
      }
      da <- da[!sapply(da, is.null)]
      if (length(da)>0) {
        Summary <- cbind(as.data.frame(da), " "=':', Summary)
      }
    }
  }
  return(Summary)
}

summary.rvsummary_numeric <- function (object, ...)
{
  x <- object
  if (is.null(qs <- rvpar("summary.quantiles.numeric"))) {
    qs <- c(0.025, 0.25, 0.5, 0.75, 0.975)
  }
  S <- t(sims(x))
  Q <- attr(S, "quantiles")
  qa <- (Q%in%qs)
  q <- S[,qa,drop=FALSE]
  m <- rvmean(x)
  s <- rvsd(x)
  S <- data.frame(mean=m, sd=s)
  S <- cbind(S, as.data.frame(q))
  rownames(S) <- .dim.index(x)
  attr(object, "summary") <- S
  NextMethod()
}

summary.rvsummary_logical <- function (object, ...)
{
  x <- object
  S <- data.frame(mean=rvmean(x), sd=rvsd(x))
  rownames(S) <- .dim.index(x)
  attr(object, "summary") <- S
  NextMethod()
}

summary.rvsummary_integer <- function (object, ...)
{
  x <- object
  if (is.null(qs <- rvpar("summary.quantiles.integer"))) {
    qs <- c(0, 0.025, 0.25, 0.5, 0.75, 0.975, 1)
  }
  S <- t(sims(x))
  Q <- attr(S, "quantiles")
  qa <- (Q%in%qs)
  q <- S[,qa,drop=FALSE]
  names_quantiles <- dimnames(q)[[2]]
  names_quantiles[names_quantiles=="0%"] <- "min"
  names_quantiles[names_quantiles=="100%"] <- "max"
  dimnames(q)[[2]] <- names_quantiles
  m <- rvmean(x)
  s <- rvsd(x)
  S <- data.frame(mean=m, sd=s)
  S <- cbind(S, as.data.frame(q))
  rownames(S) <- .dim.index(x)
  attr(object, "summary") <- S
  NextMethod()
}



summary.rvsummary_rvfactor <- function (object, all.levels=TRUE, ...) 
{
  x <- object
  levels <- levels(x)
  llev <- length(levels)
  num.levels <- seq_along(levels)
  #
  maxlev <- if (is.null(maxlev <- rvpar("max.levels"))) { 10 } else maxlev
  too.many.levels.to.show <- ((!all.levels) && (llev>maxlev))
  last.lev.no <- llev
  proportions <- t(sims(x))
  if (too.many.levels.to.show) {
    P1 <- proportions[,1:(maxlev-1),drop=FALSE]
    P2 <- proportions[,last.lev.no,drop=FALSE]
    omit_levels <- (!seq_along(levels) %in% c(1:(maxlev-1), last.lev.no))
    rest <- rowSums(proportions[,omit_levels,drop=FALSE])
    M <- cbind(P1, "*"=rest, P2)
    colnames(M) <- c(levels[1:(maxlev-1)], "*", levels[last.lev.no])
  } else {
    M <- proportions
    colnames(M) <- levels
  }
  S <- as.data.frame(M)
  rownames(S) <- .dim.index(x)
  attr(object, "summary") <- S
  NextMethod()
}

