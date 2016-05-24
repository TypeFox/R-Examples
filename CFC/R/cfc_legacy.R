#### SPLIT UP THE CODE FOR TIME BASIS AND QUANTILE BASIS: IT REDUCES CONFUSION AND IT ALSO CREATES MORE FUNCTIONS AND "FELESH"!!

cfc.tbasis <- function(p1, p2, unity.tol = 1e-6,
  diff.tol = 1e-2, diff.tol.policy = c("mean","all"),
  check = TRUE) {
  # checks for p1,p2: 1) same dimensions, 2) between 0.0 and 1.0, 3) non-increasing with time, 4) start at 1.0, 5) check for large steps
  
  diff.tol.policy <- match.arg(diff.tol.policy)
  
  if (is.null(dim(p1))) {
    nt <- length(p1)
    if (length(p2) != nt) stop("p1 and p2 have unequal lengths")
    nother <- 1
    dim.p1 <- c(nt, 1)
  } else {
    dim.p1 <- dim(p1)
    if (!identical(dim.p1, dim(p2))) stop("p1 and p2 dimensions do not match")
    nt <- dim.p1[1]
    nother <- prod(dim.p1)/nt
  }
  if (check) if (any(p1<0.0 | p1>1.0 | p2<0.0 | p2>1.0)) stop("out-of-range probabilities")

  p1.2d <- array(p1, dim = c(nt, nother))
  p2.2d <- array(p2, dim = c(nt, nother))
  if (check) if (any(abs(p1.2d[1,] - 1.0) > unity.tol)) stop("p1 probabilities must start at 1.0")
  if (check) if (any(abs(p2.2d[1,] - 1.0) > unity.tol)) stop("p2 probabilities must start at 1.0")


  seq.left <- 1:(nt-1)
  seq.right <- seq.left+1
  
  dp1 <- apply(p1.2d, 2, diff)
  if (check) {
    if (any(dp1>0.0)) stop("increasing probabilities with time detected for p1")
    if (diff.tol.policy == "mean") {
      if (mean(dp1) < -1.0*diff.tol) stop("average change in p1 exceeds threshold")
    } else if (diff.tol.policy == "all") {
      if (any(dp1 < -1.0*diff.tol)) stop("one or more changes in p1 exceed threshold")
    }
  }
  dci1 <- -0.5 * (p2.2d[seq.left,] + p2.2d[seq.right,]) * dp1
  ci1 <- rbind(0, apply(dci1, 2, cumsum))
  ci1 <- array(ci1, dim = dim.p1)
  
  dp2 <- apply(p2.2d, 2, diff)
  if (check) {
    if (any(dp2>0.0)) stop("increasing probabilities with time detected for p2")
    if (diff.tol.policy == "mean") {
      if (mean(dp2) < -1.0*diff.tol) stop("average change in p2 exceeds threshold")
    } else if (diff.tol.policy == "all") {
      if (any(dp2 < -1.0*diff.tol)) stop("one or more changes in p2 exceed threshold")
    }
  }
  dci2 <- -0.5 * (p1.2d[seq.left,] + p1.2d[seq.right,]) * dp2
  ci2 <- rbind(0, apply(dci2, 2, cumsum))
  ci2 <- array(ci2, dim = dim.p1)

  if (is.null(dim(p1))) {
    ci1 <- drop(ci1)
    ci2 <- drop(ci2)
  }
  
  if (is.null(dim(p1))) {
    ret <- cbind(ci1, ci2, p1*p2)
    colnames(ret) <- c("ci1", "ci2", "efp")
  } else {
    ret <- list(ci1 = ci1, ci2 = ci2, efp = p1*p2)
  }
  class(ret) <- c("cfc.tbasis", class(ret))
  return (ret)
}

summary.cfc.tbasis <- function(object,
  MARGIN = if (class(object)[2] == "matrix") NULL else 1, ...) {
  if (class(object)[2] == "matrix") {
    class(object)[1] <- "summary.cfc.tbasis"
    attr(object, "popavg") <- FALSE
    return (object)
  }
  MARGIN <- as.integer(MARGIN)
  if (!(1 %in% MARGIN)) stop("time dimension cannot be aggregated")
  if (identical(MARGIN, 1:length(dim(object$ci1)))) stop("cannot keep all dimensions during aggregation")
  #cat("MARGIN:", MARGIN, "\n")
  ci1 <- apply(object$ci1, MARGIN = MARGIN, mean)
  ci2 <- apply(object$ci2, MARGIN = MARGIN, mean)
  efp <- apply(object$efp, MARGIN = MARGIN, mean)
  if (is.null(dim(ci1))) {
    ret <- cbind(ci1, ci2, efp)
    colnames(ret) <- c("ci1", "ci2", "efp")
  } else {
    ret <- list(ci1 = ci1, ci2 = ci2, efp = efp)
  }
  class(ret) <- c("summary.cfc.tbasis", class(ret))
  return (invisible(ret))
}

plot.summary.cfc.tbasis <- function(x, t = 1, ci = 0.95, ...) {
  if (class(x)[2] == "matrix") {
    nt <- dim(x)[1]
    if (length(t) == 1) t <- (0:(nt-1))*t
    else if (length(t) != nt) stop("bad length for t vector")
    plot(t, x[,"efp"], type = "l", ylim = c(0.0,1.0)
         , xlab = "Time", ylab = "Probability", col = "black")
    lines(t, x[,"ci1"], col = "red")
    lines(t, x[,"ci2"], col = "green")
    legend("topright", legend = c("Event-Free", "CI - Cause 1", "CI - Cause 2")
           , col = c("black", "red", "green"), lty = rep(1,3))
  } else {
    dims <- dim(x$ci1)
    nt <- dims[1]
    if (length(t) == 1) t <- (0:(nt-1))*t
    else if (length(t) != nt) stop("bad length for t vector")
    nother <- prod(dims)/nt
    ci1.2d <- array(x$ci1, dim = c(nt, nother))
    ci2.2d <- array(x$ci2, dim = c(nt, nother))
    efp.2d <- array(x$efp, dim = c(nt, nother))
    qvec <- c(0.5*(1-ci), 0.5, 0.5*(1+ci))
    efp.q <- t(apply(efp.2d, 1, quantile, probs = qvec))
    ci1.q <- t(apply(ci1.2d, 1, quantile, probs = qvec))
    ci2.q <- t(apply(ci2.2d, 1, quantile, probs = qvec))
    plot(t, efp.q[,2], type = "l", ylim = c(0.0, 1.0)
         , xlab = "Time", ylab = "Population Average", col = "black")
    lines(t, efp.q[,1], col = "black", lty = 2)
    lines(t, efp.q[,3], col = "black", lty = 2)
    lines(t, ci1.q[,2], col = "red")
    lines(t, ci1.q[,1], col = "red", lty = 2)
    lines(t, ci1.q[,3], col = "red", lty = 2)
    lines(t, ci2.q[,2], col = "green")
    lines(t, ci2.q[,1], col = "green", lty = 2)
    lines(t, ci2.q[,3], col = "green", lty = 2)
    legend("topright", legend = c("Event-Free", "CI - Cause 1", "CI - Cause 2")
           , col = c("black", "red", "green"), lty = rep(1,3))
    return (invisible(list(efp = efp.q, ci1 = ci1.q, ci2 = ci2.q)))
  }
}

cfc.pbasis <- function(t1, t2, probs, unity.tol = 1e-6,
  diff.tol = 1e-2, diff.tol.policy = c("all", "mean")) {
  # TODO: consider allowing unsorted vectors; sort and then check for validity
  diff.tol.policy <- match.arg(diff.tol.policy)

  if (abs(probs[1] - 1.0) > unity.tol) stop("probability vector must start at 1.0")
  if (any(diff(probs) >= 0.0)) stop("probabilities must be decreasing with time")
  if (diff.tol.policy == "all") {
    if (any(diff(probs) < -1.0*diff.tol)) stop("one or more changes in probs exceed threshold")
  } else if (diff.tol.policy == "mean") {
    if (mean(diff(probs)) < -1.0*diff.tol) stop("average change in probs exceeds threshold")
  }

  if (is.null(dim(t1))) {
    nt <- length(t1)
    if (!is.null(dim(t2)) || length(t2) != nt) stop("t1 and t2 dimensions do not match")
    nother <- 1
    dim.t1 <- c(nt, 1)
  } else {
    dim.t1 <- dim(t1)
    if (!identical(dim.t1, dim(t2))) stop("t1 and t2 dimensions do not match")
    nt <- dim.t1[1]
    nother <- prod(dim.t1)/nt
  }
  t1.2d <- array(t1, dim = c(nt, nother))
  t2.2d <- array(t2, dim = c(nt, nother))
  dt1 <- apply(t1.2d, 2, diff)
  dt2 <- apply(t2.2d, 2, diff)

  if (any(dt1 <= 0.0)) stop("non-increasing times detected in t1")
  if (any(dt2 <= 0.0)) stop("non-increasing times detected in t2")

  ret <- lapply(1:nother, function(n) {
    ta <- t1.2d[,n]
    tb <- t2.2d[,n]
    tmax <- min(max(ta), max(tb))
    tcomb <- sort(unique(c(ta, tb)))
    tcomb <- tcomb[which(tcomb <= tmax)]
    pa <- approx(ta, probs, tcomb)$y
    pb <- approx(tb, probs, tcomb)$y
    rettmp <- cbind(tcomb, cfc.tbasis(pa, pb, check = FALSE))
    colnames(rettmp) <- c("time", "ci1", "ci2", "efp")
    return (rettmp)
  })
  if (nother == 1) ret <- ret[[1]]
  class(ret) <- c("cfc.pbasis", class(ret))
  return (ret)
}

summary.cfc.pbasis <- function(object, ...) {
  if (class(object)[2] == "matrix") {
    class(object)[1] <- "summary.cfc.pbasis"
    attr(object, "popavg") <- FALSE
    return (object)
  }
  tmax <- min(sapply(object, function(x) max(x[,"time"])))
  tvec <- unique(sort(unlist(sapply(object, function(x) x[,"time"]))))
  tvec <- tvec[tvec < tmax]
  ci1 <- rowMeans(sapply(object, function(x) {
    approx(x[,"time"], x[,"ci1"], tvec)$y
  }))
  ci2 <- rowMeans(sapply(object, function(x) {
    approx(x[,"time"], x[,"ci2"], tvec)$y
  }))
  efp <- 1 - (ci1 + ci2)
  
  ret <- cbind(tvec, ci1, ci2, efp)
  colnames(ret) <- c("time", "ci1", "ci2", "efp")
  attr(ret, "popavg") <- TRUE
  class(ret) <- c("summary.cfc.pbasis", class(ret))
  
  return (invisible(ret))
}

plot.summary.cfc.pbasis <- function(x, ...) {
  popavg <- attr(x, "popavg")
  ylim <- c(0.0, 1.0)
  ylab <- if (popavg) "Population Average" else "Probability"
  plot(x[,"time"], x[,"efp"], type = "l", col = "black"
       , xlab = "Time", ylab = ylab, ylim = ylim)
  lines(x[,"time"], x[,"ci1"], col = "red")
  lines(x[,"time"], x[,"ci2"], col = "green")
  legend("topright", col = c("black", "red", "green")
         , legend = c("Event-Free", "CI - Cause 1", "CI - Cause 2")
         , lty = rep(1,3))
  return (invisible(NULL))
}

# print.summary.cfc.pbasis <- function(x, ...) {
#   if (attr(x, "popavg")) cat("Population averages:\n")
#   nprint <- 6
#   print(head(x, nprint))
#   if (nrow(x) > nprint) cat("(", nrow(x)-nprint, " more rows ...)\n", sep = "")
#   return (invisible(NULL))
# }




