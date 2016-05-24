
## Interface to C version of LCS provided by Dominik Reusser

LCS <- function(a, b){
    stopifnot(is.character(a), is.character(b))
    if(any(is.na(c(a, b)))){
       out <- list(a = a, b = b, LLCS = NA, LCS = NA, QSI = NA, va = NA, vb = NA)
    } else {
       out <- .Call("lcs", as.character(a), as.character(b), max(nchar(c(a, b))), PACKAGE="qualV")
    }
    invisible(out)
}


#LCS <- function (a, b) {
#  m <- length(a)
#  n <- length(b)
#  if (m == 0 || n == 0) stop ("vector of length zero")
#
#  # creates a table
#  M <- matrix(nrow = m + 1, ncol = n + 1)
#  M[, 1] <- 0
#  M[1, ] <- 0
#
#  # fills the table
#  for (i in 2:(m + 1)) {
#    for (j in 2:(n + 1)) {
#      if (a[i - 1] == b[j - 1]) { M[i, j] <- M[i - 1, j - 1] + 1 }
#      else { M[i, j] <- max(c(M[i, j - 1], M[i - 1, j])) }
#    }
#  }
#
#  # length of the longest common subsequence
#  LLCS <- M[m + 1, n + 1]
#
#  # determines one possible longest common subsequence
#  # by means of "trace-back" by the table filled out
#  i <- m + 1; j <- n + 1
#  LCS <- va <- vb <- NULL
#  while (i > 1 & j > 1) {
#    if (M[i, j] == M[i - 1, j - 1] + 1 & a[i - 1] == b[j - 1]) {
#      LCS <- c(a[i - 1], LCS)
#      va <- c(i - 1, va); vb <- c(j - 1, vb)
#      i <- i - 1; j <- j - 1
#    }
#    else {
#      if (M[i - 1, j] > M[i, j - 1]) { i <- i - 1 }
#      else { j <- j - 1 }
#    }
#  }
#
#  # calculates the quality similarity index
#  QSI <- round(LLCS / max(m, n), digits = 2)
#  invisible(list(a = a, b = b, LLCS = LLCS, LCS = LCS, QSI = QSI , va = va, vb = vb))
#}

qvalLCS <- function (o, p,
                     o.t = seq(0, 1, length.out = length(o)),
                     p.t = seq(0, 1, length.out = length(p)),
                     smooth = c("none", "both", "obs", "sim"),
                     feature = c("f.slope", "f.curve", "f.steep", "f.level")
                     ) {
  use <- match.arg(smooth, c("none", "both", "obs", "sim"))
  feature <- match.arg(feature, c("f.slope", "f.curve", "f.steep", "f.level"))
  f <- get(feature)

  if (length(o) <= 0 || length(o.t) <= 0) stop ("observed object of length zero")
  if (length(p) <= 0 || length(p.t) <= 0) stop ("simulated object of length zero")

  obs <- data.frame(x = o.t, y = o)
  sim <- data.frame(x = p.t, y = p)

  if (use == "both" || use == "obs") {
    o.bw  <- dpill(o.t, o)
    n     <- max(tail(o.t, n = 1) - o.t[1] + 1, length(o.t))
    if (is.na(o.bw) == FALSE) {
      obs <- ksmooth(o.t, o, kernel = "normal", bandwidth = o.bw, n.points = n)
      obs <- data.frame(x = obs$x, y = obs$y)
    }
  }
  if (use == "both" || use == "sim") {
    s.bw <- dpill(p.t, p)
    n    <- max(tail(p.t, n = 1) - p.t[1] + 1, length(p.t))
    if(is.na(s.bw) == FALSE) {
      sim <- ksmooth(p.t, p, kernel = "normal", bandwidth = s.bw, n.points = n)
      sim <- data.frame(x = sim$x, y = sim$y)
    }
  }

  # correct time code
  obs <- obs[na.omit(match(sim$x, obs$x)) , ]
  sim <- sim[na.omit(match(obs$x, sim$x)) , ]

  # determine features
  obsf <- f(obs$x, obs$y)
  simf <- f(sim$x, sim$y)

  lcs  <- LCS(obsf, simf)

  erg <- structure(
           list(smooth = use,
                feature = feature,
                o = data.frame(x = o.t, y = o),
                p = data.frame(x = p.t, y = p),
                obs = obs,
                sim = sim,
                obsf = obsf,
                simf = simf,
                lcs = lcs
               ),
               class = "qvalLCS")
  erg
}

print.qvalLCS <- function (x, ...) {
  erg <- unclass(x[c("lcs")])
  print(paste("QSI:", erg$lcs$QSI))
}

summary.qvalLCS <- function (object, ...) {
  print(unclass(object))
}

plot.feature <- function(data.f, lsc.v) {
  pf <- numeric(length(lsc.v))
  pf[1] <- 1
  for (i in 2:length(lsc.v)) {
    if (data.f[lsc.v[i]] == data.f[lsc.v[i-1]]) {
      pf[i] <- pf[i-1]
    } else {
        pf[i] <- pf[i-1] + 1
    }
  }
  pf
}

plot.qvalLCS <- function (x, y = NULL, ...,
                          xlim = range(c(x$obs$x, x$sim$x)),
                          ylim = range(c(x$obs$y, x$sim$y)),
                          xlab = "time",
                          ylab = " ",
                          col.obs = "black",
                          col.pred = "red",
                          plot.title = paste("LLCS =", x$lcs$LLCS, ", QSI =", x$lcs$QSI),
                          legend = TRUE
                         ) {
  leg.lty <- NULL
  leg.pch <- NULL
  leg.col <- NULL
  leg.name <- c("measurement")

  ca <- plot.feature(x$obsf, x$lcs$va)
  cb <- plot.feature(x$simf, x$lcs$vb)

  plot(x$sim$x, x$sim$y, xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab,
       type = "l", ..., col = col.pred, lwd = 2)
  lines(x$obs$x, x$obs$y, lwd = 2, col = col.obs)
  title(main = plot.title)

  if (x$smooth == "both" || x$smooth == "obs") {
    points(x$o$x, x$o$y, lwd=2, pch=2, col=col.obs)
    leg.lty <- c(leg.lty, 0)
    leg.pch <- c(leg.pch, 2)
    leg.col <- c(leg.col, col.obs)
    leg.name <- c(leg.name, "smoothed measurement")
  }
  leg.lty <- c(leg.lty, 1, 0)
  leg.pch <- c(leg.pch, NA, 21)
  leg.col <- c(leg.col, col.obs, col.obs)
  leg.name <- c(leg.name, "LCS segments of measurement", "simulation")
  if (x$smooth == "both" || x$smooth == "sim") {
    points(x$p$x, x$p$y, lwd = 2, pch = 2, col = col.pred)
    leg.lty <- c(leg.lty, 0)
    leg.pch <- c(leg.pch, 2)
    leg.col <- c(leg.col, col.pred)
    leg.name <- c(leg.name, "smoothed simulation")
  }
  leg.lty <- c(leg.lty, 1, 0)
  leg.pch <- c(leg.pch, NA, 21)
  leg.col <- c(leg.col, col.pred, col.pred)
  leg.name <- c(leg.name, "LCS segments of simulation")

  points(x$obs$x[x$lcs$va], x$obs$y[x$lcs$va], pch = 21, bg = 1, col = ca)
  points(x$sim$x[x$lcs$vb], x$sim$y[x$lcs$vb], pch = 21, bg = 2, col = cb)

  if (legend == TRUE) {
    usr <- par("usr")
    legend(usr[1],usr[4], leg.name, lty = leg.lty,
           pch = leg.pch, lwd = c(2), col=  leg.col)
  }
}
