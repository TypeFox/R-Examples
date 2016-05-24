transBeta <- function (x, p, interval = c(0, 1),
                       inv = FALSE, pmin = -3, pmax = 3, p0 = c(0, 0)) {
  x <- (x - interval[1]) / (interval[2] - interval[1])
  erg <- if (inv) qbeta(x, exp(p[1]), exp(p[2])) else pbeta(x, exp(p[1]), exp(p[2]))
  erg * (interval[2] - interval[1]) + interval[1]
}

transSimplex <- function (x, p, interval = c(0, 1),
                          inv = FALSE, pmin = -2, pmax = 2, p0 = c(0, 0, 0, 0, 0)) {
  t <- seq(interval[1], interval[2], length.out = length(p) + 2)
  y <- c(0, cumsum(exp(c(p, 0)))) / sum(exp(c(p, 0))) * (interval[2] - interval[1]) + interval[1]
  if (inv)
    fun <- approxfun(y, t, rule = 2)
  else
    fun <- approxfun(t, y, rule = 2)
  fun(x)
}

Bernstein <- function (x, N) {
  n <- length(x)
  X <- x %o% rep(1, N + 1)
  I <- rep(1, n) %o% 0 : N
  (rep(1, n) %o% choose(N, 0 : N)) * ((1 - X)^(N - I)) * X^I
}

Bezier1D <- function (x, p) {
  c(Bernstein(x, length(p) - 1) %*% p)
}

parallelRoot <- function (f, y, ... , interval, lower = min(interval),
                          upper = max(interval), iter = 12, eps = 1E-10,
                          do.plot = FALSE) {
  x <- seq(lower, upper, length.out = 200)
  yn <- f(x, ...)
  return (approxfun(yn, x, rule = 2) (y))
}

transBezier <- function (x, p, interval = c(0, 1),
                         inv = FALSE, pmin = 0, pmax = 1,
                         p0 = c(0.25, 0.25, 0.75, 0.75)) {
  n <- length(p) / 2
  p1 <- c(0, p[(1:n) * 2 - 1], 1)
  p2 <- c(0, p[(1:n) * 2], 1)
  #p1 <- c(p[(1:n) * 2 - 1], 0)
  #p2 <- c(p[(1:n) * 2], 0)
  #p1 <- c(0, cumsum(exp(p1)) / sum(exp(p1)))
  #p2 <- c(0, cumsum(exp(p2)) / sum(exp(p2)))
  if (inv) { p0 <- p1 ; p1 <- p2 ; p2 <- p0 }
  t <- parallelRoot(Bezier1D, (x - interval[1]) / (interval[2] - interval[1]),
                    p = p1, do.plot = FALSE, interval = c(0, 1))
  y <- Bezier1D(t, p = p2) * (interval[2] - interval[1]) + interval[1]
  y
}

solidOptim <- function (p0, f, mi, ma, n = 100, debug = FALSE) {
  wr <- getOption("warn")
  if(!debug) options(warn = 3)
  opt <- Inf
  stp <- FALSE
  env <- environment()
  li <- lapply(1:n, function(i) {
    cat(".")
    if(i %% 20 == 0) cat("\n")
    if (stp) return (list(value = Inf, convergence = 11))
    withCallingHandlers(
      tr <- try ({
        if (i==1) {
          optim(p0, f)
        } else if (i==2) {
          optim(p0, f, method = "L-BFGS-B", lower = mi, upper = ma)
        } else {
          p <- runif(length(p0)) * (ma - mi) + mi
          if(i %% 3 == 0)
            optim(p, f, method = "Nelder-Mead")
          #else if (i==5)
          #  optim(p, f, method = "SANN", control = list(maxit = 5000))
          else
            optim(p, f, method = "L-BFGS-B", lower = mi, upper = ma)
        }
      }, silent = TRUE), interrupt = function (e) assign("stp", TRUE, envir = env)
    )
#    optim(runif(length(p0)) * (ma - mi) + mi, f, method = "L-BFGS-B", lower = mi, upper = ma) })
    if (inherits(tr, "try-error")) {
      if (debug) print(tr)
      return (list(value = Inf, convergence = 10))
    }
    if (tr$convergence > 2) tr$value <- Inf
    nopt <- tr$value
    if (is.finite(nopt) && nopt < opt) {
      assign("opt", nopt, envir = env)
      if (debug) cat("=>", opt, "\n")
    }
    tr
  })
  if (debug) {
    print(sapply(li, function (x) {
      k <- x$convergence
      if (is.null(k)) { print(x); 10 } else k
    }))
  }
#  print(sapply(li, function(x) x$message))
#  print(sapply(li, function(x) x$par))
  v <- sapply(li, function (x) x$value)
  if (debug) print(v)
  o <- order(v)
  options(warn = wr)
  erg <- li[[o[1]]]
  erg$OtherParams <- v <- sapply(li, function(x) x$par)
  erg$OtherValues <- v <- sapply(li, function(x) x$value)
  cat("\n")
  erg
}

compareME <- function (o, p,
                       o.t = seq(0, 1, length.out = length(o)),
                       p.t = seq(0, 1, length.out = length(p)),
                       ignore = c("raw", "centered", "scaled", "ordered"),
                       geometry = c("real", "logarithmic", "geometric", "ordinal"),
                       measure = c("mad", "var", "sd"),
                       type = "normalized",
                       time = "fixed", ... , col.vars = c("time", "ignore")
                      ) {
  hasNoTime <- missing(o.t) && missing(p.t) && missing(time)
  type = match.arg(type, c("dissimilarity", "normalized", "similarity",
    "reference", "formula", "name"), several.ok = TRUE)
  time = match.arg(time, c("fixed", "transform"), several.ok = TRUE)
  ignore = match.arg(ignore, several.ok = TRUE)
  geometry = match.arg(geometry, several.ok = TRUE)
  dim = c(geometry = length(geometry),
    measure = length(measure), ignore = length(ignore), time = length(time))
  dimnames = list (geometry = geometry, measure = measure, ignore = ignore,
    time = time)
  names(type) <- type
  kapply <- function (l, f) c(sapply(l, f))
  if(hasNoTime) {
    # Use fixed measures
    structure(lapply(type, function (type) {
      ftable(structure(
        kapply(time, function (time) {
          kapply(ignore, function (ignore) {
            kapply(measure, function (measure) {
              kapply(geometry, function(geometry) {
                generalME(o, p, ignore = ignore, geometry = geometry,
                          measure = measure, type = type)
              } )
            } )
          } )
        } ),
        dim = dim, dimnames = dimnames), col.vars = col.vars)
    }), class="compareME")
  } else {
    # Use transformed measures
    tixList <- list(o = o, p = p, o.t = o.t, p.t = p.t, ...)
    structure(lapply(type, function(type) {
      ftable(structure (
        kapply(time, function (time) {
          kapply(ignore, function (ignore) {
            kapply(measure, function (measure) {
              kapply(geometry, function(geometry) {
                do.call(timeTransME,
                  c(list(ignore = ignore, geometry = geometry,
                         measure = measure, type = type, time = time),
                    tixList))$totalME
              } )
            } )
          } )
        } ),
        dim = dim, dimnames = dimnames), col.vars=col.vars)
    } ), class="compareME")
  }
}

timeTransME <- function(o, p,
                        o.t = seq(0, 1, length.out = length(o)),
                        p.t = seq(0, 1, length.out = length(p)),
                        ignore = "scaled",
                        geometry = "real",
                        measure = "mad",
                        type = c("dissimilarity", "normalized",
                                 "similarity", "reference"),
                        interval = range(c(o.t, p.t)),
                        time = c("transformed", "fixed"),
                        trans = transBeta,
                        p0 = eval(formals(trans)$p0),
                        pmin = eval(formals(trans)$pmin, list(p = p0)),
                        pmax = eval(formals(trans)$pmax, list(p = p0)),
                        timeMEFactor = 0,
                        timeME = MAE,
                        timeMEtype = "normalized",
                        timeScale = 1,
                        ME = generalME(o, p, ignore, geometry, measure, type = "function"),
                        MEtype = c("dissimilarity", "normalized"),
                        trials = 100,
                        debug = FALSE) {

  if (is.character(ME)) ME <- generalME(method = ME)
  O    <- xy.coords(o, o.t)
  P    <- xy.coords(p, p.t)
  o    <- O$x
  o.t  <- O$y
  p    <- P$x
  p.t  <- P$y
  time   <- match.arg(time)
  type   <- match.arg(type)
  MEtype <- match.arg(MEtype)
  if (length(o) == 0 || length(o.t) != length(o))
    stop ("observations not specified")
  if (length(p) == 0 || length(p.t) != length(p))
    stop ("model not specified")

  # transform time scales into interval [0, 1] and approximate functions
  o.f <- approxfun(o.t, o, rule = 2)
  p.f <- approxfun(p.t, p, rule = 2)
  # define verbose error function
  cll <- match.call()
  error <- function (timep) {
    x   <- sort(c(o.t, trans(p.t, timep, interval, inv = TRUE)))
    ## thpe: bug! x may exceed interval !!!
    x <- ifelse(x > interval[2], interval[2],       # fix
           ifelse(x < interval[1], interval[1], x)  # may be necessary, too
         )
    xtr <- trans(x, timep, interval)

    timeme <- timeME(x * timeScale, xtr * timeScale, type = timeMEtype)
    timemeReference <- switch(timeMEtype,
                              dissimilarity = timeME(x * timeScale,
                                xtr * timeScale, type = "reference"),
                              normalized    = 1,
                              reference     = 1)
    yo <- o.f(x)
    yp <- p.f(xtr)
    has <- !(is.na(yo) | is.na(yp))
    x  <- x[has]
    yo <- yo[has]
    yp <- yp[has]
    me <- ME(yo, yp, type = MEtype)
    meref <- switch(MEtype,
                    dissimilarity = ME(yo, yp, type = "reference"),
                    normalized    = 1)
    structure(
      list(criterium = me + timeMEFactor * timeme,
           reference = meref, #timeMEFactor * timemeReference,
           call = cll,
           x = x, xp = xtr, yo = yo, yp = yp,
           ME = me, MEref = meref,
           timeME = timeme,
           timeMEref = timemeReference,
           timeMEFactor = timeMEFactor,
           timeScale = timeScale,
           p = timep, interval = interval,
           trans = trans
          ), class = "timeTransME")
  }
  # define simplified error function
  quickError <- function(timep) {
    x <- c(o.t, trans(p.t, timep, interval, inv = TRUE))
    ## thpe: bug?? x may exceed interval see error() above
    x <- ifelse(x > interval[2], interval[2],       # fix
           ifelse(x < interval[1], interval[1], x)  # may be necessary, too
         )
    xtr <- trans(x, timep, interval)
    ME(o.f(x), p.f(xtr), type = MEtype) +
      timeMEFactor * timeME(x * timeScale, xtr * timeScale, type = timeMEtype)
  }
  if (time == "LCS") {
    stop ("LCS not yet implemented in timetransME")
  } else if(time == "fixed") {
    erg <- error(p0)
  } else {
  # transformation
  # error function
  # optimizes the error function by means of General-purpose Optimization
    pneuE <- solidOptim(p0, quickError, pmin, pmax, n = trials, debug = debug)
    erg <- error(pneuE$par)
    erg$optim <- pneuE
  }
  erg$totalME <- switch(type,
                        dissimilarity = erg$criterium,
                        normalized = erg$criterium / erg$reference,
                        similarity = 1 - erg$criterium / erg$reference,
                        reference = erg$reference
                       )
  erg
}

print.timeTransME <- function (x, ..., digits=3) {
  oo <- options("digits")
  options("digits"=digits)
  print(unlist(x[c("totalME", "timeME")]))
  options(oo)
}

summary.timeTransME <- function(object, ...) {
  print(unclass(object))
}

plot.timeTransME <- function(x, y = NULL, ..., col.obs = "black",
                             col.pred = "green", col.map = "red",
                             sub = x$call, xlab = "t",
                             xlim = range(x$x),
                             ylim = range(c(0, x$yo, x$yp))) {
  plot(x$x, x$yo, xlim = xlim, ylim = ylim, xlab = xlab, ylab = "y",
       type = "b", ..., col = col.obs, sub = sub)
  lines(x$xp, x$yp, type = "b", col = col.pred)
  lines(x$x, x$yp, type = "b", col = col.map)
}

print.compareME <- function (x, ..., digits = 3) {
  oo <- options("digits")
  options("digits" = digits)
  print(unclass(x))
  options(oo)
}

summary.compareME <- function(object, ...) {
  print(object, ...)
}
