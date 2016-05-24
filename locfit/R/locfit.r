"locfit"<-
function(formula, data = sys.frame(sys.parent()), weights = 1, cens = 0, base = 0, subset,
  geth = FALSE, ..., lfproc = locfit.raw)
{
  Terms <- terms(formula, data = data)
  attr(Terms, "intercept") <- 0
  m <- match.call()
  m[[1]] <- as.name("model.frame")
  z <- pmatch(names(m), c("formula", "data", "weights", "cens", "base",
    "subset"))
  for(i in length(z):2)
    if(is.na(z[i])) m[[i]] <- NULL
  frm <- eval(m, sys.frame(sys.parent()))
  if (nrow(frm) < 1) stop("fewer than one row in the data")
  vnames <- as.character(attributes(Terms)$variables)[-1]
  if(attr(Terms, "response")) {
    y <- model.extract(frm, "response")
    yname <- deparse(formula[[2]])
    vnames <- vnames[-1]
  }
  else {
    y <- yname <- NULL
  }
  x <- as.matrix(frm[, vnames])
  if(!inherits(x, "lp")) {
    if(length(vnames) == dim(x)[2]) {
      dimnames(x) <- list(NULL, vnames)
    }
  }
  if(!missing(weights))
    weights <- model.extract(frm, weights)
  if(!missing(cens))
    cens <- model.extract(frm, cens)
  if(!missing(base))
    base <- model.extract(frm, base)
  ret <- lfproc(x, y, weights = weights, cens = cens, base = base, geth = geth,
    ...)
  if(geth == 0) {
    ret$terms <- Terms
    ret$call <- match.call()
    if(!is.null(yname))
      ret$yname <- yname
    ret$frame <- sys.frame(sys.parent())
  }
  ret
}

"locfit.raw"<-
function(x, y, weights = 1, cens = 0, base = 0, scale = FALSE, alpha = 0.7,
         deg = 2, kern = "tricube", kt = "sph", acri = "none",
         basis = list(NULL), deriv = numeric(0), dc = FALSE, family,
         link = "default", xlim, renorm = FALSE, ev = rbox(),
         maxk = 100, itype = "default", mint = 20, maxit = 20, debug = 0,
         geth = FALSE, sty = "none")
{
  if(inherits(x, "lp")) {
    alpha <- attr(x, "alpha")
    deg <- attr(x, "deg")
    sty <- attr(x, "style")
    acri <- attr(x, "acri")
    scale <- attr(x, "scale")
  }
  if(!is.matrix(x)) {
    vnames <- deparse(substitute(x))
    x <- matrix(x, ncol = 1)
    d <- 1
  }
  else {
    d <- ncol(x)
    if(is.null(dimnames(x)))
      vnames <- paste("x", 1:d, sep = "")
    else vnames <- dimnames(x)[[2]]
  }
  n <- nrow(x)
  if((!missing(y)) && (!is.null(y))) {
    yname <- deparse(substitute(y))
    if(missing(family))
      family <- if(is.logical(y)) "binomial" else "qgaussian"
  }
  else {
    if(missing(family))
      family <- "density"
    y <- 0
    yname <- family
  }
  if(!missing(basis)) {
    ## assign("basis", basis, 1)
    deg0 <- deg <- length(basis(matrix(0, nrow = 1, ncol = d), rep(0, d)))
  }
  if(length(deg) == 1)
      deg = c(deg, deg)

  xl <- rep(0, 2 * d)
  lset <- 0
  if(!missing(xlim)) {
    xl <- lflim(xlim, vnames, xl)
    lset <- 1
  }
  if(is.character(ev)) {
    stop("Character ev argument no longer used.")
  }
  if(is.numeric(ev)) {
    xev <- ev
    mg <- length(xev)/d
    ev <- list(type = "pres", xev = xev, mg = mg, cut = 0, ll = 0, ur = 0)
    if(mg == 0)
      stop("Invalid ev argument")
  }
  fl <- c(rep(ev$ll,length.out=d), rep(ev$ur,length.out=d))
  mi <- c(n, 0, deg, d, 0, 0, 0, 0, mint, maxit, renorm, 0, 0, 0, dc, maxk,
    debug, geth, 0, !missing(basis))
  if(any(is.na(mi)))
    print(mi)
  if(is.logical(scale))
    scale <- 1 - as.numeric(scale)
  if(length(scale) == 1)
    scale <- rep(scale, d)
  if(is.character(deriv))
    deriv <- match(deriv, vnames)
  alpha <- c(alpha, 0, 0, 0)[1:3]
  style <- pmatch(sty, c("none", "z1", "z2", "angle", "left", "right", "cpar"))
  if(length(style) == 1)
    style <- rep(style, d)
  dp <- c(alpha, ev$cut, 0, 0, 0, 0, 0, 0)
  size <- .C("guessnv",
    lw = integer(7),
    evt = as.character(c(ev$type, kt)),
    dp = as.numeric(dp),
    mi = as.integer(mi),
    nvc = integer(5),
    mg = as.integer(ev$mg), PACKAGE="locfit")
  nvc <- size$nvc
  lw <- size$lw
  z <- .C("slocfit",
    x = as.numeric(x),
    y = as.numeric(rep(y, length.out = n)),
    cens = as.numeric(rep(cens, length.out = n)),
    w = as.numeric(rep(weights, length.out = n)),
    base = as.numeric(rep(base, length.out = n)),
    lim = as.numeric(c(xl, fl)),
    mi = as.integer(size$mi),
    dp = as.numeric(size$dp),
    strings = c(kern, family, link, itype, acri, kt),
    scale = as.numeric(scale),
    xev = if(ev$type == "pres") as.numeric(xev) else numeric(d * nvc[1]),
    wdes = numeric(lw[1]),
    wtre = numeric(lw[2]),
    wpc = numeric(lw[4]),
    nvc = as.integer(size$nvc),
    iwk1 = integer(lw[3]),
    iwk2 = integer(lw[7]),
    lw = as.integer(lw),
    mg = as.integer(ev$mg),
    L = numeric(lw[5]),
    kap = numeric(lw[6]),
    deriv = as.integer(deriv),
    nd = as.integer(length(deriv)),
    sty = as.integer(style),
    basis = list(basis, lfbas), PACKAGE="locfit")
  nvc <- z$nvc
  names(nvc) <- c("nvm", "ncm", "vc", "nv", "nc")
  nvm <- nvc["nvm"]
  ncm <- nvc["ncm"]
  nv <- max(nvc["nv"], 1)
  nc <- nvc["nc"]
  if(geth == 1)
    return(matrix(z$L[1:(nv * n)], ncol = nv))
  if(geth == 2)
    return(list(const = z$kap, d = d))
  if(geth == 3)
    return(z$kap)
  dp <- z$dp
  mi <- z$mi
  names(mi) <- c("n", "p", "deg0", "deg", "d", "acri", "ker", "kt", "it",
    "mint", "mxit", "renorm", "ev", "tg", "link", "dc", "mk", "debug", "geth",
    "pc", "ubas")
  names(dp) <- c("nnalph", "fixh", "adpen", "cut", "lk", "df1", "df2", "rv",
    "swt", "rsc")
  if(geth == 4) {
    p <- mi["p"]
    return(list(residuals = z$y, var = z$wdes[n * (p + 2) + p * p + (1:n)],
      nl.df = dp["df1"] - 2))
  }
  if(geth == 6)
    return(z$L)
  if(length(deriv) > 0)
    trans <- function(x)
    x
  else trans <- switch(mi["link"] - 2,
      function(x)
      x,
      exp,
      expit,
      function(x)
      1/x,
      function(x)
      pmax(x, 0)^2,
      function(x)
      pmax(sin(x), 0)^2)
  t1 <- z$wtre
  t2 <- z$iwk1
  xev <- z$xev[1:(d * nv)]
  if(geth == 7)
    return(list(x = xev, y = trans(t1[1:nv])))
  coef <- matrix(t1[1:((3 * d + 8) * nvm)], nrow = nvm)[1:nv,  ]
  if(nv == 1)
    coef <- matrix(coef, nrow = 1)
  if(geth >= 70) {
    data <- list(x = x, y = y, cens = cens, base = base, w = weights)
    return(list(xev = matrix(xev, ncol = d, byrow = TRUE), coef = coef[, 1], sd =
      coef[, d + 2], lower = z$L[1:nv], upper = z$L[nvm + (1:nv)], trans =
      trans, d = d, vnames = vnames, kap = z$kap, data = data, mi = mi))
  }
  eva <- list(ev = ev, xev = xev, coef = coef, scale = z$scale, pc = z$wpc)
  class(eva) <- "lfeval"
  if(nc == 0) {
    cell <- list(sv = integer(0), ce = integer(0), s = integer(0), lo =
      as.integer(rep(0, nv)), hi = as.integer(rep(0, nv)))
  }
  else {
    mvc <- max(nv, nc)
    mvcm <- max(nvm, ncm)
    vc <- nvc["vc"]
    cell <- list(sv = t1[nvm * (3 * d + 8) + 1:nc], ce = t2[1:(vc * nc)], s =
      t2[vc * ncm + 1:mvc], lo = t2[vc * ncm + mvcm + 1:mvc], hi = t2[vc * ncm +
      2 * mvcm + 1:mvc])
  }
  ret <- list(eva = eva, cell = cell, terms = NULL, nvc = nvc, box = z$lim[2 *
    d + 1:(2 * d)], sty = style, deriv = deriv, mi = mi, dp = dp, trans = trans,
    critval = crit(const = c(rep(0, d), 1), d = d), vnames = vnames, yname =
    yname, call = match.call(), frame = sys.frame(sys.parent()))
  class(ret) <- "locfit"
  ret
}

"ang" <-
function(x, ...)
{
  ret <- lp(x, ..., style = "angle")
  dimnames(ret) <- list(NULL, deparse(substitute(x)))
  ret
}

"gam.lf"<-
function(x, y, w, xeval, ...)
{
  if(!missing(xeval)) {
    fit <- locfit.raw(x, y, weights = w, geth = 5, ...)
    return(predict(fit, xeval))
  }
  ret <- locfit.raw(x, y, weights = w, geth = 4, ...)
  names(ret) <- c("residuals", "var", "nl.df")
  ret
}

"gam.slist"<-
c("s", "lo", "random", "lf")

"lf"<-
function(..., alpha = 0.7, deg = 2, scale = 1, kern = "tcub", ev = rbox(), maxk
   = 100)
{
  if(!any(gam.slist == "lf"))
    warning("gam.slist does not include \"lf\" -- fit will be incorrect")
  x <- cbind(...)
  scall <- deparse(sys.call())
  attr(x, "alpha") <- alpha
  attr(x, "deg") <- deg
  attr(x, "scale") <- scale
  attr(x, "kern") <- kern
  attr(x, "ev") <- ev
  attr(x, "maxk") <- maxk
  attr(x, "call") <- substitute(gam.lf(data[[scall]], z, w, alpha = alpha, deg
     = deg, scale = scale, kern = kern, ev = ev, maxk = maxk))
  attr(x, "class") <- "smooth"
  x
}

"lfbas" <-
function(dim, indices, tt, ...)
{
  indices <- indices + 1
  # C starts at 0, S at 1
  x <- cbind(...)[indices,  ]
  res <- basis(x, tt)
  as.numeric(t(res))
}

"left"<-
function(x, ...)
{
  ret <- lp(x, ..., style = "left")
  dimnames(ret) <- list(NULL, deparse(substitute(x)))
  ret
}

"right"<-
function(x, ...)
{
  ret <- lp(x, ..., style = "right")
  dimnames(ret) <- list(NULL, deparse(substitute(x)))
  ret
}

"cpar"<-
function(x, ...)
{
  ret <- lp(x, ..., style = "cpar")
  dimnames(ret) <- list(NULL, deparse(substitute(x)))
  ret
}

"lp"<-
function(..., nn = 0, h = 0, adpen = 0, deg = 2, acri = "none", scale = FALSE,
  style = "none")
{
  x <- cbind(...)
  z <- as.list(match.call())
  z[[1]] <- z$nn <- z$h <- z$adpen <- z$deg <- z$acri <- z$scale <- z$style <-
    NULL
  dimnames(x) <- list(NULL, z)
  if(missing(nn) & missing(h) & missing(adpen))
    nn <- 0.7
  attr(x, "alpha") <- c(nn, h, adpen)
  attr(x, "deg") <- deg
  attr(x, "acri") <- acri
  attr(x, "style") <- style
  attr(x, "scale") <- scale
  class(x) <- c("lp", class(x))
  x
}


"[.lp" <- function (x, ..., drop = FALSE) {
    cl <- oldClass(x)
    oldClass(x) <- NULL
    ats <- attributes(x)
    ats$dimnames <- NULL
    ats$dim <- NULL
    ats$names <- NULL
    y <- x[..., drop = drop]
    attributes(y) <- c(attributes(y), ats)
    oldClass(y) <- cl
    y
}

"fitted.locfit"<-
function(object, data = NULL, what = "coef", cv = FALSE, studentize = FALSE,
  type = "fit", tr, ...)
{
  if(missing(data)) {
    data <- if(is.null(object$call$data)) sys.frame(sys.parent()) else eval(object$call$
        data)
  }
  if(missing(tr))
    tr <- if((what == "coef") & (type == "fit")) object$trans else function(x)
      x
  mm <- locfit.matrix(object, data = data)
  n <- object$mi["n"]
  pred <- .C("sfitted",
    x = as.numeric(mm$x),
    y = as.numeric(rep(mm$y, length.out = n)),
    w = as.numeric(rep(mm$w, length.out = n)),
    ce = as.numeric(rep(mm$ce, length.out = n)),
    ba = as.numeric(rep(mm$base, length.out = n)),
    fit = numeric(n),
    cv = as.integer(cv),
    st = as.integer(studentize),
    xev = as.numeric(object$eva$xev),
    coef = as.numeric(object$eva$coef),
    sv = as.numeric(object$cell$sv),
    ce = as.integer(c(object$cell$ce, object$cell$s, object$cell$lo, object$
      cell$hi)),
    wpc = as.numeric(object$eva$pc),
    scale = as.numeric(object$eva$scale),
    nvc = as.integer(object$nvc),
    mi = as.integer(object$mi),
    dp = as.numeric(object$dp),
    mg = as.integer(object$eva$ev$mg),
    deriv = as.integer(object$deriv),
    nd = as.integer(length(object$deriv)),
    sty = as.integer(object$sty),
    what = as.character(c(what, type)),
    basis = list(eval(object$call$basis)), PACKAGE="locfit")
  tr(pred$fit)
}

"formula.locfit"<-
function(x, ...)
x$call$formula

"predict.locfit"<-
function(object, newdata = NULL, where = "fitp", se.fit = FALSE, band = "none",
  what = "coef", ...)
{
  if((se.fit) && (band == "none"))
    band <- "global"
  for(i in 1:length(what)) {
    pred <- preplot.locfit(object, newdata, where = where, band = band, what =
      what[i], ...)
    fit <- pred$trans(pred$fit)
    if(i == 1)
      res <- fit
    else res <- cbind(res, fit)
  }
  if(band == "none")
    return(res)
  return(list(fit = res, se.fit = pred$se.fit, residual.scale = pred$
    residual.scale))
}

"lines.locfit"<-
function(x, m = 100, tr = x$trans, ...)
{
  newx <- lfmarg(x, m = m)[[1]]
  y <- predict(x, newx, tr = tr)
  lines(newx, y, ...)
}

"points.locfit"<-
function(x, tr, ...)
{
  d <- x$mi["d"]
  p <- x$mi["p"]
  nv <- x$nvc["nv"]
  if(d == 1) {
    if(missing(tr))
      tr <- x$trans
    x1 <- x$eva$xev
    x2 <- x$eva$coef[, 1]
    points(x1, tr(x2), ...)
  }
  if(d == 2) {
    xx <- lfknots(x, what = "x")
    points(xx[, 1], xx[, 2], ...)
  }
}

"print.locfit"<-
function(x, ...)
{
  if(!is.null(cl <- x$call)) {
    cat("Call:\n")
    dput(cl)
  }
  cat("\n")
  cat("Number of observations:         ", x$mi["n"], "\n")
  cat("Family: ", c("Density", "PP Rate", "Hazard", "Gaussian", "Logistic",
    "Poisson", "Gamma", "Geometric", "Circular", "Huber", "Robust Binomial",
    "Weibull", "Cauchy")[x$mi["tg"] %% 64], "\n")
  cat("Fitted Degrees of freedom:      ", round(x$dp["df2"], 3), "\n")
  cat("Residual scale:                 ", signif(sqrt(x$dp["rv"]), 3), "\n")
  invisible(x)
}

"residuals.locfit"<-
function(object, data = NULL, type = "deviance", ...)
{
  if(missing(data)) {
    data <- if(is.null(object$call$data)) sys.frame(sys.parent()) else eval(object$call$
        data)
  }
  fitted.locfit(object, data, ..., type = type)
}

"summary.locfit"<-
function(object, ...)
{
  mi <- object$mi
  fam <- c("Density Estimation", "Poisson process rate estimation",
    "Hazard Rate Estimation", "Local Regression", "Local Likelihood - Binomial",
    "Local Likelihood - Poisson", "Local Likelihood - Gamma",
    "Local Likelihood - Geometric", "Local Robust Regression")[mi["tg"] %% 64]
  estr <- c("Rectangular Tree", "Triangulation", "Data", "Rectangular Grid",
    "k-d tree", "k-d centres", "Cross Validation", "User-provided")[mi["ev"]]
  ret <- list(call = object$call, fam = fam, n = mi["n"], d = mi["d"], estr =
    estr, nv = object$nvc["nv"], deg = mi["deg"], dp = object$dp, vnames =
    object$vnames)
  class(ret) <- "summary.locfit"
  ret
}

"print.summary.locfit"<-
function(x, ...)
{
  cat("Estimation type:", x$fam, "\n")
  cat("\nCall:\n")
  print(x$call)
  cat("\nNumber of data points: ", x$n, "\n")
  cat("Independent variables: ", x$vnames, "\n")
  cat("Evaluation structure:", x$estr, "\n")
  cat("Number of evaluation points: ", x$nv, "\n")
  cat("Degree of fit: ", x$deg, "\n")
  cat("Fitted Degrees of Freedom: ", round(x$dp["df2"], 3), "\n")
  invisible(x)
}

"rbox"<-
function(cut = 0.8, type = "tree", ll = rep(0, 10), ur = rep(0, 10))
{
  if(!any(type == c("tree", "kdtree", "kdcenter", "phull")))
    stop("Invalid type argument")
  ret <- list(type = type, xev = 0, mg = 0, cut = as.numeric(cut), ll =
    as.numeric(ll), ur = as.numeric(ur))
  class(ret) <- "lf_evs"
  ret
}

"lfgrid"<-
function(mg = 10, ll = rep(0, 10), ur = rep(0, 10))
{
  if(length(mg) == 1)
    mg <- rep(mg, 10)
  ret <- list(type = "grid", xev = 0, mg = as.integer(mg), cut = 0, ll =
    as.numeric(ll), ur = as.numeric(ur))
  class(ret) <- "lf_evs"
  ret
}

"dat"<-
function(cv = FALSE)
{
  type <- if(cv) "crossval" else "data"
  ret <- list(type = type, xev = 0, mg = 0, cut = 0, ll = 0, ur = 0)
  class(ret) <- "lf_evs"
  ret
}

"xbar"<-
function()
{
  ret <- list(type = "xbar", xev = 0, mg = 0, cut = 0, ll = 0, ur = 0)
  class(ret) <- "lf_evs"
  ret
}

"none"<-
function()
{
  ret <- list(type = "none", xev = 0, mg = 0, cut = 0, ll = 0, ur = 0)
  class(ret) <- "lf_evs"
  ret
}

"plot.locfit"<-
function(x, xlim, pv, tv, m, mtv = 6, band = "none", tr = NULL, what = "coef",
  get.data = FALSE, f3d = (d == 2) && (length(tv) > 0), ...)
{
  d <- x$mi["d"]
  ev <- x$mi["ev"]
  where <- "grid"
  if(missing(pv))
    pv <- if(d == 1) 1 else c(1, 2)
  if(is.character(pv))
    pv <- match(pv, x$vnames)
  if(missing(tv))
    tv <- (1:d)[ - pv]
  if(is.character(tv))
    tv <- match(tv, x$vnames)
  vrs <- c(pv, tv)
  if(any(duplicated(vrs)))
    warning("Duplicated variables in pv, tv")
  if(any((vrs <= 0) | (vrs > d)))
    stop("Invalid variable numbers in pv, tv")
  if(missing(m))
    m <- if(d == 1) 100 else 40
  m <- rep(m, d)
  m[tv] <- mtv
  xl <- x$box
  if(!missing(xlim))
    xl <- lflim(xlim, x$vnames, xl)
  if((d != 2) & (any(ev == c(3, 7, 8))))
    pred <- preplot.locfit(x, where = "fitp", band = band, tr = tr, what = what,
      get.data = get.data, f3d = f3d)
  else {
    marg <- lfmarg(xl, m)
    pred <- preplot.locfit(x, marg, band = band, tr = tr, what = what, get.data
       = get.data, f3d = f3d)
  }
  plot(pred, pv = pv, tv = tv, ...)
}

"preplot.locfit"<-
function(object, newdata = NULL, where, tr = NULL, what = "coef", band = "none",
  get.data = FALSE, f3d = FALSE, ...)
{
  mi <- object$mi
  dim <- mi["d"]
  ev <- mi["ev"]
  nointerp <- any(ev == c(3, 7, 8))
  wh <- 1
  n <- 1
  if(is.null(newdata)) {
    if(missing(where))
      where <- if(nointerp) "fitp" else "grid"
    if(where == "grid")
      newdata <- lfmarg(object)
    if(any(where == c("fitp", "ev", "fitpoints"))) {
      where <- "fitp"
      newdata <- lfknots(object, what = "x", delete.pv = FALSE)
    }
    if(where == "data")
      newdata <- locfit.matrix(object)$x
    if(where == "vect")
      stop("you must give the vector points")
  }
  else {
    where <- "vect"
    if(is.data.frame(newdata))
      newdata <- as.matrix(model.frame(delete.response(object$terms), newdata))
    else if(is.list(newdata))
      where <- "grid"
    else newdata <- as.matrix(newdata)
  }
  if(is.null(tr)) {
    if(what == "coef")
      tr <- object$trans
    else tr <- function(x)
      x
  }
  if((nointerp) && (where == "grid") && (dim == 2)) {
    nv <- object$nvc["nv"]
    x <- object$eva$xev[2 * (1:nv) - 1]
    y <- object$eva$xev[2 * (1:nv)]
    z <- preplot.locfit.raw(object, 0, "fitp", what, band)$y
	haveAkima <- require(akima)
	if (! haveAkima) stop("The akima package is needed for the interp() function.  Please note its no-compercial-use license.")
    fhat <- akima::interp(x, y, z, newdata[[1]], newdata[[2]], ncp = 2)$z
  }
  else {
    z <- preplot.locfit.raw(object, newdata, where, what, band)
    fhat <- z$y
  }
  fhat[fhat == 0.1278433] <- NA
  band <- pmatch(band, c("none", "global", "local", "prediction"))
  if(band > 1)
    sse <- z$se
  else sse <- numeric(0)
  if(where != "grid")
    newdata <- list(xev = newdata, where = where)
  else newdata$where <- where
  data <- if(get.data) locfit.matrix(object) else list()
  if((f3d) | (dim > 3))
    dim <- 3
  ret <- list(xev = newdata, fit = fhat, se.fit = sse, residual.scale = sqrt(
    object$dp["rv"]), critval = object$critval, trans = tr, vnames = object$
    vnames, yname = object$yname, dim = as.integer(dim), data = data)
  class(ret) <- "preplot.locfit"
  ret
}

"preplot.locfit.raw"<-
function(object, newdata, where, what, band, ...)
{
  wh <- pmatch(where, c("vect", "grid", "data", "fitp"))
  switch(wh,
    {
      mg <- n <- nrow(newdata)
      xev <- newdata
    }
    ,
    {
      xev <- unlist(newdata)
      mg <- sapply(newdata, length)
      n <- prod(mg)
    }
    ,
    {
      mg <- n <- object$mi["n"]
      xev <- newdata
    }
    ,
    {
      mg <- n <- object$nvc["nv"]
      xev <- newdata
    }
    )
  .C("spreplot",
    xev = as.numeric(object$eva$xev),
    coef = as.numeric(object$eva$coef),
    sv = as.numeric(object$cell$sv),
    ce = as.integer(c(object$cell$ce, object$cell$s, object$cell$lo, object$
      cell$hi)),
    x = as.numeric(xev),
    y = numeric(n),
    se = numeric(n),
    wpc = as.numeric(object$eva$pc),
    scale = as.numeric(object$eva$scale),
    m = as.integer(mg),
    nvc = as.integer(object$nvc),
    mi = as.integer(object$mi),
    dp = as.numeric(object$dp),
    mg = as.integer(object$eva$ev$mg),
    deriv = as.integer(object$deriv),
    nd = as.integer(length(object$deriv)),
    sty = as.integer(object$sty),
    wh = as.integer(wh),
    what = c(what, band),
    bs = list(eval(object$call$basis)), PACKAGE="locfit")
}

"print.preplot.locfit"<-
function(x, ...)
{
  print(x$trans(x$fit))
  invisible(x)
}

"plot.locfit.1d"<-
function(x, add=FALSE, main="", xlab="default", ylab=x$yname, type="l",
         ylim, lty = 1, col = 1, ...) {
  y <- x$fit
  nos <- !is.na(y)
  xev <- x$xev[[1]][nos]
  y <- y[nos]
  ord <- order(xev)
  if(xlab == "default")
    xlab <- x$vnames
  tr <- x$trans
  yy <- tr(y)
  if(length(x$se.fit) > 0) {
    crit <- x$critval$crit.val
    cup <- tr((y + crit * x$se.fit))[ord]
    clo <- tr((y - crit * x$se.fit))[ord]
  }
  ndat <- 0
  if(length(x$data) > 0) {
    ndat <- nrow(x$data$x)
    xdsc <- rep(x$data$sc, length.out = ndat)
    xdyy <- rep(x$data$y, length.out = ndat)
    dok <- xdsc > 0
  }
  if(missing(ylim)) {
    if(length(x$se.fit) > 0)
      ylim <- c(min(clo), max(cup))
    else ylim <- range(yy)
    if(ndat > 0)
      ylim <- range(c(ylim, xdyy[dok]/xdsc[dok]))
  }
  if(!add) {
    plot(xev[ord], yy[ord], type = "n", xlab = xlab, ylab = ylab, main = main,
      xlim = range(x$xev[[1]]), ylim = ylim, ...)
  }
  lines(xev[ord], yy[ord], type = type, lty = lty, col = col)
  if(length(x$se.fit) > 0) {
    lines(xev[ord], cup, lty = 2)
    lines(xev[ord], clo, lty = 2)
  }
  if(ndat > 0) {
    xd <- x$data$x[dok]
    yd <- xdyy[dok]/xdsc[dok]
    cd <- rep(x$data$ce, length.out = ndat)[dok]
    if(length(x$data$y) < 2) {
      rug(xd[cd == 0])
      if(any(cd == 1))
        rug(xd[cd == 1], ticksize = 0.015)
    }
    else {
      plotbyfactor(xd, yd, cd, col = col, pch = c("o", "+"), add = TRUE)
    }
  }
  invisible(NULL)
}

"plot.locfit.2d"<-
function(x, type="contour", main, xlab, ylab, zlab=x$yname, ...)
{
  if(x$xev$where != "grid")
    stop("Can only plot from grids")
  if(missing(xlab))
    xlab <- x$vnames[1]
  if(missing(ylab))
    ylab <- x$vnames[2]
  tr <- x$trans
  m1 <- x$xev[[1]]
  m2 <- x$xev[[2]]
  y <- matrix(tr(x$fit))
  if(type == "contour")
    contour(m1, m2, matrix(y, nrow = length(m1)), ...)
  if(type == "image")
    image(m1, m2, matrix(y, nrow = length(m1)), ...)
  if((length(x$data) > 0) && any(type == c("contour", "image"))) {
    xd <- x$data$x
    ce <- rep(x$data$ce, length.out = nrow(xd))
    points(xd[ce == 0, 1], xd[ce == 0, 2], pch = "o")
    if(any(ce == 1))
      points(xd[ce == 1, 1], xd[ce == 1, 2], pch = "+")
  }
  if(type == "persp") {
    nos <- is.na(y)
    y[nos] <- min(y[!nos])
    persp(m1, m2, matrix(y, nrow = length(m1)), zlab=zlab, ...)
  }
  if(!missing(main))
    title(main = main)
  invisible(NULL)
}

"plot.locfit.3d"<-
function(x, main = "", pv, tv, type = "level", pred.lab = x$vnames, resp.lab =
  x$yname, crit = 1.96, ...)
{
  xev <- x$xev
  if(xev$where != "grid")
    stop("Can only plot from grids")
  xev$where <- NULL
  newx <- as.matrix(expand.grid(xev))
  newy <- x$trans(x$fit)
  wh <- rep("f", length(newy))
  if(length(x$data) > 0) {
    dat <- x$data
    for(i in tv) {
      m <- xev[[i]]
      dat$x[, i] <- m[1 + round((dat$x[, i] - m[1])/(m[2] - m[1]))]
    }
    newx <- rbind(newx, dat$x)
    if(is.null(dat$y))
      newy <- c(newy, rep(NA, nrow(dat$x)))
    else {
      newy <- c(newy, dat$y/dat$sc)
      newy[is.na(newy)] <- 0
    }
    wh <- c(wh, rep("d", nrow(dat$x)))
  }
  if(length(tv) == 0) {
    newdat <- data.frame(newy, newx[, pv])
    names(newdat) <- c("y", paste("pv", 1:length(pv), sep = ""))
  }
  else {
    newdat <- data.frame(newx[, tv], newx[, pv], newy)
    names(newdat) <- c(paste("tv", 1:length(tv), sep = ""), paste("pv", 1:
      length(pv), sep = ""), "y")
    for(i in 1:length(tv))
      newdat[, i] <- as.factor(signif(newdat[, i], 5))
  }
  loc.strip <- function(...)
  strip.default(..., strip.names = c(TRUE, TRUE), style = 1)
  if(length(pv) == 1) {
    clo <- cup <- numeric(0)
    if(length(x$se.fit) > 0) {
      if((!is.null(class(crit))) && (class(crit) == "kappa"))
        crit <- crit$crit.val
      cup <- x$trans((x$fit + crit * x$se.fit))
      clo <- x$trans((x$fit - crit * x$se.fit))
    }
    formula <- switch(1 + length(tv),
      y ~ pv1,
      y ~ pv1 | tv1,
      y ~ pv1 | tv1 * tv2,
      y ~ pv1 | tv1 * tv2 * tv3)
    pl <- xyplot(formula, xlab = pred.lab[pv], ylab = resp.lab, main = main,
      type = "l", cup = cup, wh = wh, panel = panel.xyplot.lf, data
       = newdat, strip = loc.strip, ...)
  }
  if(length(pv) == 2) {
    formula <- switch(1 + length(tv),
      y ~ pv1 * pv2,
      y ~ pv1 * pv2 | tv1,
      y ~ pv1 * pv2 | tv1 * tv2,
      y ~ pv1 * pv2 | tv1 * tv2 * tv3)
    if(type == "contour")
      pl <- contourplot(formula, xlab = pred.lab[pv[1]], ylab = pred.lab[pv[2]],
        main = main, data = newdat, strip = loc.strip, ...)
    if(type == "level")
      pl <- levelplot(formula, xlab = pred.lab[pv[1]], ylab = pred.lab[pv[2]],
        main = main, data = newdat, strip = loc.strip, ...)
    if((type == "persp") | (type == "wireframe"))
      pl <- wireframe(formula, xlab = pred.lab[pv[1]], ylab = pred.lab[pv[2]],
        zlab = resp.lab, data = newdat, strip = loc.strip, ...)
  }
  if(length(tv) > 0) {
    if(exists("is.R") && is.function(is.R) && is.R())
      names(pl$cond) <- pred.lab[tv]
    else names(attr(pl$glist, "endpts")) <- attr(pl$glist, "names") <- names(
        attr(pl$glist, "index")) <- pred.lab[tv]
  }
  pl
}

"panel.xyplot.lf"<-
function(x, y, subscripts, clo, cup, wh, type = "l", ...)
{
  wh <- wh[subscripts]
  panel.xyplot(x[wh == "f"], y[wh == "f"], type = type, ...)
  if(length(clo) > 0) {
    panel.xyplot(x[wh == "f"], clo[subscripts][wh == "f"], type = "l", lty = 2,
      ...)
    panel.xyplot(x[wh == "f"], cup[subscripts][wh == "f"], type = "l", lty = 2,
      ...)
  }
  if(any(wh == "d")) {
    yy <- y[wh == "d"]
    if(any(is.na(yy)))
      rug(x[wh == "d"])
    else panel.xyplot(x[wh == "d"], yy)
  }
}

"plot.preplot.locfit"<-
function(x, pv, tv, ...)
{
  if(x$dim == 1)
    plot.locfit.1d(x, ...)
  if(x$dim == 2)
    plot.locfit.2d(x, ...)
  if(x$dim >= 3)
    print(plot.locfit.3d(x, pv=pv, tv=tv, ...))
  invisible(NULL)
}

"summary.preplot.locfit"<-
function(object, ...)
object$trans(object$fit)


## Deepayan Sarkar's patched version:
"panel.locfit"<-
    function(x, y, subscripts, z,
             xyz.labs, xyz.axes, xyz.mid, xyz.minmax, xyz.range,
             col.regions, at, drape, contour, region, groups,
             ...)
{
    if(!missing(z)) {
        zs <- z[subscripts]
        fit <- locfit.raw(cbind(x, y), zs, ...)
        marg <- lfmarg(fit, m = 10)
        zp <- predict(fit, marg)
        if(!missing(contour)) {
            print("contour")
            print(range(zp))
            lattice::render.contour.trellis(marg[[1]], marg[[2]], zp, at = at)
        }
        else {
            loc.dat <-
                cbind(as.matrix(expand.grid(x = marg[[1]],
                                            y = marg[[1]])),
                      z = zp)
            lattice::render.3d.trellis(cbind(x = x, y = y, z = z[subscripts]),
                              type = "cloud",
                              xyz.labs = xyz.labs,
                              xyz.axes = xyz.axes,
                              xyz.mid = xyz.mid,
                              xyz.minmax = xyz.minmax,
                              xyz.range = xyz.range,
                              col.regions = col.regions,
                              at = at,
                              drape = drape)
        }
    }
    else {
        panel.xyplot(x, y, ...)
        args <- list(x = x, y = y, ...)
        ok <- names(formals(locfit.raw))
        llines.locfit(do.call("locfit.raw",
                              args[ok[ok %in% names(args)]]))
    }
}

llines.locfit <-
function (x, m = 100, tr = x$trans, ...)
{
    newx <- lfmarg(x, m = m)[[1]]
    y <- predict(x, newx, tr = tr)
    llines(newx, y, ...)
}

## "panel.locfit"<-
# function(x, y, subscripts, z, xyz.labs, xyz.axes, xyz.mid, xyz.minmax,
#   xyz.range, col.regions, at, drape, contour, region, groups, ...)
# {
#   if(!missing(z)) {
#     zs <- z[subscripts]
#     fit <- locfit.raw(cbind(x, y), zs, ...)
#     marg <- lfmarg(fit, m = 10)
#     zp <- predict(fit, marg)
#     if(!missing(contour)) {
#       print("contour")
#       print(range(zp))
#       render.contour.trellis(marg[[1]], marg[[2]], zp, at = at)
#     }
#     else {
#       loc.dat <- cbind(as.matrix(expand.grid(x = marg[[1]], y = marg[[1]])), z
#          = zp)
#       render.3d.trellis(cbind(x = x, y = y, z = z[subscripts]), type = "cloud",
#         xyz.labs = xyz.labs, xyz.axes = xyz.axes, xyz.mid = xyz.mid, xyz.minmax
#          = xyz.minmax, xyz.range = xyz.range, col.regions = col.regions, at =
#         at, drape = drape)
#     }
#   }
#   else {
#     panel.xyplot(x, y)
#     lines(locfit.raw(x, y, ...))
#   }
# }

"lfmarg"<-
function(xlim, m = 40)
{
  if(!is.numeric(xlim)) {
    d <- xlim$mi["d"]
    xlim <- xlim$box
  }
  else d <- length(m)
  marg <- vector("list", d)
  m <- rep(m, length.out = d)
  for(i in 1:d)
    marg[[i]] <- seq(xlim[i], xlim[i + d], length.out = m[i])
  marg
}

"lfeval"<-
function(object)
object$eva

"plot.lfeval"<-
function(x, add = FALSE, txt = FALSE, ...)
{
  if(class(x) == "locfit")
    x <- x$eva
  d <- length(x$scale)
  v <- matrix(x$xev, nrow = d)
  if(d == 1) {
    xx <- v[1,  ]
    y <- x$coef[, 1]
  }
  if(d == 2) {
    xx <- v[1,  ]
    y <- v[2,  ]
  }
  if(!add) {
    plot(xx, y, type = "n", ...)
  }
  points(xx, y, ...)
  if(txt)
    text(xx, y, (1:length(xx)) - 1)
  invisible(x)
}

"print.lfeval"<-
function(x, ...)
{
  if(class(x) == "locfit")
    x <- x$eva
  d <- length(x$scale)
  ret <- matrix(x$xev, ncol = d, byrow = TRUE)
  print(ret)
}

"lflim"<-
function(limits, nm, ret)
{
  d <- length(nm)
  if(is.numeric(limits))
    ret <- limits
  else {
    z <- match(nm, names(limits))
    for(i in 1:d)
      if(!is.na(z[i])) ret[c(i, i + d)] <- limits[[z[i]]]
  }
  as.numeric(ret)
}

"plot.eval"<-
function(x, add = FALSE, text = FALSE, ...)
{
  d <- x$mi["d"]
  v <- matrix(x$eva$xev, nrow = d)
  ev <- x$mi["ev"]
  pv <- if(any(ev == c(1, 2))) as.logical(x$cell$s) else rep(FALSE, ncol(v))
  if(!add) {
    plot(v[1,  ], v[2,  ], type = "n", xlab = x$vnames[1], ylab = x$vnames[2])
  }
  if(text)
    text(v[1,  ], v[2,  ], (1:x$nvc["nv"]) - 1)
  else {
    if(any(!pv))
      points(v[1, !pv], v[2, !pv], ...)
    if(any(pv))
      points(v[1, pv], v[2, pv], pch = "*", ...)
  }
  if(any(x$mi["ev"] == c(1, 2))) {
    zz <- .C("triterm",
      as.numeric(v),
      h = as.numeric(lfknots(x, what = "h", delete.pv = FALSE)),
      as.integer(x$cell$ce),
      lo = as.integer(x$cell$lo),
      hi = as.integer(x$cell$hi),
      as.numeric(x$eva$scale),
      as.integer(x$nvc),
      as.integer(x$mi),
      as.numeric(x$dp),
      nt = integer(1),
      term = integer(600),
      box = x$box, PACKAGE="locfit")
    ce <- zz$term + 1
  }
  else ce <- x$cell$ce + 1
  if(any(x$mi["ev"] == c(1, 5, 7))) {
    vc <- 2^d
    ce <- matrix(ce, nrow = vc)
    segments(v[1, ce[1,  ]], v[2, ce[1,  ]], v[1, ce[2,  ]], v[2, ce[2,  ]],
      ...)
    segments(v[1, ce[1,  ]], v[2, ce[1,  ]], v[1, ce[3,  ]], v[2, ce[3,  ]],
      ...)
    segments(v[1, ce[2,  ]], v[2, ce[2,  ]], v[1, ce[4,  ]], v[2, ce[4,  ]],
      ...)
    segments(v[1, ce[3,  ]], v[2, ce[3,  ]], v[1, ce[4,  ]], v[2, ce[4,  ]],
      ...)
  }
  if(any(x$mi["ev"] == c(2, 8))) {
    vc <- d + 1
    m <- matrix(ce, nrow = 3)
    segments(v[1, m[1,  ]], v[2, m[1,  ]], v[1, m[2,  ]], v[2, m[2,  ]], ...)
    segments(v[1, m[1,  ]], v[2, m[1,  ]], v[1, m[3,  ]], v[2, m[3,  ]], ...)
    segments(v[1, m[2,  ]], v[2, m[2,  ]], v[1, m[3,  ]], v[2, m[3,  ]], ...)
  }
  invisible(NULL)
}

"rv"<-
function(fit)
fit$dp["rv"]

"rv<-"<-
function(fit, value)
{
  fit$dp["rv"] <- value
  fit
}

"regband"<-
function(formula, what = c("CP", "GCV", "GKK", "RSW"), deg = 1, ...)
{
  m <- match.call()
  m$geth <- 3
  m$deg <- c(deg, 4)
  m$what <- NULL
  m$deriv <- match(what, c("CP", "GCV", "GKK", "RSW"))
  m[[1]] <- as.name("locfit")
  z <- eval(m, sys.frame(sys.parent()))
  names(z) <- what
  z[1:length(what)]
}

"kdeb"<-
function(x, h0 = 0.01 * sd, h1 = sd, meth = c("AIC", "LCV", "LSCV", "BCV",
  "SJPI", "GKK"), kern = "gauss", gf = 2.5)
{
  n <- length(x)
  sd <- sqrt(var(x))
  z <- .C("kdeb",
    x = as.numeric(x),
    mi = as.integer(n),
    band = numeric(length(meth)),
    ind = integer(n),
    h0 = as.numeric(gf * h0),
    h1 = as.numeric(gf * h1),
    meth = as.integer(match(meth, c("AIC", "LCV", "LSCV", "BCV", "SJPI", "GKK")
      )),
    nmeth = as.integer(length(meth)),
    kern = pmatch(kern, c("rect", "epan", "bisq", "tcub", "trwt", "gauss")),
    PACKAGE="locfit")
  band <- z$band
  names(band) <- meth
  band
}

"lfknots"<-
function(x, tr, what = c("x", "coef", "h", "nlx"), delete.pv = TRUE)
{
  nv <- x$nvc["nv"]
  d <- x$mi["d"]
  p <- x$mi["p"]
  z <- 0:(nv - 1)
  ret <- matrix(0, nrow = nv, ncol = 1)
  rname <- character(0)
  if(missing(tr))
    tr <- x$trans
  coef <- x$eva$coef
  for(wh in what) {
    if(wh == "x") {
      ret <- cbind(ret, matrix(x$eva$xev, ncol = d, byrow = TRUE))
      rname <- c(rname, x$vnames)
    }
    if(wh == "coef") {
      d0 <- coef[, 1]
      d0[d0 == 0.1278433] <- NA
      ret <- cbind(ret, tr(d0))
      rname <- c(rname, "mu hat")
    }
    if(wh == "f1") {
      ret <- cbind(ret, coef[, 1 + (1:d)])
      rname <- c(rname, paste("d", 1:d, sep = ""))
    }
    if(wh == "nlx") {
      ret <- cbind(ret, coef[, d + 2])
      rname <- c(rname, "||l(x)||")
    }
    if(wh == "nlx1") {
      ret <- cbind(ret, coef[, d + 2 + (1:d)])
      rname <- c(rname, paste("nlx-d", 1:d, sep = ""))
    }
    if(wh == "se") {
      ret <- cbind(ret, sqrt(x$dp["rv"]) * coef[, d + 2])
      rname <- c(rname, "StdErr")
    }
    if(wh == "infl") {
      z <- coef[, 2 * d + 3]
      ret <- cbind(ret, z * z)
      rname <- c(rname, "Influence")
    }
    if(wh == "infla") {
      ret <- cbind(ret, coef[, 2 * d + 3 + (1:d)])
      rname <- c(rname, paste("inf-d", 1:d, sep = ""))
    }
    if(wh == "lik") {
      ret <- cbind(ret, coef[, 3 * d + 3 + (1:3)])
      rname <- c(rname, c("LocLike", "fit.df", "res.df"))
    }
    if(wh == "h") {
      ret <- cbind(ret, coef[, 3 * d + 7])
      rname <- c(rname, "h")
    }
    if(wh == "deg") {
      ret <- cbind(ret, coef[, 3 * d + 8])
      rname <- c(rname, "deg")
    }
  }
  ret <- as.matrix(ret[, -1])
  if(nv == 1)
    ret <- t(ret)
  dimnames(ret) <- list(NULL, rname)
  if((delete.pv) && (any(x$mi["ev"] == c(1, 2))))
    ret <- ret[!as.logical(x$cell$s),  ]
  ret
}

"locfit.matrix"<-
function(fit, data)
{
  m <- fit$call
  n <- fit$mi["n"]
  y <- ce <- base <- 0
  w <- 1
  if(m[[1]] == "locfit.raw") {
    x <- as.matrix(eval(m$x, fit$frame))
    if(!is.null(m$y))
      y <- eval(m$y, fit$frame)
    if(!is.null(m$weights))
      w <- eval(m$weights, fit$frame)
    if(!is.null(m$cens))
      ce <- eval(m$cens, fit$frame)
    if(!is.null(m$base))
      base <- eval(m$base, fit$frame)
  }
  else {
    Terms <- terms(as.formula(m$formula))
    attr(Terms, "intercept") <- 0
    m[[1]] <- as.name("model.frame")
    z <- pmatch(names(m), c("formula", "data", "weights", "cens", "base",
      "subset"))
    for(i in length(z):2)
      if(is.na(z[i])) m[[i]] <- NULL
    frm <- eval(m, fit$frame)
    vnames <- as.character(attributes(Terms)$variables)[-1]
    if(attr(Terms, "response")) {
      y <- model.extract(frm, "response")
      vnames <- vnames[-1]
    }
    x <- as.matrix(frm[, vnames])
    if(any(names(m) == "weights"))
      w <- model.extract(frm, weights)
    if(any(names(m) == "cens"))
      ce <- model.extract(frm, "cens")
    if(any(names(m) == "base"))
      base <- model.extract(frm, base)
  }
  sc <- if(any((fit$mi["tg"] %% 64) == c(5:8, 11, 12))) w else 1
  list(x = x, y = y, w = w, sc = sc, ce = ce, base = base)
}

"expit"<-
function(x)
{
  y <- x
  ix <- (x < 0)
  y[ix] <- exp(x[ix])/(1 + exp(x[ix]))
  y[!ix] <- 1/(1 + exp( - x[!ix]))
  y
}

"plotbyfactor"<-
function(x, y, f, data, col = 1:10, pch = "O", add = FALSE, lg, xlab = deparse(
  substitute(x)), ylab = deparse(substitute(y)), log = "", ...)
{
  if(!missing(data)) {
    x <- eval(substitute(x), data)
    y <- eval(substitute(y), data)
    f <- eval(substitute(f), data)
  }
  f <- as.factor(f)
  if(!add)
    plot(x, y, type = "n", xlab = xlab, ylab = ylab, log = log, ...)
  lv <- levels(f)
  col <- rep(col, length.out = length(lv))
  pch <- rep(pch, length.out = length(lv))
  for(i in 1:length(lv)) {
    ss <- f == lv[i]
    if(any(ss))
      points(x[ss], y[ss], col = col[i], pch = pch[i])
  }
  if(!missing(lg))
    legend(lg[1], lg[2], legend = levels(f), col = col, pch = paste(pch,
      collapse = ""))
}

"hatmatrix"<-
function(formula, dc = TRUE, ...)
{
  m <- match.call()
  m$geth <- 1
  m[[1]] <- as.name("locfit")
  z <- eval(m, sys.frame(sys.parent()))
  nvc <- z[[2]]
  nvm <- nvc[1]
  nv <- nvc[4]
  matrix(z[[1]], ncol = nvm)[, 1:nv]
}

"locfit.robust"<-
function(x, y, weights, ..., iter = 3)
{
  m <- match.call()
  if((!is.numeric(x)) && (class(x) == "formula")) {
    m1 <- m[[1]]
    m[[1]] <- as.name("locfit")
    m$lfproc <- m1
    names(m)[[2]] <- "formula"
    return(eval(m, sys.frame(sys.parent())))
  }
  n <- length(y)
  lfr.wt <- rep(1, n)
  m[[1]] <- as.name("locfit.raw")
  for(i in 0:iter) {
    m$weights <- lfr.wt
    fit <- eval(m, sys.frame(sys.parent()))
    res <- residuals(fit, type = "raw")
    s <- median(abs(res))
    lfr.wt <- pmax(1 - (res/(6 * s))^2, 0)^2
  }
  fit
}

"locfit.censor"<-
function(x, y, cens, ..., iter = 3, km = FALSE)
{
  m <- match.call()
  if((!is.numeric(x)) && (class(x) == "formula")) {
    m1 <- m[[1]]
    m[[1]] <- as.name("locfit")
    m$lfproc <- m1
    names(m)[[2]] <- "formula"
    return(eval(m, sys.frame(sys.parent())))
  }
  lfc.y <- y
  cens <- as.logical(cens)
  m$cens <- m$iter <- m$km <- NULL
  m[[1]] <- as.name("locfit.raw")
  for (i in 0:iter) {
    m$y <- lfc.y
    fit <- eval(m, sys.frame(sys.parent()))
    fh <- fitted(fit)
    if(km) {
      sr <- y - fh
      lfc.y <- y + km.mrl(sr, cens)
    }
    else {
      rdf <- sum(1 - cens) - 2 * fit$dp["df1"] + fit$dp["df2"]
      sigma <- sqrt(sum((y - fh) * (lfc.y - fh))/rdf)
      sr <- (y - fh)/sigma
      lfc.y <- fh + (sigma * dnorm(sr))/pnorm( - sr)
    }
    lfc.y[!cens] <- y[!cens]
  }
  m$cens <- substitute(cens)
  m$y <- substitute(y)
  fit$call <- m
  fit
}

"km.mrl"<-
function(times, cens)
{
  n <- length(times)
  if(length(cens) != length(times))
    stop("times and cens must have equal length")
  ord <- order(times)
  times <- times[ord]
  cens <- cens[ord]
  n.alive <- n:1
  haz.km <- (1 - cens)/n.alive
  surv.km <- exp(cumsum(log(1 - haz.km[ - n])))
  int.surv <- c(diff(times) * surv.km)
  mrl.km <- c(rev(cumsum(rev(int.surv)))/surv.km, 0)
  mrl.km[!cens] <- 0
  mrl.km.ord <- numeric(n)
  mrl.km.ord[ord] <- mrl.km
  mrl.km.ord
}

"locfit.quasi"<-
function(x, y, weights, ..., iter = 3, var = abs)
{
  m <- match.call()
  if((!is.numeric(x)) && (class(x) == "formula")) {
    m1 <- m[[1]]
    m[[1]] <- as.name("locfit")
    m$lfproc <- m1
    names(m)[[2]] <- "formula"
    return(eval(m, sys.frame(sys.parent())))
  }
  n <- length(y)
  w0 <- lfq.wt <- if(missing(weights)) rep(1, n) else weights
  m[[1]] <- as.name("locfit.raw")
  for(i in 0:iter) {
    m$weights <- lfq.wt
    fit <- eval(m, sys.frame(sys.parent()))
    fh <- fitted(fit)
    lfq.wt <- w0/var(fh)
  }
  fit
}

"density.lf"<-
function(x, n=50, window="gaussian", width, from, to, cut=if(iwindow == 4)
         0.75 else 0.5, ev=lfgrid(mg=n, ll=from, ur=to), deg=0,
         family="density", link="ident", ...)
{
  if(!exists("logb"))
    logb <- log
  # for R
  x <- sort(x)
  r <- range(x)
  iwindow <- pmatch(window, c("rectangular", "triangular", "cosine", "gaussian"
    ), -1.)
  if(iwindow < 0.)
    kern <- window
  else kern <- c("rect", "tria", NA, "gauss")[iwindow]
  if(missing(width)) {
    nbar <- logb(length(x), base = 2.) + 1.
    width <- diff(r)/nbar * 0.5
  }
  if(missing(from))
    from <- r[1.] - width * cut
  if(missing(to))
    to <- r[2.] + width * cut
  if(to <= from)
    stop("Invalid from/to values")
  h <- width/2
  if(kern == "gauss")
    h <- h * 1.25
  fit <- locfit.raw(lp(x, h = h, deg = deg), ev = ev, kern = kern, link = link,
    family = family, ...)
  list(x = fit$eva$xev, y = fit$eva$coef[, 1])
}

"smooth.lf"<-
function(x, y, xev = x, direct = FALSE, ...)
{
  # just a simple smooth with (x,y) input, mu-hat output.
  # locfit.raw options are valid.
  if(missing(y)) {
    y <- x
    x <- 1:length(y)
  }
  if(direct) {
    fit <- locfit.raw(x, y, ev = xev, geth = 7, ...)
    fv <- fit$y
    xev <- fit$x
    if(is.matrix(x))
      xev <- matrix(xev, ncol = ncol(x), byrow = TRUE)
  }
  else {
    fit <- locfit.raw(x, y, ...)
    fv <- predict(fit, xev)
  }
  list(x = xev, y = fv, call = match.call())
}

"gcv"<-
function(x, ...)
{
  m <- match.call()
  if(is.numeric(x))
    m[[1]] <- as.name("locfit.raw")
  else {
    m[[1]] <- as.name("locfit")
    names(m)[2] <- "formula"
  }
  fit <- eval(m, sys.frame(sys.parent()))
  z <- fit$dp[c("lk", "df1", "df2")]
  n <- fit$mi["n"]
  z <- c(z, (-2 * n * z[1])/(n - z[2])^2)
  names(z) <- c("lik", "infl", "vari", "gcv")
  z
}

"gcvplot"<-
function(..., alpha, df = 2)
{
  m <- match.call()
  m[[1]] <- as.name("gcv")
  m$df <- NULL
  if(!is.matrix(alpha))
    alpha <- matrix(alpha, ncol = 1)
  k <- nrow(alpha)
  z <- matrix(nrow = k, ncol = 4)
  for(i in 1:k) {
    m$alpha <- alpha[i,  ]
    z[i,  ] <- eval(m, sys.frame(sys.parent()))
  }
  ret <- list(alpha = alpha, cri = "GCV", df = z[, df], values = z[, 4])
  class(ret) <- "gcvplot"
  ret
}

"plot.gcvplot"<-
function(x, xlab = "Fitted DF", ylab = x$cri, ...)
{
  plot(x$df, x$values, xlab = xlab, ylab = ylab, ...)
}

"print.gcvplot"<-
function(x, ...)
plot.gcvplot(x = x, ...)

"summary.gcvplot"<-
function(object, ...)
{
  z <- cbind(object$df, object$values)
  dimnames(z) <- list(NULL, c("df", object$cri))
  z
}

"aic"<-
function(x, ..., pen = 2)
{
  m <- match.call()
  if(is.numeric(x))
    m[[1]] <- as.name("locfit.raw")
  else {
    m[[1]] <- as.name("locfit")
    names(m)[2] <- "formula"
  }
  m$pen <- NULL
  fit <- eval(m, sys.frame(sys.parent()))
  dp <- fit$dp
  z <- dp[c("lk", "df1", "df2")]
  z <- c(z, -2 * z[1] + pen * z[2])
  names(z) <- c("lik", "infl", "vari", "aic")
  z
}

"aicplot"<-
function(..., alpha)
{
  m <- match.call()
  m[[1]] <- as.name("aic")
  if(!is.matrix(alpha))
    alpha <- matrix(alpha, ncol = 1)
  k <- nrow(alpha)
  z <- matrix(nrow = k, ncol = 4)
  for(i in 1:k) {
    m$alpha <- alpha[i,  ]
    z[i,  ] <- eval(m, sys.frame(sys.parent()))
  }
  ret <- list(alpha = alpha, cri = "AIC", df = z[, 2], values = z[, 4])
  class(ret) <- "gcvplot"
  ret
}

"cp"<-
function(x, ..., sig2 = 1)
{
  m <- match.call()
  if(is.numeric(x))
    m[[1]] <- as.name("locfit.raw")
  else {
    m[[1]] <- as.name("locfit")
    names(m)[2] <- "formula"
  }
  m$sig2 <- NULL
  fit <- eval(m, sys.frame(sys.parent()))
  z <- c(fit$dp[c("lk", "df1", "df2")], fit$mi["n"])
  z <- c(z, (-2 * z[1])/sig2 - z[4] + 2 * z[2])
  names(z) <- c("lik", "infl", "vari", "n", "cp")
  z
}

"cpplot"<-
function(..., alpha, sig2)
{
  m <- match.call()
  m[[1]] <- as.name("cp")
  m$sig2 <- NULL
  if(!is.matrix(alpha))
    alpha <- matrix(alpha, ncol = 1)
  k <- nrow(alpha)
  z <- matrix(nrow = k, ncol = 5)
  for(i in 1:k) {
    m$alpha <- alpha[i,  ]
    z[i,  ] <- eval(m, sys.frame(sys.parent()))
  }
  if(missing(sig2)) {
    s <- (1:k)[z[, 3] == max(z[, 3])][1]
    sig2 <- (-2 * z[s, 1])/(z[s, 4] - 2 * z[s, 2] + z[s, 3])
  }
  ret <- list(alpha = alpha, cri = "CP", df = z[, 3], values = (-2 * z[, 1])/
    sig2 - z[, 4] + 2 * z[, 2])
  class(ret) <- "gcvplot"
  ret
}

"lcv"<-
function(x, ...)
{
  m <- match.call()
  if(is.numeric(x))
    m[[1]] <- as.name("locfit.raw")
  else {
    m[[1]] <- as.name("locfit")
    names(m)[2] <- "formula"
  }
  fit <- eval(m, sys.frame(sys.parent()))
  z <- fit$dp[c("lk", "df1", "df2")]
  res <- residuals(fit, type = "d2", cv = TRUE)
  z <- c(z, sum(res))
  names(z) <- c("lik", "infl", "vari", "cv")
  z
}

"lcvplot"<-
function(..., alpha)
{
  m <- match.call()
  m[[1]] <- as.name("lcv")
  if(!is.matrix(alpha))
    alpha <- matrix(alpha, ncol = 1)
  k <- nrow(alpha)
  z <- matrix(nrow = k, ncol = 4)
  for(i in 1:k) {
    m$alpha <- alpha[i,  ]
    z[i,  ] <- eval(m, sys.frame(sys.parent()))
  }
  ret <- list(alpha = alpha, cri = "LCV", df = z[, 2], values = z[, 4])
  class(ret) <- "gcvplot"
  ret
}

"lscv"<-
function(x, ..., exact = FALSE)
{
  if(exact) {
    ret <- lscv.exact(x, ...)
  }
  else {
    m <- match.call()
    m$exact <- NULL
    if(is.numeric(x))
      m[[1]] <- as.name("locfit.raw")
    else {
      m[[1]] <- as.name("locfit")
      names(m)[2] <- "formula"
    }
    m$geth <- 6
    ret <- eval(m, sys.frame(sys.parent()))
  }
  ret
}

"lscv.exact"<-
function(x, h = 0)
{
  if(!is.null(attr(x, "alpha")))
    h <- attr(x, "alpha")[2]
  if(h <= 0)
    stop("lscv.exact: h must be positive.")
  ret <- .C("slscv",
    x = as.numeric(x),
    n = as.integer(length(x)),
    h = as.numeric(h),
    ret = numeric(2), PACKAGE="locfit")$ret
  ret
}

"lscvplot"<-
function(..., alpha)
{
  m <- match.call()
  m[[1]] <- as.name("lscv")
  if(!is.matrix(alpha))
    alpha <- matrix(alpha, ncol = 1)
  k <- nrow(alpha)
  z <- matrix(nrow = k, ncol = 2)
  for(i in 1:k) {
    m$alpha <- alpha[i,  ]
    z[i,  ] <- eval(m, sys.frame(sys.parent()))
  }
  ret <- list(alpha = alpha, cri = "LSCV", df = z[, 2], values = z[, 1])
  class(ret) <- "gcvplot"
  ret
}

"sjpi"<-
function(x, a)
{
  dnorms <- function(x, k)
  {
    if(k == 0)
      return(dnorm(x))
    if(k == 1)
      return( - x * dnorm(x))
    if(k == 2)
      return((x * x - 1) * dnorm(x))
    if(k == 3)
      return(x * (3 - x * x) * dnorm(x))
    if(k == 4)
      return((3 - x * x * (6 - x * x)) * dnorm(x))
    if(k == 6)
      return((-15 + x * x * (45 - x * x * (15 - x * x))) * dnorm(x))
    stop("k too large in dnorms")
  }
  alpha <- a * sqrt(2)
  n <- length(x)
  M <- outer(x, x, "-")
  s <- numeric(length(alpha))
  for(i in 1:length(alpha)) {
    s[i] <- sum(dnorms(M/alpha[i], 4))
  }
  s <- s/(n * (n - 1) * alpha^5)
  h <- (s * 2 * sqrt(pi) * n)^(-0.2)
  lambda <- diff(summary(x)[c(2, 5)])
  A <- 0.92 * lambda * n^(-1/7)
  B <- 0.912 * lambda * n^(-1/9)
  tb <-  - sum(dnorms(M/B, 6))/(n * (n - 1) * B^7)
  sa <- sum(dnorms(M/A, 4))/(n * (n - 1) * A^5)
  ah <- 1.357 * (sa/tb * h^5)^(1/7)
  cbind(h, a, ah/sqrt(2), s)
}

"scb"<-
function(x, ..., ev = lfgrid(20), simul = TRUE, type = 1)
{
    oc <- m <- match.call()
  if(is.numeric(x))
    m[[1]] <- as.name("locfit.raw")
  else {
    m[[1]] <- as.name("locfit")
    names(m)[2] <- "formula"
  }
  m$type <- m$simul <- NULL
  m$geth <- 70 + type + 10 * simul
  m$ev <- substitute(ev)
  fit <- eval(m, sys.frame(sys.parent()))
  fit$call <- oc
  class(fit) <- "scb"
  fit
}

"plot.scb"<-
function(x, add = FALSE, ...)
{
  fit <- x$trans(x$coef)
  lower <- x$trans(x$lower)
  upper <- x$trans(x$upper)
  d <- x$d
  if(d == 1)
    plot.scb.1d(x, fit, lower, upper, add, ...)
  if(d == 2)
    plot.scb.2d(x, fit = fit, lower = lower, upper = upper, ...)
  if(!any(d == c(1, 2)))
    stop("Can't plot this scb")
}

"plot.scb.1d"<-
function(x, fit, lower, upper, add = FALSE, style = "band", ...)
{
  if(style == "test") {
    lower <- lower - fit
    upper <- upper - fit
  }
  if(!add) {
    yl <- range(c(lower, fit, upper))
    plot(x$xev, fit, type = "l", ylim = yl, xlab = x$vnames[1])
  }
  lines(x$xev, lower, lty = 2)
  lines(x$xev, upper, lty = 2)
  if(is.null(x$call$deriv)) {
    dx <- x$data$x
    sc <- if(any((x$mi["tg"] %% 64) == c(5:8, 11, 12))) x$data$w else 1
    dy <- x$data$y/sc
    points(dx, dy)
  }
  if(style == "test")
    abline(h = 0, lty = 3)
}

"plot.scb.2d" <- function(x, fit, lower, upper, style = "tl", ylim, ...) {
    plot.tl <- function(x, y, z, nint = c(16, 15), v1, v2,
                        xlab=deparse(substitute(x)),
                        ylab=deparse(substitute(y)),
                        legend=FALSE, pch="", ...) {
        xl <- range(x)
        if (legend) {
            mar <- par()$mar
            if (mar[4] < 6.1)
                par(mar = c(mar[1:3], 6.1))
            on.exit(par(mar = mar))
            dlt <- diff(xl)
            xl[2] <- xl[2] + 0.02 * dlt
        }
        plot(1, 1, type = "n", xlim = xl, ylim = range(y), xlab = xlab,
             ylab = ylab, ...)
        nx <- length(x)
        ny <- length(y)
        if (missing(v)) {
          v <- seq(min(z) - 0.0001, max(z), length.out = nint + 1)
	    } else {
			nint <- length(v) - 1
		}
        ix <- rep(1:nx, ny)
        iy <- rep(1:ny, rep(nx, ny))
        r1 <- range(z[, 1])
        r2 <- range(z[, 2])
        hue <- if (missing(v1)) {
            floor((nint[1] * (z[, 1] - r1[1]))/(r1[2] - r1[1]) * 0.999999999)
        } else cut(z[, 1], v1) - 1
        sat <- if (missing(v2)) {
            floor((nint[2] * (z[, 2] - r2[1]))/(r2[2] - r2[1]) * 0.999999999)
        } else cut(z[, 2], v2) - 1
        col <- hue + nint[1] * sat + 1
        x <- c(2 * x[1] - x[2], x, 2 * x[nx] - x[nx - 1])
        y <- c(2 * y[1] - y[2], y, 2 * y[ny] - y[ny - 1])
        x <- (x[1:(nx + 1)] + x[2:(nx + 2)])/2
        y <- (y[1:(ny + 1)] + y[2:(ny + 2)])/2
        for (i in unique(col)) {
            u <- col == i
            if(pch == "") {
                xx <- rbind(x[ix[u]], x[ix[u] + 1], x[ix[u] + 1], x[ix[u]], NA)
                yy <- rbind(y[iy[u]], y[iy[u]], y[iy[u] + 1], y[iy[u] + 1], NA)
                polygon(xx, yy, col = i, border = 0)
            }
            else points(x[ix[u]], y[iy[u]], col = i, pch = pch)
        }
        if(legend) {
            yv <- seq(min(y), max(y), length = length(v))
            x1 <- max(x) + 0.02 * dlt
            x2 <- max(x) + 0.06 * dlt
            for(i in 1:nint) {
                polygon(c(x1, x2, x2, x1), rep(yv[i:(i + 1)], c(2, 2)),
                        col = i, border = 0)
            }
            axis(side = 4, at = yv, labels = v, adj = 0)
        }
    }

    if(style == "trell") {
        if(missing(ylim))
            ylim <- range(c(fit, lower, upper))
        loc.dat = data.frame(x1 = x$xev[, 1], x2 = x$xev[, 2], y = fit)
        pl <- xyplot(y ~ x1 | as.factor(x2), data = loc.dat,
                     panel = panel.xyplot.lf, clo=lower, cup=upper,
                     wh=rep("f", nrow(loc.dat)))
        plot(pl)
    }
    if(style == "tl") {
        ux <- unique(x$xev[, 1])
        uy <- unique(x$xev[, 2])
        sig <- abs(x$coef/x$sd)
        rv1 <- max(abs(fit)) * 1.0001
        v1 <- seq( - rv1, rv1, length.out = 17)
        v2 <-  - c(-1e-100, crit(const = x$kap, cov = c(0.5, 0.7, 0.8, 0.85,
                                                0.9, 0.95, 0.98, 0.99, 0.995,
                                                0.999, 0.9999))$crit.val,
                   1e+300)
        plot.tl(ux, uy, cbind(fit,  - sig), v1 = v1, v2 = v2,
                xlab = x$vnames[1], ylab = x$vnames[2])
    }
}

"print.scb"<-
function(x, ...)
{
  m <- cbind(x$xev, x$trans(x$coef), x$trans(x$lower), x$trans(x$upper))
  dimnames(m) <- list(NULL, c(x$vnames, "fit", "lower", "upper"))
  print(m)
}

"kappa0"<-
    function(formula, cov=0.95, ev=lfgrid(20), ...)
{
  if(class(formula) == "locfit") {
    m <- formula$call
  }
  else {
    m <- match.call()
    m$cov <- NULL
  }
  m$dc <- TRUE
  m$geth <- 2
  m$ev <- substitute(ev)
  m[[1]] <- as.name("locfit")
  z <- eval(m, sys.frame(sys.parent()))
  crit(const = z$const, d = z$d, cov = cov)
}

"crit"<-
function(fit, const = c(0, 1), d = 1, cov = 0.95, rdf = 0)
{
  if(!missing(fit)) {
    z <- fit$critval
    if(missing(const) & missing(d) & missing(cov))
      return(z)
    if(!missing(const))
      z$const <- const
    if(!missing(d))
      z$d <- d
    if(!missing(cov))
      z$cov <- cov
    if(!missing(rdf))
      z$rdf <- rdf
  }
  else {
    z <- list(const = const, d = d, cov = cov, rdf = rdf, crit.val = 0)
    class(z) <- "kappa"
  }
  z$crit.val <- .C("scritval",
    k0 = as.numeric(z$const),
    d = as.integer(z$d),
    cov = as.numeric(z$cov),
    m = as.integer(length(z$const)),
    rdf = as.numeric(z$rdf),
    x = numeric(1),
    k = as.integer(1), PACKAGE="locfit")$x
  z
}

"crit<-"<-
function(fit, value)
{
  if(is.numeric(value))
    fit$critval$crit.val <- value[1]
  else {
    if(class(value) != "kappa")
      stop("crit<-: value must be numeric or class kappa")
    fit$critval <- value
  }
  fit
}


"spence.15"<-
function(y)
{
  n <- length(y)
  y <- c(rep(y[1], 7), y, rep(y[n], 7))
  n <- length(y)
  k <- 3:(n - 2)
  a3 <- y[k - 1] + y[k] + y[k + 1]
  a2 <- y[k - 2] + y[k + 2]
  y1 <- y[k] + 3 * (a3 - a2)
  n <- length(y1)
  k <- 1:(n - 3)
  y2 <- y1[k] + y1[k + 1] + y1[k + 2] + y1[k + 3]
  n <- length(y2)
  k <- 1:(n - 3)
  y3 <- y2[k] + y2[k + 1] + y2[k + 2] + y2[k + 3]
  n <- length(y3)
  k <- 1:(n - 4)
  y4 <- y3[k] + y3[k + 1] + y3[k + 2] + y3[k + 3] + y3[k + 4]
  y4/320
}

"spence.21"<-
function(y)
{
  n <- length(y)
  y <- c(rep(y[1], 10), y, rep(y[n], 10))
  n <- length(y)
  k <- 4:(n - 3)
  y1 <-  - y[k - 3] + y[k - 1] + 2 * y[k] + y[k + 1] - y[k + 3]
  n <- length(y1)
  k <- 4:(n - 3)
  y2 <- y1[k - 3] + y1[k - 2] + y1[k - 1] + y1[k] + y1[k + 1] + y1[k + 2] + y1[
    k + 3]
  n <- length(y2)
  k <- 3:(n - 2)
  y3 <- y2[k - 2] + y2[k - 1] + y2[k] + y2[k + 1] + y2[k + 2]
  n <- length(y3)
  k <- 3:(n - 2)
  y4 <- y3[k - 2] + y3[k - 1] + y3[k] + y3[k + 1] + y3[k + 2]
  y4/350
}

"store"<-
function(data = FALSE, grand = FALSE)
{
  lfmod <- c("ang", "gam.lf", "gam.slist", "lf", "lfbas", "left", "right",
    "cpar", "lp")
  lfmeth <- c("fitted.locfit", "formula.locfit", "predict.locfit",
    "lines.locfit", "points.locfit", "print.locfit", "residuals.locfit",
    "summary.locfit", "print.summary.locfit")
  lfev <- c("rbox", "gr", "dat", "xbar", "none")
  lfplo <- c("plot.locfit", "preplot.locfit", "preplot.locfit.raw",
    "print.preplot.locfit", "plot.locfit.1d", "plot.locfit.2d",
    "plot.locfit.3d", "panel.xyplot.lf", "plot.preplot.locfit",
    "summary.preplot.locfit", "panel.locfit", "lfmarg")
  lffre <- c("hatmatrix", "locfit.robust", "locfit.censor", "km.mrl",
    "locfit.quasi", "density.lf", "smooth.lf")
  lfscb <- c("scb", "plot.scb", "plot.scb.1d", "plot.scb.2d", "print.scb",
    "kappa0", "crit", "crit<-", "plot.tl")
  lfgcv <- c("gcv", "gcvplot", "plot.gcvplot", "print.gcvplot",
    "summary.gcvplot", "aic", "aicplot", "cp", "cpplot", "lcv", "lcvplot",
    "lscv", "lscv.exact", "lscvplot", "sjpi")
  lfspen <- c("spence.15", "spence.21")
  lffuns <- c("locfit", "locfit.raw", lfmod, lfmeth, lfev, lfplo, "lfeval",
    "plot.lfeval", "print.lfeval", "lflim", "plot.eval", "rv", "rv<-",
    "regband", "kdeb", "lfknots", "locfit.matrix", "expit", "plotbyfactor",
    lffre, lfgcv, lfscb, lfspen, "store")
  lfdata <- c("bad", "cltest", "cltrain", "co2", "diab", "geyser", "ethanol",
    "mcyc", "morths", "border", "heart", "trimod", "insect", "iris", "spencer",
    "stamp")
  lfgrand <- c("locfit.raw", "crit", "predict.locfit", "preplot.locfit",
    "preplot.locfit.raw", "lfbas", "expit", "rv", "rv<-", "knots")
  dump(lffuns, "S/locfit.s")
  if(data)
    dump(lfdata, "S/locfit.dat")
  if(grand)
    dump(lfgrand, "src-gr/lfgrand.s")
  dump(lffuns, "R/locfit.s")
}

