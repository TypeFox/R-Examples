
fbvpot <- function(x, threshold, model = c("log", "bilog", "alog", "neglog", "negbilog", "aneglog", "ct", "hr", "amix"), likelihood = c("censored","poisson"), start, ..., sym = FALSE, cshape = cscale, cscale = FALSE, std.err = TRUE, corr = FALSE, method = "BFGS", warn.inf = TRUE) {
  call <- match.call()
  likelihood <- match.arg(likelihood)
  ft <- switch(likelihood,
    censored = fbvcpot(x = x, u = threshold, model = model, start = start, ...,
      sym = sym, cshape = cshape, cscale = cscale, std.err =
      std.err, corr = corr, method = method, warn.inf = warn.inf),
    poisson = fbvppot(x = x, u = threshold, model = model, start = start, ...,
      sym = sym, cshape = cshape, cscale = cscale, std.err =
      std.err, corr = corr, method = method, warn.inf = warn.inf))
  structure(c(ft, call = call), class = c("bvpot", "evd"))
}

fbvcpot <- function(x, u, model = c("log", "bilog", "alog", "neglog", "negbilog", "aneglog", "ct", "hr", "amix"), start, ..., sym = FALSE, cshape = cscale, cscale = FALSE, std.err = TRUE, corr = FALSE, method = "BFGS", warn.inf = TRUE) {
  model <- match.arg(model)
  if(sym && !(model %in% c("alog","aneglog","ct")))
    warning("Argument `sym' was ignored")
  switch(model,
    log = fbvclog(x = x, u = u, start = start, ..., sym = FALSE,
      cshape = cshape, cscale = cscale, std.err = std.err,
      corr = corr, method = method, warn.inf = warn.inf),
    bilog = fbvcbilog(x = x, u = u, start = start, ..., sym = FALSE,
      cshape = cshape, cscale = cscale, std.err = std.err,
      corr = corr, method = method, warn.inf = warn.inf),
    alog = fbvcalog(x = x, u = u, start = start, ..., sym = sym,
      cshape = cshape, cscale = cscale, std.err = std.err,
      corr = corr, method = method, warn.inf = warn.inf),
    neglog = fbvcneglog(x = x, u = u, start = start, ..., sym = FALSE,
      cshape = cshape, cscale = cscale, std.err = std.err,
      corr = corr, method = method, warn.inf = warn.inf),
    negbilog = fbvcnegbilog(x = x, u = u, start = start, ..., sym = FALSE,
      cshape = cshape, cscale = cscale, std.err = std.err,
      corr = corr, method = method, warn.inf = warn.inf),
    aneglog = fbvcaneglog(x = x, u = u, start = start, ..., sym = sym,
      cshape = cshape, cscale = cscale, std.err = std.err,
      corr = corr, method = method, warn.inf = warn.inf),
    ct = fbvcct(x = x, u = u, start = start, ..., sym = sym,
      cshape = cshape, cscale = cscale, std.err = std.err,
      corr = corr, method = method, warn.inf = warn.inf),
    hr = fbvchr(x = x, u = u, start = start, ..., sym = FALSE,
      cshape = cshape, cscale = cscale, std.err = std.err,
      corr = corr, method = method, warn.inf = warn.inf),
    amix = fbvcamix(x = x, u = u, start = start, ..., sym = FALSE,
      cshape = cshape, cscale = cscale, std.err = std.err,
      corr = corr, method = method, warn.inf = warn.inf))
}

fbvppot <- function(x, u, model = c("log", "bilog", "alog", "neglog", "negbilog", "aneglog", "ct", "hr", "amix"), start, ..., sym = FALSE, cshape = cscale, cscale = FALSE, std.err = TRUE, corr = FALSE, method = "BFGS", warn.inf = TRUE) {
  model <- match.arg(model)
  if(model %in% c("alog","aneglog","amix"))
    stop("This model is not appropriate for poisson likelihood")
  if(sym && (model != "ct"))
    warning("Argument `sym' was ignored")
  switch(model,
    log = fbvplog(x = x, u = u, start = start, ..., sym = FALSE,
      cshape = cshape, cscale = cscale, std.err = std.err,
      corr = corr, method = method, warn.inf = warn.inf),
    bilog = fbvpbilog(x = x, u = u, start = start, ..., sym = FALSE,
      cshape = cshape, cscale = cscale, std.err = std.err,
      corr = corr, method = method, warn.inf = warn.inf),
    neglog = fbvpneglog(x = x, u = u, start = start, ..., sym = FALSE,
      cshape = cshape, cscale = cscale, std.err = std.err,
      corr = corr, method = method, warn.inf = warn.inf),
    negbilog = fbvpnegbilog(x = x, u = u, start = start, ..., sym = FALSE,
      cshape = cshape, cscale = cscale, std.err = std.err,
      corr = corr, method = method, warn.inf = warn.inf),
    ct = fbvpct(x = x, u = u, start = start, ..., sym = sym,
      cshape = cshape, cscale = cscale, std.err = std.err,
      corr = corr, method = method, warn.inf = warn.inf),
    hr = fbvphr(x = x, u = u, start = start, ..., sym = FALSE,
      cshape = cshape, cscale = cscale, std.err = std.err,
      corr = corr, method = method, warn.inf = warn.inf))
}

### Censored Likelihood Fitting ###

fbvclog <- function(x, u, start, ..., sym = FALSE, cshape = cscale, cscale = FALSE, std.err = TRUE, corr = FALSE, method = "BFGS", warn.inf = TRUE) {

  nllbvclog <- function(scale1, shape1, scale2, shape2, dep) {
    if(cshape) shape2 <- shape1
    if(cscale) scale2 <- scale1
    .C("nllbvclog", spx$x1, spx$x2, spx$nn, spx$n, spx$thdi,
      spx$lambda, dep, scale1, shape1, scale2, shape2, dns = double(1),
      PACKAGE = "evd")$dns
  }
  param <- c("scale1", "shape1")
  if(!cscale) param <- c(param, "scale2")
  if(!cshape) param <- c(param, "shape2")
  param <- c(param, "dep")
  nmdots <- names(list(...))
  start <- bvstart.vals(x = x, start = start, nmdots = nmdots,
    param = param, method = "pot", u = u, model = "log")
  spx <- sep.bvdata(x = x, method = "cpot", u = u)  
  nm <- names(start)
  l <- length(nm)
  fixed.param <- list(...)[nmdots %in% param]
  if(any(!(param %in% c(nm,names(fixed.param)))))
    stop("unspecified parameters")
  prind <- c(TRUE, TRUE, !cscale, !cshape, TRUE)
  f <- formals(nllbvclog)[prind]
  names(f) <- param
  m <- match(nm, param)
  if(any(is.na(m))) 
    stop("`start' specifies unknown arguments")     
  formals(nllbvclog) <- c(f[m], f[-m])
  nll <- function(p, ...) nllbvclog(p, ...)
  if(l > 1) {
    body(nll) <- parse(text = paste("nllbvclog(", paste("p[", 1:l, "]", collapse=", "), ", ...)"))
  }
  start.arg <- c(list(p = unlist(start)), fixed.param)
  if(warn.inf && do.call("nll", start.arg) == 1e+06)
    warning("negative log-likelihood is infinite at starting values")
  opt <- optim(start, nll, hessian = std.err, ..., method = method)
  cmar <- c(cscale, cshape); nat <- spx$nat
  bvpost.optim(x = x, opt = opt, nm = nm, fixed.param = fixed.param,
    std.err = std.err, corr = corr, sym = sym, cmar = cmar, method = "pot",
    u = u, nat = nat, likelihood = "censored", model = "log")
}

fbvcbilog <- function(x, u, start, ..., sym = FALSE, cshape = cscale, cscale = FALSE, std.err = TRUE, corr = FALSE, method = "BFGS", warn.inf = TRUE) {

  nllbvcbilog <- function(scale1, shape1, scale2, shape2, alpha, beta) {
    if(cshape) shape2 <- shape1
    if(cscale) scale2 <- scale1
    .C("nllbvcbilog", spx$x1, spx$x2, spx$nn, spx$n, spx$thdi, spx$lambda,
      alpha, beta, scale1, shape1, scale2, shape2, dns = double(1),
      PACKAGE = "evd")$dns
  }
  param <- c("scale1", "shape1")
  if(!cscale) param <- c(param, "scale2")
  if(!cshape) param <- c(param, "shape2")
  param <- c(param, "alpha", "beta")
  nmdots <- names(list(...))
  start <- bvstart.vals(x = x, start = start, nmdots = nmdots,
    param = param, method = "pot", u = u, model = "bilog")
  spx <- sep.bvdata(x = x, method = "cpot", u = u) 
  nm <- names(start)
  l <- length(nm)
  fixed.param <- list(...)[nmdots %in% param]
  if(any(!(param %in% c(nm,names(fixed.param)))))
    stop("unspecified parameters")
  prind <- c(TRUE, TRUE, !cscale, !cshape, TRUE, TRUE)
  f <- formals(nllbvcbilog)[prind]
  names(f) <- param
  m <- match(nm, param)
  if(any(is.na(m))) 
    stop("`start' specifies unknown arguments")   
  formals(nllbvcbilog) <- c(f[m], f[-m])
  nll <- function(p, ...) nllbvcbilog(p, ...)
  if(l > 1) {
    body(nll) <- parse(text = paste("nllbvcbilog(", paste("p[", 1:l,
      "]", collapse=", "), ", ...)"))
  }
  start.arg <- c(list(p = unlist(start)), fixed.param)
  if(warn.inf && do.call("nll", start.arg) == 1e+06)
    warning("negative log-likelihood is infinite at starting values")
  opt <- optim(start, nll, hessian = std.err, ..., method = method)
  cmar <- c(cscale, cshape); nat <- spx$nat
  bvpost.optim(x = x, opt = opt, nm = nm, fixed.param = fixed.param,
    std.err = std.err, corr = corr, sym = sym, cmar = cmar, method = "pot",
    u = u, nat = nat, likelihood = "censored", model = "bilog")
}

fbvcalog <- function(x, u, start, ..., sym = FALSE, cshape = cscale, cscale = FALSE, std.err = TRUE, corr = FALSE, method = "BFGS", warn.inf = TRUE) {
  nllbvcalog <- function(scale1, shape1, scale2, shape2, asy1, asy2, dep) {
    if(sym) asy2 <- asy1
    if(cshape) shape2 <- shape1
    if(cscale) scale2 <- scale1
    .C("nllbvcalog", spx$x1, spx$x2, spx$nn, spx$n, spx$thdi, spx$lambda,
       dep, asy1, asy2, scale1, shape1, scale2, shape2, dns = double(1),
       PACKAGE = "evd")$dns
  }
  param <- c("scale1", "shape1")
  if(!cscale) param <- c(param, "scale2")
  if(!cshape) param <- c(param, "shape2")
  if(!sym) param <- c(param, "asy1", "asy2", "dep")
  else param <- c(param, "asy1", "dep")
  nmdots <- names(list(...))
  start <- bvstart.vals(x = x, start = start, nmdots = nmdots,
    param = param, method = "pot", u = u, model = "alog")
  spx <- sep.bvdata(x = x, method = "cpot", u = u)
  nm <- names(start)
  l <- length(nm)
  fixed.param <- list(...)[nmdots %in% param]
  if(any(!(param %in% c(nm,names(fixed.param)))))
    stop("unspecified parameters")
  prind <- c(TRUE, TRUE, !cscale, !cshape, TRUE, !sym, TRUE)
  f <- formals(nllbvcalog)[prind]
  names(f) <- param
  m <- match(nm, param)
  if(any(is.na(m))) 
    stop("`start' specifies unknown arguments") 
  formals(nllbvcalog) <- c(f[m], f[-m])
  nll <- function(p, ...) nllbvcalog(p, ...)
  if(l > 1) {
    body(nll) <- parse(text = paste("nllbvcalog(", paste("p[", 1:l, "]", collapse=", "), ", ...)"))
  }
  start.arg <- c(list(p = unlist(start)), fixed.param)
  if(warn.inf && do.call("nll", start.arg) == 1e+06)
    warning("negative log-likelihood is infinite at starting values")
  opt <- optim(start, nll, hessian = std.err, ..., method = method)
  cmar <- c(cscale, cshape); nat <- spx$nat
  bvpost.optim(x = x, opt = opt, nm = nm, fixed.param = fixed.param,
    std.err = std.err, corr = corr, sym = sym, cmar = cmar, method = "pot",
    u = u, nat = nat, likelihood = "censored", model = "alog")
}

fbvcneglog <- function(x, u, start, ..., sym = FALSE, cshape = cscale, cscale = FALSE, std.err = TRUE, corr = FALSE, method = "BFGS", warn.inf = TRUE) {

  nllbvcneglog <- function(scale1, shape1, scale2, shape2, dep) {
    if(cshape) shape2 <- shape1
    if(cscale) scale2 <- scale1
    .C("nllbvcneglog", spx$x1, spx$x2, spx$nn, spx$n, spx$thdi, spx$lambda,
      dep, scale1, shape1, scale2, shape2, dns = double(1),
      PACKAGE = "evd")$dns
  }
  param <- c("scale1", "shape1")
  if(!cscale) param <- c(param, "scale2")
  if(!cshape) param <- c(param, "shape2")
  param <- c(param, "dep")
  nmdots <- names(list(...))
  start <- bvstart.vals(x = x, start = start, nmdots = nmdots,
    param = param, method = "pot", u = u, model = "neglog")
  spx <- sep.bvdata(x = x, method = "cpot", u = u) 
  nm <- names(start)
  l <- length(nm)
  fixed.param <- list(...)[nmdots %in% param]
  if(any(!(param %in% c(nm,names(fixed.param)))))
    stop("unspecified parameters")
  prind <- c(TRUE, TRUE, !cscale, !cshape, TRUE)
  f <- formals(nllbvcneglog)[prind]
  names(f) <- param
  m <- match(nm, param)
  if(any(is.na(m))) 
    stop("`start' specifies unknown arguments") 
  formals(nllbvcneglog) <- c(f[m], f[-m])
  nll <- function(p, ...) nllbvcneglog(p, ...)
  if(l > 1) {
    body(nll) <- parse(text = paste("nllbvcneglog(", paste("p[", 1:l, "]", collapse=", "), ", ...)"))
  }
  start.arg <- c(list(p = unlist(start)), fixed.param)
  if(warn.inf && do.call("nll", start.arg) == 1e+06)
    warning("negative log-likelihood is infinite at starting values")
  opt <- optim(start, nll, hessian = std.err, ..., method = method)
  cmar <- c(cscale, cshape); nat <- spx$nat
  bvpost.optim(x = x, opt = opt, nm = nm, fixed.param = fixed.param,
    std.err = std.err, corr = corr, sym = sym, cmar = cmar, method = "pot",
    u = u, nat = nat, likelihood = "censored", model = "neglog")
}

fbvcnegbilog <- function(x, u, start, ..., sym = FALSE, cshape = cscale, cscale = FALSE, std.err = TRUE, corr = FALSE, method = "BFGS", warn.inf = TRUE) {

  nllbvcnegbilog <- function(scale1, shape1, scale2, shape2, alpha, beta) {
    if(cshape) shape2 <- shape1
    if(cscale) scale2 <- scale1
    .C("nllbvcnegbilog", spx$x1, spx$x2, spx$nn, spx$n, spx$thdi, spx$lambda,
       alpha, beta, scale1, shape1, scale2, shape2, dns = double(1),
       PACKAGE = "evd")$dns
  }
  param <- c("scale1", "shape1")
  if(!cscale) param <- c(param, "scale2")
  if(!cshape) param <- c(param, "shape2")
  param <- c(param, "alpha", "beta")
  nmdots <- names(list(...))
  start <- bvstart.vals(x = x, start = start, nmdots = nmdots,
    param = param, method = "pot", u = u, model = "negbilog")
  spx <- sep.bvdata(x = x, method = "cpot", u = u) 
  nm <- names(start)
  l <- length(nm)
  fixed.param <- list(...)[nmdots %in% param]
  if(any(!(param %in% c(nm,names(fixed.param)))))
    stop("unspecified parameters")
  prind <- c(TRUE, TRUE, !cscale, !cshape, TRUE, TRUE)
  f <- formals(nllbvcnegbilog)[prind]
  names(f) <- param
  m <- match(nm, param)
  if(any(is.na(m))) 
    stop("`start' specifies unknown arguments") 
  formals(nllbvcnegbilog) <- c(f[m], f[-m])
  nll <- function(p, ...) nllbvcnegbilog(p, ...)
  if(l > 1) {
    body(nll) <- parse(text = paste("nllbvcnegbilog(", paste("p[", 1:l,
      "]", collapse=", "), ", ...)"))
  }
  start.arg <- c(list(p = unlist(start)), fixed.param)
  if(warn.inf && do.call("nll", start.arg) == 1e+06)
    warning("negative log-likelihood is infinite at starting values")
  opt <- optim(start, nll, hessian = std.err, ..., method = method)
  cmar <- c(cscale, cshape); nat <- spx$nat
  bvpost.optim(x = x, opt = opt, nm = nm, fixed.param = fixed.param,
    std.err = std.err, corr = corr, sym = sym, cmar = cmar, method = "pot",
    u = u, nat = nat, likelihood = "censored", model = "negbilog")
}

fbvcaneglog <- function(x, u, start, ..., sym = FALSE, cshape = cscale, cscale = FALSE, std.err = TRUE, corr = FALSE, method = "BFGS", warn.inf = TRUE) {
  
  nllbvcaneglog <- function(scale1, shape1, scale2, shape2, asy1, asy2, dep) {
    if(sym) asy2 <- asy1
    if(cshape) shape2 <- shape1
    if(cscale) scale2 <- scale1
    .C("nllbvcaneglog", spx$x1, spx$x2, spx$nn, spx$n, spx$thdi, spx$lambda,
       dep, asy1, asy2, scale1, shape1, scale2, shape2, dns = double(1),
       PACKAGE = "evd")$dns
  }
  param <- c("scale1", "shape1")
  if(!cscale) param <- c(param, "scale2")
  if(!cshape) param <- c(param, "shape2")
  if(!sym) param <- c(param, "asy1", "asy2", "dep")
  else param <- c(param, "asy1", "dep")
  nmdots <- names(list(...))
  start <- bvstart.vals(x = x, start = start, nmdots = nmdots,
    param = param, method = "pot", u = u, model = "aneglog")
  spx <- sep.bvdata(x = x, method = "cpot", u = u) 
  nm <- names(start)
  l <- length(nm)
  fixed.param <- list(...)[nmdots %in% param]
  if(any(!(param %in% c(nm,names(fixed.param)))))
    stop("unspecified parameters")
  prind <- c(TRUE, TRUE, !cscale, !cshape, TRUE, !sym, TRUE)
  f <- formals(nllbvcaneglog)[prind]
  names(f) <- param
  m <- match(nm, param)
  if(any(is.na(m))) 
    stop("`start' specifies unknown arguments") 
  formals(nllbvcaneglog) <- c(f[m], f[-m])
  nll <- function(p, ...) nllbvcaneglog(p, ...)
  if(l > 1) {
    body(nll) <- parse(text = paste("nllbvcaneglog(", paste("p[", 1:l, "]", collapse=", "), ", ...)"))
  }
  start.arg <- c(list(p = unlist(start)), fixed.param)
  if(warn.inf && do.call("nll", start.arg) == 1e+06)
    warning("negative log-likelihood is infinite at starting values")
  opt <- optim(start, nll, hessian = std.err, ..., method = method)
  cmar <- c(cscale, cshape); nat <- spx$nat
  bvpost.optim(x = x, opt = opt, nm = nm, fixed.param = fixed.param,
    std.err = std.err, corr = corr, sym = sym, cmar = cmar, method = "pot",
    u = u, nat = nat, likelihood = "censored", model = "aneglog")
}

fbvcct <- function(x, u, start, ..., sym = FALSE, cshape = cscale, cscale = FALSE, std.err = TRUE, corr = FALSE, method = "BFGS", warn.inf = TRUE) {

  nllbvcct <- function(scale1, shape1, scale2, shape2, alpha, beta) {
    if(sym) beta <- alpha
    if(cshape) shape2 <- shape1
    if(cscale) scale2 <- scale1
    .C("nllbvcct", spx$x1, spx$x2, spx$nn, spx$n, spx$thdi, spx$lambda,
      alpha, beta, scale1, shape1, scale2, shape2, dns = double(1),
      PACKAGE = "evd")$dns
  }
  param <- c("scale1", "shape1")
  if(!cscale) param <- c(param, "scale2")
  if(!cshape) param <- c(param, "shape2")
  if(!sym) param <- c(param, "alpha", "beta")
  else param <- c(param, "alpha")
  nmdots <- names(list(...))
  start <- bvstart.vals(x = x, start = start, nmdots = nmdots,
    param = param, method = "pot", u = u, model = "ct")
  spx <- sep.bvdata(x = x, method = "cpot", u = u) 
  nm <- names(start)
  l <- length(nm)
  fixed.param <- list(...)[nmdots %in% param]
  if(any(!(param %in% c(nm,names(fixed.param)))))
    stop("unspecified parameters")
  prind <- c(TRUE, TRUE, !cscale, !cshape, TRUE, !sym)
  f <- formals(nllbvcct)[prind]
  names(f) <- param
  m <- match(nm, param)
  if(any(is.na(m))) 
    stop("`start' specifies unknown arguments") 
  formals(nllbvcct) <- c(f[m], f[-m])
  nll <- function(p, ...) nllbvcct(p, ...)
  if(l > 1) {
    body(nll) <- parse(text = paste("nllbvcct(", paste("p[", 1:l, "]", collapse=", "), ", ...)"))
  }
  start.arg <- c(list(p = unlist(start)), fixed.param)
  if(warn.inf && do.call("nll", start.arg) == 1e+06)
    warning("negative log-likelihood is infinite at starting values")
  opt <- optim(start, nll, hessian = std.err, ..., method = method)
  cmar <- c(cscale, cshape); nat <- spx$nat
  bvpost.optim(x = x, opt = opt, nm = nm, fixed.param = fixed.param,
    std.err = std.err, corr = corr, sym = sym, cmar = cmar, method = "pot",
    u = u, nat = nat, likelihood = "censored", model = "ct")
}

fbvchr <- function(x, u, start, ..., sym = FALSE, cshape = cscale, cscale = FALSE, std.err = TRUE, corr = FALSE, method = "BFGS", warn.inf = TRUE) {

  nllbvchr <- function(scale1, shape1, scale2, shape2, dep) {
    if(cshape) shape2 <- shape1
    if(cscale) scale2 <- scale1
    .C("nllbvchr", spx$x1, spx$x2, spx$nn, spx$n, spx$thdi,
      spx$lambda, dep, scale1, shape1, scale2, shape2, dns = double(1),
      PACKAGE = "evd")$dns
  }
  param <- c("scale1", "shape1")
  if(!cscale) param <- c(param, "scale2")
  if(!cshape) param <- c(param, "shape2")
  param <- c(param, "dep")
  nmdots <- names(list(...))
  start <- bvstart.vals(x = x, start = start, nmdots = nmdots,
    param = param, method = "pot", u = u, model = "hr")
  spx <- sep.bvdata(x = x, method = "cpot", u = u)  
  nm <- names(start)
  l <- length(nm)
  fixed.param <- list(...)[nmdots %in% param]
  if(any(!(param %in% c(nm,names(fixed.param)))))
    stop("unspecified parameters")
  prind <- c(TRUE, TRUE, !cscale, !cshape, TRUE)
  f <- formals(nllbvchr)[prind]
  names(f) <- param
  m <- match(nm, param)
  if(any(is.na(m))) 
    stop("`start' specifies unknown arguments")     
  formals(nllbvchr) <- c(f[m], f[-m])
  nll <- function(p, ...) nllbvchr(p, ...)
  if(l > 1) {
    body(nll) <- parse(text = paste("nllbvchr(", paste("p[", 1:l, "]", collapse=", "), ", ...)"))
  }
  start.arg <- c(list(p = unlist(start)), fixed.param)
  if(warn.inf && do.call("nll", start.arg) == 1e+06)
    warning("negative log-likelihood is infinite at starting values")
  opt <- optim(start, nll, hessian = std.err, ..., method = method)
  cmar <- c(cscale, cshape); nat <- spx$nat
  bvpost.optim(x = x, opt = opt, nm = nm, fixed.param = fixed.param,
    std.err = std.err, corr = corr, sym = sym, cmar = cmar, method = "pot",
    u = u, nat = nat, likelihood = "censored", model = "hr")
}

fbvcamix <- function(x, u, start, ..., sym = FALSE, cshape = cscale, cscale = FALSE, std.err = TRUE, corr = FALSE, method = "BFGS", warn.inf = TRUE) {

  nllbvcamix <- function(scale1, shape1, scale2, shape2, alpha, beta) {
    if(cshape) shape2 <- shape1
    if(cscale) scale2 <- scale1
    .C("nllbvcamix", spx$x1, spx$x2, spx$nn, spx$n, spx$thdi, spx$lambda,
      alpha, beta, scale1, shape1, scale2, shape2, dns = double(1),
      PACKAGE = "evd")$dns
  }
  param <- c("scale1", "shape1")
  if(!cscale) param <- c(param, "scale2")
  if(!cshape) param <- c(param, "shape2")
  param <- c(param, "alpha", "beta")
  nmdots <- names(list(...))
  start <- bvstart.vals(x = x, start = start, nmdots = nmdots,
    param = param, method = "pot", u = u, model = "amix")
  spx <- sep.bvdata(x = x, method = "cpot", u = u) 
  nm <- names(start)
  l <- length(nm)
  fixed.param <- list(...)[nmdots %in% param]
  if(any(!(param %in% c(nm,names(fixed.param)))))
    stop("unspecified parameters")
  prind <- c(TRUE, TRUE, !cscale, !cshape, TRUE, TRUE)
  f <- formals(nllbvcamix)[prind]
  names(f) <- param
  m <- match(nm, param)
  if(any(is.na(m))) 
    stop("`start' specifies unknown arguments") 
  formals(nllbvcamix) <- c(f[m], f[-m])
  nll <- function(p, ...) nllbvcamix(p, ...)
  if(l > 1) {
    body(nll) <- parse(text = paste("nllbvcamix(", paste("p[", 1:l, "]", collapse=", "), ", ...)"))
  }
  start.arg <- c(list(p = unlist(start)), fixed.param)
  if(warn.inf && do.call("nll", start.arg) == 1e+06)
    warning("negative log-likelihood is infinite at starting values")
  opt <- optim(start, nll, hessian = std.err, ..., method = method)
  cmar <- c(cscale, cshape); nat <- spx$nat
  bvpost.optim(x = x, opt = opt, nm = nm, fixed.param = fixed.param,
    std.err = std.err, corr = corr, sym = sym, cmar = cmar, method = "pot",
    u = u, nat = nat, likelihood = "censored", model = "amix")
}

### Poisson Likelihood Fitting ###

fbvplog <- function(x, u, start, ..., sym = FALSE, cshape = cscale, cscale = FALSE, std.err = TRUE, corr = FALSE, method = "BFGS", warn.inf = TRUE) {

  nllbvplog <- function(scale1, shape1, scale2, shape2, dep) {
    if(cshape) shape2 <- shape1
    if(cscale) scale2 <- scale1
    .C("nllbvplog", spx$x1, spx$x2, spx$nn, spx$thdi, spx$r1, spx$r2,
      spx$lambda, dep, scale1, shape1, scale2, shape2, dns = double(1),
      PACKAGE = "evd")$dns
  }
  param <- c("scale1", "shape1")
  if(!cscale) param <- c(param, "scale2")
  if(!cshape) param <- c(param, "shape2")
  param <- c(param, "dep")
  nmdots <- names(list(...))
  start <- bvstart.vals(x = x, start = start, nmdots = nmdots,
    param = param, method = "pot", u = u, model = "log")
  spx <- sep.bvdata(x = x, method = "ppot", u = u) 
  nm <- names(start)
  l <- length(nm)
  fixed.param <- list(...)[nmdots %in% param]
  if(any(!(param %in% c(nm,names(fixed.param)))))
    stop("unspecified parameters")
  prind <- c(TRUE, TRUE, !cscale, !cshape, TRUE)
  f <- formals(nllbvplog)[prind]
  names(f) <- param
  m <- match(nm, param)
  if(any(is.na(m))) 
    stop("`start' specifies unknown arguments")   
  formals(nllbvplog) <- c(f[m], f[-m])
  nll <- function(p, ...) nllbvplog(p, ...)
  if(l > 1) {
    body(nll) <- parse(text = paste("nllbvplog(", paste("p[", 1:l, "]", collapse=", "), ", ...)"))
  }
  start.arg <- c(list(p = unlist(start)), fixed.param)
  if(warn.inf && do.call("nll", start.arg) == 1e+06)
    warning("negative log-likelihood is infinite at starting values")
  opt <- optim(start, nll, hessian = std.err, ..., method = method)
  cmar <- c(cscale, cshape); nat <- spx$nat
  bvpost.optim(x = x, opt = opt, nm = nm, fixed.param = fixed.param,
    std.err = std.err, corr = corr, sym = sym, cmar = cmar, method = "pot",
    u = u, nat = nat, likelihood = "poisson", model = "log")
}

fbvpneglog <- function(x, u, start, ..., sym = FALSE, cshape = cscale, cscale = FALSE, std.err = TRUE, corr = FALSE, method = "BFGS", warn.inf = TRUE) {

  nllbvpneglog <- function(scale1, shape1, scale2, shape2, dep) {
    if(cshape) shape2 <- shape1
    if(cscale) scale2 <- scale1
    .C("nllbvpneglog", spx$x1, spx$x2, spx$nn, spx$thdi, spx$r1,
      spx$r2, spx$lambda, dep, scale1, shape1, scale2, shape2,
      dns = double(1), PACKAGE = "evd")$dns    
  }
  param <- c("scale1", "shape1")
  if(!cscale) param <- c(param, "scale2")
  if(!cshape) param <- c(param, "shape2")
  param <- c(param, "dep")
  nmdots <- names(list(...))
  start <- bvstart.vals(x = x, start = start, nmdots = nmdots,
    param = param, method = "pot", u = u, model = "neglog")
  spx <- sep.bvdata(x = x, method = "ppot", u = u)
  nm <- names(start)
  l <- length(nm)
  fixed.param <- list(...)[nmdots %in% param]
  if(any(!(param %in% c(nm,names(fixed.param)))))
    stop("unspecified parameters")
  prind <- c(TRUE, TRUE, !cscale, !cshape, TRUE)
  f <- formals(nllbvpneglog)[prind]
  names(f) <- param
  m <- match(nm, param)
  if(any(is.na(m))) 
    stop("`start' specifies unknown arguments") 
  formals(nllbvpneglog) <- c(f[m], f[-m])
  nll <- function(p, ...) nllbvpneglog(p, ...)
  if(l > 1) {
    body(nll) <- parse(text = paste("nllbvpneglog(", paste("p[", 1:l, "]", collapse=", "), ", ...)"))
  }
  start.arg <- c(list(p = unlist(start)), fixed.param)
  if(warn.inf && do.call("nll", start.arg) == 1e+06)
    warning("negative log-likelihood is infinite at starting values")
  opt <- optim(start, nll, hessian = std.err, ..., method = method)
  cmar <- c(cscale, cshape); nat <- spx$nat
  bvpost.optim(x = x, opt = opt, nm = nm, fixed.param = fixed.param,
    std.err = std.err, corr = corr, sym = sym, cmar = cmar, method = "pot",
    u = u, nat = nat, likelihood = "poisson", model = "neglog")
}

fbvpct <- function(x, u, start, ..., sym = FALSE, cshape = cscale, cscale = FALSE, std.err = TRUE, corr = FALSE, method = "BFGS", warn.inf = TRUE) {

  nllbvpct <- function(scale1, shape1, scale2, shape2, alpha, beta) {
    if(sym) beta <- alpha
    if(cshape) shape2 <- shape1
    if(cscale) scale2 <- scale1
    .C("nllbvpct", spx$x1, spx$x2, spx$nn, spx$thdi, spx$r1, spx$r2,
      spx$lambda, alpha, beta, scale1, shape1, scale2, shape2,
      dns = double(1), PACKAGE = "evd")$dns
  }
  param <- c("scale1", "shape1")
  if(!cscale) param <- c(param, "scale2")
  if(!cshape) param <- c(param, "shape2")
  if(!sym) param <- c(param, "alpha", "beta")
  else param <- c(param, "alpha")
  nmdots <- names(list(...))
  start <- bvstart.vals(x = x, start = start, nmdots = nmdots,
    param = param, method = "pot", u = u, model = "ct")
  spx <- sep.bvdata(x = x, method = "ppot", u = u)
  nm <- names(start)
  l <- length(nm)
  fixed.param <- list(...)[nmdots %in% param]
  if(any(!(param %in% c(nm,names(fixed.param)))))
    stop("unspecified parameters")
  prind <- c(TRUE, TRUE, !cscale, !cshape, TRUE, !sym)
  f <- formals(nllbvpct)[prind]
  names(f) <- param
  m <- match(nm, param)
  if(any(is.na(m))) 
    stop("`start' specifies unknown arguments") 
  formals(nllbvpct) <- c(f[m], f[-m])
  nll <- function(p, ...) nllbvpct(p, ...)
  if(l > 1) {
    body(nll) <- parse(text = paste("nllbvpct(", paste("p[", 1:l, "]", collapse=", "), ", ...)"))
  }
  start.arg <- c(list(p = unlist(start)), fixed.param)
  if(warn.inf && do.call("nll", start.arg) == 1e+06)
    warning("negative log-likelihood is infinite at starting values")
  opt <- optim(start, nll, hessian = std.err, ..., method = method)
  cmar <- c(cscale, cshape); nat <- spx$nat
  bvpost.optim(x = x, opt = opt, nm = nm, fixed.param = fixed.param,
    std.err = std.err, corr = corr, sym = sym, cmar = cmar, method = "pot",
    u = u, nat = nat, likelihood = "poisson", model = "ct")
}

fbvpbilog <- function(x, u, start, ..., sym = FALSE, cshape = cscale, cscale = FALSE, std.err = TRUE, corr = FALSE, method = "BFGS", warn.inf = TRUE) {

  nllbvpbilog <- function(scale1, shape1, scale2, shape2, alpha, beta) {
    if(cshape) shape2 <- shape1
    if(cscale) scale2 <- scale1
    .C("nllbvpbilog", spx$x1, spx$x2, spx$nn, spx$thdi, spx$r1,
      spx$r2, spx$lambda, alpha, beta, scale1, shape1, scale2, shape2,
      dns = double(1), PACKAGE = "evd")$dns    
  }
  param <- c("scale1", "shape1")
  if(!cscale) param <- c(param, "scale2")
  if(!cshape) param <- c(param, "shape2")
  param <- c(param, "alpha", "beta")
  nmdots <- names(list(...))
  start <- bvstart.vals(x = x, start = start, nmdots = nmdots,
    param = param, method = "pot", u = u, model = "bilog")
  spx <- sep.bvdata(x = x, method = "ppot", u = u) 
  nm <- names(start)
  l <- length(nm)
  fixed.param <- list(...)[nmdots %in% param]
  if(any(!(param %in% c(nm,names(fixed.param)))))
    stop("unspecified parameters")
  prind <- c(TRUE, TRUE, !cscale, !cshape, TRUE, TRUE)
  f <- formals(nllbvpbilog)[prind]
  names(f) <- param
  m <- match(nm, param)
  if(any(is.na(m))) 
    stop("`start' specifies unknown arguments") 
  formals(nllbvpbilog) <- c(f[m], f[-m])
  nll <- function(p, ...) nllbvpbilog(p, ...)
  if(l > 1) {
    body(nll) <- parse(text = paste("nllbvpbilog(", paste("p[", 1:l,
      "]", collapse=", "), ", ...)"))
  }
  start.arg <- c(list(p = unlist(start)), fixed.param)
  if(warn.inf && do.call("nll", start.arg) == 1e+06)
    warning("negative log-likelihood is infinite at starting values")
  opt <- optim(start, nll, hessian = std.err, ..., method = method)
  cmar <- c(cscale, cshape); nat <- spx$nat
  bvpost.optim(x = x, opt = opt, nm = nm, fixed.param = fixed.param,
    std.err = std.err, corr = corr, sym = sym, cmar = cmar, method = "pot",
    u = u, nat = nat, likelihood = "poisson", model = "bilog")
}

fbvpnegbilog <- function(x, u, start, ..., sym = FALSE, cshape = cscale, cscale = FALSE, std.err = TRUE, corr = FALSE, method = "BFGS", warn.inf = TRUE) {

  nllbvpnegbilog <- function(scale1, shape1, scale2, shape2, alpha, beta) {
    if(cshape) shape2 <- shape1
    if(cscale) scale2 <- scale1
    .C("nllbvpnegbilog", spx$x1, spx$x2, spx$nn, spx$thdi, spx$r1,
      spx$r2, spx$lambda, alpha, beta, scale1, shape1, scale2, shape2,
      dns = double(1), PACKAGE = "evd")$dns 
  }
  param <- c("scale1", "shape1")
  if(!cscale) param <- c(param, "scale2")
  if(!cshape) param <- c(param, "shape2")
  param <- c(param, "alpha", "beta")
  nmdots <- names(list(...))
  start <- bvstart.vals(x = x, start = start, nmdots = nmdots,
    param = param, method = "pot", u = u, model = "negbilog")
  spx <- sep.bvdata(x = x, method = "ppot", u = u) 
  nm <- names(start)
  l <- length(nm)
  fixed.param <- list(...)[nmdots %in% param]
  if(any(!(param %in% c(nm,names(fixed.param)))))
    stop("unspecified parameters")
  prind <- c(TRUE, TRUE, !cscale, !cshape, TRUE, TRUE)
  f <- formals(nllbvpnegbilog)[prind]
  names(f) <- param
  m <- match(nm, param)
  if(any(is.na(m))) 
    stop("`start' specifies unknown arguments") 
  formals(nllbvpnegbilog) <- c(f[m], f[-m])
  nll <- function(p, ...) nllbvpnegbilog(p, ...)
  if(l > 1) {
    body(nll) <- parse(text = paste("nllbvpnegbilog(", paste("p[", 1:l,
      "]", collapse=", "), ", ...)"))
  }
  start.arg <- c(list(p = unlist(start)), fixed.param)
  if(warn.inf && do.call("nll", start.arg) == 1e+06)
    warning("negative log-likelihood is infinite at starting values")
  opt <- optim(start, nll, hessian = std.err, ..., method = method)
  cmar <- c(cscale, cshape); nat <- spx$nat
  bvpost.optim(x = x, opt = opt, nm = nm, fixed.param = fixed.param,
    std.err = std.err, corr = corr, sym = sym, cmar = cmar, method = "pot",
    u = u, nat = nat, likelihood = "poisson", model = "negbilog")
}

fbvphr <- function(x, u, start, ..., sym = FALSE, cshape = cscale, cscale = FALSE, std.err = TRUE, corr = FALSE, method = "BFGS", warn.inf = TRUE) {

  nllbvphr <- function(scale1, shape1, scale2, shape2, dep) {
    if(cshape) shape2 <- shape1
    if(cscale) scale2 <- scale1
    .C("nllbvphr", spx$x1, spx$x2, spx$nn, spx$thdi, spx$r1, spx$r2,
      spx$lambda, dep, scale1, shape1, scale2, shape2, dns = double(1),
      PACKAGE = "evd")$dns
  }
  param <- c("scale1", "shape1")
  if(!cscale) param <- c(param, "scale2")
  if(!cshape) param <- c(param, "shape2")
  param <- c(param, "dep")
  nmdots <- names(list(...))
  start <- bvstart.vals(x = x, start = start, nmdots = nmdots,
    param = param, method = "pot", u = u, model = "hr")
  spx <- sep.bvdata(x = x, method = "ppot", u = u) 
  nm <- names(start)
  l <- length(nm)
  fixed.param <- list(...)[nmdots %in% param]
  if(any(!(param %in% c(nm,names(fixed.param)))))
    stop("unspecified parameters")
  prind <- c(TRUE, TRUE, !cscale, !cshape, TRUE)
  f <- formals(nllbvphr)[prind]
  names(f) <- param
  m <- match(nm, param)
  if(any(is.na(m))) 
    stop("`start' specifies unknown arguments")
  formals(nllbvphr) <- c(f[m], f[-m])
  nll <- function(p, ...) nllbvphr(p, ...)
  if(l > 1) {
    body(nll) <- parse(text = paste("nllbvphr(", paste("p[", 1:l, "]", collapse=", "), ", ...)"))
  }
  start.arg <- c(list(p = unlist(start)), fixed.param)
  if(warn.inf && do.call("nll", start.arg) == 1e+06)
    warning("negative log-likelihood is infinite at starting values")
  opt <- optim(start, nll, hessian = std.err, ..., method = method)
  cmar <- c(cscale, cshape); nat <- spx$nat
  bvpost.optim(x = x, opt = opt, nm = nm, fixed.param = fixed.param,
    std.err = std.err, corr = corr, sym = sym, cmar = cmar, method = "pot",
    u = u, nat = nat, likelihood = "poisson", model = "hr")
}

### Method Function ###

"print.bvpot" <-  function(x, digits = max(3, getOption("digits") - 3), ...) 
{
    cat("\nCall:", deparse(x$call), "\n")
		cat("Likelihood:", x$likelihood, "\n")
    cat("Deviance:", deviance(x), "\n")
    cat("AIC:", AIC(x), "\n")
    if(!is.null(x$dep.summary)) cat("Dependence:", x$dep.summary, "\n")
    
    cat("\nThreshold:", round(x$threshold, digits), "\n")
    cat("Marginal Number Above:", x$nat[1:2], "\n")
    cat("Marginal Proportion Above:", round(x$nat[1:2]/x$n, digits), "\n")
    cat("Number Above:", x$nat[3], "\n")
    cat("Proportion Above:", round(x$nat[3]/x$n, digits), "\n")
    
    cat("\nEstimates\n")
    print.default(format(fitted(x), digits = digits), print.gap = 2, 
        quote = FALSE)
    if(!is.null(std.errors(x))) {
      cat("\nStandard Errors\n")
      print.default(format(std.errors(x), digits = digits),
          print.gap = 2, quote = FALSE)
    }
    if(!is.null(x$corr)) {
      cat("\nCorrelations\n")
      print.default(format(x$corr, digits = digits), print.gap = 2, 
          quote = FALSE)
    }
    cat("\nOptimization Information\n")
    cat("  Convergence:", x$convergence, "\n")
    cat("  Function Evaluations:", x$counts["function"], "\n")
    if(!is.na(x$counts["gradient"]))
        cat("  Gradient Evaluations:", x$counts["gradient"], "\n")
    if(!is.null(x$message)) cat("  Message:", x$message, "\n")
    cat("\n")
    invisible(x)
}



