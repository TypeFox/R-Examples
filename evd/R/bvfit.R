
"fbvevd" <-
function(x, model = c("log", "alog", "hr", "neglog", "aneglog", "bilog", "negbilog", "ct", "amix"), start, ..., sym = FALSE, nsloc1 = NULL, nsloc2 = NULL, cshape = cscale, cscale = cloc, cloc = FALSE, std.err = TRUE, corr = FALSE, method = "BFGS", warn.inf = TRUE)
{
  call <- match.call()
  model <- match.arg(model)
  if(!is.matrix(x) && !is.data.frame(x))
    stop("`x' must be a matrix or data frame")
  if(ncol(x) != 2) {
    if(ncol(x) != 3) stop("`x' has incorrect number of columns")
    if(!is.logical(x[,3])) stop("third column of `x' must be logical")
  }
  if(sym && !(model %in% c("alog","aneglog","ct")))
    warning("Argument `sym' was ignored")
  
  ft <- switch(model,
    log = fbvlog(x=x, start=start, ..., sym=FALSE, nsloc1=nsloc1,
      nsloc2=nsloc2, cshape=cshape, cscale=cscale, cloc=cloc, std.err=std.err,
      corr=corr, method=method, warn.inf=warn.inf),
    alog = fbvalog(x=x, start=start, ..., sym=sym, nsloc1=nsloc1,
      nsloc2=nsloc2, cshape=cshape, cscale=cscale, cloc=cloc,
      std.err=std.err, corr=corr, method=method, warn.inf=warn.inf),
    hr = fbvhr(x=x, start=start, ..., sym=FALSE, nsloc1=nsloc1, nsloc2=
      nsloc2, cshape=cshape, cscale=cscale, cloc=cloc, std.err=std.err,
      corr=corr, method=method, warn.inf=warn.inf),
    neglog = fbvneglog(x=x, start=start, ..., sym=FALSE, nsloc1=nsloc1,
      nsloc2=nsloc2, cshape=cshape, cscale=cscale, cloc=cloc, std.err=std.err,
      corr=corr, method=method, warn.inf=warn.inf),
    aneglog = fbvaneglog(x=x, start=start, ..., sym=sym, nsloc1=nsloc1,
      nsloc2=nsloc2, cshape=cshape, cscale=cscale, cloc=cloc, std.err=std.err,
      corr=corr, method=method, warn.inf=warn.inf),
    bilog = fbvbilog(x=x, start=start, ..., sym=FALSE, nsloc1=nsloc1,
      nsloc2=nsloc2, cshape=cshape, cscale=cscale, cloc=cloc, std.err=std.err,
      corr=corr, method=method, warn.inf=warn.inf),
    negbilog = fbvnegbilog(x=x, start=start, ..., sym=FALSE, nsloc1=nsloc1,
      nsloc2=nsloc2, cshape=cshape, cscale=cscale, cloc=cloc, std.err=std.err,
      corr= corr, method=method, warn.inf=warn.inf),
    ct = fbvct(x=x, start=start, ..., sym=sym, nsloc1=nsloc1,
      nsloc2=nsloc2, cshape=cshape, cscale=cscale, cloc=cloc,
      std.err=std.err, corr=corr, method=method, warn.inf=warn.inf),
    amix = fbvamix(x=x, start=start, ..., sym=FALSE, nsloc1=nsloc1,
      nsloc2=nsloc2, cshape=cshape, cscale=cscale, cloc=cloc,
      std.err=std.err, corr=corr, method=method, warn.inf=warn.inf))

  structure(c(ft, call = call), class = c("bvevd","evd"))
}

"fbvlog"<- 
function(x, start, ..., sym = FALSE, nsloc1 = NULL, nsloc2 = NULL, cshape = cscale, cscale = cloc, cloc = FALSE, std.err = TRUE, corr = FALSE, method = "BFGS", warn.inf = TRUE)
{
    nlbvlog <- function(loc1, scale1, shape1, loc2, scale2, shape2, dep)
    {
        if(cshape) shape2 <- shape1
        if(cscale) scale2 <- scale1
        
        if(any(c(scale1,scale2) < 0.01) || dep < 0.1 || dep > 1)
            return(1e6)
        if(!is.null(nsloc1)) {
            ns <- numeric(length(loc.param1))
            for(i in 1:length(ns))
                ns[i] <- get(loc.param1[i])
            loc1 <- drop(nslocmat1 %*% ns)
        }
        else loc1 <- rep(loc1, length.out = nrow(x))
        if(cloc) loc2 <- loc1 else {
          if(!is.null(nsloc2)) {
            ns <- numeric(length(loc.param2))
            for(i in 1:length(ns))
                ns[i] <- get(loc.param2[i])
            loc2 <- drop(nslocmat2 %*% ns)
          }
          else loc2 <- rep(loc2, length.out = nrow(x))
        }
        if(spx$n.m1)
            m1l <- .C("nlgev", spx$x.m1, spx$n.m1, loc1[spx$na == 2], scale1,
                      shape1, dns = double(1), PACKAGE = "evd")$dns
        else m1l <- 0
        if(spx$n.m2)
            m2l <- .C("nlgev", spx$x.m2, spx$n.m2, loc2[spx$na == 1], scale2,
                      shape2, dns = double(1), PACKAGE = "evd")$dns
        else m2l <- 0   
        bvl <- .C("nlbvlog", spx$x1, spx$x2, spx$n, spx$si, dep,
                  loc1[spx$na == 0], scale1, shape1, loc2[spx$na == 0],
                  scale2, shape2, cfalse, dns = double(1), PACKAGE = "evd")$dns
        if(any(is.nan(c(m1l,m2l,bvl)))) {
            warning("NaN returned in likelihood")
            return(1e6)
        }
        if(any(c(m1l,m2l,bvl) == 1e6)) return(1e6)
        else return(m1l + m2l + bvl)
    }
    if(cloc && !identical(nsloc1, nsloc2))
      stop("nsloc1 and nsloc2 must be identical")
    if(!is.null(nsloc1)) {
        if(is.vector(nsloc1)) nsloc1 <- data.frame(trend = nsloc1)
        if(nrow(nsloc1) != nrow(x))
          stop("`nsloc1' and data are not compatible")
        nslocmat1 <- cbind(1, as.matrix(nsloc1))
    }
    if(!is.null(nsloc2)) {
        if(is.vector(nsloc2)) nsloc2 <- data.frame(trend = nsloc2)
        if(nrow(nsloc2) != nrow(x))
          stop("`nsloc2' and data are not compatible")
        nslocmat2 <- cbind(1,as.matrix(nsloc2))
    }
    loc.param1 <- paste("loc1", c("",names(nsloc1)), sep="")
    loc.param2 <- paste("loc2", c("",names(nsloc2)), sep="")
    param <- c(loc.param1, "scale1", "shape1")
    if(!cloc) param <- c(param, loc.param2) else loc.param2 <- NULL
    if(!cscale) param <- c(param, "scale2")
    if(!cshape) param <- c(param, "shape2")
    param <- c(param, "dep")
    nmdots <- names(list(...))
    start <- bvstart.vals(x = x, start = start, nmdots = nmdots, param =
      param, nsloc1 = nsloc1, nsloc2 = nsloc2, model = "log")
    spx <- sep.bvdata(x = x)
    cfalse <- as.integer(0)
    nm <- names(start)
    l <- length(nm)
    fixed.param <- list(...)[nmdots %in% param]
    if(any(!(param %in% c(nm,names(fixed.param)))))
        stop("unspecified parameters")
    prind <- (5:7)[c(!cscale, !cshape, TRUE)]
    f <- c(as.list(numeric(length(loc.param1))), formals(nlbvlog)[2:3],
           as.list(numeric(length(loc.param2))), formals(nlbvlog)[prind])
    names(f) <- param
    m <- match(nm, param)
    if(any(is.na(m))) 
        stop("`start' specifies unknown arguments")    
    formals(nlbvlog) <- c(f[m], f[-m])
    nllh <- function(p, ...) nlbvlog(p, ...)
    if(l > 1)
        body(nllh) <- parse(text = paste("nlbvlog(", paste("p[",1:l,
            "]", collapse = ", "), ", ...)"))
    start.arg <- c(list(p = unlist(start)), fixed.param)
    if(warn.inf && do.call("nllh", start.arg) == 1e6)
        warning("negative log-likelihood is infinite at starting values")
    opt <- optim(start, nllh, hessian = TRUE, ..., method = method)
	if(is.null(names(opt$par))) names(opt$par) <- nm
    cmar <- c(cloc, cscale, cshape)
    bvpost.optim(x = x, opt = opt, nm = nm, fixed.param = fixed.param,
      std.err = std.err, corr = corr, sym = sym, cmar = cmar,
      nsloc1 = nsloc1, nsloc2 = nsloc2, model = "log")
}

"fbvalog"<- 
function(x, start, ..., sym = FALSE, nsloc1 = NULL, nsloc2 = NULL, cshape = cscale, cscale = cloc, cloc = FALSE, std.err = TRUE, corr = FALSE, method = "BFGS", warn.inf = TRUE)
{   
    nlbvalog <- function(loc1, scale1, shape1, loc2, scale2, shape2,
                          asy1, asy2, dep)
    {
        if(sym) asy2 <- asy1
        if(cshape) shape2 <- shape1
        if(cscale) scale2 <- scale1
          
        if(any(c(scale1,scale2) < 0.01) || any(c(dep,asy1,asy2) > 1) ||
           any(c(asy1,asy2) < 0.001) || dep < 0.1)
            return(1e6)
        if(!is.null(nsloc1)) {
            ns <- numeric(length(loc.param1))
            for(i in 1:length(ns))
                ns[i] <- get(loc.param1[i])
            loc1 <- drop(nslocmat1 %*% ns)
        }
        else loc1 <- rep(loc1, length.out = nrow(x))
        if(cloc) loc2 <- loc1 else {
          if(!is.null(nsloc2)) {
            ns <- numeric(length(loc.param2))
            for(i in 1:length(ns))
              ns[i] <- get(loc.param2[i])
            loc2 <- drop(nslocmat2 %*% ns)
          }
          else loc2 <- rep(loc2, length.out = nrow(x))
        }
        if(spx$n.m1)
            m1l <- .C("nlgev", spx$x.m1, spx$n.m1, loc1[spx$na == 2], scale1,
                shape1, dns = double(1), PACKAGE = "evd")$dns
        else m1l <- 0
        if(spx$n.m2)
            m2l <- .C("nlgev", spx$x.m2, spx$n.m2, loc2[spx$na == 1], scale2,
                shape2, dns = double(1), PACKAGE = "evd")$dns
        else m2l <- 0
        bvl <- .C("nlbvalog", spx$x1, spx$x2, spx$n, spx$si, dep, asy1, asy2,
                loc1[spx$na == 0], scale1, shape1, loc2[spx$na == 0],
                scale2, shape2, cfalse, dns = double(1), PACKAGE = "evd")$dns
        if(any(is.nan(c(m1l,m2l,bvl)))) {
            warning("NaN returned in likelihood")
            return(1e6)
        }
        if(any(c(m1l,m2l,bvl) == 1e6)) return(1e6)
        else return(m1l + m2l + bvl)  
    }
    if(cloc && !identical(nsloc1, nsloc2))
      stop("nsloc1 and nsloc2 must be identical")
    if(!is.null(nsloc1)) {
        if(is.vector(nsloc1)) nsloc1 <- data.frame(trend = nsloc1)
        if(nrow(nsloc1) != nrow(x))
          stop("`nsloc1' and data are not compatible")
        nslocmat1 <- cbind(1,as.matrix(nsloc1))
    }
    if(!is.null(nsloc2)) {
        if(is.vector(nsloc2)) nsloc2 <- data.frame(trend = nsloc2)
        if(nrow(nsloc2) != nrow(x))
          stop("`nsloc2' and data are not compatible")
        nslocmat2 <- cbind(1,as.matrix(nsloc2))
    }
    loc.param1 <- paste("loc1", c("",names(nsloc1)), sep="")
    loc.param2 <- paste("loc2", c("",names(nsloc2)), sep="")
    param <- c(loc.param1, "scale1", "shape1")
    if(!cloc) param <- c(param, loc.param2) else loc.param2 <- NULL
    if(!cscale) param <- c(param, "scale2")
    if(!cshape) param <- c(param, "shape2")
    if(!sym) param <- c(param, "asy1", "asy2", "dep")
    else param <- c(param, "asy1", "dep")
    nmdots <- names(list(...))
    start <- bvstart.vals(x = x, start = start, nmdots = nmdots, param =
      param, nsloc1 = nsloc1, nsloc2 = nsloc2, model = "alog")
    spx <- sep.bvdata(x = x)
    cfalse <- as.integer(0)
    nm <- names(start)
    l <- length(nm)
    fixed.param <- list(...)[nmdots %in% param]
    if(any(!(param %in% c(nm,names(fixed.param)))))
        stop("unspecified parameters")
    prind <- (5:9)[c(!cscale, !cshape, TRUE, !sym, TRUE)]
    f <- c(as.list(numeric(length(loc.param1))), formals(nlbvalog)[2:3],
      as.list(numeric(length(loc.param2))), formals(nlbvalog)[prind])
    names(f) <- param
    m <- match(nm, param)
    if(any(is.na(m))) 
        stop("`start' specifies unknown arguments")
    formals(nlbvalog) <- c(f[m], f[-m])
    nllh <- function(p, ...) nlbvalog(p, ...)
    if(l > 1) body(nllh) <- parse(text = paste("nlbvalog(", paste("p[",1:l,"]",
      collapse = ", "), ", ...)"))
    start.arg <- c(list(p = unlist(start)), fixed.param)
    if(warn.inf && do.call("nllh", start.arg) == 1e6)
        warning("negative log-likelihood is infinite at starting values")
    opt <- optim(start, nllh, hessian = TRUE, ..., method = method)
	if(is.null(names(opt$par))) names(opt$par) <- nm
    cmar <- c(cloc, cscale, cshape)
    bvpost.optim(x = x, opt = opt, nm = nm, fixed.param = fixed.param,
      std.err = std.err, corr = corr, sym = sym, cmar = cmar,
      nsloc1 = nsloc1, nsloc2 = nsloc2, model = "alog")
}

"fbvhr"<- 
function(x, start, ..., sym = FALSE, nsloc1 = NULL, nsloc2 = NULL, cshape = cscale, cscale = cloc, cloc = FALSE, std.err = TRUE, corr = FALSE, method = "BFGS", warn.inf = TRUE)
{
    nlbvhr <- function(loc1, scale1, shape1, loc2, scale2, shape2,
                       dep)
    {
        if(cshape) shape2 <- shape1
        if(cscale) scale2 <- scale1
        
        if(any(c(scale1,scale2) < 0.01) || dep < 0.2 || dep > 10)
            return(1e6)
        if(!is.null(nsloc1)) {
            ns <- numeric(length(loc.param1))
            for(i in 1:length(ns))
                ns[i] <- get(loc.param1[i])
            loc1 <- drop(nslocmat1 %*% ns)
        }
        else loc1 <- rep(loc1, length.out = nrow(x))
        if(cloc) loc2 <- loc1 else {
          if(!is.null(nsloc2)) {
            ns <- numeric(length(loc.param2))
            for(i in 1:length(ns))
                ns[i] <- get(loc.param2[i])
            loc2 <- drop(nslocmat2 %*% ns)
          }
          else loc2 <- rep(loc2, length.out = nrow(x))
        }
        if(spx$n.m1)
            m1l <- .C("nlgev", spx$x.m1, spx$n.m1, loc1[spx$na == 2], scale1,
                shape1, dns = double(1), PACKAGE = "evd")$dns
        else m1l <- 0
        if(spx$n.m2)
            m2l <- .C("nlgev", spx$x.m2, spx$n.m2, loc2[spx$na == 1], scale2,
                shape2, dns = double(1), PACKAGE = "evd")$dns
        else m2l <- 0
        bvl <- .C("nlbvhr", spx$x1, spx$x2, spx$n, spx$si, dep,
                loc1[spx$na == 0], scale1, shape1, loc2[spx$na == 0], scale2,
                shape2, cfalse, dns = double(1), PACKAGE = "evd")$dns
        if(any(is.nan(c(m1l,m2l,bvl)))) {
            warning("NaN returned in likelihood")
            return(1e6)
        }
        if(any(c(m1l,m2l,bvl) == 1e6)) return(1e6)
        else return(m1l + m2l + bvl)
    }
    if(cloc && !identical(nsloc1, nsloc2))
      stop("nsloc1 and nsloc2 must be identical")
    if(!is.null(nsloc1)) {
        if(is.vector(nsloc1)) nsloc1 <- data.frame(trend = nsloc1)
        if(nrow(nsloc1) != nrow(x))
          stop("`nsloc1' and data are not compatible") 
        nslocmat1 <- cbind(1,as.matrix(nsloc1))
    }
    if(!is.null(nsloc2)) {
        if(is.vector(nsloc2)) nsloc2 <- data.frame(trend = nsloc2)
        if(nrow(nsloc2) != nrow(x))
          stop("`nsloc2' and data are not compatible")
        nslocmat2 <- cbind(1,as.matrix(nsloc2))
    }
    loc.param1 <- paste("loc1", c("",names(nsloc1)), sep="")
    loc.param2 <- paste("loc2", c("",names(nsloc2)), sep="")
    param <- c(loc.param1, "scale1", "shape1")
    if(!cloc) param <- c(param, loc.param2) else loc.param2 <- NULL
    if(!cscale) param <- c(param, "scale2")
    if(!cshape) param <- c(param, "shape2")
    param <- c(param, "dep")
    nmdots <- names(list(...))
    start <- bvstart.vals(x = x, start = start, nmdots = nmdots, param =
      param, nsloc1 = nsloc1, nsloc2 = nsloc2, model = "hr")
    spx <- sep.bvdata(x = x)
    cfalse <- as.integer(0)
    nm <- names(start)
    l <- length(nm)
    fixed.param <- list(...)[nmdots %in% param]
    if(any(!(param %in% c(nm,names(fixed.param)))))
        stop("unspecified parameters")
    prind <- (5:7)[c(!cscale, !cshape, TRUE)]
    f <- c(as.list(numeric(length(loc.param1))), formals(nlbvhr)[2:3],
           as.list(numeric(length(loc.param2))), formals(nlbvhr)[prind])
    names(f) <- param
    m <- match(nm, param)
    if(any(is.na(m))) 
        stop("`start' specifies unknown arguments")    
    formals(nlbvhr) <- c(f[m], f[-m])
    nllh <- function(p, ...) nlbvhr(p, ...)
    if(l > 1)
        body(nllh) <- parse(text = paste("nlbvhr(", paste("p[",1:l,
            "]", collapse = ", "), ", ...)"))
    start.arg <- c(list(p = unlist(start)), fixed.param)
    if(warn.inf && do.call("nllh", start.arg) == 1e6)
        warning("negative log-likelihood is infinite at starting values")
    opt <- optim(start, nllh, hessian = TRUE, ..., method = method)
	if(is.null(names(opt$par))) names(opt$par) <- nm
    cmar <- c(cloc, cscale, cshape)
    bvpost.optim(x = x, opt = opt, nm = nm, fixed.param = fixed.param,
      std.err = std.err, corr = corr, sym = sym, cmar = cmar,
      nsloc1 = nsloc1, nsloc2 = nsloc2, model = "hr")
}

"fbvneglog"<- 
function(x, start, ..., sym = FALSE, nsloc1 = NULL, nsloc2 = NULL, cshape = cscale, cscale = cloc, cloc = FALSE, std.err = TRUE, corr = FALSE, method = "BFGS", warn.inf = TRUE)
{
    nlbvneglog <- function(loc1, scale1, shape1, loc2, scale2, shape2,
                           dep)
    {
        if(cshape) shape2 <- shape1
        if(cscale) scale2 <- scale1
        
        if(any(c(scale1,scale2) < 0.01) || dep < 0.05 || dep > 5)
            return(1e6)
        if(!is.null(nsloc1)) {
            ns <- numeric(length(loc.param1))
            for(i in 1:length(ns))
                ns[i] <- get(loc.param1[i])
            loc1 <- drop(nslocmat1 %*% ns)
        }
        else loc1 <- rep(loc1, length.out = nrow(x))
        if(cloc) loc2 <- loc1 else {
          if(!is.null(nsloc2)) {
            ns <- numeric(length(loc.param2))
            for(i in 1:length(ns))
                ns[i] <- get(loc.param2[i])
            loc2 <- drop(nslocmat2 %*% ns)
          }
          else loc2 <- rep(loc2, length.out = nrow(x))
        }
        if(spx$n.m1)
            m1l <- .C("nlgev", spx$x.m1, spx$n.m1, loc1[spx$na == 2], scale1,
                shape1, dns = double(1), PACKAGE = "evd")$dns
        else m1l <- 0
        if(spx$n.m2)
            m2l <- .C("nlgev", spx$x.m2, spx$n.m2, loc2[spx$na == 1], scale2,
                shape2, dns = double(1), PACKAGE = "evd")$dns
        else m2l <- 0
        bvl <- .C("nlbvneglog", spx$x1, spx$x2, spx$n, spx$si, dep,
                  loc1[spx$na == 0], scale1, shape1, loc2[spx$na == 0],
                  scale2, shape2, cfalse, dns = double(1), PACKAGE = "evd")$dns
        if(any(is.nan(c(m1l,m2l,bvl)))) {
            warning("NaN returned in likelihood")
            return(1e6)
        }
        if(any(c(m1l,m2l,bvl) == 1e6)) return(1e6)
        else return(m1l + m2l + bvl)
    }
    if(cloc && !identical(nsloc1, nsloc2))
      stop("nsloc1 and nsloc2 must be identical")
    if(!is.null(nsloc1)) {
        if(is.vector(nsloc1)) nsloc1 <- data.frame(trend = nsloc1)
        if(nrow(nsloc1) != nrow(x))
          stop("`nsloc1' and data are not compatible")
        nslocmat1 <- cbind(1,as.matrix(nsloc1))
    }
    if(!is.null(nsloc2)) {
        if(is.vector(nsloc2)) nsloc2 <- data.frame(trend = nsloc2)
        if(nrow(nsloc2) != nrow(x))
          stop("`nsloc2' and data are not compatible")
        nslocmat2 <- cbind(1,as.matrix(nsloc2))
    }
    loc.param1 <- paste("loc1", c("",names(nsloc1)), sep="")
    loc.param2 <- paste("loc2", c("",names(nsloc2)), sep="")
    param <- c(loc.param1, "scale1", "shape1")
    if(!cloc) param <- c(param, loc.param2) else loc.param2 <- NULL
    if(!cscale) param <- c(param, "scale2")
    if(!cshape) param <- c(param, "shape2")
    param <- c(param, "dep")
    nmdots <- names(list(...))
    start <- bvstart.vals(x = x, start = start, nmdots = nmdots, param =
      param, nsloc1 = nsloc1, nsloc2 = nsloc2, model = "neglog")
    spx <- sep.bvdata(x = x)
    cfalse <- as.integer(0)
    nm <- names(start)
    l <- length(nm)
    fixed.param <- list(...)[nmdots %in% param]
    if(any(!(param %in% c(nm,names(fixed.param)))))
        stop("unspecified parameters")
    prind <- (5:7)[c(!cscale, !cshape, TRUE)]
    f <- c(as.list(numeric(length(loc.param1))), formals(nlbvneglog)[2:3],
           as.list(numeric(length(loc.param2))), formals(nlbvneglog)[prind])
    names(f) <- param   
    m <- match(nm, param)
    if(any(is.na(m))) 
        stop("`start' specifies unknown arguments")    
    formals(nlbvneglog) <- c(f[m], f[-m])
    nllh <- function(p, ...) nlbvneglog(p, ...)
    if(l > 1)
        body(nllh) <- parse(text = paste("nlbvneglog(", paste("p[",1:l,
            "]", collapse = ", "), ", ...)"))
    start.arg <- c(list(p = unlist(start)), fixed.param)
    if(warn.inf && do.call("nllh", start.arg) == 1e6)
        warning("negative log-likelihood is infinite at starting values")
    opt <- optim(start, nllh, hessian = TRUE, ..., method = method)
	if(is.null(names(opt$par))) names(opt$par) <- nm
    cmar <- c(cloc, cscale, cshape)
    bvpost.optim(x = x, opt = opt, nm = nm, fixed.param = fixed.param,
      std.err = std.err, corr = corr, sym = sym, cmar = cmar,
      nsloc1 = nsloc1, nsloc2 = nsloc2, model = "neglog")
}

"fbvaneglog"<- 
function(x, start, ..., sym = FALSE, nsloc1 = NULL, nsloc2 = NULL, cshape = cscale, cscale = cloc, cloc = FALSE, std.err = TRUE, corr = FALSE, method = "BFGS", warn.inf = TRUE)
{
    nlbvaneglog <- function(loc1, scale1, shape1, loc2, scale2, shape2,
                            asy1, asy2, dep)
    {
        if(sym) asy2 <- asy1
        if(cshape) shape2 <- shape1
        if(cscale) scale2 <- scale1
        
        if(any(c(scale1,scale2) < 0.01) || any(c(asy1,asy2) > 1) ||
           any(c(asy1,asy2) < 0.001) || dep < 0.05 || dep > 5)
            return(1e6)
        if(!is.null(nsloc1)) {
            ns <- numeric(length(loc.param1))
            for(i in 1:length(ns))
                ns[i] <- get(loc.param1[i])
            loc1 <- drop(nslocmat1 %*% ns)
        }
        else loc1 <- rep(loc1, length.out = nrow(x))
        if(cloc) loc2 <- loc1 else {
          if(!is.null(nsloc2)) {
            ns <- numeric(length(loc.param2))
            for(i in 1:length(ns))
              ns[i] <- get(loc.param2[i])
            loc2 <- drop(nslocmat2 %*% ns)
          }
          else loc2 <- rep(loc2, length.out = nrow(x))
        }
        if(spx$n.m1)
            m1l <- .C("nlgev", spx$x.m1, spx$n.m1, loc1[spx$na == 2],
                scale1, shape1, dns = double(1), PACKAGE = "evd")$dns
        else m1l <- 0
        if(spx$n.m2)
            m2l <- .C("nlgev", spx$x.m2, spx$n.m2, loc2[spx$na == 1],
                scale2, shape2, dns = double(1), PACKAGE = "evd")$dns
        else m2l <- 0
        bvl <- .C("nlbvaneglog", spx$x1, spx$x2, spx$n, spx$si, dep, asy1,
                asy2, loc1[spx$na == 0], scale1, shape1, loc2[spx$na == 0],
                scale2, shape2, cfalse, dns = double(1), PACKAGE = "evd")$dns
        if(any(is.nan(c(m1l,m2l,bvl)))) {
            warning("NaN returned in likelihood")
            return(1e6)
        }
        if(any(c(m1l,m2l,bvl) == 1e6)) return(1e6)
        else return(m1l + m2l + bvl)
    }
    if(cloc && !identical(nsloc1, nsloc2))
      stop("nsloc1 and nsloc2 must be identical")
    if(!is.null(nsloc1)) {
        if(is.vector(nsloc1)) nsloc1 <- data.frame(trend = nsloc1)
        if(nrow(nsloc1) != nrow(x))
          stop("`nsloc1' and data are not compatible")
        nslocmat1 <- cbind(1,as.matrix(nsloc1))
    }
    if(!is.null(nsloc2)) {
        if(is.vector(nsloc2)) nsloc2 <- data.frame(trend = nsloc2)
        if(nrow(nsloc2) != nrow(x))
          stop("`nsloc2' and data are not compatible")
        nslocmat2 <- cbind(1,as.matrix(nsloc2))
    }
    loc.param1 <- paste("loc1", c("",names(nsloc1)), sep="")
    loc.param2 <- paste("loc2", c("",names(nsloc2)), sep="")
    param <- c(loc.param1, "scale1", "shape1")
    if(!cloc) param <- c(param, loc.param2) else loc.param2 <- NULL
    if(!cscale) param <- c(param, "scale2")
    if(!cshape) param <- c(param, "shape2")
    if(!sym) param <- c(param, "asy1", "asy2", "dep")
    else param <- c(param, "asy1", "dep")
    nmdots <- names(list(...))
    start <- bvstart.vals(x = x, start = start, nmdots = nmdots, param =
      param, nsloc1 = nsloc1, nsloc2 = nsloc2, model = "aneglog")
    spx <- sep.bvdata(x = x)
    cfalse <- as.integer(0)
    nm <- names(start)
    l <- length(nm)
    fixed.param <- list(...)[nmdots %in% param]
    if(any(!(param %in% c(nm,names(fixed.param)))))
        stop("unspecified parameters")
    prind <- (5:9)[c(!cscale, !cshape, TRUE, !sym, TRUE)]
    f <- c(as.list(numeric(length(loc.param1))), formals(nlbvaneglog)[2:3],
      as.list(numeric(length(loc.param2))), formals(nlbvaneglog)[prind])
    names(f) <- param
    m <- match(nm, param)
    if(any(is.na(m))) 
        stop("`start' specifies unknown arguments")
    formals(nlbvaneglog) <- c(f[m], f[-m])
    nllh <- function(p, ...) nlbvaneglog(p, ...)
    if(l > 1) body(nllh) <- parse(text = paste("nlbvaneglog(",
      paste("p[",1:l,"]", collapse = ", "), ", ...)"))
    start.arg <- c(list(p = unlist(start)), fixed.param)
    if(warn.inf && do.call("nllh", start.arg) == 1e6)
        warning("negative log-likelihood is infinite at starting values")
    opt <- optim(start, nllh, hessian = TRUE, ..., method = method)
	if(is.null(names(opt$par))) names(opt$par) <- nm
    cmar <- c(cloc, cscale, cshape)
    bvpost.optim(x = x, opt = opt, nm = nm, fixed.param = fixed.param,
      std.err = std.err, corr = corr, sym = sym, cmar = cmar,
      nsloc1 = nsloc1, nsloc2 = nsloc2, model = "aneglog")
}

"fbvbilog"<- 
function(x, start, ..., sym = FALSE, nsloc1 = NULL, nsloc2 = NULL, cshape = cscale, cscale = cloc, cloc = FALSE, std.err = TRUE, corr = FALSE, method = "BFGS", warn.inf = TRUE)
{
    nlbvbilog <- function(loc1, scale1, shape1, loc2, scale2, shape2,
                          alpha, beta)
    {
        if(cshape) shape2 <- shape1
        if(cscale) scale2 <- scale1
        
        if(any(c(scale1,scale2) < 0.01) || any(c(alpha,beta) < 0.1) ||
           any(c(alpha,beta) > 0.999))
            return(1e6)
        if(!is.null(nsloc1)) {
            ns <- numeric(length(loc.param1))
            for(i in 1:length(ns))
                ns[i] <- get(loc.param1[i])
            loc1 <- drop(nslocmat1 %*% ns)
        }
        else loc1 <- rep(loc1, length.out = nrow(x))
        if(cloc) loc2 <- loc1 else {
          if(!is.null(nsloc2)) {
            ns <- numeric(length(loc.param2))
            for(i in 1:length(ns))
                ns[i] <- get(loc.param2[i])
            loc2 <- drop(nslocmat2 %*% ns)
          }
          else loc2 <- rep(loc2, length.out = nrow(x))
        }
        if(spx$n.m1)
            m1l <- .C("nlgev", spx$x.m1,  spx$n.m1, loc1[spx$na == 2],
                scale1, shape1, dns = double(1), PACKAGE = "evd")$dns
        else m1l <- 0
        if(spx$n.m2)
            m2l <- .C("nlgev", spx$x.m2,  spx$n.m2, loc2[spx$na == 1],
                scale2, shape2, dns = double(1), PACKAGE = "evd")$dns
        else m2l <- 0
        bvl <- .C("nlbvbilog", spx$x1, spx$x2, spx$n, spx$si, alpha, beta,
                loc1[spx$na == 0], scale1, shape1, loc2[spx$na == 0],
                scale2, shape2, cfalse, dns = double(1), PACKAGE = "evd")$dns
        if(any(is.nan(c(m1l,m2l,bvl)))) {
            warning("NaN returned in likelihood")
            return(1e6)
        }
        if(any(c(m1l,m2l,bvl) == 1e6)) return(1e6)
        else return(m1l + m2l + bvl)
    }
    if(cloc && !identical(nsloc1, nsloc2))
      stop("nsloc1 and nsloc2 must be identical")
    if(!is.null(nsloc1)) {
        if(is.vector(nsloc1)) nsloc1 <- data.frame(trend = nsloc1)
        if(nrow(nsloc1) != nrow(x))
          stop("`nsloc1' and data are not compatible")
        nslocmat1 <- cbind(1,as.matrix(nsloc1))
    }
    if(!is.null(nsloc2)) {
        if(is.vector(nsloc2)) nsloc2 <- data.frame(trend = nsloc2)
        if(nrow(nsloc2) != nrow(x))
          stop("`nsloc2' and data are not compatible")
        nslocmat2 <- cbind(1,as.matrix(nsloc2))
    }
    loc.param1 <- paste("loc1", c("",names(nsloc1)), sep="")
    loc.param2 <- paste("loc2", c("",names(nsloc2)), sep="")
    param <- c(loc.param1, "scale1", "shape1")
    if(!cloc) param <- c(param, loc.param2) else loc.param2 <- NULL
    if(!cscale) param <- c(param, "scale2")
    if(!cshape) param <- c(param, "shape2")
    param <- c(param, "alpha", "beta")
    nmdots <- names(list(...))
    start <- bvstart.vals(x = x, start = start, nmdots = nmdots, param =
      param, nsloc1 = nsloc1, nsloc2 = nsloc2, model = "bilog")
    spx <- sep.bvdata(x = x)
    cfalse <- as.integer(0)
    nm <- names(start)
    l <- length(nm)
    fixed.param <- list(...)[nmdots %in% param]
    if(any(!(param %in% c(nm,names(fixed.param)))))
        stop("unspecified parameters")
    prind <- (5:8)[c(!cscale, !cshape, TRUE, TRUE)]
    f <- c(as.list(numeric(length(loc.param1))), formals(nlbvbilog)[2:3],
           as.list(numeric(length(loc.param2))), formals(nlbvbilog)[prind])
    names(f) <- param
    m <- match(nm, param)
    if(any(is.na(m))) 
        stop("`start' specifies unknown arguments")    
    formals(nlbvbilog) <- c(f[m], f[-m])
    nllh <- function(p, ...) nlbvbilog(p, ...)
    if(l > 1)
        body(nllh) <- parse(text = paste("nlbvbilog(", paste("p[",1:l,
            "]", collapse = ", "), ", ...)"))
    start.arg <- c(list(p = unlist(start)), fixed.param)
    if(warn.inf && do.call("nllh", start.arg) == 1e6)
        warning("negative log-likelihood is infinite at starting values")
    opt <- optim(start, nllh, hessian = TRUE, ..., method = method)
	if(is.null(names(opt$par))) names(opt$par) <- nm
    cmar <- c(cloc, cscale, cshape)
    bvpost.optim(x = x, opt = opt, nm = nm, fixed.param = fixed.param,
      std.err = std.err, corr = corr, sym = sym, cmar = cmar,
      nsloc1 = nsloc1, nsloc2 = nsloc2, model = "bilog")
}

"fbvnegbilog"<- 
function(x, start, ..., sym = FALSE, nsloc1 = NULL, nsloc2 = NULL, cshape = cscale, cscale = cloc, cloc = FALSE, std.err = TRUE, corr = FALSE, method = "BFGS", warn.inf = TRUE)
{
    nlbvnegbilog <- function(loc1, scale1, shape1, loc2, scale2, shape2,
                             alpha, beta)
    {
        if(cshape) shape2 <- shape1
        if(cscale) scale2 <- scale1
        
        if(any(c(scale1,scale2) < 0.01) || any(c(alpha,beta) < 0.1) ||
           any(c(alpha,beta) > 20))
            return(1e6)
        if(!is.null(nsloc1)) {
            ns <- numeric(length(loc.param1))
            for(i in 1:length(ns))
                ns[i] <- get(loc.param1[i])
            loc1 <- drop(nslocmat1 %*% ns)
        }
        else loc1 <- rep(loc1, length.out = nrow(x))
        if(cloc) loc2 <- loc1 else {
          if(!is.null(nsloc2)) {
            ns <- numeric(length(loc.param2))
            for(i in 1:length(ns))
                ns[i] <- get(loc.param2[i])
            loc2 <- drop(nslocmat2 %*% ns)
          }
          else loc2 <- rep(loc2, length.out = nrow(x))
        }
        if(spx$n.m1)
            m1l <- .C("nlgev", spx$x.m1, spx$n.m1, loc1[spx$na == 2],
                scale1, shape1, dns = double(1), PACKAGE = "evd")$dns
        else m1l <- 0
        if(spx$n.m2)
            m2l <- .C("nlgev", spx$x.m2, spx$n.m2, loc2[spx$na == 1],
                scale2, shape2, dns = double(1), PACKAGE = "evd")$dns
        else m2l <- 0
        bvl <- .C("nlbvnegbilog", spx$x1, spx$x2, spx$n, spx$si, alpha, beta,
                  loc1[spx$na == 0], scale1, shape1, loc2[spx$na == 0],
                  scale2, shape2, cfalse, dns = double(1), PACKAGE = "evd")$dns
        if(any(is.nan(c(m1l,m2l,bvl)))) {
            warning("NaN returned in likelihood")
            return(1e6)
        }
        if(any(c(m1l,m2l,bvl) == 1e6)) return(1e6)
        else return(m1l + m2l + bvl)
    }
    if(cloc && !identical(nsloc1, nsloc2))
      stop("nsloc1 and nsloc2 must be identical")
    if(!is.null(nsloc1)) {
        if(is.vector(nsloc1)) nsloc1 <- data.frame(trend = nsloc1)
        if(nrow(nsloc1) != nrow(x))
          stop("`nsloc1' and data are not compatible")
        nslocmat1 <- cbind(1,as.matrix(nsloc1))
    }
    if(!is.null(nsloc2)) {
        if(is.vector(nsloc2)) nsloc2 <- data.frame(trend = nsloc2)
        if(nrow(nsloc2) != nrow(x))
          stop("`nsloc2' and data are not compatible")
        nslocmat2 <- cbind(1,as.matrix(nsloc2))
    }
    loc.param1 <- paste("loc1", c("",names(nsloc1)), sep="")
    loc.param2 <- paste("loc2", c("",names(nsloc2)), sep="")
    param <- c(loc.param1, "scale1", "shape1")
    if(!cloc) param <- c(param, loc.param2) else loc.param2 <- NULL
    if(!cscale) param <- c(param, "scale2")
    if(!cshape) param <- c(param, "shape2")
    param <- c(param, "alpha", "beta")
    nmdots <- names(list(...))
    start <- bvstart.vals(x = x, start = start, nmdots = nmdots, param =
      param, nsloc1 = nsloc1, nsloc2 = nsloc2, model = "negbilog") 
    spx <- sep.bvdata(x = x)
    cfalse <- as.integer(0)
    nm <- names(start)
    l <- length(nm)
    fixed.param <- list(...)[nmdots %in% param]
    if(any(!(param %in% c(nm,names(fixed.param)))))
        stop("unspecified parameters")
    prind <- (5:8)[c(!cscale, !cshape, TRUE, TRUE)]
    f <- c(as.list(numeric(length(loc.param1))), formals(nlbvnegbilog)[2:3],
           as.list(numeric(length(loc.param2))), formals(nlbvnegbilog)[prind])
    names(f) <- param
    m <- match(nm, param)
    if(any(is.na(m))) 
        stop("`start' specifies unknown arguments")    
    formals(nlbvnegbilog) <- c(f[m], f[-m])
    nllh <- function(p, ...) nlbvnegbilog(p, ...)
    if(l > 1)
        body(nllh) <- parse(text = paste("nlbvnegbilog(", paste("p[",1:l,
            "]", collapse = ", "), ", ...)"))
    start.arg <- c(list(p = unlist(start)), fixed.param)
    if(warn.inf && do.call("nllh", start.arg) == 1e6)
        warning("negative log-likelihood is infinite at starting values") 
    opt <- optim(start, nllh, hessian = TRUE, ..., method = method)
	if(is.null(names(opt$par))) names(opt$par) <- nm
    cmar <- c(cloc, cscale, cshape)
    bvpost.optim(x = x, opt = opt, nm = nm, fixed.param = fixed.param,
      std.err = std.err, corr = corr, sym = sym, cmar = cmar,
      nsloc1 = nsloc1, nsloc2 = nsloc2, model = "negbilog")
}

"fbvct"<- 
function(x, start, ..., sym = FALSE, nsloc1 = NULL, nsloc2 = NULL, cshape = cscale, cscale = cloc, cloc = FALSE, std.err = TRUE, corr = FALSE, method = "BFGS", warn.inf = TRUE)
{
    nlbvct <- function(loc1, scale1, shape1, loc2, scale2, shape2,
                       alpha, beta)
    {
        if(sym) beta <- alpha
        if(cshape) shape2 <- shape1
        if(cscale) scale2 <- scale1
        
        if(any(c(scale1,scale2) < 0.01) || any(c(alpha,beta) < 0.001) ||
           any(c(alpha,beta) > 30))
            return(1e6)
        if(!is.null(nsloc1)) {
            ns <- numeric(length(loc.param1))
            for(i in 1:length(ns))
                ns[i] <- get(loc.param1[i])
            loc1 <- drop(nslocmat1 %*% ns)
        }
        else loc1 <- rep(loc1, length.out = nrow(x))
        if(cloc) loc2 <- loc1 else {
          if(!is.null(nsloc2)) {
            ns <- numeric(length(loc.param2))
            for(i in 1:length(ns))
              ns[i] <- get(loc.param2[i])
            loc2 <- drop(nslocmat2 %*% ns)
          }
          else loc2 <- rep(loc2, length.out = nrow(x))
        }
        if(spx$n.m1)
            m1l <- .C("nlgev", spx$x.m1, spx$n.m1, loc1[spx$na == 2],
                scale1, shape1, dns = double(1), PACKAGE = "evd")$dns
        else m1l <- 0
        if(spx$n.m2)
            m2l <- .C("nlgev", spx$x.m2, spx$n.m2, loc2[spx$na == 1],
                scale2, shape2, dns = double(1), PACKAGE = "evd")$dns
        else m2l <- 0
        bvl <- .C("nlbvct", spx$x1, spx$x2, spx$n, spx$si, alpha, beta,
                loc1[spx$na == 0], scale1, shape1, loc2[spx$na == 0],
                scale2, shape2, cfalse, dns = double(1), PACKAGE = "evd")$dns
        if(any(is.nan(c(m1l,m2l,bvl)))) {
            warning("NaN returned in likelihood")
            return(1e6)
        }
        if(any(c(m1l,m2l,bvl) == 1e6)) return(1e6)
        else return(m1l + m2l + bvl)
    }
    if(cloc && !identical(nsloc1, nsloc2))
      stop("nsloc1 and nsloc2 must be identical")
    if(!is.null(nsloc1)) {
        if(is.vector(nsloc1)) nsloc1 <- data.frame(trend = nsloc1)
        if(nrow(nsloc1) != nrow(x))
          stop("`nsloc1' and data are not compatible")
        nslocmat1 <- cbind(1,as.matrix(nsloc1))
    }
    if(!is.null(nsloc2)) {
        if(is.vector(nsloc2)) nsloc2 <- data.frame(trend = nsloc2)
        if(nrow(nsloc2) != nrow(x))
          stop("`nsloc2' and data are not compatible")
        nslocmat2 <- cbind(1,as.matrix(nsloc2))
    }
    loc.param1 <- paste("loc1", c("",names(nsloc1)), sep="")
    loc.param2 <- paste("loc2", c("",names(nsloc2)), sep="")
    param <- c(loc.param1, "scale1", "shape1")
    if(!cloc) param <- c(param, loc.param2) else loc.param2 <- NULL
    if(!cscale) param <- c(param, "scale2")
    if(!cshape) param <- c(param, "shape2")
    if(!sym) param <- c(param, "alpha", "beta")
    else param <- c(param, "alpha")
    nmdots <- names(list(...))
    start <- bvstart.vals(x = x, start = start, nmdots = nmdots, param =
      param, nsloc1 = nsloc1, nsloc2 = nsloc2, model = "ct")
    spx <- sep.bvdata(x = x)
    cfalse <- as.integer(0)
    nm <- names(start)
    l <- length(nm)
    fixed.param <- list(...)[nmdots %in% param]
    if(any(!(param %in% c(nm,names(fixed.param)))))
        stop("unspecified parameters")
    prind <- (5:8)[c(!cscale, !cshape, TRUE, !sym)]
    f <- c(as.list(numeric(length(loc.param1))), formals(nlbvct)[2:3],
      as.list(numeric(length(loc.param2))), formals(nlbvct)[prind])
    names(f) <- param
    m <- match(nm, param)
    if(any(is.na(m))) 
        stop("`start' specifies unknown arguments")
    formals(nlbvct) <- c(f[m], f[-m])
    nllh <- function(p, ...) nlbvct(p, ...)
    if(l > 1) body(nllh) <- parse(text = paste("nlbvct(",
      paste("p[",1:l,"]", collapse = ", "), ", ...)"))
    start.arg <- c(list(p = unlist(start)), fixed.param)
    if(warn.inf && do.call("nllh", start.arg) == 1e6)
        warning("negative log-likelihood is infinite at starting values")
    opt <- optim(start, nllh, hessian = TRUE, ..., method = method)
	if(is.null(names(opt$par))) names(opt$par) <- nm
    cmar <- c(cloc, cscale, cshape)
    bvpost.optim(x = x, opt = opt, nm = nm, fixed.param = fixed.param,
      std.err = std.err, corr = corr, sym = sym, cmar = cmar,
      nsloc1 = nsloc1, nsloc2 = nsloc2, model = "ct")
}

"fbvamix"<- 
function(x, start, ..., sym = FALSE, nsloc1 = NULL, nsloc2 = NULL, cshape = cscale, cscale = cloc, cloc = FALSE, std.err = TRUE, corr = FALSE, method = "BFGS", warn.inf = TRUE)
{
    nlbvamix <- function(loc1, scale1, shape1, loc2, scale2, shape2,
                       alpha, beta)
    {
        if(cshape) shape2 <- shape1
        if(cscale) scale2 <- scale1
        
        if(any(c(scale1,scale2) < 0.01)) return(1e6)
        if(alpha < 0 || (alpha + 3*beta) < 0) return(1e6)
        if((alpha + beta) > 1 || (alpha + 2*beta) > 1) return(1e6) 
        
        if(!is.null(nsloc1)) {
            ns <- numeric(length(loc.param1))
            for(i in 1:length(ns))
                ns[i] <- get(loc.param1[i])
            loc1 <- drop(nslocmat1 %*% ns)
        }
        else loc1 <- rep(loc1, length.out = nrow(x))
        if(cloc) loc2 <- loc1 else {
          if(!is.null(nsloc2)) {
            ns <- numeric(length(loc.param2))
            for(i in 1:length(ns))
              ns[i] <- get(loc.param2[i])
            loc2 <- drop(nslocmat2 %*% ns)
          }
          else loc2 <- rep(loc2, length.out = nrow(x))
        }
        if(spx$n.m1)
            m1l <- .C("nlgev", spx$x.m1, spx$n.m1, loc1[spx$na == 2],
                scale1, shape1, dns = double(1), PACKAGE = "evd")$dns
        else m1l <- 0
        if(spx$n.m2)
            m2l <- .C("nlgev", spx$x.m2, spx$n.m2, loc2[spx$na == 1],
                scale2, shape2, dns = double(1), PACKAGE = "evd")$dns
        else m2l <- 0
        bvl <- .C("nlbvamix", spx$x1, spx$x2, spx$n, spx$si, alpha, beta,
                loc1[spx$na == 0], scale1, shape1, loc2[spx$na == 0],
                scale2, shape2, cfalse, dns = double(1),
                PACKAGE = "evd")$dns
        if(any(is.nan(c(m1l,m2l,bvl)))) {
            warning("NaN returned in likelihood")
            return(1e6)
        }
        if(any(c(m1l,m2l,bvl) == 1e6)) return(1e6)
        else return(m1l + m2l + bvl)
    }
    if(cloc && !identical(nsloc1, nsloc2))
      stop("nsloc1 and nsloc2 must be identical")
    if(!is.null(nsloc1)) {
        if(is.vector(nsloc1)) nsloc1 <- data.frame(trend = nsloc1)
        if(nrow(nsloc1) != nrow(x))
          stop("`nsloc1' and data are not compatible")
        nslocmat1 <- cbind(1,as.matrix(nsloc1))
    }
    if(!is.null(nsloc2)) {
        if(is.vector(nsloc2)) nsloc2 <- data.frame(trend = nsloc2)
        if(nrow(nsloc2) != nrow(x))
          stop("`nsloc2' and data are not compatible")
        nslocmat2 <- cbind(1,as.matrix(nsloc2))
    }
    loc.param1 <- paste("loc1", c("",names(nsloc1)), sep="")
    loc.param2 <- paste("loc2", c("",names(nsloc2)), sep="")
    param <- c(loc.param1, "scale1", "shape1")
    if(!cloc) param <- c(param, loc.param2) else loc.param2 <- NULL
    if(!cscale) param <- c(param, "scale2")
    if(!cshape) param <- c(param, "shape2")
    param <- c(param, "alpha", "beta")
    nmdots <- names(list(...))
    start <- bvstart.vals(x = x, start = start, nmdots = nmdots, param =
      param, nsloc1 = nsloc1, nsloc2 = nsloc2, model = "amix")
    spx <- sep.bvdata(x = x)
    cfalse <- as.integer(0)
    nm <- names(start)
    l <- length(nm)
    fixed.param <- list(...)[nmdots %in% param]
    if(any(!(param %in% c(nm,names(fixed.param)))))
        stop("unspecified parameters")
    prind <- (5:8)[c(!cscale, !cshape, TRUE, TRUE)]
    f <- c(as.list(numeric(length(loc.param1))), formals(nlbvamix)[2:3],
      as.list(numeric(length(loc.param2))), formals(nlbvamix)[prind])
    names(f) <- param
    m <- match(nm, param)
    if(any(is.na(m))) 
        stop("`start' specifies unknown arguments")
    formals(nlbvamix) <- c(f[m], f[-m])
    nllh <- function(p, ...) nlbvamix(p, ...)
    if(l > 1) body(nllh) <- parse(text = paste("nlbvamix(",
      paste("p[",1:l,"]", collapse = ", "), ", ...)"))
    start.arg <- c(list(p = unlist(start)), fixed.param)
    if(warn.inf && do.call("nllh", start.arg) == 1e6)
        warning("negative log-likelihood is infinite at starting values")
    opt <- optim(start, nllh, hessian = TRUE, ..., method = method)
	if(is.null(names(opt$par))) names(opt$par) <- nm
    cmar <- c(cloc, cscale, cshape)
    bvpost.optim(x = x, opt = opt, nm = nm, fixed.param = fixed.param,
      std.err = std.err, corr = corr, sym = sym, cmar = cmar,
      nsloc1 = nsloc1, nsloc2 = nsloc2, model = "amix")
}

### Method Function ###

"print.bvevd" <-  function(x, digits = max(3, getOption("digits") - 3), ...) 
{
    cat("\nCall:", deparse(x$call), "\n")
    cat("Deviance:", deviance(x), "\n")
    cat("AIC:", AIC(x), "\n")
    if(!is.null(x$dep.summary)) cat("Dependence:", x$dep.summary, "\n")
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

### Ancillary Functions ###

"bvstart.vals" <-
# Calculate Starting Values For Bivariate Models
function(x, start, nmdots, param, method = c("evd","pot"), nsloc1 = NULL,
         nsloc2 = NULL, u = NULL, model)
{
  method <- match.arg(method)
  if(method == "evd") {
    loc.param1 <- paste("loc1", c("",names(nsloc1)), sep="")
    loc.param2 <- paste("loc2", c("",names(nsloc2)), sep="")
  }
  if(method == "pot") loc.param1 <- loc.param2 <- NULL
  if(missing(start)) {
    start <- as.list(numeric(length(param)))
    names(start) <- param
    if(method == "evd") {
      loc.param1 <- paste("loc1", c("",names(nsloc1)), sep="")
      loc.param2 <- paste("loc2", c("",names(nsloc2)), sep="")
      st1 <- fitted(fgev(x[,1], nsloc = nsloc1, std.err = FALSE))
      st2 <- fitted(fgev(x[,2], nsloc = nsloc2, std.err = FALSE))
      st1 <- as.list(st1); st2 <- as.list(st2)
    }
    if(method == "pot") {
      st1 <- fitted(fpot(x[,1], u[1], std.err = FALSE))
      st2 <- fitted(fpot(x[,2], u[2], std.err = FALSE))
      st1 <- as.list(st1); st2 <- as.list(st2)
    }
    start[c(loc.param1, "scale1", "shape1")] <- st1
    tmp2 <- loc.param2
    if("scale2" %in% param) tmp2 <- c(tmp2, "scale2")
    if("shape2" %in% param) tmp2 <- c(tmp2, "shape2")
    tmp <- sub("2", "", tmp2)
    start[tmp2] <- st2[tmp]
    if(model == "log") start[["dep"]] <- 0.75
    if(model == "alog") {
        start[["asy1"]] <- 0.75
        if("asy2" %in% param) start[["asy2"]] <- 0.75 
        start[["dep"]] <- 0.65
    }
    if(model == "hr") start[["dep"]] <- 1
    if(model == "neglog") start[["dep"]] <- 0.6
    if(model == "aneglog") {
        start[["asy1"]] <- 0.75
        if("asy2" %in% param) start[["asy2"]] <- 0.75 
        start[["dep"]] <- 0.8
    }
    if(model == "bilog") start[["alpha"]] <- start[["beta"]] <- 0.75
    if(model == "negbilog") start[["alpha"]] <- start[["beta"]] <- 1/0.6
    if(model == "ct") {
      start[["alpha"]] <- 0.6
      if("beta" %in% param) start[["beta"]] <- 0.6 
    }
    if(model == "amix") {
      start[["alpha"]] <- 0.75
      start[["beta"]] <- 0 
    }
    start <- start[!(param %in% nmdots)]
  }
  if(any(!is.na(match(names(start),c("mar1","mar2","asy"))))) {
    if(("mar1" %in% names(start)) && (length(start$mar1) !=
      (2+length(loc.param1)))) stop("mar1 in `start' has incorrect length")
    if(("mar2" %in% names(start)) && (length(start$mar2) !=
      (2+length(loc.param2)))) stop("mar2 in `start' has incorrect length")
    if(("asy" %in% names(start)) && (length(start$asy) != 2))
      stop("asy in `start' should have length two")
    start <- unlist(start)
    names(start)[grep("mar1",names(start))] <- c(loc.param1,"scale1","shape1")
    names(start)[grep("mar2",names(start))] <- c(loc.param2,"scale2","shape2")
    start <- as.list(start)
  }
  if(!is.list(start)) 
    stop("`start' must be a named list")
  start
}

"sep.bvdata" <-
# Separate Bivariate Data For Bivariate Models
function(x, method = c("evd","cpot","ppot"), u = NULL)
{
    method <- match.arg(method)
    if(method == "evd") {
      na <- rowSums(cbind(is.na(x[,1]), 2*is.na(x[,2])))
      if(!any(na == 0))
        stop("`x' must have at least one complete observation")
      x.m1 <- as.double(x[na == 2, 1])
      n.m1 <- as.integer(length(x.m1))
      x.m2 <- as.double(x[na == 1, 2])
      n.m2 <- as.integer(length(x.m2))
      x.full <- x[na == 0, , drop = FALSE]
      x1 <- as.double(x.full[,1])
      x2 <- as.double(x.full[,2])      
      n <- as.integer(nrow(x.full))
      if(ncol(x) == 3) {
        si <- x.full[,3]
        si[is.na(si)] <- 2
      }
      else si <- rep(2, n)
      si <- as.integer(si)
      spx <- list(x.m1 = x.m1, n.m1 = n.m1, x.m2 = x.m2, n.m2 = n.m2,
        x1 = x1, x2 = x2, n = n, si = si, na = na)
    } else {
      x1 <- x[,1]
      x2 <- x[,2]
      n <- length(x1)
      r1 <- r2 <- NULL
      iau1 <- (x1 > u[1]) & !is.na(x1)
      iau2 <- (x2 > u[2]) & !is.na(x2)
      nat <- c(sum(iau1), sum(iau2), sum(iau1 & iau2))
      lambda1 <- sum(iau1) / (n + 1)
      lambda2 <- sum(iau2) / (n + 1)
      lambda <- c(lambda1, lambda2)
      if(method == "ppot") {
        x1[is.na(x1)] <- mean(x1[!iau1], na.rm = TRUE)
        x2[is.na(x2)] <- mean(x2[!iau2], na.rm = TRUE)
        r1 <- 1 - rank(x1) / (n + 1)
        r2 <- 1 - rank(x2) / (n + 1)
        r1[iau1] <- lambda1
        r2[iau2] <- lambda2
      }
      x1 <- x1 - u[1]
      x2 <- x2 - u[2]
      x1[!iau1] <- 0
      x2[!iau2] <- 0
      i0 <- iau1 | iau2
      x1 <- x1[i0] ; x2 <- x2[i0]
      if(method == "ppot") {
       r1 <- r1[i0] ; r2 <- r2[i0]
      }
      nn <- length(x1)
      thdi <- as.logical(x1) + 2*as.logical(x2)
      spx <- list(x1 = x1, x2 = x2, nn = nn, n = n, thdi = thdi, lambda =
        lambda, r1 = r1, r2 = r2, nat = nat)
    }
    spx
}

"bvpost.optim" <-
# Post-optimization Processing
function(x, opt, nm, fixed.param, std.err, corr, sym, cmar, method = c("evd","pot"), nsloc1 = NULL, nsloc2 = NULL, u = NULL, nat = NULL, likelihood = NULL, model)
{
    method <- match.arg(method)
    if(opt$convergence != 0) {
      warning(paste("optimization for", model, "may not have succeeded"), call. = FALSE)
      if(opt$convergence == 1) opt$convergence <- "iteration limit reached"
    }
    else opt$convergence <- "successful"
    if(std.err) {
        tol <- .Machine$double.eps^0.5
        var.cov <- qr(opt$hessian, tol = tol)
        if(var.cov$rank != ncol(var.cov$qr)) 
            stop(paste("observed information matrix for", model,
                       "is singular; use std.err = FALSE"))
        var.cov <- solve(var.cov, tol = tol)
        std.err <- diag(var.cov)
        if(any(std.err <= 0))
            stop(paste("observed information matrix for", model,
                       "is singular; use std.err = FALSE"))
        std.err <- sqrt(std.err)
        names(std.err) <- nm
        if(corr) {
            .mat <- diag(1/std.err, nrow = length(std.err))
            corr <- structure(.mat %*% var.cov %*% .mat, dimnames = list(nm,nm))
            diag(corr) <- rep(1, length(std.err))
        }
        else corr <- NULL
    }
    else std.err <- var.cov <- corr <- NULL
    fixed <- unlist(fixed.param)
    param <- c(opt$par, fixed)
    fixed2 <- NULL
    if(method == "evd") {
      if(cmar[1]) {
        loc.param1 <- paste("loc1", c("",names(nsloc1)), sep="")
        fixed2 <- c(fixed2, param[loc.param1])
      }
      if(cmar[2]) fixed2 <- c(fixed2, param["scale1"])
      if(cmar[3]) fixed2 <- c(fixed2, param["shape1"])
    }
    if(method == "pot") {
      if(cmar[1]) fixed2 <- c(fixed2, param["scale1"])
      if(cmar[2]) fixed2 <- c(fixed2, param["shape1"])
    }
    if(sym) {
      if(model %in% c("alog","aneglog")) fixed2 <- c(fixed2, param["asy1"])
      if(model == "ct") fixed2 <- c(fixed2, param["alpha"])
    }
    if(!is.null(fixed2)) {
      names(fixed2) <- sub("1", "2", names(fixed2))
      names(fixed2) <- sub("alpha", "beta", names(fixed2))
    }
    param <- c(param, fixed2)
    # Transform to stationarity
    x2 <- x
    if(!is.null(nsloc1)) {
        trend <- param[paste("loc1", names(nsloc1), sep="")]
        trend <- drop(as.matrix(nsloc1) %*% trend)
        x2[,1] <- x[,1] - trend
    }
    if(!is.null(nsloc2)) {
        trend <- param[paste("loc2", names(nsloc2), sep="")]
        trend <- drop(as.matrix(nsloc2) %*% trend)
        x2[,2] <- x[,2] - trend
    }
    # End transform
    # Dependence chi
    if(model %in% c("log", "hr", "neglog")) {
      dep <- param["dep"]
      dep.sum <- 2*(1 - abvevd(dep = dep, model = model))
    }
    if(model %in% c("alog", "aneglog")) {
      dep <- param["dep"]
      asy <- param[c("asy1", "asy2")]
      dep.sum <- 2*(1 - abvevd(dep = dep, asy = asy, model = model))
    }
    if(model %in% c("bilog", "negbilog", "ct", "amix")) {
      alpha <- param["alpha"]
      beta <- param["beta"]
      dep.sum <- 2*(1-abvevd(alpha = alpha, beta = beta, model = model))
    }
    # End dependence chi
    out <- list(estimate = opt$par, std.err = std.err, fixed = fixed,
    fixed2 = fixed2, param = param, deviance = 2*opt$value,
    dep.summary = dep.sum, corr = corr, var.cov = var.cov, convergence =
    opt$convergence, counts = opt$counts, message = opt$message, data = x)
    if(method == "evd")
      out <- c(out, list(tdata = x2, nsloc1 = nsloc1, nsloc2 = nsloc2))
    if(method == "pot")
      out <- c(out, list(threshold = u, nat = nat, likelihood = likelihood))
    c(out, list(n = nrow(x), sym = sym, cmar = cmar, model = model))
}





