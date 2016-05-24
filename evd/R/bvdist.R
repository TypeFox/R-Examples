
"rbvevd" <-
function(n, dep, asy = c(1,1), alpha, beta, model = c("log", "alog",
    "hr", "neglog", "aneglog", "bilog", "negbilog", "ct", "amix"),
    mar1 = c(0,1,0), mar2 = mar1)
{
  model <- match.arg(model)
  m1 <- c("bilog", "negbilog", "ct", "amix")
  m2 <- c(m1, "log", "hr", "neglog")
  m3 <- c("log", "alog", "hr", "neglog", "aneglog")
  if((model %in% m1) && !missing(dep))
    warning("ignoring `dep' argument")
  if((model %in% m2) && !missing(asy))
    warning("ignoring `asy' argument")
  if((model %in% m3) && !missing(alpha))
    warning("ignoring `alpha' argument")
  if((model %in% m3) && !missing(beta))
    warning("ignoring `beta' argument")
    
  switch(model,
    log = rbvlog(n = n, dep = dep, mar1 = mar1, mar2 = mar2),
    alog = rbvalog(n = n, dep = dep, asy = asy, mar1 = mar1, mar2 = mar2),
    hr = rbvhr(n = n, dep = dep, mar1 = mar1, mar2 = mar2),
    neglog = rbvneglog(n = n, dep = dep, mar1 = mar1, mar2 = mar2),
    aneglog = rbvaneglog(n = n, dep = dep, asy = asy, mar1 = mar1,
      mar2 = mar2),
    bilog = rbvbilog(n = n, alpha = alpha, beta = beta, mar1 = mar1,
      mar2 = mar2),
    negbilog = rbvnegbilog(n = n, alpha = alpha, beta = beta, mar1 = mar1,
      mar2 = mar2),
    ct = rbvct(n = n, alpha = alpha, beta = beta, mar1 = mar1, mar2 = mar2),
    amix = rbvamix(n = n, alpha = alpha, beta = beta, mar1 = mar1,
      mar2 = mar2)) 
}

"rbvlog"<-
# Uses Algorithm 1.1 in Stephenson(2003)
function(n, dep, mar1 = c(0,1,0), mar2 = mar1)
{
    if(length(dep) != 1 || mode(dep) != "numeric" || dep <= 0 ||
        dep > 1) stop("invalid argument for `dep'")
    sim <- .C("rbvlog_shi",
               as.integer(n), as.double(dep), sim = double(2*n),
               PACKAGE = "evd")$sim
    sim <- matrix(sim, nrow = n, ncol = 2, byrow = TRUE)
    mtransform(1/sim, list(mar1, mar2), inv = TRUE, drp = TRUE)
}

"rbvalog"<-
# Uses Algorithm 1.2 in Stephenson(2003)
function(n, dep, asy = c(1,1), mar1 = c(0,1,0), mar2 = mar1)
{
    if(length(dep) != 1 || mode(dep) != "numeric" || dep <= 0 ||
       dep > 1) stop("invalid argument for `dep'")
    if(length(asy) != 2 || mode(asy) != "numeric" || min(asy) < 0 ||
       max(asy) > 1) stop("invalid argument for `asy'")
    if(dep == 1 || any(asy == 0)) {
        asy <- c(0,0)
        dep <- 1
    }
    sim <- .C("rbvalog_shi",
              as.integer(n), as.double(dep), as.double(asy),
              sim = double(2*n), PACKAGE = "evd")$sim
    sim <- matrix(sim, nrow = n, ncol = 2, byrow = TRUE)
    mtransform(1/sim, list(mar1, mar2), inv = TRUE, drp = TRUE)
}

"rbvhr" <-
function(n, dep, mar1 = c(0,1,0), mar2 = mar1)
{
    if(length(dep) != 1 || mode(dep) != "numeric" || dep <= 0)
        stop("invalid argument for `dep'")
    sim <- .C("rbvhr",
        as.integer(n), as.double(dep), sim = runif(2*n),
        PACKAGE = "evd")$sim
    sim <- matrix(sim, nrow = n, ncol = 2, byrow = TRUE)
    mtransform(-log(sim), list(mar1, mar2), inv = TRUE, drp = TRUE)
}

"rbvneglog"<- 
function(n, dep, mar1 = c(0,1,0), mar2 = mar1)
{
    if(length(dep) != 1 || mode(dep) != "numeric" || dep <= 0)
        stop("invalid argument for `dep'")
    sim <- .C("rbvneglog",
        as.integer(n), as.double(dep), sim = runif(2*n),
        PACKAGE = "evd")$sim
    sim <- matrix(sim, nrow = n, ncol = 2, byrow = TRUE)
    mtransform(-log(sim), list(mar1, mar2), inv = TRUE, drp = TRUE)
}

"rbvaneglog"<- 
function(n, dep, asy = c(1,1), mar1 = c(0,1,0), mar2 = mar1)
{
    if(length(dep) != 1 || mode(dep) != "numeric" || dep <= 0)
        stop("invalid argument for `dep'")
    if(length(asy) != 2 || mode(asy) != "numeric" || min(asy) < 0 ||
       max(asy) > 1) stop("invalid argument for `asy'")
    sim <- .C("rbvaneglog",
        as.integer(n), as.double(dep), as.double(asy), sim = runif(2*n),
        PACKAGE = "evd")$sim
    sim <- matrix(sim, nrow = n, ncol = 2, byrow = TRUE)
    mtransform(-log(sim), list(mar1, mar2), inv = TRUE, drp = TRUE)
}

"rbvbilog"<- 
function(n, alpha, beta, mar1 = c(0,1,0), mar2 = mar1)
{
    if(length(alpha) != 1 || mode(alpha) != "numeric")
        stop("invalid argument for `alpha'")
    if(length(beta) != 1 || mode(beta) != "numeric")
        stop("invalid argument for `beta'")
    if(any(c(alpha,beta) <= 0) || any(c(alpha,beta) >= 1))
        stop("`alpha' and `beta' must be in the open interval (0,1)")
    sim <- .C("rbvbilog",
        as.integer(n), as.double(alpha), as.double(beta), sim = runif(2*n),
        PACKAGE = "evd")$sim
    sim <- matrix(sim, nrow = n, ncol = 2, byrow = TRUE)
    mtransform(-log(sim), list(mar1, mar2), inv = TRUE, drp = TRUE)
}

"rbvnegbilog"<- 
function(n, alpha, beta, mar1 = c(0,1,0), mar2 = mar1)
{
    if(length(alpha) != 1 || mode(alpha) != "numeric")
        stop("invalid argument for `alpha'")
    if(length(beta) != 1 || mode(beta) != "numeric")
        stop("invalid argument for `beta'")
    if(any(c(alpha,beta) <= 0))
        stop("`alpha' and `beta' must be non-negative")
    sim <- .C("rbvnegbilog",
        as.integer(n), as.double(alpha), as.double(beta), sim = runif(2*n),
        PACKAGE = "evd")$sim
    sim <- matrix(sim, nrow = n, ncol = 2, byrow = TRUE)
    mtransform(-log(sim), list(mar1, mar2), inv = TRUE, drp = TRUE)
}

"rbvct" <-
function(n, alpha, beta, mar1 = c(0,1,0), mar2 = mar1)
{
    if(length(alpha) != 1 || mode(alpha) != "numeric")
        stop("invalid argument for `alpha'")
    if(length(beta) != 1 || mode(beta) != "numeric")
        stop("invalid argument for `beta'")
    if(any(c(alpha,beta) <= 0))
        stop("`alpha' and `beta' must be non-negative")
    sim <- .C("rbvct",
        as.integer(n), as.double(alpha), as.double(beta), sim = runif(2*n),
        PACKAGE = "evd")$sim
    sim <- matrix(sim, nrow = n, ncol = 2, byrow = TRUE)
    mtransform(-log(sim), list(mar1, mar2), inv = TRUE, drp = TRUE)
}

"rbvamix" <-
function(n, alpha, beta, mar1 = c(0,1,0), mar2 = mar1)
{
    if(length(alpha) != 1 || mode(alpha) != "numeric")
        stop("invalid argument for `alpha'")
    if(length(beta) != 1 || mode(beta) != "numeric")
        stop("invalid argument for `beta'")
    if(alpha < 0)
        stop("`alpha' must be non-negative")
    if((alpha + beta) > 1)
        stop("`alpha' + `beta' cannot be greater than one")
    if((alpha + 2*beta) > 1)
        stop("`alpha' + `2*beta' cannot be greater than one")
    if((alpha + 3*beta) < 0)
        stop("`alpha' + `3*beta' must be non-negative")
    sim <- .C("rbvamix",
        as.integer(n), as.double(alpha), as.double(beta),
        sim = runif(2*n), PACKAGE = "evd")$sim
    sim <- matrix(sim, nrow = n, ncol = 2, byrow = TRUE)
    mtransform(-log(sim), list(mar1, mar2), inv = TRUE, drp = TRUE)
}

"pbvevd" <-
function(q, dep, asy = c(1,1), alpha, beta, model = c("log", "alog",
    "hr", "neglog", "aneglog", "bilog", "negbilog", "ct", "amix"),
    mar1 = c(0,1,0), mar2 = mar1, lower.tail = TRUE)
{
  model <- match.arg(model)
  m1 <- c("bilog", "negbilog", "ct", "amix")
  m2 <- c(m1, "log", "hr", "neglog")
  m3 <- c("log", "alog", "hr", "neglog", "aneglog")
  if((model %in% m1) && !missing(dep))
    warning("ignoring `dep' argument")
  if((model %in% m2) && !missing(asy))
    warning("ignoring `asy' argument")
  if((model %in% m3) && !missing(alpha))
    warning("ignoring `alpha' argument")
  if((model %in% m3) && !missing(beta))
    warning("ignoring `beta' argument")
    
  switch(model,
    log = pbvlog(q = q, dep = dep, mar1 = mar1, mar2 = mar2,
      lower.tail = lower.tail),
    alog = pbvalog(q = q, dep = dep, asy = asy, mar1 = mar1, mar2 = mar2,
      lower.tail = lower.tail),
    hr = pbvhr(q = q, dep = dep, mar1 = mar1, mar2 = mar2,
      lower.tail = lower.tail),
    neglog = pbvneglog(q = q, dep = dep, mar1 = mar1, mar2 = mar2,
      lower.tail = lower.tail),
    aneglog = pbvaneglog(q = q, dep = dep, asy = asy, mar1 = mar1,
      mar2 = mar2, lower.tail = lower.tail),
    bilog = pbvbilog(q = q, alpha = alpha, beta = beta, mar1 = mar1,
      mar2 = mar2, lower.tail = lower.tail),
    negbilog = pbvnegbilog(q = q, alpha = alpha, beta = beta, mar1 = mar1,
      mar2 = mar2, lower.tail = lower.tail),
    ct = pbvct(q = q, alpha = alpha, beta = beta, mar1 = mar1,
      mar2 = mar2, lower.tail = lower.tail),
    amix = pbvamix(q = q, alpha = alpha, beta = beta, mar1 = mar1,
      mar2 = mar2, lower.tail = lower.tail)) 
}

"pbvlog"<- 
function(q, dep, mar1 = c(0,1,0), mar2 = mar1, lower.tail = TRUE)
{
    if(length(dep) != 1 || mode(dep) != "numeric" || dep <= 0 ||
        dep > 1) stop("invalid argument for `dep'")
    if(is.null(dim(q))) dim(q) <- c(1,2)
    q <- mtransform(q, list(mar1, mar2))
    v <- apply(q^(1/dep),1,sum)^dep
    pp <- exp(-v)
    if(!lower.tail) {
      pp <- 1 - pgev(-log(q[,1])) - pgev(-log(q[,2])) + pp
    }
    pp
}

"pbvalog"<- 
function(q, dep, asy = c(1,1), mar1 = c(0,1,0), mar2 = mar1, lower.tail = TRUE)
{
    if(length(dep) != 1 || mode(dep) != "numeric" || dep <= 0 ||
        dep > 1) stop("invalid argument for `dep'")
    if(length(asy) != 2 || mode(asy) != "numeric" || min(asy) < 0 ||
       max(asy) > 1) stop("invalid argument for `asy'")
    if(is.null(dim(q))) dim(q) <- c(1,2)
    q <- mtransform(q, list(mar1, mar2))
    asy <- rep(asy,rep(nrow(q),2))
    v <- apply((asy*q)^(1/dep),1,sum)^dep + apply((1-asy)*q,1,sum)
    pp <- exp(-v)
    if(!lower.tail) {
      pp <- 1 - pgev(-log(q[,1])) - pgev(-log(q[,2])) + pp
    }
    pp
}

"pbvhr" <-
function(q, dep, mar1 = c(0,1,0), mar2 = mar1, lower.tail = TRUE)
{
    if(length(dep) != 1 || mode(dep) != "numeric" || dep <= 0)
        stop("invalid argument for `dep'")
    if(is.null(dim(q))) dim(q) <- c(1,2)
    q <- mtransform(q, list(mar1, mar2))
    fn <- function(x1,x2) x1*pnorm(1/dep + dep * log(x1/x2) / 2)
    v <- fn(q[,1],q[,2]) + fn(q[,2],q[,1])
    pp <- exp(-v)
    if(!lower.tail) {
      pp <- 1 - pgev(-log(q[,1])) - pgev(-log(q[,2])) + pp
    }
    pp
}

"pbvneglog"<- 
function(q, dep, mar1 = c(0,1,0), mar2 = mar1, lower.tail = TRUE)
{
    if(length(dep) != 1 || mode(dep) != "numeric" || dep <= 0)
        stop("invalid argument for `dep'")
    if(is.null(dim(q))) dim(q) <- c(1,2)
    q <- mtransform(q, list(mar1, mar2))
    v <- apply(q,1,sum) - apply(q^(-dep),1,sum)^(-1/dep)
    pp <- exp(-v)
    if(!lower.tail) {
      pp <- 1 - pgev(-log(q[,1])) - pgev(-log(q[,2])) + pp
    }
    pp
}

"pbvaneglog"<- 
function(q, dep, asy = c(1,1), mar1 = c(0,1,0), mar2 = mar1, lower.tail = TRUE)
{
    if(length(dep) != 1 || mode(dep) != "numeric" || dep <= 0)
        stop("invalid argument for `dep'")
    if(length(asy) != 2 || mode(asy) != "numeric" || min(asy) < 0 ||
       max(asy) > 1) stop("invalid argument for `asy'")
    if(is.null(dim(q))) dim(q) <- c(1,2)
    q <- mtransform(q, list(mar1, mar2))
    asy <- rep(asy,rep(nrow(q),2))
    v <- apply(q,1,sum) - apply((asy*q)^(-dep),1,sum)^(-1/dep)
    pp <- exp(-v)
    if(!lower.tail) {
      pp <- 1 - pgev(-log(q[,1])) - pgev(-log(q[,2])) + pp
    }
    pp
}

"pbvbilog"<- 
function(q, alpha, beta, mar1 = c(0,1,0), mar2 = mar1, lower.tail = TRUE)
{
    if(length(alpha) != 1 || mode(alpha) != "numeric")
        stop("invalid argument for `alpha'")
    if(length(beta) != 1 || mode(beta) != "numeric")
        stop("invalid argument for `beta'")
    if(any(c(alpha,beta) <= 0) || any(c(alpha,beta) >= 1))
        stop("`alpha' and `beta' must be in the open interval (0,1)")
    if(is.null(dim(q))) dim(q) <- c(1,2)
    q <- mtransform(q, list(mar1, mar2))
    gma <- numeric(nrow(q))
    for(i in 1:nrow(q)) {
        gmafn <- function(x)
            (1-alpha) * q[i,1] * (1-x)^beta - (1-beta) * q[i,2] * x^alpha
        if(any(is.na(q[i,]))) gma[i] <- NA
        else if(any(is.infinite(q[i,]))) gma[i] <- 0.5
        else if(q[i,1] == 0) gma[i] <- 0
        else if(q[i,2] == 0) gma[i] <- 1
        else gma[i] <- uniroot(gmafn, lower = 0, upper = 1,
                               tol = .Machine$double.eps^0.5)$root
    }
    v <- q[,1] * gma^(1-alpha) + q[,2] * (1 - gma)^(1-beta)
    pp <- exp(-v)
    if(!lower.tail) {
      pp <- 1 - pgev(-log(q[,1])) - pgev(-log(q[,2])) + pp
    }
    pp
}

"pbvnegbilog"<- 
function(q, alpha, beta, mar1 = c(0,1,0), mar2 = mar1, lower.tail = TRUE)
{
    if(length(alpha) != 1 || mode(alpha) != "numeric")
        stop("invalid argument for `alpha'")
    if(length(beta) != 1 || mode(beta) != "numeric")
        stop("invalid argument for `beta'")
    if(any(c(alpha,beta) <= 0))
        stop("`alpha' and `beta' must be non-negative")
    if(is.null(dim(q))) dim(q) <- c(1,2)
    q <- mtransform(q, list(mar1, mar2))
    gma <- numeric(nrow(q))
    for(i in 1:nrow(q)) {
        gmafn <- function(x)
            (1+alpha) * q[i,1] * x^alpha - (1+beta) * q[i,2] * (1-x)^beta
        if(any(is.na(q[i,]))) gma[i] <- NA
        else if(any(is.infinite(q[i,]))) gma[i] <- Inf
        else if(q[i,1] == 0) gma[i] <- 1
        else if(q[i,2] == 0) gma[i] <- 0
        else gma[i] <- uniroot(gmafn, lower = 0, upper = 1,
                               tol = .Machine$double.eps^0.5)$root
    }
    v <- q[,1] + q[,2] - q[,1] * gma^(1+alpha) - q[,2] * (1 - gma)^(1+beta)
    v[is.infinite(gma)] <- Inf
    pp <- exp(-v)
    if(!lower.tail) {
      pp <- 1 - pgev(-log(q[,1])) - pgev(-log(q[,2])) + pp
    }
    pp
}

"pbvct" <-
function(q, alpha, beta, mar1 = c(0,1,0), mar2 = mar1, lower.tail = TRUE)
{
    if(length(alpha) != 1 || mode(alpha) != "numeric")
        stop("invalid argument for `alpha'")
    if(length(beta) != 1 || mode(beta) != "numeric")
        stop("invalid argument for `beta'")
    if(any(c(alpha,beta) <= 0))
        stop("`alpha' and `beta' must be non-negative")
    if(is.null(dim(q))) dim(q) <- c(1,2)
    q <- mtransform(q, list(mar1, mar2))
    u <- (alpha * q[,2]) / (alpha * q[,2] + beta * q[,1])  
    v <- q[,2] * pbeta(u, shape1 = alpha, shape2 = beta + 1) +
      q[,1] * pbeta(u, shape1 = alpha + 1, shape2 = beta, lower.tail = FALSE)
    v[is.infinite(q[,1]) || is.infinite(q[,2])] <- Inf
    v[(q[,1] + q[,2]) == 0] <- 0
    pp <- exp(-v)
    if(!lower.tail) {
      pp <- 1 - pgev(-log(q[,1])) - pgev(-log(q[,2])) + pp
    }
    pp
}

"pbvamix"<- 
function(q, alpha, beta, mar1 = c(0,1,0), mar2 = mar1, lower.tail = TRUE)
{
    if(length(alpha) != 1 || mode(alpha) != "numeric")
        stop("invalid argument for `alpha'")
    if(length(beta) != 1 || mode(beta) != "numeric")
        stop("invalid argument for `beta'")
    if(alpha < 0)
        stop("`alpha' must be non-negative")
    if((alpha + beta) > 1)
        stop("`alpha' + `beta' cannot be greater than one")
    if((alpha + 2*beta) > 1)
        stop("`alpha' + `2*beta' cannot be greater than one")
    if((alpha + 3*beta) < 0)
        stop("`alpha' + `3*beta' must be non-negative")
    if(is.null(dim(q))) dim(q) <- c(1,2)
    q <- mtransform(q, list(mar1, mar2))
    qsum <- apply(q, 1, sum)  
    v <- qsum - (alpha + beta) * q[,1] + alpha * (q[,1]^2)/qsum +
      beta * (q[,1]^3)/(qsum^2)
    v[is.infinite(q[,1]) || is.infinite(q[,2])] <- Inf
    v[(q[,1] + q[,2]) == 0] <- 0
    pp <- exp(-v)
    if(!lower.tail) {
      pp <- 1 - pgev(-log(q[,1])) - pgev(-log(q[,2])) + pp
    }
    pp
}

"abvevd" <-
function(x = 0.5, dep, asy = c(1,1), alpha, beta, model = c("log", "alog",
    "hr", "neglog", "aneglog", "bilog", "negbilog", "ct", "amix"),
     rev = FALSE, plot = FALSE, add = FALSE, lty = 1, lwd = 1, col = 1,
     blty = 3, blwd = 1, xlim = c(0,1), ylim = c(0.5,1), xlab = "t",
     ylab = "A(t)", ...)
{
  if(any(x < 0, na.rm = TRUE) || any(x > 1, na.rm = TRUE))
    stop("invalid argument for `x'")
  model <- match.arg(model)
  m1 <- c("bilog", "negbilog", "ct", "amix")
  m2 <- c(m1, "log", "hr", "neglog")
  m3 <- c("log", "alog", "hr", "neglog", "aneglog")
  if((model %in% m1) && !missing(dep))
    warning("ignoring `dep' argument")
  if((model %in% m2) && !missing(asy))
    warning("ignoring `asy' argument")
  if((model %in% m3) && !missing(alpha))
    warning("ignoring `alpha' argument")
  if((model %in% m3) && !missing(beta))
    warning("ignoring `beta' argument")

  if(rev && (model %in% c("aneglog", "alog"))) asy <- asy[2:1]
  if(rev && (model %in% c("bilog", "negbilog", "ct"))) {
    tmpalpha <- alpha
    alpha <- beta
    beta <- tmpalpha
  }
  if(rev && (model == "amix")) {
    tmpalpha <- alpha
    alpha <- alpha + 3*beta
    beta <- -beta
  }
     
  switch(model,
    log = abvlog(x = x, dep = dep, plot = plot, add = add, lty = lty,
      lwd = lwd, col = col, blty = blty, blwd = blwd, xlim = xlim, ylim = ylim,
      xlab = xlab, ylab = ylab, ...),
    alog = abvalog(x = x, dep = dep, asy = asy, plot = plot, add = add,
      lty = lty, lwd = lwd, col = col, blty = blty, blwd = blwd, xlim = xlim,
      ylim = ylim, xlab = xlab, ylab = ylab, ...),
    hr = abvhr(x = x, dep = dep, plot = plot, add = add, lty = lty,
      lwd = lwd, col = col, blty = blty, blwd = blwd, xlim = xlim, ylim = ylim,
      xlab = xlab, ylab = ylab, ...),
    neglog = abvneglog(x = x, dep = dep, plot = plot, add = add, lty = lty,
      lwd = lwd, col = col, blty = blty, blwd = blwd, xlim = xlim, ylim = ylim,
      xlab = xlab, ylab = ylab, ...),
    aneglog = abvaneglog(x = x, dep = dep, asy = asy, plot = plot, add = add,
      lty = lty, lwd = lwd, col = col, blty = blty, blwd = blwd, xlim = xlim,
      ylim = ylim, xlab = xlab, ylab = ylab, ...),
    bilog = abvbilog(x = x, alpha = alpha, beta = beta, plot = plot,
      add = add, lty = lty, lwd = lwd, col = col, blty = blty, blwd = blwd,
      xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, ...),
    negbilog = abvnegbilog(x = x, alpha = alpha, beta = beta, plot = plot,
      add = add, lty = lty, lwd = lwd, col = col, blty = blty, blwd = blwd,
      xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, ...),
    ct = abvct(x = x, alpha = alpha, beta = beta, plot = plot, add = add,
      lty = lty, lwd = lwd, col = col, blty = blty, blwd = blwd, xlim = xlim,
      ylim = ylim, xlab = xlab, ylab = ylab, ...),
    amix = abvamix(x = x, alpha = alpha, beta = beta, plot = plot, add = add,
      lty = lty, lwd = lwd, col = col, blty = blty, blwd = blwd, xlim = xlim,
      ylim = ylim, xlab = xlab, ylab = ylab, ...)) 
}

"abvlog"<- 
function(x = 0.5, dep, plot = FALSE, add = FALSE, lty = 1, lwd = 1, col = 1,
         blty = 3, blwd = 1, xlim = c(0,1), ylim = c(0.5,1), xlab = "",
         ylab = "", ...)
{
    if(length(dep) != 1 || mode(dep) != "numeric" || dep <= 0 ||
        dep > 1) stop("invalid argument for `dep'")
    if(plot || add) x <- seq(0, 1, length = 100)
    idep <- 1/dep
    a <- (x^idep + (1-x)^idep)^dep
    if(plot || add) {
      if(!add)  { 
        plot(x, a, type="n", xlab = xlab, ylab = ylab,
             xlim = xlim, ylim = ylim, ...) 
        polygon(c(0, 0.5, 1), c(1, 0.5, 1), lty = blty, lwd = blwd)  
      }
      lines(x, a, lty = lty, lwd = lwd, col = col) 
      return(invisible(list(x = x, y = a)))
    }
    a
}

"abvalog"<- 
function(x = 0.5, dep, asy = c(1,1), plot = FALSE, add = FALSE,
         lty = 1, lwd = 1, col = 1, blty = 3, blwd = 1, xlim = c(0,1),
         ylim = c(0.5,1), xlab = "", ylab = "", ...)
{
    if(length(dep) != 1 || mode(dep) != "numeric" || dep <= 0 ||
        dep > 1) stop("invalid argument for `dep'")
    if(length(asy) != 2 || mode(asy) != "numeric" || min(asy) < 0 ||
       max(asy) > 1) stop("invalid argument for `asy'")
    if(plot || add) x <- seq(0, 1, length = 100)
    idep <- 1/dep
    a <- ((asy[1]*x)^idep + (asy[2]*(1-x))^idep)^dep +
        (1-asy[1])*x + (1-asy[2])*(1-x)    
    if(plot || add) {
      if(!add)  { 
        plot(x, a, type="n", xlab = xlab, ylab = ylab,
             xlim = xlim, ylim = ylim, ...) 
        polygon(c(0, 0.5, 1), c(1, 0.5, 1), lty = blty, lwd = blwd)  
      }
      lines(x, a, lty = lty, lwd = lwd, col = col) 
      return(invisible(list(x = x, y = a)))
    }
    a
}

"abvhr" <-
function(x = 0.5, dep, plot = FALSE, add = FALSE, lty = 1, lwd = 1, col = 1,
         blty = 3, blwd = 1, xlim = c(0,1), ylim = c(0.5,1), xlab = "",
         ylab = "", ...)
{
    if(length(dep) != 1 || mode(dep) != "numeric" || dep <= 0)
        stop("invalid argument for `dep'")
    if(plot || add) x <- seq(0, 1, length = 100)
    fn <- function(z) z*pnorm(1/dep + dep * log(z/(1-z)) / 2)
    a <- fn(x) + fn(1-x)
    if(plot || add) {
      if(!add)  { 
        plot(x, a, type="n", xlab = xlab, ylab = ylab,
             xlim = xlim, ylim = ylim, ...) 
        polygon(c(0, 0.5, 1), c(1, 0.5, 1), lty = blty, lwd = blwd)  
      }
      lines(x, a, lty = lty, lwd = lwd, col = col) 
      return(invisible(list(x = x, y = a)))
    }
    a
}

"abvneglog"<- 
function(x = 0.5, dep, plot = FALSE, add = FALSE, lty = 1, lwd = 1, col = 1,
         blty = 3, blwd = 1, xlim = c(0,1), ylim = c(0.5,1), xlab = "",
         ylab = "", ...)
{ 
    if(length(dep) != 1 || mode(dep) != "numeric" || dep <= 0)
        stop("invalid argument for `dep'")
    if(plot || add) x <- seq(0, 1, length = 100)
    a <- 1 - (x^(-dep) + (1-x)^(-dep))^(-1/dep)
    if(plot || add) {
      if(!add)  { 
        plot(x, a, type="n", xlab = xlab, ylab = ylab,
             xlim = xlim, ylim = ylim, ...) 
        polygon(c(0, 0.5, 1), c(1, 0.5, 1), lty = blty, lwd = blwd)  
      }
      lines(x, a, lty = lty, lwd = lwd, col = col) 
      return(invisible(list(x = x, y = a)))
    }
    a
}

"abvaneglog"<- 
function(x = 0.5, dep, asy = c(1,1), plot = FALSE, add = FALSE,
         lty = 1, lwd = 1, col = 1, blty = 3, blwd = 1, xlim = c(0,1),
         ylim = c(0.5,1), xlab = "", ylab = "", ...)
{
    if(length(dep) != 1 || mode(dep) != "numeric" || dep <= 0)
        stop("invalid argument for `dep'")
    if(plot || add) x <- seq(0, 1, length = 100)
    if(length(asy) != 2 || mode(asy) != "numeric" || min(asy) < 0 ||
       max(asy) > 1) stop("invalid argument for `asy'")
    a <- 1 - ((asy[1]*x)^(-dep) + (asy[2]*(1-x))^(-dep))^(-1/dep)
    if(plot || add) {
      if(!add)  { 
        plot(x, a, type="n", xlab = xlab, ylab = ylab,
             xlim = xlim, ylim = ylim, ...) 
        polygon(c(0, 0.5, 1), c(1, 0.5, 1), lty = blty, lwd = blwd)  
      }
      lines(x, a, lty = lty, lwd = lwd, col = col) 
      return(invisible(list(x = x, y = a)))
    }
    a
}

"abvbilog"<- 
function(x = 0.5, alpha, beta, plot = FALSE, add = FALSE,
         lty = 1, lwd = 1, col = 1, blty = 3, blwd = 1, xlim = c(0,1),
         ylim = c(0.5,1), xlab = "", ylab = "", ...)
{
    if(length(alpha) != 1 || mode(alpha) != "numeric")
        stop("invalid argument for `alpha'")
    if(length(beta) != 1 || mode(beta) != "numeric")
        stop("invalid argument for `beta'")
    if(any(c(alpha,beta) <= 0) || any(c(alpha,beta) >= 1))
        stop("`alpha' and `beta' must be in the open interval (0,1)")
    if(plot || add) x <- seq(0, 1, length = 100)
    gma <- numeric(length(x))
    for(i in 1:length(x)) {
        gmafn <- function(z)
            (1-alpha) * x[i] * (1-z)^beta - (1-beta) * (1-x[i]) * z^alpha
        if(is.na(x[i])) gma[i] <- NA
        else if(x[i] == 0) gma[i] <- 0
        else if(x[i] == 1) gma[i] <- 1
        else gma[i] <- uniroot(gmafn, lower = 0, upper = 1,
                               tol = .Machine$double.eps^0.5)$root
    }
    a <- x * gma^(1-alpha) + (1-x) * (1 - gma)^(1-beta)
    if(plot || add) {
      if(!add)  { 
        plot(x, a, type="n", xlab = xlab, ylab = ylab,
             xlim = xlim, ylim = ylim, ...) 
        polygon(c(0, 0.5, 1), c(1, 0.5, 1), lty = blty, lwd = blwd)  
      }
      lines(x, a, lty = lty, lwd = lwd, col = col) 
      return(invisible(list(x = x, y = a)))
    }
    a
}

"abvnegbilog"<- 
function(x = 0.5, alpha, beta, plot = FALSE, add = FALSE,
         lty = 1, lwd = 1, col = 1, blty = 3, blwd = 1, xlim = c(0,1),
         ylim = c(0.5,1), xlab = "", ylab = "", ...)
{
    if(length(alpha) != 1 || mode(alpha) != "numeric")
        stop("invalid argument for `alpha'")
    if(length(beta) != 1 || mode(beta) != "numeric")
        stop("invalid argument for `beta'")
    if(any(c(alpha,beta) <= 0))
        stop("`alpha' and `beta' must be non-negative")
    if(plot || add) x <- seq(0, 1, length = 100)
    gma <- numeric(length(x))
    for(i in 1:length(x)) {
        gmafn <- function(z)
            (1+alpha) * x[i] * z^alpha - (1+beta) * (1-x[i]) * (1-z)^beta
        if(is.na(x[i])) gma[i] <- NA
        else if(x[i] == 0) gma[i] <- 1
        else if(x[i] == 1) gma[i] <- 0
        else gma[i] <- uniroot(gmafn, lower = 0, upper = 1,
                               tol = .Machine$double.eps^0.5)$root
    }
    a <- 1 - x * gma^(1+alpha) - (1-x) * (1 - gma)^(1+beta)
    if(plot || add) {
      if(!add)  { 
        plot(x, a, type="n", xlab = xlab, ylab = ylab,
             xlim = xlim, ylim = ylim, ...) 
        polygon(c(0, 0.5, 1), c(1, 0.5, 1), lty = blty, lwd = blwd)  
      }
      lines(x, a, lty = lty, lwd = lwd, col = col) 
      return(invisible(list(x = x, y = a)))
    }
    a
}

"abvct" <-
function(x = 0.5, alpha, beta, plot = FALSE, add = FALSE,
         lty = 1, lwd = 1, col = 1, blty = 3, blwd = 1, xlim = c(0,1),
         ylim = c(0.5,1), xlab = "", ylab = "", ...)
{
    if(length(alpha) != 1 || mode(alpha) != "numeric")
        stop("invalid argument for `alpha'")
    if(length(beta) != 1 || mode(beta) != "numeric")
        stop("invalid argument for `beta'")
    if(any(c(alpha,beta) <= 0))
        stop("`alpha' and `beta' must be non-negative")
    if(plot || add) x <- seq(0, 1, length = 100)
    u <- (alpha * (1-x)) / (alpha * (1-x) + beta * x)
    a <- (1-x) * pbeta(u, shape1 = alpha, shape2 = beta + 1) +
      x * pbeta(u, shape1 = alpha + 1, shape2 = beta, lower.tail = FALSE)
    if(plot || add) {
      if(!add)  { 
        plot(x, a, type="n", xlab = xlab, ylab = ylab,
             xlim = xlim, ylim = ylim, ...) 
        polygon(c(0, 0.5, 1), c(1, 0.5, 1), lty = blty, lwd = blwd)  
      }
      lines(x, a, lty = lty, lwd = lwd, col = col) 
      return(invisible(list(x = x, y = a)))
    }
    a
}

"abvamix" <-
function(x = 0.5, alpha, beta, plot = FALSE, add = FALSE,
         lty = 1, lwd = 1, col = 1, blty = 3, blwd = 1, xlim = c(0,1),
         ylim = c(0.5,1), xlab = "", ylab = "", ...)
{
    if(length(alpha) != 1 || mode(alpha) != "numeric")
        stop("invalid argument for `alpha'")
    if(length(beta) != 1 || mode(beta) != "numeric")
        stop("invalid argument for `beta'")
    if(alpha < 0)
        stop("`alpha' must be non-negative")
    if((alpha + beta) > 1)
        stop("`alpha' + `beta' cannot be greater than one")
    if((alpha + 2*beta) > 1)
        stop("`alpha' + `2*beta' cannot be greater than one")
    if((alpha + 3*beta) < 0)
        stop("`alpha' + `3*beta' must be non-negative")
    if(plot || add) x <- seq(0, 1, length = 100)
    a <- 1 - (alpha + beta) * x + alpha * (x^2) + beta * (x^3)
    if(plot || add) {
      if(!add)  { 
        plot(x, a, type="n", xlab = xlab, ylab = ylab,
             xlim = xlim, ylim = ylim, ...) 
        polygon(c(0, 0.5, 1), c(1, 0.5, 1), lty = blty, lwd = blwd)  
      }
      lines(x, a, lty = lty, lwd = lwd, col = col) 
      return(invisible(list(x = x, y = a)))
    }
    a
}

"hbvevd" <-
function(x = 0.5, dep, asy = c(1,1), alpha, beta, model = c("log", "alog",
    "hr", "neglog", "aneglog", "bilog", "negbilog", "ct", "amix"),
     half = FALSE, plot = FALSE, add = FALSE, lty = 1, ...)
{
  if(any(x < 0, na.rm = TRUE) || any(x > 1, na.rm = TRUE))
    stop("invalid argument for `x'")
  model <- match.arg(model)
  m1 <- c("bilog", "negbilog", "ct", "amix")
  m2 <- c(m1, "log", "hr", "neglog")
  m3 <- c("log", "alog", "hr", "neglog", "aneglog")
  if((model %in% m1) && !missing(dep))
    warning("ignoring `dep' argument")
  if((model %in% m2) && !missing(asy))
    warning("ignoring `asy' argument")
  if((model %in% m3) && !missing(alpha))
    warning("ignoring `alpha' argument")
  if((model %in% m3) && !missing(beta))
    warning("ignoring `beta' argument")
     
  switch(model,
    log = hbvlog(x = x, dep = dep, plot = plot, add = add, half = half,
      lty = lty, ...),
    alog = hbvalog(x = x, dep = dep, asy = asy, plot = plot, add = add,
      half = half, lty = lty, ...),
    hr = hbvhr(x = x, dep = dep, plot = plot, add = add, half = half,
      lty = lty, ...),
    neglog = hbvneglog(x = x, dep = dep, plot = plot, add = add,
      half = half, lty = lty, ...),
    aneglog = hbvaneglog(x = x, dep = dep, asy = asy, plot = plot, add = add,
      half = half, lty = lty, ...),
    bilog = hbvbilog(x = x, alpha = alpha, beta = beta, plot = plot,
      add = add, half = half, lty = lty, ...),
    negbilog = hbvnegbilog(x = x, alpha = alpha, beta = beta, plot = plot,
      add = add, half = half, lty = lty, ...),
    ct = hbvct(x = x, alpha = alpha, beta = beta, plot = plot, add = add,
      half = half, lty = lty, ...),
    amix = hbvamix(x = x, alpha = alpha, beta = beta, plot = plot, add = add,
      half = half, lty = lty, ...)) 
}

"hbvlog"<- 
function(x = 0.5, dep, plot = FALSE, add = FALSE, half = FALSE, lty = 1, xlab = "t",
         ylab = "h(t)", xlim = c(0,1), ylim = c(0, max(h, na.rm = TRUE)), ...)
{
    if(length(dep) != 1 || mode(dep) != "numeric" || dep <= 0 ||
        dep > 1) stop("invalid argument for `dep'")
    if(plot || add) x <- seq(0, 1, length = 100)
    idep <- 1/dep
    h <- (idep - 1) * (x * (1-x))^(-1-idep) * (x^(-idep) + (1-x)^(-idep))^(dep-2)
    if(half) h <- h/2
    if(plot || add) {
      if(!add) {
        plot(x, h, type = "l", xlab = xlab, ylab = ylab, xlim = xlim,
          ylim = ylim, lty = lty, ...) 
      }
      lines(x, h, lty = lty) 
      return(invisible(list(x = x, y = h)))
    }
    h
}

"hbvalog"<- 
function(x = 0.5, dep, asy = c(1,1), plot = FALSE, add = FALSE,
         half = FALSE, lty = 1, xlab = "t", ylab = "h(t)", xlim = c(0,1),
         ylim = c(0, max(h, na.rm = TRUE)), ...)
{
    if(length(dep) != 1 || mode(dep) != "numeric" || dep <= 0 ||
        dep > 1) stop("invalid argument for `dep'")
    if(length(asy) != 2 || mode(asy) != "numeric" || min(asy) < 0 ||
       max(asy) > 1) stop("invalid argument for `asy'")
    if(plot || add) x <- seq(0, 1, length = 100)
    idep <- 1/dep
    h <- (idep - 1) * (asy[1] * asy[2])^idep * (x * (1-x))^(-1-idep) *
      ((asy[1]/x)^idep + (asy[2]/(1-x))^idep)^(dep-2)
    if(half) h <- h/2
    if(plot || add) {
      if(!add) {
        plot(x, h, type = "l", xlab = xlab, ylab = ylab, xlim = xlim,
          ylim = ylim, lty = lty, ...) 
      }
      lines(x, h, lty = lty) 
      return(invisible(list(x = x, y = h)))
    }
    h
}

"hbvhr" <-
function(x = 0.5, dep, plot = FALSE, add = FALSE, half = FALSE, lty = 1, xlab = "t",
         ylab = "h(t)", xlim = c(0,1), ylim = c(0, max(h, na.rm = TRUE)), ...)
{
    if(length(dep) != 1 || mode(dep) != "numeric" || dep <= 0)
        stop("invalid argument for `dep'")
    if(plot || add) x <- seq(0, 1, length = 100)
    h <- dep * dnorm(1/dep + dep * log(x/(1-x)) / 2)
    h <- h / (2 * x * (1-x)^2)
    if(half) h <- h/2
    if(plot || add) {
      if(!add) {
        plot(x, h, type = "l", xlab = xlab, ylab = ylab, xlim = xlim,
          ylim = ylim, lty = lty, ...) 
      }
      lines(x, h, lty = lty) 
      return(invisible(list(x = x, y = h)))
    }
    h
}

"hbvneglog"<- 
function(x = 0.5, dep, plot = FALSE, add = FALSE, half = FALSE, lty = 1, xlab = "t",
         ylab = "h(t)", xlim = c(0,1), ylim = c(0, max(h, na.rm = TRUE)), ...)
{ 
    if(length(dep) != 1 || mode(dep) != "numeric" || dep <= 0)
        stop("invalid argument for `dep'")
    if(plot || add) x <- seq(0, 1, length = 100)
    h <- (1 + dep) * (x * (1-x))^(dep-1) * (x^dep + (1-x)^dep)^(-1/dep-2)
    if(half) h <- h/2
    if(plot || add) {
      if(!add) {
        plot(x, h, type = "l", xlab = xlab, ylab = ylab, xlim = xlim,
          ylim = ylim, lty = lty, ...) 
      }
      lines(x, h, lty = lty) 
      return(invisible(list(x = x, y = h)))
    }
    h
}

"hbvaneglog"<- 
function(x = 0.5, dep, asy = c(1,1), plot = FALSE, add = FALSE,
         half = FALSE, lty = 1, xlab = "t", ylab = "h(t)", xlim = c(0,1),
         ylim = c(0, max(h, na.rm = TRUE)), ...)
{
    if(length(dep) != 1 || mode(dep) != "numeric" || dep <= 0)
        stop("invalid argument for `dep'")
    if(plot || add) x <- seq(0, 1, length = 100)
    if(length(asy) != 2 || mode(asy) != "numeric" || min(asy) < 0 ||
       max(asy) > 1) stop("invalid argument for `asy'")
    h <- (1 + dep) * (asy[1] * asy[2])^(-dep) * (x * (1-x))^(dep-1) *
      ((x/asy[1])^dep + ((1-x)/asy[2])^dep)^(-1/dep-2)
    if(half) h <- h/2
    if(plot || add) {
      if(!add) {
        plot(x, h, type = "l", xlab = xlab, ylab = ylab, xlim = xlim,
          ylim = ylim, lty = lty, ...) 
      }
      lines(x, h, lty = lty) 
      return(invisible(list(x = x, y = h)))
    }
    h
}

"hbvbilog"<- 
function(x = 0.5, alpha, beta, plot = FALSE, add = FALSE,
         half = FALSE, lty = 1, xlab = "t", ylab = "h(t)", xlim = c(0,1),
         ylim = c(0, max(h, na.rm = TRUE)), ...)
{
    if(length(alpha) != 1 || mode(alpha) != "numeric")
        stop("invalid argument for `alpha'")
    if(length(beta) != 1 || mode(beta) != "numeric")
        stop("invalid argument for `beta'")
    if(any(c(alpha,beta) <= 0) || any(c(alpha,beta) >= 1))
        stop("`alpha' and `beta' must be in the open interval (0,1)")
    if(plot || add) x <- seq(0, 1, length = 100)
    gma <- numeric(length(x))
    for(i in 1:length(x)) {
        gmafn <- function(z)
            (1-alpha) * (1-x[i]) * (1-z)^beta - (1-beta) * x[i] * z^alpha
        if(is.na(x[i])) gma[i] <- NA
        else if(x[i] == 0) gma[i] <- 0
        else if(x[i] == 1) gma[i] <- 1
        else gma[i] <- uniroot(gmafn, lower = 0, upper = 1,
                               tol = .Machine$double.eps^0.5)$root
    }
    a <- x * gma^(1-alpha) + (1-x) * (1 - gma)^(1-beta)
    h <- exp(log(1-alpha) + log(beta) + (beta - 1)*log(1-gma) + log(1-x)) +
      exp(log(1-beta) + log(alpha) + (alpha - 1)*log(gma) + log(x))
    h <- exp(log(1-alpha) + log(1-beta) - log(x * (1-x)) - log(h)) 
    if(half) h <- h/2
    if(plot || add) {
      if(!add) {
        plot(x, h, type = "l", xlab = xlab, ylab = ylab, xlim = xlim,
          ylim = ylim, lty = lty, ...) 
      }
      lines(x, h, lty = lty) 
      return(invisible(list(x = x, y = h)))
    }
    h
}

"hbvnegbilog"<- 
function(x = 0.5, alpha, beta, plot = FALSE, add = FALSE,
         half = FALSE, lty = 1, xlab = "t", ylab = "h(t)", xlim = c(0,1),
         ylim = c(0, max(h, na.rm = TRUE)), ...)
{
    if(length(alpha) != 1 || mode(alpha) != "numeric")
        stop("invalid argument for `alpha'")
    if(length(beta) != 1 || mode(beta) != "numeric")
        stop("invalid argument for `beta'")
    if(any(c(alpha,beta) <= 0))
        stop("`alpha' and `beta' must be non-negative")
    if(plot || add) x <- seq(0, 1, length = 100)
    gma <- numeric(length(x))
    for(i in 1:length(x)) {
        gmafn <- function(z)
            (1+alpha) * (1-x[i]) * z^alpha - (1+beta) * x[i] * (1-z)^beta
        if(is.na(x[i])) gma[i] <- NA
        else if(x[i] == 0) gma[i] <- 1
        else if(x[i] == 1) gma[i] <- 0
        else gma[i] <- uniroot(gmafn, lower = 0, upper = 1,
                               tol = .Machine$double.eps^0.5)$root
    }
    h <- exp(log(1+alpha) + log(alpha) + (alpha - 1)*log(gma) + log(1-x)) +
      exp(log(1+beta) + log(beta) + (beta - 1)*log(1-gma) + log(x))
    h <- exp(log(1+alpha) + log(1+beta) + alpha * log(gma) + beta * log(1-gma) -
             log(x * (1-x)) - log(h))
    if(half) h <- h/2
    if(plot || add) {
      if(!add) {
        plot(x, h, type = "l", xlab = xlab, ylab = ylab, xlim = xlim,
          ylim = ylim, lty = lty, ...) 
      }
      lines(x, h, lty = lty) 
      return(invisible(list(x = x, y = h)))
    }
    h
}

"hbvct" <-
function(x = 0.5, alpha, beta, plot = FALSE, add = FALSE, half = FALSE, lty = 1,
         xlab = "t", ylab = "h(t)", xlim = c(0,1),
         ylim = c(0, max(h, na.rm = TRUE)), ...)
{
    if(length(alpha) != 1 || mode(alpha) != "numeric")
        stop("invalid argument for `alpha'")
    if(length(beta) != 1 || mode(beta) != "numeric")
        stop("invalid argument for `beta'")
    if(any(c(alpha,beta) <= 0))
        stop("`alpha' and `beta' must be non-negative")
    if(plot || add) x <- seq(0, 1, length = 100)
    u <- (alpha * x) / (alpha * x + beta * (1-x))
    c1 <- alpha * beta / (alpha + beta + 1)
    h <- dbeta(u, shape1 = alpha + 1, shape2 = beta + 1) /
          (alpha * x^2 * (1-x) + beta * x * (1-x)^2) * c1
    if(half) h <- h/2
    if(plot || add) {
      if(!add) {
        plot(x, h, type = "l", xlab = xlab, ylab = ylab, xlim = xlim,
          ylim = ylim, lty = lty, ...) 
      }
      lines(x, h, lty = lty) 
      return(invisible(list(x = x, y = h)))
    }
    h
}

"hbvamix" <-
function(x = 0.5, alpha, beta, plot = FALSE, add = FALSE, half = FALSE, lty = 1,
         xlab = "t", ylab = "h(t)", xlim = c(0,1),
         ylim = c(0, max(h, na.rm = TRUE)), ...)
{
    if(length(alpha) != 1 || mode(alpha) != "numeric")
        stop("invalid argument for `alpha'")
    if(length(beta) != 1 || mode(beta) != "numeric")
        stop("invalid argument for `beta'")
    if(alpha < 0)
        stop("`alpha' must be non-negative")
    if((alpha + beta) > 1)
        stop("`alpha' + `beta' cannot be greater than one")
    if((alpha + 2*beta) > 1)
        stop("`alpha' + `2*beta' cannot be greater than one")
    if((alpha + 3*beta) < 0)
        stop("`alpha' + `3*beta' must be non-negative")
    if(plot || add) x <- seq(0, 1, length = 100)
    h <- 2 * alpha + 6 * beta * (1-x)
    if(half) h <- h/2
    if(plot || add) {
      if(!add) {
        plot(x, h, type = "l", xlab = xlab, ylab = ylab, xlim = xlim,
          ylim = ylim, lty = lty, ...) 
      }
      lines(x, h, lty = lty) 
      return(invisible(list(x = x, y = h)))
    }
    h
}

"dbvevd" <-
function(x, dep, asy = c(1,1), alpha, beta, model = c("log", "alog",
    "hr", "neglog", "aneglog", "bilog", "negbilog", "ct", "amix"),
    mar1 = c(0,1,0), mar2 = mar1, log = FALSE)
{
  model <- match.arg(model)
  m1 <- c("bilog", "negbilog", "ct", "amix")
  m2 <- c(m1, "log", "hr", "neglog")
  m3 <- c("log", "alog", "hr", "neglog", "aneglog")
  if((model %in% m1) && !missing(dep))
    warning("ignoring `dep' argument")
  if((model %in% m2) && !missing(asy))
    warning("ignoring `asy' argument")
  if((model %in% m3) && !missing(alpha))
    warning("ignoring `alpha' argument")
  if((model %in% m3) && !missing(beta))
    warning("ignoring `beta' argument")
    
  switch(model,
    log = dbvlog(x = x, dep = dep, mar1 = mar1, mar2 = mar2, log = log),
    alog = dbvalog(x = x, dep = dep, asy = asy, mar1 = mar1,
      mar2 = mar2, log = log),
    hr = dbvhr(x = x, dep = dep, mar1 = mar1, mar2 = mar2, log = log),
    neglog = dbvneglog(x = x, dep = dep, mar1 = mar1, mar2 = mar2, log = log),
    aneglog = dbvaneglog(x = x, dep = dep, asy = asy, mar1 = mar1,
      mar2 = mar2, log = log),
    bilog = dbvbilog(x = x, alpha = alpha, beta = beta, mar1 = mar1,
      mar2 = mar2, log = log),
    negbilog = dbvnegbilog(x = x, alpha = alpha, beta = beta, mar1 = mar1,
      mar2 = mar2, log = log),
    ct = dbvct(x = x, alpha = alpha, beta = beta, mar1 = mar1,
      mar2 = mar2, log = log),
    amix = dbvamix(x = x, alpha = alpha, beta = beta, mar1 = mar1,
      mar2 = mar2, log = log)) 
}

"dbvlog"<- 
function(x, dep, mar1 = c(0,1,0), mar2 = mar1, log = FALSE)
{
    if(length(dep) != 1 || mode(dep) != "numeric" || dep <= 0 ||
        dep > 1) stop("invalid argument for `dep'")
    if(is.null(dim(x))) dim(x) <- c(1,2)
    mar1 <- matrix(t(mar1), nrow = nrow(x), ncol = 3, byrow = TRUE)
    mar2 <- matrix(t(mar2), nrow = nrow(x), ncol = 3, byrow = TRUE)
    d <- numeric(nrow(x))
    x <- mtransform(x, list(mar1, mar2))
    ext <- apply(x,1,function(z) any(z %in% c(0,Inf)))
    d[ext] <- -Inf    
    if(any(!ext)) {
        x <- x[!ext, ,drop=FALSE]
        mar1 <- mar1[!ext, ,drop=FALSE]
        mar2 <- mar2[!ext, ,drop=FALSE]
        idep <- 1/dep
        z <- apply(x^idep,1,sum)^dep
        lx <- log(x)
        .expr1 <- (idep+mar1[,3])*lx[,1] + (idep+mar2[,3])*lx[,2] -
            log(mar1[,2]*mar2[,2])
        d[!ext] <- .expr1 + (1-2*idep)*log(z) + log(idep-1+z) - z
    }
    if(!log) d <- exp(d)
    d
}

"dbvalog"<- 
function(x, dep, asy = c(1,1), mar1 = c(0,1,0), mar2 = mar1, log = FALSE)
{
    if(length(dep) != 1 || mode(dep) != "numeric" || dep <= 0 ||
        dep > 1) stop("invalid argument for `dep'")
    if(length(asy) != 2 || mode(asy) != "numeric" || min(asy) < 0 ||
       max(asy) > 1) stop("invalid argument for `asy'")
    if(is.null(dim(x))) dim(x) <- c(1,2)
    mar1 <- matrix(t(mar1), nrow = nrow(x), ncol = 3, byrow = TRUE)
    mar2 <- matrix(t(mar2), nrow = nrow(x), ncol = 3, byrow = TRUE)
    d <- numeric(nrow(x))
    x <- mtransform(x, list(mar1, mar2))
    ext <- apply(x,1,function(z) any(z %in% c(0,Inf)))
    d[ext] <- -Inf
    if(any(!ext)) {
        x <- x[!ext, ,drop=FALSE]
        mar1 <- mar1[!ext, ,drop=FALSE]
        mar2 <- mar2[!ext, ,drop=FALSE]
        asy <- matrix(asy, ncol = 2, nrow = nrow(x), byrow = TRUE)
        idep <- 1/dep
        z <- apply((asy*x)^idep,1,sum)^dep
        v <- z + apply((1-asy)*x,1,sum)
        f1asy <- (idep)*log(asy)
        f2asy <- log(1-asy)
        lx <- log(x)
        fx <- (idep-1)*lx
        jac <- (1+mar1[,3])*lx[,1] + (1+mar2[,3])*lx[,2] -
            log(mar1[,2]*mar2[,2])
        .expr1 <- apply(f2asy,1,sum)
        .expr2 <- f2asy[,1] + f1asy[,2] + fx[,2]
        .expr3 <- f2asy[,2] + f1asy[,1] + fx[,1]
        .expr4 <- (1-idep)*log(z) + log(exp(.expr2)+exp(.expr3))
        .expr5 <- apply(cbind(f1asy,fx),1,sum) + (1-2*idep)*log(z) +
            log(idep-1+z)
        d[!ext] <- log(exp(.expr1)+exp(.expr4)+exp(.expr5))-v+jac
    }
    if(!log) d <- exp(d)
    d
}

"dbvhr" <-
function(x, dep, mar1 = c(0,1,0), mar2 = mar1, log = FALSE)
{
    if(length(dep) != 1 || mode(dep) != "numeric" || dep <= 0)
        stop("invalid argument for `dep'")
    if(is.null(dim(x))) dim(x) <- c(1,2)
    mar1 <- matrix(t(mar1), nrow = nrow(x), ncol = 3, byrow = TRUE)
    mar2 <- matrix(t(mar2), nrow = nrow(x), ncol = 3, byrow = TRUE)
    d <- numeric(nrow(x))
    x <- mtransform(x, list(mar1, mar2))
    ext <- apply(x,1,function(z) any(z %in% c(0,Inf)))
    d[ext] <- -Inf
    if(any(!ext)) {
        x <- x[!ext, ,drop=FALSE]
        mar1 <- mar1[!ext, ,drop=FALSE]
        mar2 <- mar2[!ext, ,drop=FALSE]
        fn <- function(x1, x2, nm = pnorm) x1 *
            nm(1/dep + dep * log(x1/x2) / 2)
        v <- fn(x[,1], x[,2]) + fn(x[,2], x[,1])
        .expr1 <- fn(x[,1], x[,2]) * fn(x[,2], x[,1]) +
            dep * fn(x[,1], x[,2], nm = dnorm) / 2
        lx <- log(x)
        jac <- mar1[,3]*lx[,1] + mar2[,3]*lx[,2] - log(mar1[,2]*mar2[,2])
        d[!ext] <- log(.expr1)+jac-v
    }
    if(!log) d <- exp(d)
    d    
}

"dbvneglog"<- 
function(x, dep, mar1 = c(0,1,0), mar2 = mar1, log = FALSE)
{
    if(length(dep) != 1 || mode(dep) != "numeric" || dep <= 0)
        stop("invalid argument for `dep'")
    if(is.null(dim(x))) dim(x) <- c(1,2)
    mar1 <- matrix(t(mar1), nrow = nrow(x), ncol = 3, byrow = TRUE)
    mar2 <- matrix(t(mar2), nrow = nrow(x), ncol = 3, byrow = TRUE)
    d <- numeric(nrow(x))
    x <- mtransform(x, list(mar1, mar2))
    ext <- apply(x,1,function(z) any(z %in% c(0,Inf)))
    d[ext] <- -Inf
    if(any(!ext)) {
        x <- x[!ext, ,drop=FALSE]
        mar1 <- mar1[!ext, ,drop=FALSE]
        mar2 <- mar2[!ext, ,drop=FALSE]
        idep <- 1/dep
        z <- apply(x^(-dep),1,sum)^(-idep)
        v <- apply(x,1,sum) - z
        lx <- log(x)
        fx <- (-dep-1)*lx
        jac <- (1+mar1[,3])*lx[,1] + (1+mar2[,3])*lx[,2] -
            log(mar1[,2]*mar2[,2])
        .expr1 <- (1+dep)*log(z) + log(exp(fx[,1])+exp(fx[,2]))
        .expr2 <- fx[,1] + fx[,2] + (1+2*dep)*log(z) + log(1+dep+z)
        d[!ext] <- log(1-exp(.expr1)+exp(.expr2))-v+jac
    }
    if(!log) d <- exp(d)
    d    
}

"dbvaneglog"<- 
function(x, dep, asy = c(1,1), mar1 = c(0,1,0), mar2 = mar1, log = FALSE)
{
    if(length(dep) != 1 || mode(dep) != "numeric" || dep <= 0)
        stop("invalid argument for `dep'")
    if(length(asy) != 2 || mode(asy) != "numeric" || min(asy) < 0 ||
       max(asy) > 1) stop("invalid argument for `asy'")
    if(is.null(dim(x))) dim(x) <- c(1,2)
    mar1 <- matrix(t(mar1), nrow = nrow(x), ncol = 3, byrow = TRUE)
    mar2 <- matrix(t(mar2), nrow = nrow(x), ncol = 3, byrow = TRUE)
    d <- numeric(nrow(x))
    x <- mtransform(x, list(mar1, mar2))
    ext <- apply(x,1,function(z) any(z %in% c(0,Inf)))
    d[ext] <- -Inf
    if(any(!ext)) {
        x <- x[!ext, ,drop=FALSE]
        mar1 <- mar1[!ext, ,drop=FALSE]
        mar2 <- mar2[!ext, ,drop=FALSE]
        asy <- matrix(asy, ncol = 2, nrow = nrow(x), byrow = TRUE)
        idep <- 1/dep
        z <- apply((asy*x)^(-dep),1,sum)^(-idep)
        v <- apply(x,1,sum) - z
        fasy <- (-dep)*log(asy)
        lx <- log(x)
        fx <- (-dep-1)*lx
        jac <- (1+mar1[,3])*lx[,1] + (1+mar2[,3])*lx[,2] -
            log(mar1[,2]*mar2[,2])
        .expr1 <- fasy[,1] + fx[,1]
        .expr2 <- fasy[,2] + fx[,2]
        .expr3 <- (1+dep)*log(z) + log(exp(.expr1)+exp(.expr2))
        .expr4 <- apply(cbind(fasy,fx),1,sum) + (1+2*dep)*log(z) + log(1+dep+z)
        d[!ext] <- log(1-exp(.expr3)+exp(.expr4))-v+jac
    }
    if(!log) d <- exp(d)
    d
}

"dbvbilog"<- 
function(x, alpha, beta, mar1 = c(0,1,0), mar2 = mar1, log = FALSE)
{
    if(length(alpha) != 1 || mode(alpha) != "numeric")
        stop("invalid argument for `alpha'")
    if(length(beta) != 1 || mode(beta) != "numeric")
        stop("invalid argument for `beta'")
    if(any(c(alpha,beta) <= 0) || any(c(alpha,beta) >= 1))
        stop("`alpha' and `beta' must be in the open interval (0,1)")
    if(is.null(dim(x))) dim(x) <- c(1,2)
    mar1 <- matrix(t(mar1), nrow = nrow(x), ncol = 3, byrow = TRUE)
    mar2 <- matrix(t(mar2), nrow = nrow(x), ncol = 3, byrow = TRUE)
    d <- numeric(nrow(x))
    x <- mtransform(x, list(mar1, mar2))
    ext <- apply(x,1,function(z) any(z %in% c(0,Inf)))
    d[ext] <- -Inf
    if(any(!ext)) {
        x <- x[!ext, ,drop=FALSE]
        mar1 <- mar1[!ext, ,drop=FALSE]
        mar2 <- mar2[!ext, ,drop=FALSE]
        gma <- numeric(nrow(x))
        for(i in 1:nrow(x)) {
            gmafn <- function(z)
                (1-alpha) * x[i,1] * (1-z)^beta - (1-beta) * x[i,2] * z^alpha
            if(any(is.na(x[i,]))) gma[i] <- NA
            else gma[i] <- uniroot(gmafn, lower = 0, upper = 1,
                               tol = .Machine$double.eps^0.5)$root
        }
        v <- x[,1] * gma^(1-alpha) + x[,2] * (1 - gma)^(1-beta)
        lx <- log(x)
        jac <- (1+mar1[,3])*lx[,1] + (1+mar2[,3])*lx[,2] -
            log(mar1[,2]*mar2[,2])
        .expr1 <- exp((1-alpha)*log(gma) + (1-beta)*log(1-gma))
        .expr2 <- exp(log(1-alpha) + log(beta) + (beta - 1)*log(1-gma) +
            lx[,1]) + exp(log(1-beta) + log(alpha) + (alpha - 1)*log(gma) +
            lx[,2])
        d[!ext] <- log(.expr1 + (1-alpha)*(1-beta)/.expr2) - v + jac
    }
    if(!log) d <- exp(d)
    d   
}

"dbvnegbilog"<- 
function(x, alpha, beta, mar1 = c(0,1,0), mar2 = mar1, log = FALSE)
{
    if(length(alpha) != 1 || mode(alpha) != "numeric")
        stop("invalid argument for `alpha'")
    if(length(beta) != 1 || mode(beta) != "numeric")
        stop("invalid argument for `beta'")
    if(any(c(alpha,beta) <= 0))
        stop("`alpha' and `beta' must be non-negative")
    if(is.null(dim(x))) dim(x) <- c(1,2)
    mar1 <- matrix(t(mar1), nrow = nrow(x), ncol = 3, byrow = TRUE)
    mar2 <- matrix(t(mar2), nrow = nrow(x), ncol = 3, byrow = TRUE)
    d <- numeric(nrow(x))
    x <- mtransform(x, list(mar1, mar2))
    ext <- apply(x,1,function(z) any(z %in% c(0,Inf)))
    d[ext] <- -Inf
    if(any(!ext)) {
        x <- x[!ext, ,drop=FALSE]
        mar1 <- mar1[!ext, ,drop=FALSE]
        mar2 <- mar2[!ext, ,drop=FALSE]
        gma <- numeric(nrow(x))
        for(i in 1:nrow(x)) {
            gmafn <- function(z)
                (1+alpha) * x[i,1] * z^alpha - (1+beta) * x[i,2] * (1-z)^beta
            if(any(is.na(x[i,]))) gma[i] <- NA
            else gma[i] <- uniroot(gmafn, lower = 0, upper = 1,
                               tol = .Machine$double.eps^0.5)$root
        }
        v <- x[,1] + x[,2] - x[,1] * gma^(1+alpha) - x[,2] * (1 - gma)^(1+beta)
        lx <- log(x)
        jac <- (1+mar1[,3])*lx[,1] + (1+mar2[,3])*lx[,2] -
            log(mar1[,2]*mar2[,2])
        .expr1 <- (1-gma^(1+alpha)) * (1 - (1-gma)^(1+beta))
        .expr2 <- exp(log(1+alpha) + log(1+beta) + alpha*log(gma) +
                      beta*log(1-gma))
        .expr3 <- exp(log(1+alpha) + log(alpha) + (alpha - 1)*log(gma) +
            lx[,1]) + exp(log(1+beta) + log(beta) + (beta - 1)*log(1-gma) +
            lx[,2])
        d[!ext] <- log(.expr1 + .expr2/.expr3) - v + jac
    }
    if(!log) d <- exp(d)
    d   
}

"dbvct"<- 
function(x, alpha, beta, mar1 = c(0,1,0), mar2 = mar1, log = FALSE)
{
    if(length(alpha) != 1 || mode(alpha) != "numeric")
        stop("invalid argument for `alpha'")
    if(length(beta) != 1 || mode(beta) != "numeric")
        stop("invalid argument for `beta'")
    if(any(c(alpha,beta) <= 0))
        stop("`alpha' and `beta' must be non-negative")
    if(is.null(dim(x))) dim(x) <- c(1,2)
    mar1 <- matrix(t(mar1), nrow = nrow(x), ncol = 3, byrow = TRUE)
    mar2 <- matrix(t(mar2), nrow = nrow(x), ncol = 3, byrow = TRUE)
    d <- numeric(nrow(x))
    x <- mtransform(x, list(mar1, mar2))
    ext <- apply(x,1,function(z) any(z %in% c(0,Inf)))
    d[ext] <- -Inf
    if(any(!ext)) {
        x <- x[!ext, ,drop=FALSE]
        mar1 <- mar1[!ext, ,drop=FALSE]
        mar2 <- mar2[!ext, ,drop=FALSE]
        u <- (alpha * x[,2]) / (alpha * x[,2] + beta * x[,1]) 
        v <- x[,2] * pbeta(u, shape1 = alpha, shape2 = beta + 1) + x[,1] *
          pbeta(u, shape1 = alpha + 1, shape2 = beta, lower.tail = FALSE)
        lx <- log(x)
        jac <- (1+mar1[,3])*lx[,1] + (1+mar2[,3])*lx[,2] -
            log(mar1[,2]*mar2[,2])
        .c1 <- alpha * beta / (alpha + beta + 1)
        .expr1 <- pbeta(u, shape1 = alpha, shape2 = beta + 1) *
          pbeta(u, shape1 = alpha + 1, shape2 = beta, lower.tail = FALSE)
        .expr2 <- dbeta(u, shape1 = alpha + 1, shape2 = beta + 1) /
          (alpha * x[,2] + beta * x[,1]) 
        d[!ext] <- log(.expr1 + .c1 * .expr2) - v + jac
    }
    if(!log) d <- exp(d)
    d   
}

"dbvamix"<- 
function(x, alpha, beta, mar1 = c(0,1,0), mar2 = mar1, log = FALSE)
{
    if(length(alpha) != 1 || mode(alpha) != "numeric")
        stop("invalid argument for `alpha'")
    if(length(beta) != 1 || mode(beta) != "numeric")
        stop("invalid argument for `beta'")
    if(alpha < 0)
        stop("`alpha' must be non-negative")
    if((alpha + beta) > 1)
        stop("`alpha' + `beta' cannot be greater than one")
    if((alpha + 2*beta) > 1)
        stop("`alpha' + `2*beta' cannot be greater than one")
    if((alpha + 3*beta) < 0)
        stop("`alpha' + `3*beta' must be non-negative")
    if(is.null(dim(x))) dim(x) <- c(1,2)
    mar1 <- matrix(t(mar1), nrow = nrow(x), ncol = 3, byrow = TRUE)
    mar2 <- matrix(t(mar2), nrow = nrow(x), ncol = 3, byrow = TRUE)
    d <- numeric(nrow(x))
    x <- mtransform(x, list(mar1, mar2))
    ext <- apply(x,1,function(z) any(z %in% c(0,Inf)))
    d[ext] <- -Inf
    if(any(!ext)) {
        x <- x[!ext, ,drop=FALSE]
        mar1 <- mar1[!ext, ,drop=FALSE]
        mar2 <- mar2[!ext, ,drop=FALSE]

        xsum <- apply(x, 1, sum)  
        v <- xsum - (alpha + beta) * x[,1] + alpha * (x[,1]^2)/xsum +
          beta * (x[,1]^3)/(xsum^2)
        lx <- log(x)
        jac <- (1+mar1[,3])*lx[,1] + (1+mar2[,3])*lx[,2] -
            log(mar1[,2]*mar2[,2])

        x1a <- x[,1]/xsum; x2a <- x[,2]/xsum
        v1 <- 1 - alpha * (x2a)^2 - beta * (3 * x2a^2 - 2 * x2a^3)
        v2 <- 1 - alpha * (x1a)^2 - 2 * beta * x1a^3
        v12 <- (-2 * alpha * x1a * x2a - 6 * beta * x1a^2 * x2a) / xsum
        d[!ext] <- log(v1 * v2 - v12) - v + jac
    }
    if(!log) d <- exp(d)
    d   
}




