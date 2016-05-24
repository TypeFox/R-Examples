
"rmvevd" <-
function(n, dep, asy, model = c("log", "alog"), d = 2, mar = c(0,1,0))
{
  model <- match.arg(model)
  if(model == "log" && !missing(asy))
    warning("ignoring `asy' argument")
    
  switch(model,
    log = rmvlog(n = n, dep = dep, d = d, mar = mar),
    alog = rmvalog(n = n, dep = dep, asy = asy, d = d, mar = mar)) 
}

"rmvlog"<-
# Uses Algorithm 2.1 in Stephenson(2003)
function(n, dep, d = 2, mar = c(0,1,0))
{
    if(length(dep) != 1 || mode(dep) != "numeric" || dep <= 0 ||
       dep > 1) stop("invalid argument for `dep'")
    sim <- .C("rmvlog_tawn",
              as.integer(n), as.integer(d), as.double(dep),
              sim = double(d*n), PACKAGE = "evd")$sim
    mtransform(matrix(1/sim, ncol=d, byrow=TRUE), mar, inv = TRUE, drp = TRUE)
}

"rmvalog"<-
# Uses Algorithm 2.2 in Stephenson(2003)
function(n, dep, asy, d = 2, mar = c(0,1,0))
{
    nb <- 2^d-1
    dep <- rep(dep, length.out = nb-d)
    asy <- mvalog.check(asy, dep, d = d)
    dep <- c(rep(1,d), dep)
    sim <- .C("rmvalog_tawn",
              as.integer(n), as.integer(d), as.integer(nb), as.double(dep),
              as.double(t(asy)), sim = double(n*d), PACKAGE = "evd")$sim
    mtransform(matrix(1/sim, ncol=d, byrow=TRUE), mar, inv = TRUE, drp = TRUE)
}

"pmvevd" <-
function(q, dep, asy, model = c("log", "alog"), d = 2,
         mar = c(0,1,0), lower.tail = TRUE)
{
  model <- match.arg(model)
  if(model == "log" && !missing(asy))
    warning("ignoring `asy' argument")
    
  switch(model,
    log = pmvlog(q = q, dep = dep, d = d, mar = mar,
      lower.tail = lower.tail),
    alog = pmvalog(q = q, dep = dep, asy = asy, d = d, mar = mar,
      lower.tail = lower.tail)) 
}

"pmvlog"<- 
function(q, dep, d = 2, mar = c(0,1,0), lower.tail = TRUE)
{
    if(lower.tail) {
      if(length(dep) != 1 || mode(dep) != "numeric" || dep <= 0 ||
        dep > 1) stop("invalid argument for `dep'")
      if(is.null(dim(q))) dim(q) <- c(1,d)
      if(ncol(q) != d) stop("`q' and `d' are not compatible")
      q <- mtransform(q, mar)
      pp <- exp(-apply(q^(1/dep),1,sum)^dep)
    } else {
      pp <- numeric(1)
      ss <- c(list(numeric(0)), subsets(d))
      ssl <- d - sapply(ss, length)
      for(i in 1:(2^d)) {
        tmpq <- q
        tmpq[, ss[[i]]] <- Inf
        pp <- (-1)^ssl[i] * Recall(tmpq, dep, d, mar) + pp
      }
    }
    pp
}

"pmvalog"<-
function(q, dep, asy, d = 2, mar = c(0,1,0), lower.tail = TRUE)
{
    if(lower.tail) {
      nb <- 2^d-1
      dep <- rep(dep, length.out = nb-d)
      asy <- mvalog.check(asy, dep, d = d)
      dep <- c(rep(1,d), dep)
      if(is.null(dim(q))) dim(q) <- c(1,d)
      if(ncol(q) != d) stop("`q' and `d' are not compatible")
      q <- mtransform(q, mar)
      inner <- function(par)
        apply((rep(par[1:d], rep(nrow(q),d))*q)^(1/par[d+1]), 1, sum)^par[d+1]
      comps <- apply(cbind(asy,dep),1,inner)
      if(is.null(dim(comps))) dim(comps) <- c(1,nb)
      pp <- exp(-apply(comps,1,sum))
    } else {
      pp <- numeric(1)
      ss <- c(list(numeric(0)), subsets(d))
      ssl <- d - sapply(ss, length)
      for(i in 1:(2^d)) {
        tmpq <- q
        tmpq[, ss[[i]]] <- Inf
        pp <- (-1)^ssl[i] * Recall(tmpq, dep, asy, d, mar) + pp
      }
    }
    pp
}

"dmvevd" <-
function(x, dep, asy, model = c("log", "alog"), d = 2,
         mar = c(0,1,0), log = FALSE)
{
  model <- match.arg(model)
  if(model == "log" && !missing(asy))
    warning("ignoring `asy' argument")

  switch(model,
    log = dmvlog(x = x, dep = dep, d = d, mar = mar, log = log),
    alog = dmvalog(x = x, dep = dep, asy = asy, d = d, mar = mar, log = log))
}   

"dmvlog"<- 
function(x, dep, d = 2, mar = c(0,1,0), log = FALSE)
{
    if(length(dep) != 1 || mode(dep) != "numeric" || dep <= 0 ||
        dep > 1) stop("invalid argument for `dep'")
    if(is.null(dim(x))) dim(x) <- c(1,d)
    if(ncol(x) != d) stop("`x' and `d' are not compatible")
    if(is.list(mar)) {
      if(length(mar) != d) stop("`mar' and `d' are not compatible")
      for(i in 1:d) mar[[i]] <- matrix(t(mar[[i]]), nrow = nrow(x),
        ncol = 3, byrow = TRUE)
    }
    else mar <- matrix(t(mar), nrow = nrow(x), ncol = 3, byrow = TRUE)
    dns <- numeric(nrow(x))
    x <- mtransform(x, mar)
    ext <- apply(x, 1, function(z) any(z %in% c(0,Inf)))
    dns[ext] <- -Inf
    idep <- 1/dep
    cf <- matrix(0, nrow = d, ncol = d)
    diag(cf) <- dep^(1:d - 1)
    cf[,1] <- exp(lgamma(1:d - dep) - lgamma(1:d) - lgamma(1 - dep))
    if(d >= 3) {
      for(i in 3:d) {
        for(j in 2:(d-1))
          cf[i,j] <- ((i - 1 - j*dep) * cf[i-1,j] + dep*(j-1)*cf[i-1,j-1])/
            (i - 1)
      }
    }
    cf <- log(cf[d,]) - lgamma(1:d)
    if(any(!ext)) {
      x <- x[!ext, ,drop=FALSE]
      if(is.list(mar)) {
        marscl <- marshp <- matrix(NA, nrow = nrow(x), ncol = d)
        for(i in 1:d) {
          mar[[i]] <- mar[[i]][!ext, ,drop=FALSE]
          marscl[,i] <- mar[[i]][,2]
          marshp[,i] <- mar[[i]][,3]
        }
      }  
      else {
        mar <- mar[!ext, ,drop=FALSE]
        marscl <- mar[,2]
        marshp <- mar[,3]
      }
      z <- rowSums(x^(idep))^dep
      lx <- log(x)
      .expr1 <- rowSums(lx * (idep + marshp) - log(marscl))
      dns[!ext] <- .expr1 + (1 - d*idep) * log(z) - z
      lz <- matrix(log(z), nrow = d, ncol = nrow(x), byrow = TRUE)
      lz <- (0:(d-1)) * lz + cf
      dns[!ext] <- dns[!ext] + log(colSums(exp(lz))) +
        lgamma(d) - (d-1)*log(dep)
    }
    if(!log) dns <- exp(dns)
    dns
}

"dmvalog"<- 
function(x,  dep, asy, d = 2, mar = c(0,1,0), log = FALSE)
{
    nb <- 2^d-1
    dep <- rep(dep, length.out = nb-d)
    ss <- mvalog.check(asy, dep, d = d, ss = TRUE)
    asy <- ss$asy ; ss <- ss$ss
    dep <- c(rep(1,d), dep)
    if(is.null(dim(x))) dim(x) <- c(1,d)
    if(ncol(x) != d) stop("`x' and `d' are not compatible")
    nn <- nrow(x)
    if(is.list(mar)) {
      if(length(mar) != d) stop("`mar' and `d' are not compatible")
      for(i in 1:d) mar[[i]] <- matrix(t(mar[[i]]), nrow = nn,
        ncol = 3, byrow = TRUE)
    }
    else mar <- matrix(t(mar), nrow = nn, ncol = 3, byrow = TRUE)
    x <- mtransform(x, mar)
    ext <- apply(x, 1, function(z) any(z %in% c(0,Inf)))
    dns <- numeric(nn)
    dns[ext] <- -Inf

    if(any(!ext)) {
      x <- x[!ext, ,drop=FALSE]
      if(is.list(mar)) {
        marscl <- marshp <- matrix(NA, nrow = nn, ncol = d)
        for(i in 1:d) {
          mar[[i]] <- mar[[i]][!ext, ,drop=FALSE]
          marscl[,i] <- mar[[i]][,2]
          marshp[,i] <- mar[[i]][,3]
        }
      }  
      else {
        mar <- mar[!ext, ,drop=FALSE]
        marscl <- mar[,2]
        marshp <- mar[,3]
      }
      cfinit <- function(dep) {
	    if(dep==1) return(diag(d))
        cf <- matrix(0, nrow = d, ncol = d)
        diag(cf) <- dep^(1:d - 1)
        cf[,1] <- exp(lgamma(1:d - dep) - lgamma(1:d) - lgamma(1 - dep))
        if(d >= 3) {
          for(i in 3:d) {
            for(j in 2:(d-1))
              cf[i,j] <- ((i - 1 - j*dep) * cf[i-1,j] + dep * (j-1) *
                cf[i-1,j-1])/(i - 1)
          }
        }
        cf
      }
      cfs <- paste("cf", 1:nb, sep = "")
      for(i in 1:nb) assign(cfs[i], cfinit(dep[i]))
      qfn <- function(lz, dep, p, cf) {
        if(p == 1) return(rep(0, nn))
        cf <- log(cf[p,1:p]) - lgamma(1:p)
        lz <- matrix(lz, nrow = p, ncol = nn, byrow = TRUE)
        lz <- ((1-p):0) * lz + cf
        log(colSums(exp(lz))) + lgamma(p) - (p-1)*log(dep) 
      }
      indm <- matrix(NA, nrow = 2^(d-1), ncol = d)
      for(i in 1:d) indm[,i] <- which(sapply(ss, function(x) i %in% x))
      indm <- as.matrix(do.call("expand.grid", as.data.frame(indm)))
      inner <- function(par) {
        int <- (rep(par[1:d], each = nn) * x)^(1/par[d+1])
        rowSums(int)^par[d+1]
      }
      zm <- log(apply(cbind(asy,dep), 1, inner))
      if(is.null(dim(zm))) dim(zm) <- c(1,nb)
      vv <- rowSums(exp(zm))
      lx <- log(x)
      mexpr <- rowSums(lx * marshp - log(marscl))
      
      tot1 <- matrix(0, nrow = nn, ncol = 2^(d*(d-1)))
      tot2 <- matrix(0, nrow = nn, ncol = d)
      for(i in 1:(2^(d*(d-1)))) {
        indmi <- indm[i,]
        pp <- tabulate(indmi, 2^d-1)[indmi]
        if(all(asy[indmi + (2^d-1)*(0:(d-1))] != 0)) {
          for(j in 1:d) 
            tot2[,j] <- 1/dep[indmi[j]] * (log(asy[indmi[j],j]) + lx[,j]) +
              (1 - 1/dep[indmi[j]]) * zm[,indmi[j]] + qfn(zm[,indmi[j]],
              dep[indmi[j]], pp[j], get(cfs[indmi[j]])) / pp[j]
          tot1[,i] <- rowSums(tot2)
        }
        else tot1[,i] <- -Inf
      }
      dns[!ext] <- log(rowSums(exp(tot1))) - vv + mexpr
    }
    if(!log) dns <- exp(dns)
    dns
}

"amvevd" <-
function(x = rep(1/d,d), dep, asy, model = c("log", "alog"), d = 3, plot =
    FALSE, col = heat.colors(12), blty = 0, grid = if(blty) 150 else 50,
    lower = 1/3, ord = 1:3, lab = as.character(1:3), lcex = 1)
{
  model <- match.arg(model)
  if(model == "log" && !missing(asy))
    warning("ignoring `asy' argument") 

  if(!plot) {
    if(is.vector(x)) x <- as.matrix(t(x))
    if(!is.matrix(x) || ncol(x) != d)
      stop("`x' must be a vector/matrix with `d' elements/columns")
    if(any(x < 0, na.rm = TRUE))
      stop("`x' must be non-negative")
    rs <- rowSums(x)
    if(any(rs <= 0, na.rm = TRUE))
      stop("row(s) of `x' must have a positive sum")
    if(max(abs(rs[!is.na(rs)] - 1)) > 1e-6)
      warning("row(s) of `x' will be rescaled")
    x <- x/rs
  }
  if(plot) {
    if(d == 2) stop("use abvnonpar for bivariate plots")
    if(d >= 4) stop("cannot plot in high dimensions")
  }
     
  switch(model,
    log = amvlog(x = x, dep = dep, d = d, plot = plot, col = col, blty =
      blty, grid = grid, lower = lower, ord = ord, lab = lab, lcex = lcex),
    alog = amvalog(x = x, dep = dep, d = d, asy = asy, plot = plot, col =
      col, blty = blty, grid = grid, lower = lower, ord = ord, lab = lab,
      lcex = lcex))
}

"amvlog"<- 
function(x, dep, d = 3, plot = FALSE, col = heat.colors(12), blty = 0,
    grid = if(blty) 150 else 50, lower = 1/3, ord = 1:3, lab =
    as.character(1:3), lcex = 1)
{
    if(length(dep) != 1 || mode(dep) != "numeric" || dep <= 0 ||
      dep > 1) stop("invalid argument for `dep'")
    depfn <- function(x, dep) rowSums(x^(1/dep))^dep
    if(plot) {
      mz <- tvdepfn(depfn = depfn, col = col, blty = blty, grid = grid,
        lower = lower, ord = ord, lab = lab, lcex = lcex, dep = dep)
      return(invisible(mz))
    }
    depfn(x = x, dep = dep)   
}

"amvalog"<- 
function(x, dep, asy, d = 3, plot = FALSE, col = heat.colors(12), blty = 0,
    grid = if(blty) 150 else 50, lower = 1/3, ord = 1:3, lab =
    as.character(1:3), lcex = 1)
{
    dep <- rep(dep, length.out = 2^d-d-1)
    asy <- mvalog.check(asy, dep, d)
    
    depfn <- function(x, dep, asy)
    {
      d <- ncol(x)
      dep <- c(rep(1,d), dep)
      tot <- matrix(0, nrow = nrow(x), ncol = 2^d-1)
      idep <- 1/dep
      x <- t(x)
      for(k in 1:(2^d-1))
        tot[,k] <- colSums((asy[k,] * x)^idep[k])^dep[k]
      rowSums(tot)
    }
    if(plot) {
      mz <- tvdepfn(depfn = depfn, col = col, blty = blty, grid = grid,
        lower = lower, ord = ord, lab = lab, lcex = lcex, dep = dep,
        asy = asy)
      return(invisible(mz))
    }
    depfn(x = x, dep = dep, asy = asy)   
}

### Ancillary Functions ###

"tvdepfn" <-
# Plots Dependence Function On Simplex S3
function(depfn, col = heat.colors(12), blty = 0, grid =
    if(blty) 150 else 50, lower = 1/3, ord = 1:3, lab = as.character(1:3),
    lcex = 1, ...)
{
    oldpar <- par(pty="s")
    on.exit(par(oldpar))

    s3 <- sqrt(1/3)
    x <- seq(-s3, s3, len = 2*grid)
    y <- seq(0, 1, len = 2*grid)
    a <- function(x,y) {
        w1 <- y
        w2 <- (x/s3 + 1 - y)/2
        vals <- numeric(length(w1))
        tpts <- cbind(w1,w2,w3=1-w1-w2)[, ord]
        nv <- apply(tpts, 1, function(b) any(b < 0))
        is.na(vals) <- nv
        vals[!nv] <- depfn(tpts[!nv,], ...)
        vals
    }
    z <- outer(x,y,a)
    mz <- min(z, na.rm = TRUE)
    mxz <- max(z, na.rm = TRUE)
    if(mz < 1/3 || mxz > 1 || any(is.nan(z)))
      warning("parameters may be too near edge of parameter space")
    else if(mz < lower)
      warning("`lower' is greater than calculated values")
    plot(c(-s3, s3), c(0.5-s3, 0.5+s3), type="n", axes=FALSE, xlab="", ylab="")
    if(!is.null(lab)) {
      lab <- lab[ord]
      eps <- 0.025
      text(0, 1+eps, lab[1], adj = c(0.5, 0), cex = lcex)
      text(s3, -eps, lab[2], adj = c(1, 1), cex = lcex)
      text(-s3, -eps, lab[3], adj = c(0, 1), cex = lcex)
    }
    image(x, y, z, zlim = c(lower,1), xlim = c(-s3,s3), ylim =
          c(0.5-s3, 0.5+s3), col = col, add = TRUE)
    polygon(c(-s3, s3, 0), c(0, 0, 1), lty = blty)
    invisible(mz)
}

"mvalog.check" <-
# Checks And Transforms Arguments For Asymmetric Logistic
function(asy, dep, d, ss = FALSE)
{
    if(mode(dep) != "numeric" || any(dep <= 0) || any(dep > 1))
      stop("invalid argument for `dep'")
    nb <- 2^d-1
    if(mode(asy) != "list" || length(asy) != nb)
      stop(paste("`asy' should be a list of length", nb))

    tasy <- function(theta, b) {
        trans <- matrix(0, nrow=nb, ncol=d)
        for(i in 1:nb) trans[i,(1:d %in% b[[i]])] <- theta[[i]]
        trans
    }
    
    b <- subsets(d)
    if(any(sapply(asy, length) != sapply(b, length)))
      stop("`asy' is not of the correct form")
    asy <- tasy(asy, b)
    if(!is.matrix(asy) || mode(asy) != "numeric")
        stop("`asy' is not of the correct form")
    if(min(asy) < 0 || max(asy) > 1)
       stop("`asy' must contain parameters in [0,1]")
    if(any(apply(asy,2,sum) != 1) || any(asy[c(rep(FALSE,d),dep==1),] != 0) ||
       any(apply(asy[-(1:d),,drop=FALSE],1,function(x) sum(x!=0)) == 1))
       stop("`asy' does not satisfy the appropriate constraints")
    if(ss) return(list(asy = asy, ss = b))
    asy
}

"subsets" <-
# Lists All Subsets Of 1:d
function(d) {
    x <- 1:d
    k <- NULL
    for(m in x) 
        k <- rbind(cbind(TRUE, k), cbind(FALSE, k))
    pset <- apply(k, 1, function(z) x[z])
    pset[sort.list(unlist(lapply(pset,length)))[-1]] 
}
  
