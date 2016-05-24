som.init <- function(data, xdim, ydim, init="linear") {
  ## make sure xdim * ydim > 1
  if (xdim == 1 && ydim == 1) stop("Need at least two map cells.")
  INIT <- c("sample", "random", "linear")
  init.type <- pmatch(init, INIT)
  if (is.na(init.type))
    stop("Init only supports `sample', `random', and `linear'")

  SampleInit <- function(xdim, ydim) {
    ind <- sample(1:dim(data)[1], size=xdim*ydim)
    data[ind,]
  }

  RandomInit <- function(xdim, ydim) {
    ##uniformly random in (min, max) for each dimension
    ans <- matrix(NA, xdim * ydim, dim(data)[2])
    mi <- apply(data, 2, min)
    ma <- apply(data, 2, max)
    for (i in 1:(xdim*ydim)) {
      ans[i,] <- mi + (ma - mi) * runif(ncol(ans))
    }
    ans
    ##matrix(rnorm(xdim * ydim * dim(data)[2]), xdim*ydim)
  }

  LinearInit <- function(xdim, ydim) {
    ## get the first two principle components
    pcm <- prcomp(data)
    pc <- pcm$rotation[,1:2]
    sd <- pcm$sdev[1:2]
    mn <- apply(data, 2, mean)
    ans <- matrix(NA, xdim * ydim, dim(data)[2])
    ## give the 1st pc to the bigger dimension
    if (xdim >= ydim) {
      xtick <- sd[1] * pc[,1]
      ytick <- sd[2] * pc[,2]
    }else {
      xtick <- sd[2] * pc[,2]
      ytick <- sd[1] * pc[,1]
    }
    if (xdim == 1) xis <- rep(0, xdim)
    else xis <- seq(-2, 2, length=xdim)
    if (ydim == 1) yis <- rep(0, ydim)
    else yis <- seq(-2, 2, length=ydim)
    for (i in 1:(xdim*ydim)) {
      ##xf <- 4 * (i - 1) %% xdim / (xdim - 1) - 2
      ##yf <- 4 * (i - 1) %/% xdim / (ydim - 1) - 2
      xi <- (i - 1) %% xdim + 1
      yi <- (i - 1) %/% xdim + 1
      ans[i, ] <- mn + xis[xi] * xtick + yis[yi] * ytick
    }
    ans
  }
  
  if (init.type == 1) code <- SampleInit(xdim, ydim)
  if (init.type == 2) code <- RandomInit(xdim, ydim)
  if (init.type == 3) code <- LinearInit(xdim, ydim)
  code
}

som <- function(data, xdim, ydim,
                init = "linear",
                alpha = NULL, alphaType = "inverse",
                neigh = "gaussian", topol = "rect",
                radius = NULL, rlen = NULL, err.radius=1, inv.alp.c=NULL) {
  code <-  som.init(data, xdim, ydim, init)

  ALPHA <- c("linear", "inverse")
  alpha.type <- pmatch(alphaType, ALPHA)
  if (is.na(alpha.type)) stop("AlphaType only supports `linear' and `inverse'.")

  NEIGH <- c("bubble", "gaussian")
  neigh.type <- pmatch(neigh, NEIGH)
  if (is.na(neigh.type)) stop("Neigh only supports `bubble' and `gaussian'")

  TOPOL <- c("rect", "hexa")
  topol.type <- pmatch(topol, TOPOL)
  if (is.na(topol.type)) stop("Topol only supports `rect'.")

  if (is.null(alpha)) alpha <- c(0.05, 0.02)
  else if (length(alpha) == 1) alpha <- c(alpha, alpha/2)
  
  if (is.null(radius)) radius <- c(min(xdim, ydim), min(3, min(xdim, ydim)))
  else if (length(radius) == 1) radius <- c(radius, max(2, radius))

  if (is.null(rlen)) rlen <- c(dim(data)[1] * 2, dim(data)[1] * 10)
  else if (length(rlen) == 1) rlen <- c(rlen, rlen*10)

  if (is.null(inv.alp.c)) inv.alp.c <- rlen/100;

  ## for first round training 
  paramv1 <- list(alpha=alpha.type, neigh=neigh.type, topol=topol.type,
                   alpha0=alpha[1], radius0=radius[1], rlen=rlen[1],
                   err.radius=err.radius, xdim=xdim, ydim=ydim,
                   inv.alp.c=inv.alp.c[1])

  ## for second round training
  paramv2 <- list(alpha=alpha.type, neigh=neigh.type, topol=topol.type,
                   alpha0=alpha[2], radius0=radius[2], rlen=rlen[2],
                   err.radius=err.radius, xdim=xdim, ydim=ydim,
                   inv.alp.c=inv.alp.c[2])

  ## sparamv <- list(alpha=alpha.type, neigh=neigh.type, topol=topol.type)
  foo <- .Call("som_bat", as.matrix(data), code, paramv1, paramv2, PACKAGE="som")

  if (!is.null(names(data))) names(foo$code) <- names(data)
  foo$visual <- as.data.frame(foo$visual)
  names(foo$visual) <-  c("x", "y", "qerror")
  foo$qerror <- foo$qerror
  
  ##foo$data <- deparse(substitute(data))

  foo <- list(data=data, code=foo$code, visual=foo$visual,
              qerror=foo$qerror, init=init,
              alpha=ALPHA[alpha.type],
              neigh=NEIGH[neigh.type],
              topol=TOPOL[topol.type],
              alpha0=alpha, radius0=radius, rlen=rlen,
              xdim=xdim, ydim=ydim,
              err.radius=err.radius,
              inv.alp.c=inv.alp.c)
  foo$code.sum <- somsum(foo)
  ## foo$som.par <- som.par
  class(foo) <- "som"
  foo
}

som.train <- function(data, code, xdim, ydim,
                      ##init = "linear",
                      alpha = NULL, alphaType = "inverse",
                      neigh = "gaussian", topol = "rect",
                      radius = NULL, rlen = NULL, err.radius=1,
                      inv.alp.c=NULL) {
  ALPHA <- c("linear", "inverse")
  alpha.type <- pmatch(alphaType, ALPHA)
  if (is.na(alpha.type)) stop("AlphaType only supports `linear' and `inverse'.")

  NEIGH <- c("bubble", "gaussian")
  neigh.type <- pmatch(neigh, NEIGH)
  if (is.na(neigh.type)) stop("Neigh only supports `bubble' and `gaussian'")

  TOPOL <- c("rect", "hexa")
  topol.type <- pmatch(topol, TOPOL)
  if (is.na(topol.type)) stop("Topol only supports `rect'.")

  if (is.null(alpha)) alpha <- 0.05
  
  if (is.null(radius)) radius <- min(xdim, ydim)
  
  if (is.null(rlen)) rlen <- dim(data)[1] * 2

  if (is.null(inv.alp.c)) inv.alp.c <- rlen/100;

  paramv <- list(alpha=alpha.type, neigh=neigh.type, topol=topol.type,
                 alpha0=alpha[1], radius0=radius[1], rlen=rlen[1],
                 err.radius=err.radius, xdim=xdim, ydim=ydim,
                 inv.alp.c=inv.alp.c[1])
  foo <- .Call("som", as.matrix(data), code, paramv, PACKAGE="som")
  if (!is.null(names(data))) names(foo$code) <- names(data)
  foo$visual <- as.data.frame(foo$visual)
  names(foo$visual) <-  c("x", "y", "qerror")
  ##foo$data <- deparse(substitute(data))
  foo <- list(data=data, code=foo$code, visual=foo$visual,
              qerror=foo$qerror, init=NULL,
              alpha=ALPHA[alpha.type],
              neigh=NEIGH[neigh.type],
              topol=TOPOL[topol.type],
              alpha0=paramv$alpha0, radius0=paramv$radius0, rlen=paramv$rlen,
              xdim=paramv$xdim, ydim=paramv$ydim,
              err.radius=paramv$err.radius,
              inv.alp.c=paramv$inv.alp.c)
  foo$code.sum <- somsum(foo)
  ## foo$som.par <- som.par
  class(foo) <- "som"
  foo
}

som.par <- function(obj) {
  ALPHA <- c("linear", "inverse")
  NEIGH <- c("bubble", "gaussian")
  TOPOL <- c("rect", "hexa")
  alpha.type <- pmatch(obj$alpha, ALPHA)
  neigh.type <- pmatch(obj$neigh, NEIGH)
  topol.type <- pmatch(obj$topol, TOPOL)
  
  paramv <- list(alpha=alpha.type, neigh=neigh.type, topol=topol.type,
                 alpha0=obj$alpha0, radius0=obj$radius0, rlen=obj$rlen,
                 err.radius=obj$err.radius, xdim=obj$xdim, ydim=obj$ydim,
                 inv.alp.c=obj$inv.alp.c[1])
  paramv
}

som.update <- function(obj, alpha = NULL, radius = NULL,
                          rlen = NULL, err.radius = 1, inv.alp.c = NULL) {
  paramv <- som.par(obj)
  if (is.null(alpha)) paramv$alpha0 <- paramv$alpha0 / 2
  if (is.null(radius)) paramv$radius0 <- min(3, min(obj$xdim, obj$ydim))
  if (is.null(rlen)) paramv$rlen <- 5 * paramv$rlen
  if (is.null(inv.alp.c)) paramv$inv.alp.c <- paramv$rlen / 100
  paramv$err.radius <- err.radius
  
  foo <- .Call("som", as.matrix(obj$data), obj$code, paramv, PACKAGE="som")
  if (!is.null(names(data))) names(foo$code) <- names(data)
  foo$visual <- as.data.frame(foo$visual)
  names(foo$visual) <-  c("x", "y", "qerror")
  foo$qerror <- foo$qerror
  ##foo$data <- deparse(substitute(data))

  foo <- list(data=obj$data, code=foo$code, visual=foo$visual,
              qerror=foo$qerror, init=obj$init,
              alpha=obj$alpha, neigh=obj$neigh, topol=obj$topol,
              alpha0=paramv$alpha0, radius0=paramv$radius0, rlen=paramv$rlen,
              xdim=paramv$xdim, ydim=paramv$ydim,
              err.radius=paramv$err.radius,
              inv.alp.c=paramv$inv.alp.c)
  foo$code.sum <- somsum(foo)
  ## foo$som.par <- som.par
  class(foo) <- "som"
  foo
}
  
som.project <- function(obj, newdat) {
  obj$data <- newdat
  foo <- som.update(obj, rlen = 0)
  foo$visual
}

qerror <- function(obj, err.radius = 1) {
  if (err.radius == 1) obj$qerror
  else {
    ## som.update(obj, rlen=0, err.radius=err.radius)$qerror
    paramv <- som.par(obj)
    paramv$err.radius <- err.radius
    .Call("som", as.matrix(obj$data), obj$code, paramv, PACKAGE="som")$qerror
  }
}

summary.som <- function(object, ...) {
  cat("Initialization: ", object$init, "\n")
  cat("Topology: ", object$topol, "\n")
  cat("Neighborhood type: ", object$neigh, "\n")
  cat("Learning rate type: ", object$alpha, "\n")
  cat("Initial learning rate parameter: ", object$alpha0, "\n")
  cat("Initial radius of training area: ", object$radius0, "\n")
  cat("Average quantization error: ", mean(object$visual$qerror), "\n")
  cat("Average distortion measure: ", object$qerror,
        "with error radius: ", object$err.radius, "\n")
}

print.som <- function(x, ...) {
  summary(x)
  cat("The code book is:\n")
  print(cbind(x$code.sum, x$code))
}

somsum <- function(obj) { 
  xdim <- obj$xdim 
  ydim <- obj$ydim 
  m <- nrow(obj$code) 
  x <- (1:m - 1) %% xdim 
  y <- (1:m - 1) %/% xdim 
  f <- function(ii) { 
    x <- (ii - 1) %% xdim 
    y <- (ii - 1) %/% xdim 
    ind <- obj$visual$x == x & obj$visual$y == y 
    n <- length(ind[ind]) 
    n 
  } 
  nobs <- sapply(1:m, f) 
  data.frame(x, y, nobs) 
} 

filtering <- function(x, lt=20, ut=16000, mmr=3, mmd=200) {
  if (!is.matrix(x)) x <- as.matrix(x)
  ## floor and ceiling
  n <- nrow(x)
  x[x < lt] <- lt
  x[x > ut] <- ut
  tmp1 <- apply(x, 1, max)
  tmp2 <- apply(x, 1, min)
  wch <- (1:n)[ (tmp1/tmp2 > mmr) & (tmp1 - tmp2 > mmd)]
  x[wch,]
}

normalize <- function(x, byrow=TRUE) {
  if (is.vector(x))
    scale(x)
  else if (is.matrix(x) || is.data.frame(x)) {
    if (byrow) t(apply(x, 1, scale))
    else apply(x, 2, scale)
    }
  else stop("The object to be normalized must be vector, matrix, or dataframe.\n")
}

## The following are to be orphaned, not clean code

inrange <- function (x, xlim) {
  (x > xlim[1] & x < xlim[2])
}
  
ciplot <- function(x, y, se, n=1, d, ywindow) {
  ind <-  (inrange(y+n*se, ywindow) & inrange(y-n*se, ywindow))
  if (any(ind)) {
    x <- x[ind];y <- y[ind];se <- se[ind]
    segments(x, y+n*se, x, y-n*se)
    segments(x-d, y+n*se, x+d, y+n*se)
    segments(x-d, y-n*se, x+d, y-n*se)
  }
}

plotcell <- function(x, y, dat, code, n, sdbar=1, ylim, yadj) { 
##  yadj <- 0.1
  text(x+1/2, y+(1-yadj/2), paste("n=", n, sep=""))
  if (!is.data.frame(dat)) dat <- as.data.frame(dat)
  ylen <- diff(ylim)
  ## n <- nrow(dat)
  mm <- code
  l <- length(code)
  if (n > 1) {
##      l <- ncol(dat)
##      mm <- sapply(dat, mean)
    ss <- sapply(dat, sd)/sqrt(n)
  }
  else {
##      mm <- dat[,1]
##      l <- length(mm)
##      print(mm)
    ss <- rep(0, l)
  }
  lines(x+((1:l)-1/2)/l, y+1+(ylim[1] + mm)*(1-yadj)/ylen)
  if (sdbar > 0 && n > 1)
    ciplot(x+((1:l)-1/2)/l, y+1+(ylim[1] + mm)*(1-yadj)/ylen,
           ss, n=sdbar, 1/100, ywindow=c(y, y+1-yadj))
}

somgrids <- function(xdim, ydim, color,
                     yadj=0.1, hexa, ntik, ylim) {
  if (color) color <- rainbow(ydim*xdim, start=.7, end=.1)
  else color <- rep("gray", xdim*ydim)
  for (i in 0:(ydim-1)) {
    if (hexa) d <- (i %% 2)/2
    else d <- 0
    lines(c(0+d,xdim+d), c(i,i))
    for (j in 0:xdim) {
      segments(j+d, i, j+d, i+1)
      if (j == xdim) break
      rect(j+d, i+1-yadj, j+1+d, i+1, col=color[j*ydim+i+1])
    }
    lines(c(0+d,xdim+d), c(i+1,i+1))
    if (i %% 2 == 1) axis(2, seq(i, i+1-yadj, length=ntik), seq(ylim[1], ylim[2], length=ntik))
    else axis(4, seq(i, i+1-yadj, length=ntik), seq(ylim[1], ylim[2], length=ntik))
  }
}

plot.som <- function(x, sdbar=1, ylim=c(-3, 3), color=TRUE, ntik=3, yadj=0.1,
                     xlab="", ylab="", ...) {
 if (class(x) != "som" ) stop("The funciton must apply to a som object.\n")
 hexa <- (x$topol == "hexa")
 if (hexa) d <- 1/2
 else d <- 0
 xdim <- x$xdim; ydim <- x$ydim
 plot(c(0,xdim+d), c(0,ydim), xlim=c(0,xdim+d), ylim=c(0, ydim),
      type="n", xlab=xlab, ylab=ylab, axes=FALSE, ...)
 axis(1, 0:xdim, 0:xdim)
 
 somgrids(xdim, ydim, color=color, yadj=yadj, hexa=hexa, ylim=ylim, ntik=ntik)
 
 for (i in 0:(ydim-1)) {
   if (hexa) d <- (i %% 2)/2
   else d <- 0
   for (j in 0:(xdim-1)) {
     ind <- x$visual$x==j & x$visual$y==i
     n <- length(ind[ind])
     plotcell(j+d, i,
              x$data[ind, ], x$code[i*xdim + j+1,], n,
              sdbar=sdbar, ylim=ylim, yadj=yadj)
   }
 }
}
