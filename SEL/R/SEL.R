SEL <- function(x, alpha, bounds = c(0,1), d = 4, inknts = x, N,
                gamma, Delta, fitbnds = c(1e-8, 10)*diff(bounds),
                pistar = NULL,
                constr = c("none", "unimodal", "decreasing", "increasing"),
                mode){
  nord <- d+1
  # Check for invalid arguments
  if(length(x) != length(alpha))
    stop("x and alpha need to be of the same length.")
  if(any(alpha <= 0) | any(alpha >= 1))
    stop("alpha values need to be in (0,1)")
  if(any(x>=bounds[2]) | any(x<=bounds[1]))
    stop("'x' needs to lie within bounds specified via 'bounds'")
  if(any(order(x) != order(alpha))){
    stop("invalid alpha selected for x (quantiles need to be increasing with increasing x).")
  }
  constr <- match.arg(constr)
  if(missing(mode) & constr == "unimodal")
    stop("need value of mode, when unimodal shape constraint is used.")
  if(nord <= 1)
    stop("degree d needs to be larger than 0.")
  if(length(bounds) != 2)
    stop("'bounds' needs to be numeric of length 2.")
  if(bounds[2] <= bounds[1])
    stop("'bounds[1]' needs to be smaller than 'bounds[2]'.")
  if(!missing(mode)){
    if((mode < bounds[1])|(mode > bounds[2])){
      stop("'mode' needs to be inside bounds.")
    }
  }
  if(any(alpha >= 1)|any(alpha <= 0))
    stop("values in 'alpha' need to be in (0,1)")
  if(!missing(N)){
    if(N==0) inknts <- NULL
    else inknts <- (2*(1:N)-1)/(2*N)*(diff(bounds))+bounds[1]
  }
  constr <- match.arg(constr)
  ord <- order(x)
  x <- x[ord]
  alpha <- alpha[ord]
  knts <- c(rep(bounds[1],nord), inknts, rep(bounds[2],nord))
  X <- splineDesign(x, knots = knts, ord = nord)
  if(!is.null(inknts)){
    if(nord != 2){
      P <- get.PmatBS(inknts, bounds, nord)
    } else {
      P <- get.PmatPC(inknts, bounds)
    }
  } else {
    P <- get.PmatBern(nord, bounds)
  }
  Ab <- get.Ab(length(inknts)+nord, constr, knts, length(inknts), d, mode)
  A <- Ab$A;bvec <- Ab$b

  if(missing(Delta)){
    xstand <- (x-bounds[1])/diff(bounds)
    mdelta2 <- mean((xstand-alpha)^2)
    Delta <- sqrt(mdelta2)/2
  }
  if(!is.null(pistar)){
    dplus <- get.dplus(inknts, bounds, nord, pistar)
  } else {
    dplus <- 0
  }
  if(missing(gamma)){
    gamma <- min.abs(Delta, X, alpha, P, A,
                     bvec, fitbnds[1], fitbnds[2], dplus)
  }
  sol <- fit.SEL(X, alpha, P, A, bvec, gamma, dplus)
  res <- list()
  res$constr <- constr
  res$inknts <- inknts
  res$nord <- nord
  res$bounds <- bounds
  res$gamma <- gamma
  res$xalpha <- list(x=x, alpha=alpha)
  res$coefs <- sol
  res$Omega <- P
  res$pistar <- pistar
  res$dplus <- dplus  
  class(res) <- "SEL"
  res
}

prod.Bspline <- function(x, nord, knots, i , j){
  X <- splineDesign(x, knots = knots, ord = nord, derivs=rep(1, length(x)))
  X[,i]*X[,j]
}

int.funct <- function(knots, nord, tol, i, j){
  diff <- j-i
  integrate(prod.Bspline, knots[i+diff], knots[j+nord-diff], nord=nord,
            knots = knots, rel.tol = tol, i=i, j=j)$value
}

get.PmatBS <- function(inknts, bounds, nord){
  knots <- c(rep(bounds[1], nord), inknts, rep(bounds[2], nord))
  tol <- .Machine$double.eps^0.5/max(diff(range(knots)), 1)
  n <- nord + length(inknts)
  P <- matrix(0, nrow=n, ncol=n)
  for(i in 1:n){
    for(j in i:min(i+nord-1,n)){
      P[i,j] <- int.funct(knots, nord, tol, i, j)
      if(i != j)
        P[j,i] <- P[i,j]
    }
  }
  P
}

itgrl <- function(alpha1, alpha2, n){
  # integral over product of beta dens.
  # B(alpha1,n-alpha1)*B(alpha2, n-alpha2)
  beta(alpha1+alpha2-1, n-alpha1+n-alpha2-1)/
    (beta(alpha1, n-alpha1)*beta(alpha2, n-alpha2))
}

get.PmatBern <- function(dim, bounds){
  P <- outer(1:(dim-1), 1:(dim-1), "itgrl", n=dim)
  P <- P/diff(bounds)
  D <- -1*diff(diag(dim))
  P <- t(D)%*%P%*%D
  P
}

get.PmatPC <- function(inknts, bounds){
  knots <- c(bounds[1], inknts, bounds[2])
  P <- diag(1/diff(knots))
  dim <- length(inknts)+2
  D <- -1*diff(diag(dim))
  P <- t(D)%*%P%*%D
  P
}

prod.densBsp <- function(x, nord, knots, dens, i){
  X <- splineDesign(x, knots = knots, ord = nord, derivs=rep(1, length(x)))
  resp <- do.call(dens, list(x))
  resp*X[,i]
}

int.funct2 <- function(knots, nord, tol, dens, i){
  integrate(prod.densBsp, knots[i], knots[i+nord], nord=nord,
            dens = dens, knots = knots, rel.tol = tol, i=i)$value
}

get.dplus <- function(inknts, bounds, nord, dens){
  knots <- c(rep(bounds[1], nord), inknts, rep(bounds[2], nord))
  tol <- .Machine$double.eps^0.5/max(diff(range(knots)), 1)
  n <- nord + length(inknts)
  dplus <- numeric(n)
  for(i in 1:(length(knots)-nord)){
    dplus[i] <- int.funct2(knots, nord, tol, dens, i)
  }
  dplus
}

get.Ab <- function(dim, constr, knts, N, d, mode){
  r1 <- c(1,rep(0,dim-1))
  r2 <- c(rep(0,dim-1),1)
  B <- diff(diag(dim))
  A <- rbind(r1,r2,B)
  bvec <- c(0, 1, rep(0,dim-1))
  if(constr != "none"){
    A2 <- get.Conv(knts, dim, constr, N, d, mode)
    A <- rbind(A,A2)
    bvec <- c(bvec, rep(0, nrow(A2)))
  }
  list(A = t(A), b = bvec)
}

get.Conv <- function(knts, dim, constr, N, d, mode){
  knts <- knts[-c(1, length(knts))]
  difs <- diff(knts, lag=d)
  res <- matrix(0, nrow=dim-2, ncol=dim)
  for(i in 1:(dim-2)){
    res[i,i:(i+2)] <- c(1/difs[i], -1/difs[i]-1/difs[i+1], 1/difs[i+1])
  }
  res <- res*(-1)^(constr == "decreasing")
  if(constr == "unimodal"){
    rsk <- rep(0.0, N+d)
    kntave <- .C("knotave", as.double(knts), as.integer(length(rsk)),
              as.integer(d-1), res=as.double(rsk))$res
    ind <- which.min(abs(mode-kntave))
    inds <- c(rep(1, ind-1), rep(-1, length(rsk)-ind))
    res <- res*inds
  }
  res
}

fit.SEL <- function(N, alpha, P, A, b, gamma, dplus = 0){
  D <- crossprod(N)+gamma*P
  d <- c(t(alpha)%*%N) + gamma*dplus
  solve.QP(Dmat=D, dvec=d, bvec=b, Amat=A, meq=2)$solution
}

min.abs <- function(Delta, N, alpha, P, A, b, lb, ub, dplus = 0){
  min.foo <- function(gamma, Delta, N, alpha, P, A, b, dplus){
    ff <- fit.SEL(N, alpha, P, A, b, gamma, dplus)
    mean((alpha-N%*%ff)^2) - Delta^2
  }
  gam <- try(uniroot(min.foo, c(lb, ub), Delta=Delta, N=N, alpha=alpha, P=P,
                     A=A, b=b, dplus = dplus))
  if(is.list(gam))
    gam <- gam$root
  else{
    stop("Cannot find gamma matching Delta.\n 
         Try more inner knots (or a smaller fitbnds[1]).")
  }
}

plot.SEL <- function(x, ..., type = c("density", "cdf"),
                     deriv, points = TRUE, n = 101, xlab = "",
                     ylab = "", ylim){
  xx <- seq(x$bounds[1], x$bounds[2], length = n)
  if(missing(deriv)){
    type <- match.arg(type)
    if(length(type) > 2)
      type <- type[1]
    deriv <- as.numeric(type == "density")
    if(missing(ylab)){
      ylab <- type
    }
  }
  knts <- c(rep(x$bounds[1],x$nord), x$inknts, rep(x$bounds[2],x$nord))
  X <- splineDesign(xx, knots = knts, ord = x$nord,
                    derivs = rep(deriv, n))
  if(deriv == 1 & missing(ylim)){
    yl <- c(0, max(X%*%x$coefs))
  } else if(missing(ylim)){
    yl <- NULL
  } else {
    yl <- ylim
  }
  plot(xx, X%*%x$coefs, type = "l", ylim = yl, ylab=ylab, xlab=xlab, ...)
  if(points & !deriv)
    points(x$xalpha$x, x$xalpha$alpha)
}

comparePlot <- function (..., type = c("density", "cdf"), deriv, points = TRUE, 
    superpose = FALSE, n = 101, xlab = "", ylab, addArgs = NULL)
{
    if (missing(deriv)) {
        type <- match.arg(type)
        if (length(type) > 2) 
            type <- type[1]
        deriv <- as.numeric(type == "density")
    }
    if(missing(ylab)){
      ylab <- type
    }
    objs <- list(...)
    nobj <- length(objs)
    nams <- as.character(match.call())[2:(nobj + 1)]
    x0 <- objs[[1]]$bounds[1]
    x1 <- objs[[1]]$bounds[2]
    xseq <- seq(x0, x1, length = n)
    xvec <- rep(xseq, nobj)
    xa <- objs[[1]]$xalpha
    plotvec <- namvec <- numeric(n * nobj)
    for (i in 1:nobj) {
        if (class(objs[[i]]) != "SEL")
          stop("need objects of class SEL")
        if ((objs[[i]]$bounds[1] != x0) | (objs[[i]]$bounds[2] != x1)) {
            stop("bounds of all SEL-objects need to be the same")
        }
        if (!identical(xa, objs[[i]]$xalpha)) {
            stop("Specified x,alpha statements should be the same for all objects")
        }
        nord <- objs[[i]]$nord
        inknts <- objs[[i]]$inknts
        knts <- c(rep(x0, nord), inknts, rep(x1, nord))
        X <- splineDesign(xseq, knots = knts, ord = nord, derivs = rep(deriv, n))
        ind <- (((i - 1) * n + 1):(i * n))
        plotvec[ind] <- as.numeric(X %*% objs[[i]]$coefs)
        namvec[ind] <- rep(nams[i], n)
    }
    namvec <- as.factor(namvec)
    plotdf <- data.frame(xvec, plotvec, namvec)
    if (deriv == 1) {
        yl <- c(-0.1, 1.1) * max(plotvec)
    }
    else {
        yl <- c(-0.1, 1.1)
    }
    if(superpose){
      spL <- trellis.par.get("superpose.line")
      spL$lty <- rep(spL$lty, nobj%/%length(spL$lty) + 1)[1:nobj]
      spL$lwd <- rep(spL$lwd, nobj%/%length(spL$lwd) + 1)[1:nobj]
      spL$col <- rep(spL$col, nobj%/%length(spL$col) + 1)[1:nobj]
      callList <-list(plotvec ~ xvec, data = plotdf, subscripts = TRUE,
             groups = plotdf$namvec, panel.data = list(f=deriv, xa=xa),
          panel = function(x, y, subscripts, groups, panel.data) {
           panel.superpose(x, y, subscripts, groups, type = "l")
           if(!panel.data$f & points){
             panel.xyplot(panel.data$xa$x, panel.data$xa$alpha, 
                          type="p", col = 1)
           }
          }, xlab = xlab, 
          ylab = ylab, key = list(lines = spL, transparent = TRUE, 
          text = list(levels(plotdf$namvec), cex = 0.9), columns = ifelse(nobj < 
          5, nobj, min(4,ceiling(nobj/min(ceiling(nobj/4),3))))))
      callList <- c(callList, addArgs)
      do.call("xyplot", callList)

    } else {
       callList <- list(plotvec ~ xvec | namvec, data = plotdf, ylim = yl, 
          as.table = T, xlab = xlab, ylab = ylab, panel.data = list(f = deriv, 
              xa = xa), panel = function(x, y, panel.data) {
              panel.xyplot(x, y, type = "l")
              if (!panel.data$f & points) {
                  panel.xyplot(panel.data$xa$x, panel.data$xa$alpha)
              }
          })
       callList <- c(callList, addArgs)
       do.call("xyplot", callList)
    }
}

print.SEL <- function(x, ...){
  knts <- c(rep(x$bounds[1],x$nord), x$inknts, rep(x$bounds[2],x$nord))
  X <- splineDesign(x$xalpha$x, knots = knts, ord = x$nord)
  fit <- x$coefs
  rss <- mean((x$xalpha$alpha-X%*%fit)^2)
  ent <- t(fit)%*%x$Omega%*%fit
  cat("Fitted SEL object\n\n")
  cat("Bounds: ", x$bounds, "\n")
  if(!is.null(x$inknts)){
    knts <- x$inknts
  } else {
    knts <- "none"
  }
  cat("Inner Knots:", knts, "\n")
  cat("Degree: ", x$nord-1, "/ Dimension: ", dim(x$Omega)[1], "\n")
  cat("Shape Constraint:", x$constr, "\n")
  cat("gamma: ", signif(x$gamma, 4), "/ Delta", signif(sqrt(rss), 4), "\n")
  if(is.null(x$pistar)){
    cat("Brier Neg-Entropy: ", ent, "\n")
  } else {
    foo <- function(y, x) x$pistar(y)^2
    div <- try(integrate(foo, x$bounds[1], x$bounds[2], x=x)$value)
    if(class(div) == "try-error"){
      cat("Brier Neg-Entropy: ", ent, "\n")      
    } else {
      cat("Brier divergence to pistar: ", div+ent-2*x$dplus%*%fit, "\n")
    }
  }
}

summary.SEL <-  function(object,...){
  oldClass(object) <- "summary.SEL"
  print(object, ...)
}

print.summary.SEL <- function(x, ...){
  res <- numeric(7)
  res[c(1,5)] <- x$bounds
  oldClass(x) <- "SEL"
  res[2:4] <- quantSEL(c(0.25, 0.5, 0.75), x)
  foo1 <- function(y, x){
    y*predict(x, newdata=y)
  }
  foo2 <- function(y, x){
    y^2*predict(x, newdata=y)
  }
  res[6] <- integrate(foo1, res[1], res[5], x=x)$value
  res[7] <- integrate(foo2, res[1], res[5], x=x)$value
  res[7] <- res[7] - res[6]^2
  names(res) <- c("Min.","1st Qu.","Median",
                 "3rd Qu.","Max.","Mean", "Var.")
  res
}

predict.SEL <- function(object, newdata = seq(object$bounds[1],
                   object$bounds[2], length = 101),
                   type = c("density", "cdf"), deriv, ...){

  if(missing(deriv)){
    type <- match.arg(type)
    if(length(type) > 2)
      type <- type[1]
    deriv <- as.numeric(type == "density")
  }

  # evaluate distribution at newdata
  knts <- c(rep(object$bounds[1],object$nord), object$inknts,
            rep(object$bounds[2],object$nord))
  X <- splineDesign(newdata, knots = knts, ord = object$nord, 
              deriv = rep(deriv, length(newdata)), outer.ok = TRUE)
  as.numeric(X%*%object$coefs)
}

coef.SEL <- function(object, ...){
  object$coefs
}

quantSEL <- function(q, object, nPoints = 1000){
  # q - vector of quantiles
  # object - SEL object
  # nPoints - number of equally spaced points on which
  #           to evaluate cdf
  if(any(q < 0) | any(q > 1))
    stop("q values need to be in [0,1]")
  if(class(object) != "SEL")
    stop("Need object of class SEL in 'object'")
  bnds <- object$bounds
  sq <- seq(bnds[1], bnds[2], length = nPoints)
  pred <- predict(object, newdata = sq, deriv = 0)
  quantfunc <- splinefun(pred, sq, method = "monoH.FC", ties = mean)
  quantfunc(q)
}

rvSEL <- function(n, object, nPoints = 1000){
  # n - number of simulated variates
  # object - SEL object
  # nPoints - number of equally spaced points on which
  #           to evaluate cdf
  u <- runif(n)
  quantSEL(u, object = object, nPoints = nPoints)
}

# update weights analytically for Bernstein polynomial
updWeights <- function(F, k, n){
  p <- diff(F)
  p[p<=0] <- 1e-20
  p <- p/sum(p)
  ks <- length(p)
  j <- 1:ks
  postp <- log(p)+lfactorial(ks)+lfactorial(k+j-1)+lfactorial(n+ks-k-j)
  postp <- postp-lfactorial(j-1)-lfactorial(ks-j)-lfactorial(n+ks-1)
  postp <- exp(postp-mean(postp))
  postp <- postp/sum(postp)
  postp
}

# posterior for binomial data
getbinPost <- function(x, object, k, n, type = c("density", "cdf"),
                       rel.tol = .Machine$double.eps^0.5){
  if(class(object) != "SEL")
    stop("Need object of class SEL in 'object'")
  bnds <- object$bounds
  if(bnds[1] != 0 | bnds[2] != 1)
    stop("support of the SEL object needs to be [0,1]")
  type <- match.arg(type)
  if(is.null(object$inknts)){ # Bernstein polynomial
    F <- object$coefs
    ks <- length(F)-1
    probs <- updWeights(F, k, n)
    if(type == "density"){
      func <- dbeta
    } else {
      func <- pbeta
    }
    mat <- apply(data.frame(1:ks), 1, function(y)
                 func(x, k+y, n+ks-k-y+1))
    out <- c(mat%*%probs)
  } else {
    fu <- function(x, k, n, object){
      dbinom(k,n,x)*predict(object, x)
    }
    nc <- integrate(fu, 0, 1, k=k, n=n,
                    object=object, rel.tol = rel.tol)$value
    if(type == "density"){
      out <- fu(x, k, n, object)/nc
    } else {
      i <- 1:10000
      grd <- (2*i-1)/(2*10000)
      vals <- cumsum(fu(grd, k, n, object))/(10000*nc)
      predfoo <- splinefun(c(0,grd,1), c(0,vals,1),
                   method = "monoH.FC", ties = mean)
      out <- predfoo(x)
    }
  }
  out
}

#knotave <- function(knots, d){
#  ln <- length(knots)
#  out <- numeric(ln-(d+1))
#  out[1] <- sum(knots[2:(d+1)])
#  for(i in 2:(ln-(d+1))){
#    out[i] <- out[i-1] - knots[i] + knots[i+d]
#  }
#  out/d
#}

knotave <- function(knots, d){
  ln <- length(knots)
  rsk <- numeric(ln-d-1)
  .C("knotave", as.double(knots), as.integer(length(rsk)),
               as.integer(d), res=as.double(rsk))$res
}
