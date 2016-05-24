make.all.branches.t.dtlik <- function(cache, control,
                                      initial.conditions.base) {
  control <- check.control.ode(control)
  branches <- make.branches.dtlik(cache$info, control)
  initial.conditions <-
    make.initial.conditions.t(cache$info, initial.conditions.base)
  function(pars, intermediates, preset=NULL)
    all.branches.matrix(pars, cache, initial.conditions,
                          branches, preset)
}

make.initial.conditions.t <- function(info, initial.conditions) {
  tm <- info$tm
  function(init, pars, t, idx) {
    tm$set(pars)
    initial.conditions(init, tm$get(t), t, idx)
  }
}

## TODO: Still cannnot prevent ROOT.EQUI being used here.
make.rootfunc.t <- function(cache, rootfunc) {
  pars.t <- cache$info$tm$get
  t.root <- cache$depth[cache$root]
  function(ans, pars, ...)             # pars here is ignored...
    rootfunc(ans, pars.t(t.root), ...) # ...because tm version used.
}

## TODO: We could lose the t.range argument here (and in time machine)
## if we did not test that the spline data was OK.  That seems risky,
## but it's a pain of a thing to carry around, especially if the
## splines already bail nicely if we're out of range.
update.info.t <- function(info,
                          functions, t.range, nonnegative=TRUE,
                          truncate=FALSE, with.q=FALSE,
                          spline.data=NULL) {
  k.for.tm <- if ( with.q ) info$k else 0
  ## TODO: This is converging on the function testing I already have.
  ## Remove the duplication at some point.
  if ( length(functions) == 1 )
    functions <- rep(functions, info$np)
  if ( is.null(names(functions)) )
    names(functions) <- info$argnames
  tm <- make.time.machine(functions, t.range,
                          nonnegative, truncate, k.for.tm,
                          spline.data)
  info$time.varying <- TRUE
  info$tm <- tm
  info$derivs.base <- info$derivs
  info$derivs <- make.derivs.t(info$derivs.base, tm)
  info$argnames <- attr(tm, "argnames")
  info$name.pretty <- sprintf("%s (time-varying)", info$name.pretty)
  info$np <- length(attr(tm, "argnames"))
  info
}

update.cache.t <- function(cache,
                           functions, nonnegative=TRUE, truncate=FALSE,
                           with.q=FALSE, spline.data=NULL) {
  t.range <- range(0, cache$depth[cache$root])
  cache$info <- update.info.t(cache$info,
                              functions, t.range, nonnegative,
                              truncate, with.q, spline.data)
  cache
}

make.derivs.t <- function(derivs, tm) {
  ret <- function(t, y, pars)
    derivs(t, y, tm$get(t))
  attr(ret, "set") <- function(pars) tm$set(pars)
  ret
}

## Making the output useful.
predict.dtlik.t <- function(object, p, t, nt=101, v=NULL, thin=10,
                            alpha=1/20, everything=FALSE, ...) {
  cache <- get.cache(object)
  tm <- cache$info$tm
  
  if ( inherits(p, "fit.mle") || inherits(p, "mcmcsamples") )
    ## TODO: Improve the coef.mcmcsamples to allow full, here, then
    ## use full.  Possibly add a 'lik' function in that can do the
    ## resolution below?
    p <- stats::coef(p)
  if ( missing(t) )
    t <- seq(min(cache$depth), max(cache$depth), length.out=nt)
  if ( is.null(v) )
    v <- tm$names
  is.matrix <- !is.null(dim(p)) && nrow(p) > 1
  ## Thin the chain to speed things up?
  if ( is.matrix && thin > 1 )
    p <- p[seq(1, nrow(p), by=thin),,drop=FALSE]
  if ( inherits(object, "constrained") ) {
    if ( is.matrix && ncol(p) == length(argnames(object)) )
      p <- t(apply(p, 1, object, pars.only=TRUE))
    else if ( !is.matrix && length(p) == length(argnames(object)) )
      p <- object(p, pars.only=TRUE)
  }

  if ( is.matrix ) {
    get1 <- function(p, t, i) {
      tm$set(p)
      tm$get(t)[i]
    }
    if ( everything ) {
      np <- if ( is.matrix ) nrow(p) else 1
      ret <- lapply(t, function(ti)
                    t(apply(p, 1, get1, ti, i=match(v, tm$names))))
      tmp <- array(unlist(ret), c(np, length(v), length(t)))
      tmp <- aperm(tmp, c(2, 3, 1))
      dimnames(tmp) <- list(v, NULL, NULL)
      ret <- vector("list", length(v))
      names(ret) <- v
      for ( i in v )
        ret[[i]] <- tmp[i,,]
    } else {
      ret <- lapply(match(v, tm$names), function(i)
                    average.over.mcmc(p, get1, t, i=i))
      names(ret) <- v
    }
  } else {
    tm$set(p)
    ret <- t(sapply(t, tm$getv))
    if ( length(v) != ncol(ret) )
      stop("I am genuinely surprised")
    colnames(ret) <- v
  }
  list(t=t, y=ret)
}

plot.dtlik.t <- function(x, p, xlab="Time", ylab="Parameter",
                         lty=1, lwd=1, v=NULL, col=NULL,
                         nt=101, legend.pos=NULL, thin=10, ...) {
  xy <- predict(x, p, v=v, nt=nt, thin=thin)
  is.matrix <- !is.null(dim(p)) && nrow(p) > 1
  xlim <- rev(range(xy$t))  
  if ( is.matrix ) {
    v <- names(xy$y)
    if ( is.null(col) )
      col <- seq_along(v)
    fill <- add.alpha(col, .5)
    names(col) <- names(fill) <- v
    ylim <- range(lapply(xy$y, range))
    plot(NA, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, las=1)
    tt <- c(xy$t, rev(xy$t))
    for ( i in rev(v) )
      graphics::polygon(tt, c(xy$y[[i]][,"lower"], rev(xy$y[[i]][,"upper"])),
                        col=fill[i], border=NA)
    for ( i in rev(v) )
      graphics::lines(xy$t, xy$y[[i]][,"mean"], col=col[i])
  } else {
    v <- colnames(xy$y)
    if ( is.null(col) )
      col <- seq_along(v)
    matplot(xy$t, xy$y, type="l", xlim=xlim, las=1, lty=lty,
            col=col, xlab=xlab, ylab=ylab, lwd=lwd, ...)
  }
  if ( !is.null(legend.pos) )
    legend(legend.pos, v, col=col, lty=lty, lwd=lwd, bty="n")
  invisible(v)
}

## Given a function that takes arguments
##   f(p, x, ...)
## where p is the parameter vector and 'x' is the domain position.
## Additional arguments are passed in.
average.over.mcmc <- function(p, f, xx, ..., alpha=1/20) {
  g <- function(xi) {
    y <- apply(p, 1, f, xi, ...)
    c(mean(y), hdr(y, 1-alpha))
  }
  if ( is.data.frame(p) )
    p <- as.matrix(p)
  ret <- t(sapply(xx, g))
  colnames(ret) <- c("mean", "lower", "upper")
  ret
}
