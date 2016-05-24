.onAttach <- function(libname,pkgname) {
  if(!havefftw()) {
    packageStartupMessage("*** ",pkgname,": FFTW not used.\n*** You should install it from http://fftw.org\n*** or check if your OS-distribution provides it, and recompile.")
  }
}


# Chebyshev transformation.  I.e. coefficients for given function values in the knots.

# The Chebyshev knots of order n on an interval
chebknots1 <- function(n, interval=NULL) {
  kn <- cos(pi*((1:n)-0.5)/n)
  if((n %% 2) == 1) kn[[(n+1)/2]] = 0
  if(is.null(interval)) return(kn)
  kn*diff(interval)/2 + mean(interval)
}

chebknots <- function(dims, intervals=NULL) {
  if(is.null(intervals)) {
    res <- lapply(dims,chebknots1)
  } else {
    if(length(dims) == 1 && is.numeric(intervals)) intervals=list(intervals)
    if(!is.list(intervals)) stop("intervals should be a list")
    res <- mapply(chebknots1,dims,intervals,SIMPLIFY=FALSE)
  }

  res
}

# evaluate a function on a Chebyshev grid
evalongrid <- function(fun,dims,intervals=NULL,...,grid=NULL) {
# do a call to stuff which doesn't really expand the grid
  if(is.numeric(grid)) grid <- list(grid)
  if(is.null(grid)) grid <- chebknots(dims,intervals)
  mf <- match.fun(fun)
  .Call(C_evalongrid,function(x) mf(x,...), grid)
#  structure(apply(expand.grid(chebknots(dims,intervals)),1,fun,...), dim=dims)
}


# Chebyshev coefficients for x, which may be an array
chebcoef <- function(val, dct=FALSE) {
  structure(.Call(C_chebcoef,as.array(val),dct),dimnames=dimnames(val))
}

chebeval <- function(x,coef,intervals=NULL) {
  if(is.null(intervals)) return(.Call(C_evalcheb,coef,x))
  # map into intervals
  .Call(C_evalcheb,coef,mapply(function(x,i) 2*(x[[1]]-mean(i))/diff(i),x,intervals))
}

# return a function which is a Chebyshev interpolation
chebappx <- function(val,intervals=NULL) {
  if(is.null(dim(val))) {
    # allow for one-dimensional
    dim(val) <- length(val)
  }
   # allow for vector, e.g. intervals=c(0,1), put it inside list
  if(is.numeric(intervals) && length(intervals) == 2) intervals <- list(intervals)

  cf <- chebcoef(val)

  if(is.null(intervals)) {
    # it's [-1,1] intervals, so drop transformation
    fun <- structure(function(x) .Call(C_evalcheb,cf,x),arity=length(dim(val)))
    rm(val)
  } else {
    # it's intervals, create mapping into [-1,1]
    if(!is.list(intervals)) stop("intervals should be a list")
    if(any(sapply(intervals,length) != 2)) stop("interval elements should have length 2")
    if(length(intervals) != length(dim(val))) stop("values should have the same dimension as intervals",
               length(intervals),length(dim(val)))
    ispan <- sapply(intervals,function(x) 2/diff(x))
    mid <- sapply(intervals,function(x) mean(x))
    imap <- cmpfun(function(x) (x-mid)*ispan)
    fun <- structure(function(x) .Call(C_evalcheb,cf,imap(x)),arity=length(dim(val)),domain=intervals)
    rm(val)
  }
  fun
}

# interpolate a function
chebappxf <- function(fun,dims,intervals=NULL,...) {
  chebappx(evalongrid(fun,dims,intervals,...),intervals)
}

# interpolate on a non-Chebyshev grid. This is useful if you for some reason
# do not have the function values on a Chebyshev-grid, but on some other grid.  It comes at a cost,
# The interpolation may not be very good compared to the Chebyshev-one.

# val are the function values on a grid, an array of appropriate dimension
# in the order of expand.grid()

# if grid is unspecified, it is assumed that it is on a Chebyshev grid in [-1,1]
# If grid is specified, it is a list of vectors. The length of the list
# is the dimension of the grid. Each vector contains grid-points in increasing or decreasing order.
# val is assumed to be the function values on expand.grid(grid)
# The user-grid is mapped to Chebshev knots in [-1,1] by splinefun in pkg stats



chebappxg <- function(val,grid=NULL,mapdim=NULL) {
  # grid is a list of grid points. val is the values as in expand.grid(grid)
  # if grid is null it is assumed to be a chebyshev grid. The dimensions
  # must be present in val
  if(is.null(grid)) return(chebappx(val))
  if(is.null(dim(val))) dim(val) <- length(val)
  if(!is.list(grid) && length(grid) == length(val)) grid <- list(grid)
  if(prod(sapply(grid,length)) != length(val)) stop('grid size must match data length')
  dim(val) <- sapply(grid,length)

  # ok, grid is something like list(c(...),c(...),c(...))
  # create piecewise linear functions which maps grid-points to chebyshev grids


  intervals <- lapply(grid,function(x) c(min(x),max(x)))
  # create monotone splines with splinefun, method monoH.FC
  gridmaps <- mapply(splinefun, grid, chebknots(dim(val)), MoreArgs=list(method='monoH.FC'))
#  gridmaps <- mapply(polyh, chebknots(dim(val)), grid, MoreArgs=list(k=1))
  gridmap <- cmpfun(function(x) mapply(function(gm,x) gm(x),gridmaps,x))
  ch <- chebappx(val)
  structure(function(x) ch(gridmap(x)),arity=length(grid),domain=intervals,grid=grid)
}

chebappxgf <- function(fun, grid, ..., mapdim=NULL) {

  if(!is.list(grid)) grid <- list(grid)
  chebappxg(evalongrid(fun, ..., grid=grid),grid,mapdim)
}

# we can actually find the grid-maps for uniform grids.
# the Chebyshev knots are cos(pi*(j+0.5)/n) for j=0..n-1 These should
# map into the n grid points. These have distance 2/(n-1), starting in -1, ending in 1
# so they are -1 + 2*j/(n-1). After some manipulation, the function is:

ugm <- function(x,n) sin(0.5*pi*x*(1-n)/n)

ucappx <- function(val, intervals=NULL) {
  if(is.null(dim(val))) dim(val) <- length(val)
  dims <- dim(val)
  ch <- chebappx(val)
  if(is.null(intervals)) {
    gridmap <- function(x) mapply(function(xi,d) ugm(xi,d),x,dims)
  } else {
    # precompute interval mid points and inverse lengths
    md <- lapply(intervals,mean)
    ispan <- lapply(intervals, function(i) 2/diff(i))
    gridmap <- function(x) mapply(function(xi,mid,is,d) ugm(is*(xi-mid),d),x,md,ispan,dims)
  }
  return(structure(function(x) ch(gridmap(x)), arity=length(dims)))
}

ucappxf <- function(fun, dims, intervals=NULL,...) {
  if(is.null(intervals))
    return(ucappx(evalongrid(fun,...,grid=lapply(dims,function(d) seq(-1,1,length.out=d)))))
  if(is.numeric(intervals) && length(intervals) == 2) intervals <- list(intervals)
  return(
    ucappx(evalongrid(fun,...,
                      grid=mapply(function(d,i) seq(min(i),max(i),length.out=d),
                        dims, intervals,SIMPLIFY=FALSE)),
           intervals))
}

mlappx <- function(val, grid, ...) {
  if(is.numeric(grid)) grid <- list(grid)
  if(is.function(val)) val <- evalongrid(val,grid=grid,...)
  gl <- prod(sapply(grid,length))
  if(length(val) != gl)
    stop("length of values ",length(val)," do not match size of grid ",gl)
  function(x) .Call(C_evalmlip,grid,as.numeric(val),as.numeric(x))
}

havefftw <- function() .Call(C_havefftw)


polyh <- function(val, knots, k=2, ...) {
# Linear polyharmonic splines. Centres are columns in matrix knots. Function values in val.
# Quite slow for ncol(knots) > 3000 or so.  k=2 yields thin-plate splines.
# There exist faster evaluation methods for dimensions <= 4
# k < 0 yields Gaussian kernel splines, with sigma^2 = -1/k
# we compute r^2, so powers etc. are adjusted for that case
# Hmm, I should take a look at http://dx.doi.org/10.1016/j.jat.2012.11.008 for the unit ball
# perhaps some normalization should be added?
  if(is.null(dim(knots))) dim(knots) <- c(1,length(knots))
  if(is.function(val)) val <- apply(knots,2,val,...)
  N <- ncol(knots)
  M <- nrow(knots)
  ki <- k/2
  if(k < 0)
    phi <- local(cmpfun(function(r2) exp(k*r2)),list(k=k))
  else if(k %% 2 == 1) {
    phi <- local(cmpfun(function(r2) r2^ki), list(ki=ki))
  } else {
    ki <- as.integer(ki)-1L # trick to handle r2=0. Works because 0^0 = 1 in R
    phi <- local(cmpfun(function(r2) 0.5*r2^ki * log(r2^r2)),list(ki=ki))
  }

  sqnm <- apply(knots,2,function(x) sum(x^2))
  # would it be faster to apply phi only on lower tri?  Without crossprod?
  # and either fill in the upper tri, or tailor a solver? I think
  # evaluation of phi on the matrix is fast compared to creating the matrix,
  # though I haven't measured it. Could've done it in C, probably somewhat faster.
  # use abs when we call phi. The argument can't be negative, but numerically it can
  A <- phi(abs(-2*crossprod(knots) + sqnm + rep(sqnm,each=N)))
  A[diag(A)] <- 1
  V <- rbind(1,knots)
  mat <- cbind(rbind(A,V),rbind(t(V),matrix(0,M+1,M+1)))
  rhs <- c(val,rep(0,M+1))
  wv <- try(solve(mat,rhs))
  if(inherits(wv,'try-error')) {
    warning('Failed to fit exactly, fallback to least squares fit')
    wv <- lm.fit(mat,rhs)$coefficients
    wv[is.na(wv)] <- 0
  }
  w <- wv[1:N]
  v <- wv[(N+1):length(wv)]
  local(cmpfun(function(x) {
    if(is.vector(x) && length(x) == M) {
      sum(w*phi(abs(-2*crossprod(x,knots) + sqnm + sum(x^2)))) + sum(v*c(1,x))
    } else {
      if(is.null(dim(x))) dim(x) <- c(M,length(x)/M)
      if(nrow(x) != M) stop('spline was built for dimension ',M,' not ',nrow(x))
      apply(x,2,function(xx) sum(w*phi(abs(-2*crossprod(xx,knots) + sqnm + sum(xx^2)))) + sum(v*c(1,xx)))
    }
  }), list(w=w,v=v,knots=knots,phi=phi,sqnm=sqnm,M=M))
}


