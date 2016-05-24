##  Part of R package BDA
##  Copyright (C) 2009-2010 Bin Wang
##
##  Unlimited use and distribution (see LICENCE).

.nkde <- function(x,x0,h,truncate=FALSE)
{
    gpoints <- x0
    n = length(x)
    M <- length(gpoints)
    a <- gpoints[1L]
    b <- gpoints[M]
    if(missing(h))
      h <-   (1/(4*pi))^(1/10)*(243/(35*n))^(1/5)*sqrt(var(x))*15^(1/5)
    if(truncate) trun=1 else trun=0
    gcounts <- .Fortran(.F_linbin, as.double(x), as.integer(n),
             as.double(a), as.double(b), as.integer(M),
                        as.integer(trun), x=double(M))$x
    ## Compute kernel weights
    delta  <- (b - a)/(h * (M-1L))
    L <- min(floor(4./delta), M)
    if (L == 0)
      warning("Binning grid too coarse for current (small) bandwidth: consider increasing 'gridsize'")

    lvec <- 0L:L
    kappa <-  dnorm(lvec*delta)/(n*h)

    ## Now combine weight and counts to obtain estimate
    ## we need P >= 2L+1L, M: L <= M.
    P <- 2^(ceiling(log(M+L+1L)/log(2)))
    kappa <- c(kappa, rep(0, P-2L*L-1L), rev(kappa[-1L]))
    tot <- sum(kappa) * (b-a)/(M-1L) * n # should have total weight one
    gcounts <- c(gcounts, rep(0L, P-M))
    kappa <- fft(kappa/tot)
    gcounts <- fft(gcounts)
    list(x = gpoints, y = (Re(fft(kappa*gcounts, TRUE))/P)[1L:M])
}


.smkde <- function(x, bandwidth, from, to, gridsize=512L){
    stopifnot(gridsize > 10)
    nbin <- length(x$mids)
  
    F <- x$counts; X <- x$mids; n <- sum(F)
    A <- x$breaks[1:nbin]; B <- x$breaks[-1]
    x0 <- runif(n,rep(A,F),rep(B,F))

    ## to approximate the initial estimate of f(x) by simulation
    h = ifelse(missing(bandwidth), bw.nrd(x0), bandwidth)
    a <- ifelse(missing(from), min(x0), from)
    b <- ifelse(missing(to), max(x0), to)
    range.x <- c(a, b)
    
    N = length(X)
    B =  (B - A)/2; A = -B  # A and B are now the half binwdiths

    ## set grid to evaluate the pdf/cdf
    M <- 2^(ceiling(log2(gridsize)))
    gpoints <- seq(range.x[1L], range.x[2L], length = M)
 
    f0 = .nkde(x0,gpoints,h)$y;

    iter=500;
    out <- .Fortran(.F_iterfx, fx=as.double(f0),
                    as.double(gpoints), as.integer(M),
                    as.double(X), as.double(F),
                    as.integer(N), as.double(B),
                    h = as.double(h),iter=as.integer(iter))
    if(out$iter>=iter) warning("Fixed point not found!")
    y = out$fx
    sele1 = is.na(y) | !is.finite(out$fx)
    y[sele1] = 0.0

    structure(list(y = y, x = gpoints,
                   conf.level = NULL,
                   type="smoothKDE",
                   pars = out$h,
                   ucb=NULL, lcb=NULL,
                   call = match.call()
                   ),
              class = 'histosmooth')
}


