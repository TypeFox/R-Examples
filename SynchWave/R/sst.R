cwt_fw <- function(x, type, nv, dt, opt) {

    # From Synchrosqueezing Toolbox authored by Eugene Brevdo (http://www.math.princeton.edu/~ebrevdo/)
    # R Port by Dongik Jang (dongik.s.jang@gmail.com)

    if (is.null(opt) ) opt <- list()
    
    # options: symmeric, replicate, circular
    if ( is.null(opt$padtype) ) opt$padtype <- "symmetric"
    if ( is.null(opt$rpadded) ) opt$rpadded <- 0
    
    n <- length(x)
    
    # Pad x first
    tmp <- padsignal(x, opt$padtype)
    x <- tmp$xpad
    N <- tmp$nup 
    n1 <- tmp$n1
    n2 <- tmp$n2

    # Choosing more than this means the wavelet window becomes too short
    noct <- log2(N)-1
    if ( (noct <= 0) | (noct%%1 != 0) ) 
        stop("Data length is too short or is not power of 2. \n")
      
    if( (nv <= 0) | (nv%%1 != 0) ) 
        stop("Data length is less than 0 or number of voices (nv) is integer. \n")
    
    if( dt <= 0 ) 
        stop("The sampling period (dt) is less than 0. \n")
    
    if( any(is.na(x)) ) 
        stop("Data x has missing. \n")
    
    na <- noct*nv
    asc <- (2^(1/nv))^(1:na)
    Wx <- matrix(0, na, N)
    
    dWx <- Wx
    opt$dt <- dt
   
    xh <- fft(x)
    
    #  for each octave
    for (ai in 1:na) {
        a <- asc[ai]
        
        tmp <- wfilth(type, N, a, opt)
        psih <- tmp$psih
        dpsih <- tmp$dpsih
    
        dxcpsi = ifftshift( fft(dpsih * xh, inverse = TRUE)/length(xh) )
        dWx[ai, ] = dxcpsi

        xcpsi = ifftshift(fft(psih * xh, inverse = TRUE)/length(x) )
        Wx[ai, ] = xcpsi
    }
    
    # Shorten W to proper size (remove padding)
    if ( !as.logical(opt$rpadded) ) {
        Wx = Wx[ , (n1+1):(n1+n) ]
        dWx = dWx[ , (n1+1):(n1+n) ]
    }
    
    # Output a for graphing purposes, scale by dt
    asc = asc * dt
    list( Wx=Wx, asc=asc, dWx=dWx )
}

p2up <- function(n) {
  
    # From Synchrosqueezing Toolbox authored by Eugene Brevdo (http://www.math.princeton.edu/~ebrevdo/)
    # R Port by Dongik Jang (dongik.s.jang@gmail.com)
  
    # Calculates next power of 2, and left/right padding to center the
    # original n locations.
    #
    # Input:
    #   n: non-dyadic integer
    # Output:
    #   up: next power of 2
    #   n1: length on left
    #   n2: length on right
    #
    eps <- .Machine$double.eps
    up <- 2^( 1 + round(log2(n + eps)) )
    n1 <- floor((up - n)/2)
    n2 <- n1
    if ( (2*n1 + n) %% 2 == 1 )
        n2 <- n1 + 1
    
    list(up=up, n1=n1, n2=n2)
}

padsignal <- function(x, padtype) {

    # From Synchrosqueezing Toolbox authored by Eugene Brevdo (http://www.math.princeton.edu/~ebrevdo/)
    # R Port by Dongik Jang (dongik.s.jang@gmail.com)
  
    tmp <- p2up(length(x))
    
    xl <- padarray(x, tmp$n1, padtype, "pre")
    xr <- padarray(x, tmp$n2, padtype, "post")

    xpad = c(xl[1:tmp$n1], x, xr[(length(xr) - tmp$n2 + 1):length(xr)])
    list(xpad=xpad, nup=tmp$up, n1=tmp$n1, n2=tmp$n2)
}

padarray <- function(x, padsize, method=c("circular", "replicate", "symmetric"), direction=c("pre", "post")) {

    # This function is adapted from Matlab.
  
    n <- length(x)
    
    if (method[1] == "replicate") {
        out <- switch(direction[1], pre=c(rep(x[1], padsize), x), post=c(x, rep(x[n], padsize)))  
    } else if (method[1] == "symmetric") {
        out <- switch(direction[1], pre=c(x[padsize:1], x), post=c(x, x[n:(n-padsize+1)]))
    } else
        out <- switch(direction[1], pre=c(x[(n-padsize+1):n], x), post=c(x, x[1:padsize]))

    out
}

wfilth <- function(type, N, a=1, opt=NULL) {

    # From Synchrosqueezing Toolbox authored by Eugene Brevdo (http://www.math.princeton.edu/~ebrevdo/)
    # R Port by Dongik Jang (dongik.s.jang@gmail.com)
  
    # Outputs the FFT of the wavelet of family 'type' with parameters
    # in 'opt', of length N at scale a: (psi(-t/a))^.
    #
    # Note that the output is made so that the inverse fft of the
    # result is zero-centered in time.  This is important for
    # convolving with the derivative(dpsih).  To get the correct
    # output, perform an ifftshift.  That is,
    #   psi = ifftshift(ifft(psih));
    #   xfilt = ifftshift(ifft(fft(x) .* psih));
    #
    # Inputs:
    #   type: wavelet type (see help wfiltfn)
    #   N: number of samples to calculate
    #   a: wavelet scale parameter (default = 1)
    #   opt: wavelet options (see help wfiltfn)
    #     opt.dt: delta t (sampling period, default = 1)
    #             important for properly scaling dpsih
    #
    # Outputs:
    #   psih: wavelet sampling in frequency domain (for use in fft)
    #   dpsih: derivative of same wavelet, sampled in frequency domain (for fft)
    #   xi: associated fourier domain frequencies of the samples.

    if ( (log2(N) %% 1) != 0 )
        stop("N is not power of 2. /n")
    
    if ( is.null(a) ) a <- 1
    if ( is.null(opt) ) opt <- list()
    if ( is.null(opt$dt) ) opt$dt <- 1
    
    # When w_k = 1/N * [0, 1, ..., N]
    # The continuous-time frequency variables are xi_k = 2*pi*w_k,
    k <- 0:(N-1)
    xi <- rep(0, N)
    
    # Oppenheim et al., Eq. (10.1), (10.4)
    #   V[k] = V(e^{jw}) @ w=2*pi*k/N, k=0,1,...,N-1
    #   V(e^{jw}) = 1/T sum_{r=-inf}^inf V_c(j (w+2 pi r)/T)
    # where T is the sampling interval and V_c is the FT of v(t)
    #
    # Oppenheim et al., S. 10.2.2, P. 1
    #   Effects of spectral sampling, k=0..N/2, N/2+1..N
    #   are observing xi_k = 2*pi*k/(N*T) for k=0..N/2,
    #                        2*pi*(N-k)/(N*T) k=N/2+1,..,N
    
    xi[ 1:(N/2+1) ] <- 2 * pi/N * (0:(N/2))
    xi[ (N/2+2):N ] <- 2 * pi/N * ((-N/2+1):(-1))
    psihfn <- wfiltfn(type, opt)
    psih <- psihfn(a*xi)
    
    #Sometimes bump gives a NaN when it means 0
    if (type =="bump")
        psih[is.na(psih)] <- 0
    
    # (1/sqrt(a)*psi(n/a))^ = a/sqrt(a) * psi^(n) = sqrt(a) * psi^(n)
    # Mallat, Wavelet Tour 3rd ed, S. 4.3.3
    #
    # Also normalize appropriately
    psih <- psih * sqrt(a) / sqrt(2*pi)
    
    # Center around zero in the time domain
    psih <- psih * (-1)^k
    
    #Calculate the fourier transform of the derivative
    
    dpsih <- (1i * xi / opt$dt) * psih
        
    list(psih=psih, dpsih=dpsih, xi=xi)
}

wfiltfn <- function(type=c("gauss", "mhat", "cmhat", "morlet", "shannon", "hshannon", "hhat", "hhhat", "bump"), opt=NULL) {

    # From Synchrosqueezing Toolbox authored by Eugene Brevdo (http://www.math.princeton.edu/~ebrevdo/)
    # R Port by Dongik Jang (dongik.s.jang@gmail.com)
  
    types <- c("gauss", "mhat", "cmhat", "morlet", "shannon", "hshannon", "hhat", "hhhat", "bump")
    if( !any(types %in% type[1]) )
        stop( paste("Unknown wavelet type: ", type[1], "\n"))
    
    if(type[1] == "gauss") {
        # Similar to morlet, but can control bandwidth
        # can be used with synsq for large enough mu/s ratio
        s <- ifelse(is.null(opt$s), 1/6, opt$s)
        mu <- ifelse(is.null(opt$mu), 2, opt$mu)
        psihfngauss <- function(w, mu, s) {
             2^(-(w-mu)^2/(2*s^2))
        }
        out <- function(w) do.call("psihfngauss", list(w=w, mu=mu, s=s))
    } else if(type[1] == "mhat") { # mexican hat
        s <- ifelse(is.null(opt$s), 1, opt$s)
        psihfnmhat <- function(w, s) {
            -sqrt(8)*s^(5/2)*pi^(1/4)/sqrt(3)*w^2*exp(-s^2*w^2/2)
        }
        out <- function(w) do.call("psihfnmhat", list(w=w, s=s))
    } else if(type[1] == "cmhat") {
        # complex mexican hat: hilbert analytic function of sombrero
        # can be used with synsq
        s <- ifelse(is.null(opt$s), 1, opt$s)
        mu <- ifelse(is.null(opt$mu), 1, opt$mu)
        
        psihfnshiftcmhat <- function(w, s) {
            2*sqrt(2/3)*pi^(-1/4)*s^(5/2)*w^2*exp(-s^2*w^2/2)*as.numeric(w>=0)
        }
        psihfncmhat <- function(w, mu, s) {
            psihfnshiftcmhat(w-mu, s)
        }
        out <- function(w) do.call("psihfncmhat", list(w=w, mu=mu, s=s))  
    } else if(type[1] == "morlet") {
        # can be used with synsq for large enough s (e.g. >5)
        mu <- ifelse(is.null(opt$mu), 2*pi, opt$mu)
        
        cs <- (1+exp(-mu^2) - 2*exp(-3/4*mu^2))^(-1/2)
        ks <- exp(-(1/2)*mu^2)
        psihfnmorlet <- function(w, mu, cs, ks) {
             cs*pi^(-1/4)*( exp(-(1/2)*(mu-w)^2) - ks*exp(-(1/2)*w^2) )
        }
        out <- function(w) do.call("psihfnmorlet", list(w=w, mu=mu, cs=cs, ks=ks))
    } else if(type[1] == "shannon") {
        psihfnshannon <- function(w) {
             exp(-1i*w/2) * as.numeric( (abs(w)>=pi) & (abs(w)<=2*pi))
        }
        out <- function(w) do.call("psihfnshannon", list(w=w))
    } else if(type[1] == "hshannon") {
        # hilbert analytic function of shannon transform
        # time decay is too slow to be of any use in synsq transform

        mu <- ifelse(is.null(opt$mu), 0, opt$mu)
        
        psihfnshifthshannon <- function(w) {
            exp(-1i*w/2) * as.numeric((w>=pi) & (w<=2*pi)) * ( 1 + sign(w) )
        }
        psihfnhshannon <- function(w, mu) {
            psihfnshifthshannon(w-mu)
        }
        out <- function(w) do.call("psihfnhshannon", list(w=w, mu=mu))   
    } else if(type[1] == "hhat") {
        # hermitian hat
        psihfnhhat <- function(w) {
             2/sqrt(5) * pi^(-1/4) * w * (1+w) * exp(-(1/2)*(w^2))
        }
        out <- function(w) do.call("psihfnhhat", list(w=w))
    } else if(type[1] == "hhhat") {
        # % hilbert analytic function of hermitian hat
        # can be used with synsq
        mu <- ifelse(is.null(opt$mu), 5, opt$mu)
        
        psihfnshifthhhat <- function(w) {
            2/sqrt(5) * pi^(-1/4) * ( w*(1+w)*exp(-(1/2)*(w^2)) ) * (1+sign(w))
        }
        psihfnhhhat <- function(w, mu) {
            psihfnshifthhhat(w-mu)
        }
        out <- function(w) do.call("psihfnhhhat", list(w=w, mu=mu))       
    } else if(type[1] == "bump") {
        # complex mexican hat: hilbert analytic function of sombrero
        # can be used with synsq

        s <- ifelse(is.null(opt$s), 1, opt$s)
        mu <- ifelse(is.null(opt$mu), 5, opt$mu)
         
        psihfnorig <- function(w) {
            exp(-1/(1-w^2)) * ( abs(w)<1 )
        }
        psihfnbump <- function(w, mu, s) {
            psihfnorig((w-mu)/s)
        }
        out <- function(w) do.call("psihfnbump", list(w=w, mu=mu, s=s))      
    }
    
    out
}

fftshift <- function(x) {

    # This function is adapted from Matlab.
  
    if (!is.vector(x))
        stop("Argument 'x' must be a vector.")

    m <- length(x)
    p <- ceiling(m/2)
    x[c((p+1):m, 1:p)]
}

ifftshift <- function(x) {

    # This function is adapted from Matlab.
  
    if (!is.vector(x))
        stop("Argument 'x' must be a vector.")

    m <- length(x)
    p <- floor(m/2)
    x[c((p+1):m, 1:p)]
}

cwt_iw <- function(Wx, type, opt=NULL) {

    # From Synchrosqueezing Toolbox authored by Eugene Brevdo (http://www.math.princeton.edu/~ebrevdo/)
    # R Port by Dongik Jang (dongik.s.jang@gmail.com)
  
    if (is.null(opt) ) opt <- list()
    
    na <- nrow(Wx)
    n <- ncol(Wx)
    
    tmp <- p2up(n)
    N <- tmp$up
    n1 <- tmp$n1
    n2 <- tmp$n2
    Wxp <- matrix(0, na, N)
    
    #TODO - do we want to pad the wavelet representation here?  Or not
    # cut off Wx in cwt_fw and then use that to reconstruct?    
    Wxp[, (n1+1):(n1+n)] <- Wx
    Wx <-  Wxp
    rm(Wxp)

    # Following the same value in cwt_fw
    noct <- log2(N)-1
    nv <- na/noct
    asc <- (2^(1/nv))^(1:na)
    
    if ( noct%%1 != 0 )
        stop("Data length is too short or is not power of 2. \n")
    
    if( (nv <= 0) | (nv%%1 != 0) )
        stop("The number of voices (nv) is not positive integer. \n")
    
    # Find the admissibility coefficient Cpsi
    Cpsi <- cwt_adm(type, opt)
        
    x <- rep(0, N)

    for ( ai in 1:na) {
        a <- asc[ai]
        Wxa <- Wx[ai, ]
        psih <- wfilth(type, N, a, opt)$psih
        
        # Convolution theorem here
        Wxah <- fft(Wxa)
        xah <- Wxah * psih
        xa <- ifftshift( fft(xah, inverse = TRUE)/length(xah) )

        x <- x + xa/a  
    }
    
    # Take real part and normalize by log_e(a)/Cpsi
    x <- log( 2^(1/nv) )/Cpsi * Re(x)
    
    #Keep the unpadded part
    x <- x[(n1+1):(n1+n)]
    
    x
}

cwt_adm <- function(type, opt=NULL) {

    # From Synchrosqueezing Toolbox authored by Eugene Brevdo (http://www.math.princeton.edu/~ebrevdo/)
    # R Port by Dongik Jang (dongik.s.jang@gmail.com)
  
    # Calculate cwt admissibility constant int(|f(w)|^2/w, w=0..inf) as
    # per Eq. (4.67) of [1].
    #
    # 1. Mallat, S., Wavelet Tour of Signal Processing 3rd ed.
    
    if (is.null(opt) ) opt <- list()
    
    if ( type == "sombrero" ) {
        s <- ifelse(is.null(opt$s), 1, opt$s)
        Cpsi <- (4/3)*s*sqrt(pi)    
    } else if ( type == "shannon" ) {
        Cpsi <- log(2)
    } else {
        psihfn <- wfiltfn(type, opt)
            
        f <- function(x, psihfn) {
            xx <- x/(1-x)
            x2x <- xx^2
            val <- ( Conj( psihfn(x2x) ) * psihfn(x2x) )/x2x
            val <-  2*xx* val / (1 - x)^2
        }
            
        Cpsi <- quadgk(f, 0, 1, psihfn = psihfn)    
    }
    
    # Normalize
    Cpsi <- Cpsi / (4*pi)
    
    Cpsi
}

quadgk <- function(f, a, b, tol = .Machine$double.eps^0.5, ...) {

    # This function is from pracma packake 1.4.5 authored by Hans W Borchers (hwborchers at googlemail.com).
    #
    # Adaptive Gauss-Kronrod Quadrature
    # Adaptive version of the (7, 15)-point Gauss-Kronrod quadrature formula,
    # where in each recursion the error is taken as the difference between these
    # two estimated integrals.
    #
    # Input:
    # f : integrand as function, may have singularities at the endpoints.
    # a, b : endpoints of the integration interval.
    # tol : relative tolerence.
    # ... : Additional parameters to be passed to the function f
    #
    # Output:
    # Value of the integration. The relative error should be of the same
    # order of magnitude as the relative tolerance (or much smaller).
    
    stopifnot(is.numeric(a), length(a) == 1,
              is.numeric(b), length(b) == 1)
    eps <- .Machine$double.eps

    fun <- match.fun(f)
    f <- function(x) fun(x, ...)
    
    if (a == b)     return(0)
    else if (a > b) return(-1 * quadgk(f, b, a, tol = tol))

    # Nodes and weights for Gauss-Kronrod (7, 15)
    n15 <- c(-0.9914553711208126, -0.9491079123427585, -0.8648644233597691,
             -0.7415311855993944, -0.5860872354676911, -0.4058451513773972,
             -0.2077849550078985,  0.0,                 0.2077849550078985,
              0.4058451513773972,  0.5860872354676911,  0.7415311855993944,
              0.8648644233597691,  0.9491079123427585,  0.9914553711208126)
    n7  <- c(-0.9491079123427585, -0.7415311855993944, -0.4058451513773972,
               0.0,
               0.4058451513773972, 0.7415311855993944,  0.9491079123427585)
    w15 <- c(0.02293532201052922, 0.06309209262997855,  0.1047900103222502,
             0.1406532597155259,  0.1690047266392679,   0.1903505780647854,
             0.2044329400752989,  0.2094821410847278,   0.2044329400752989,
             0.1903505780647854,  0.1690047266392679,   0.1406532597155259,
             0.1047900103222502,  0.06309209262997855,  0.02293532201052922)
    w7  <- c(0.1294849661688697,  0.2797053914892767,   0.3818300505051189,
             0.4179591836734694,
             0.3818300505051189,  0.2797053914892767,   0.1294849661688697)

    .gkadpt <- function(f, a, b, tol = tol) {
        # use nodes and weights from the environment
        x15 <- 0.5 * ((b - a) * n15 + b + a)
        x7  <- 0.5 * ((b - a) * n7  + b + a)
        Q7  <- sum(w7  * f(x7))  * (b-a)/2
        Q15 <- sum(w15 * f(x15)) * (b-a)/2

        if (!is.finite(Q7) || !is.finite(Q15)) {
            warning("Infinite or NA function value encountered.")
            return(Q15)
        } else if (abs(Q15 - Q7) < tol) {
            return(Q15)
        } else if (abs(b-a) < 16*eps) {
            warning("Minimum step size reached; singularity possible.")
            return(Q2)
        } # else

        Q2 <- .gkadpt(f, (a+b)/2, b, tol = tol)
        Q1 <- .gkadpt(f, a, (a+b)/2, tol = tol)

        return(Q1 + Q2)
    }

    # start the recursive procedure
    .gkadpt(f, a, b, tol = tol)
}
                    
est_riskshrink_thresh <- function (Wx, nv) {
  
    # From Synchrosqueezing Toolbox authored by Eugene Brevdo (http://www.math.princeton.edu/~ebrevdo/)
    # R Port by Dongik Jang (dongik.s.jang@gmail.com)
  
    # Estimate the RiskShrink hard thresholding level.
    #
    # Implements Defn. 1 of Sec. 2.4 in [1], using the suggested
    # noise estimator from the discussion "Estimating the noise level"
    # in that same section.
    #
    # 1. Donoho, D.L.; I.M. Johnstone (1994), "Ideal spatial adaptation by
    #    wavelet shrinkage," Biometrika, vol 81, pp. 425–455.
    #
    # Inputs:
    #  Wx: wavelet transform of a signal, see help cwt_fw
    #  opt: options structure used for forward wavelet transform.
    #
    # Output:
    #   the RiskShrink hard threshold estimate
  
    # TODO, error if opt has no nv, or no opt.
    
    na <- nrow(Wx)
    n <- ncol(Wx)

    Wx_fine <- abs( Wx[ 1:nv, ] )
    
    madval <- median(abs(Wx_fine - median(Wx_fine)))
    
    #out <- sqrt( 2*log(n) ) * mad(Wx_fine) * 1.4826
    out <- sqrt( 2*log(n) ) * madval * 1.4826
}

synsq_cwt_fw <- function(tt, x, nv=16, opt=NULL) {
  
    # From Synchrosqueezing Toolbox authored by Eugene Brevdo (http://www.math.princeton.edu/~ebrevdo/)
    # R Port by Dongik Jang (dongik.s.jang@gmail.com)
  
    if ( is.null(opt) )
        opt <- list()
        
    # Choose some default values
    
    # Wavelet type
    if ( is.null(opt$type) ) opt$type <- "morlet"
    # Display bugging information?
    if ( is.null(opt$disp) ) opt$disp <- 0
    # Direct or numerical differntiation? (default: direct)
    if ( is.null(opt$dtype) ) opt$dtype <- 1
    
    # Calculate sampling period, assuming regular spacing
    dt = tt[2] - tt[1]
    
    # Check for uniformity of spacing
    if ( any( diff(tt, difference=2)/(tt[length(tt)] - tt[1]) > 1e-5) ) 
        stop("time vector tt is not uniformly sampled. \n")
    
    # Calculate the wavelet transform - padded via symmetrization

    N <- length(x)
    tmp <- p2up(N)
    Nup <- tmp$up
    n1 <- tmp$n1
    n2 <- tmp$n2

    if ( as.logical(opt$dtype) ){
        # Calculate derivative directly in the wavelet domain before taking
        # wavelet transform.
        
        # Don't return padded data'
        opt$rpadded <- 0
        
        tmp <- cwt_fw(x, opt$type, nv, dt, opt)
        Wx <- tmp$Wx
        asc <- tmp$asc
        dWx <- tmp$dWx

        # Approximate instantaneous frequency
        w <- cwt_freq_direct(Wx, dWx, opt)
        rm(dWx)
    } else {
        # Calculate derivative numerically after calculating wavelet
        # transorm.  This requires less memory and is more accurate at
        # small a.
        
        # Return padded data
        opt$rpadded <- 1
        
        tmp <- cwt_fw(x, opt$type, nv, dt, opt) 
        Wx <- tmp$Wx
        asc <- tmp$asc
        Wx <- Wx[ , (n1-4+1):(n1+N+4) ]
        
        # Approximate instantaneous frequency
        w <- cwt_freq(Wx, dt, opt)
    }
    
    # Calculate the synchrosqueezed frequency decomposition
    tmp <- synsq_cwt_squeeze(Wx, w, tt, nv, opt)
    
    if ( !as.logical(opt$dtype) ) {
        Wx <- Wx[ , (4+1):(ncol(Wx) - 4) ]
        w <- w[ , (4+1):(ncol(w) - 4) ]
        tmp$Tx <- tmp$Tx[ , (4+1):(ncol(tmp$Tx) - 4) ]
    }
    
    list(Tx=tmp$Tx, fs=tmp$fs, Wx=Wx, asc=asc, w=w)
}

cwt_freq_direct <- function (Wx, dWx, opt) {

    # From Synchrosqueezing Toolbox authored by Eugene Brevdo (http://www.math.princeton.edu/~ebrevdo/)
    # R Port by Dongik Jang (dongik.s.jang@gmail.com)
  
    # Calculate the demodulated FM-frequency at each (scale,time) pair:
    #   w(a,b) = Im( (1/2pi) * d/db [ Wx(a,b) ] / Wx(a,b) )
    # Uses direct differentiation by calculating dWx/db in frequency
    # domain (the secondary output of cwt_fw, see help cwt_fw)
    #
    # This is the analytic implementation of Eq. (7) of [1].
    #
    # 1. E. Brevdo, N.S. Fučkar, G. Thakur, and H-T. Wu, "The
    # Synchrosqueezing algorithm: a robust analysis tool for signals
    # with time-varying spectrum," 2011.
    #
    #
    # Input:
    #  Wx: wavelet transform of x (see help cwt_fw)
    #  dWx: samples of time derivative of wavelet transform of x (see help cwt_fw)
    #  opt: options struct,
    #    opt.gamma: wavelet threshold (default: sqrt(machine epsilon))
    #
    # Output:
    #  w: demodulated FM-estimates, size(w) = size(Wx)
    eps <- .Machine$double.eps
    
    if (is.null(opt) ) opt <- list()
    # epsilon from Daubechies, H-T Wu, et al.
    # gamma from Brevdo, H-T Wu, et al.
    if ( is.null(opt$gamma) )  opt$gamma <- sqrt(eps)
    
    # Calculate inst. freq for each ai, normalize by (2*pi) for
    # true freq.
    
    w <- Im( (dWx / Wx) / (2*pi) )
    w[ abs(Wx) < opt$gamma ]  <- NA
    w
}

cwt_freq <- function (Wx, dt, opt) {
  
    # From Synchrosqueezing Toolbox authored by Eugene Brevdo (http://www.math.princeton.edu/~ebrevdo/)
    # R Port by Dongik Jang (dongik.s.jang@gmail.com)
  
    # Calculate the demodulated FM-frequency at each (scale,time) pair:
    #   w(a,b) = Im( (1/2pi) * d/db [ Wx(a,b) ] / Wx(a,b) )
    # Uses numerical differentiation (1st, 2nd, or 4th order).
    #
    # This is a numerical differentiation implementation of Eq. (7) of
    # [1].
    #
    # 1. E. Brevdo, N.S. Fučkar, G. Thakur, and H-T. Wu, "The
    # Synchrosqueezing algorithm: a robust analysis tool for signals
    # with time-varying spectrum," 2011.
    
    # Input:
    #  Wx: wavelet transform of x (see help cwt_fw)
    #  dt: delta t, the sampling period (e.g. t(2)-t(1))
    #  opt: options struct,
    #    opt.dorder: differences order (values: 1, 2, 4.  default = 4)
    #    opt.gamma: wavelet threshold (default: sqrt(machine epsilon))
    #
    # Output:
    #  w: demodulated FM-estimates, size(w) = size(Wx)
    eps <- .Machine$double.eps
    
    if (is.null(opt) ) opt <- list()
    
    # Order of differentiation (1, 2, or 4)
    if( is.null(opt$dorder) ) opt$dorder <- 4
    
    # epsilon from Daubechies, H-T Wu, et al.
    # gamma from Brevdo, H-T Wu, et al.
    if ( is.null(opt$gamma) )  opt$gamma <- sqrt(eps)
    
    # Section below replaced by the following mex code
    w <- diff_so(Wx, dt, opt$dorder)

    # switch opt.dorder
    #   case 1
    #     w = [Wx(:, 2:end) - Wx(:, 1:end-1), ...
    #          Wx(:, 1)-Wx(:, end)];
    #     w = w / dt;
    #   case 2
    #     % Append for differencing
    #     Wxr = [Wx(:, end-1:end) Wx Wx(:, 1:2)];
    #     % Calculate 2nd order forward difference
    #     w = -Wxr(:, 5:end) + 4*Wxr(:, 4:end-1) - 3*Wxr(:, 3:end-2);
    #     w = w / (2*dt);
    #   case 4
    #     % Centered difference with fourth order error
    #     % Append for differencing
    #     Wxr = [Wx(:, end-1:end) Wx Wx(:, 1:2)];

    #     % Calculate 4th order central difference
    #     w = -Wxr(:, 5:end);
    #     w = w + 8*Wxr(:, 4:end-1);
    #     w = w - 8*Wxr(:, 2:end-3);
    #     w = w + Wxr(:, 1:end-4);
    #     w = w / (12*dt);
    #   otherwise
    #     error('Differentiation order %d not supported', opt.dorder);
    # end
    
    w[ abs(Wx) < opt$gamma ] <- NA
    
    # Calculate inst. freq for each ai, normalize by (2*pi) for
    # true freq.
    
    w <- Re(-1i * w / Wx) / (2*pi)
    
    w
}

diff_so <- function(Wx, dt, dorder=c(1,2,4)) {
    na <- nrow(Wx)
    N <- ncol(Wx)
    if( any( c(1,2,4) %in% dorder[1] ) ) {
        out <- .Fortran("diff_w", Wx=as.complex(Wx), na=as.integer(na), N=as.integer(N), 
                        dt=as.double(dt), dorder=as.integer(dorder[1]),
                        out = as.complex(nrow(Wx)*ncol(Wx)) )$out
        dim(out) <- dim(Wx)
    } else {
        stop(paste("Differentiation order", dorder[1],  "not supported. \n") )
    }
    out
}

synsq_cwt_squeeze <- function(Wx, w, tt, nv, opt) {
  
    # From Synchrosqueezing Toolbox authored by Eugene Brevdo (http://www.math.princeton.edu/~ebrevdo/)
    # R Port by Dongik Jang (dongik.s.jang@gmail.com)
  
    # Calculates synchrosqueezed transform of f on a logarithmic
    # scale.  Used internally by synsq_cwt_fw.
    #
    # Input:
    #   Wx: wavelet transform of x
    #   w: estimate of frequency at locations in Wx (see synsq_cwt_fw code)
    #   tt: time vector
    #   nv: number of voices
    #   opt: options struct, not currently used
    #
    # Output:
    #   Tx: synchrosqueezed output
    #   fs: associated frequencies
    #
    # Note the multiplicative correction term f in synsq_cwt_squeeze_mex (and in
    # the matlab equivalent code commented out), required due to the fact that
    # the squeezing integral of Eq. (2.7), in, [1], is taken w.r.t. dlog(a).
    # This correction term needs to be included as a factor of Eq. (2.3), which
    # we implement here.
    #
    # A more detailed explanation is available in Sec. III of [2].
    # Specifically, this is an implementation of Sec. IIIC, Alg. 1.
    # Note the constant multiplier log(2)/nv has been moved to the
    # inverse of the normalization constant, as calculated in synsq_adm.m
    #
    # 1. I. Daubechies, J. Lu, H.T. Wu, "Synchrosqueezed Wavelet Transforms: a
    # tool for empirical mode decomposition", 2010.
    #
    # 2. E. Brevdo, N.S. Fučkar, G. Thakur, and H-T. Wu, "The
    # Synchrosqueezing algorithm: a robust analysis tool for signals
    # with time-varying spectrum," 2011.
    eps <- .Machine$double.eps
    dt <- tt[2] - tt[1]
    dT <- tt[length(tt)] - tt[1]
    
    # Maximum measurable frequency of data
    #fM = 1/(4*dt); % wavelet limit - tested
    fM <- 1/(2*dt) # standard
    
    # Minimum measurable frequency, due to wavelet
    fm <- 1/dT # really
    #fm = 1/(2*dT); % standard
    
    na <- nrow(Wx)
    N <- ncol(Wx)
    asc <-  (2^(1/nv))^(1:na)
    das <- c(1, diff(asc))
    lfm <- log2(fm)
    lfM <- log2(fM);
    fs <- 2^seq(lfm, lfM, ,na)
    # dfs = c(fs[1], diff(fs))
    
    # Harmonics of diff. frequencies but same
    # magnitude have same |Tx|
    dfs <- rep(1, length(fs))
    
    #The rest of the computation is performed efficiently in MEX
    frobnorm <- sqrt(sum(Mod(Wx)^2)/2)

    #if ( norm(Wx, "F") < eps ){
    if ( frobnorm < eps ) {
        Tx <- matrix(0, nrow(Wx), ncol(Wx))
    } else {
        out <- .Fortran("synsq_cwt_squeeze",  Wx=as.complex(Wx), na=as.integer(na), N=as.integer(N), w=as.double(w),
                        sc=as.double(asc), fs=as.double(fs), dfs=as.double(dfs),  lfm=as.double(lfm), lfM=as.double(lfM),
                        out=complex(nrow(Wx)*ncol(Wx)), NAOK=TRUE)$out
        dim(out) <- dim(Wx)
        Tx <- c(1/nv) * out
    }
    list( Tx=Tx, fs=fs)
}

synsq_cwt_iw <- function(Tx, fs, opt=NULL) {
  
    # From Synchrosqueezing Toolbox authored by Eugene Brevdo (http://www.math.princeton.edu/~ebrevdo/)
    # R Port by Dongik Jang (dongik.s.jang@gmail.com)
  
    if ( !is.matrix(Tx) ) Tx <- as.matrix(Tx)
    
    if ( is.null(opt$type) ) opt$type <- "morlet"

    # Find the admissibility coefficient Cpsi
    Css <- synsq_adm(opt$type, opt)

    # Integration
    # Due to linear discretization of integral in log(fs), this becomes
    # a simple normalized sum  
    x <- Re( 1/Css * c(apply(Tx, 2, sum)) )
    x
}

synsq_adm <- function(type, opt=NULL) {

    # From Synchrosqueezing Toolbox authored by Eugene Brevdo (http://www.math.princeton.edu/~ebrevdo/)
    # R Port by Dongik Jang (dongik.s.jang@gmail.com)
  
    # Calculate the Synchrosqueezing admissibility constant, the term
    # R_\psi in Eq. 3 of [1].  Note, here we multiply R_\psi by the
    # inverse of log(2)/nv (found in Alg. 1 of that paper).
    #
    # 1. E. Brevdo, N.S. Fu��kar, G. Thakur, and H-T. Wu, "The
    # Synchrosqueezing algorithm: a robust analysis tool for signals
    # with time-varying spectrum," 2011.
    #
    # Uses numerical integration.
    #
    # Input:
    #   type: type of wavelet (see help wfiltfn)
    #   opt: options structure (wavelet parameters, see help wfiltfn)
    #
    # Output:
    #   Css: proportional to 2*int(conj(f(w))/w, w=0..inf)
    
    # switch type
    # case 'sombrero',
    #   if ~isfield(opt,'s'), s = 1; else s = opt.s; end
    #   Cpsi = (4/3)*s*sqrt(pi);
    # case 'shannon',
    #   Cpsi = log(2);
    # otherwise
  
    psihfn <- wfiltfn(type, opt)
    f <- function(x, psihfn) {
        xx <- x/(1-x)
        x2x <- xx^2
        val <- Conj(psihfn(x2x))/x2x
        val <-  2*xx* val / (1 - x)^2
    }
    #Css <- integrate(f, lower = 0, upper = 20, psihfn = psihfn, rel.tol = .Machine$double.eps^0.4)$value
    Css <- quadgk(f, 0, 1, psihfn = psihfn)
    
    # Normalization constant, due to logarithmic scaling in wavelet
    # transform
    Css <- Css / (sqrt(2*pi)*2*log(2))
    Css
}

synsq_filter_pass <- function(Tx, fs, fm, fM) {

    # From Synchrosqueezing Toolbox authored by Eugene Brevdo (http://www.math.princeton.edu/~ebrevdo/)
    # R Port by Dongik Jang (dongik.s.jang@gmail.com)
  
    if (length(fm) == 1) fm <- rep(fm, ncol(Tx))
    if (length(fM) == 1) fM <- rep(fM, ncol(Tx))

    if (length(fm) != ncol(Tx))
        stop("fm does not match t-length of Tx. \n")

    if (length(fs) != nrow(Tx))
        stop("fs does not match Tx. \n")

    fmi <- rep(0, length(fm))
    fMi <- rep(0, length(fm))
    Txf <- matrix(0, nrow(Tx), ncol(Tx))
    for (ti in 1:length(fm)) {
        fmi[ti] <- which(fs >= fm[ti])[1]
        fMi[ti] <- rev(which(fs <= fM[ti]))[1]
        Txf[fmi[ti]:fMi[ti], ti] <- Tx[fmi[ti]:fMi[ti], ti]
    }

    list(Txf=Txf, fmi=fmi, fMi=fMi)
}


curve_ext_multi <- function(Tx, fs, nc, lambda = 1e3, clwin = 4) {

    # From Synchrosqueezing Toolbox authored by Eugene Brevdo (http://www.math.princeton.edu/~ebrevdo/)
    # R Port by Dongik Jang (dongik.s.jang@gmail.com)
  
    eps <- .Machine$double.eps
    na <- nrow(Tx)
    N <- ncol(Tx)
    Cs <- matrix(0, N, nc)
    Es <- rep(0, nc)
    
    for (ni in 1:nc) {
        tmp <- curve_ext(Tx, fs, lambda)
        Cs[ , ni] <- tmp$FreqOut
        Es[ni] <- tmp$EnergyOut
    
        # Remove this curve from the representation
        # Max/min frequencies for each time step
        
        for (cli in 0:clwin) {
            fb <- Cs[, ni] - cli
            fe <- Cs[, ni] + cli
            fb[fb < 1] <- 1
            fe[fe > na] <- na
 
            Tx[ (0:(N-1))*na + fb ] <-  sqrt(eps)
            if ( cli > 0 )
                Tx[(0:(N-1))*na + fe] <- sqrt(eps)
        }
    }
    list(Cs=Cs, Es=Es)
}

curve_ext <- function(Tx, fs, lambda=0) {
 
    N <- ncol(Tx)
    na <- nrow(Tx)
    
    out <- .Fortran("curve_ext", as.complex(Tx), as.integer(na), as.integer(N), as.double(fs), as.double(lambda),
                          FreqOut=double(N), EnergyOut=double(1))[c("FreqOut", "EnergyOut")]
    out
}

curve_ext_recon <- function(Tx, fs, Cs, opt=NULL, clwin=4) {
  
    # From Synchrosqueezing Toolbox authored by Eugene Brevdo (http://www.math.princeton.edu/~ebrevdo/)
    # R Port by Dongik Jang (dongik.s.jang@gmail.com)
  
    Nc <- ncol(Cs)
    na <- nrow(Tx)
    N <- ncol(Tx)
    
    if ( na != length(fs) )
        stop("The number of row of Tx is the same as length of fs. \n")
    
    # Reconstruct the signal for each curve
    xrs <- matrix(0, N, Nc)
    for ( ni in 1:Nc ) {
        cmask <- matrix(0, na, N)
        cmask[ (0:(N-1)) * na + Cs[, ni]] <- 1
        cmask <- imdilate(cmask, matrix(1, clwin, 1)) 
        xrs[, ni] <- synsq_cwt_iw(Tx*cmask, fs, opt)
    }   
    xrs
}
                
imdilate <- function(IM, SE) {

    # This function is adapted from Matlab.
    # dilate image
  
    if ( !is.matrix(IM) ) IM <- as.matrix(IM) 
    if ( !is.matrix(SE) ) SE <- as.matrix(SE) 
    m <- nrow(IM)
    n <- ncol(IM)
    ms <- nrow(SE)
    ns <- ncol(SE)
    
    out <- .Fortran("imdilate", as.integer(IM), as.integer(m), as.integer(n), 
                     as.integer(SE), as.integer(ms), as.integer(ns),
                     out=integer(m*n))$out
    dim(out) <- c(m, n)
    out          
}

                    
