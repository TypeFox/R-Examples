#' @keywords internal
spm_hrf <- function(TR, P) {
    ## R port of spm8 functions needed for deconvolution
    ## spm_hrf - Returns a hemodynamic response function as a difference of gammas (canonical)
    ##
    ## USAGE:
    ##   hrf = spm_hrf(TR,[p])
    ##
    ##INPUT:
    ## TR   - scan repetition time (seconds)
    ## P    - 1x7 vector of parameters for the response function (two gamma functions)
    ##                                                 (default value, in seconds)
    ## P[1] - delay of response (relative to onset)    (6)
    ## P[2] - delay of undershoot (relative to onset)  (16)
    ## P[3] - dispersion of response                   (1)
    ## P[4] - dispersion of undershoot                 (1)
    ## P[5] - ratio of response to undershoot          (6)
    ## P[6] - onset (seconds)                          (0)
    ## P[7] - length of kernel (seconds)               (32)
    ##
    ## OUTPUT:
    ##  $hrf    - hemodynamic response function
    ##  $p      - parameters of the response function
    ##_______________________________________________________________________
    ## Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

    ## Karl Friston
    ## $Id: spm_hrf.m 387 2005-12-17 18:31:23Z klaas $

    ## default parameters
    ##-----------------------------------------------------------------------
    fMRI_T = 16 #microtime resolution is 1/16th of TR

    p = c(6, 16, 1, 1, 6, 0, 32)

    if (!missing(P)) {
        p[1:length(P)] = P
    }

    ## modelled hemodynamic response function - {mixture of Gammas}
    ##-----------------------------------------------------------------------
    dt    = TR/fMRI_T
    u     = 0:(p[7]/dt) - p[6]/dt #sampling grid of HRF in microtime units (e.g., 0:256 for default params and 2s TR)

    ## MH: eliminated use of spm Gpdf in favor of built-in gamma PDF in R. Checked that this yields identical
    ## results, although the third parameter of spm_Gpdf is really rate, not shape.
    hrf   = dgamma(u, shape=p[1]/p[3], rate=dt/p[3]) - dgamma(u, shape=p[2]/p[4], rate=dt/p[4])/p[5]

    hrf   = hrf[0:(p[7]/TR)*fMRI_T + 1] #subsample hrf back onto TR grid
    hrf   = hrf/sum(hrf) #Unit normalize hrf

    #return hrf  explicitly as column vector because other functions cbind additional basis functions
    return(list(hrf=matrix(hrf, ncol=1), p=p))
}

spm_get_bf <- function(xBF) {
    ## fills in basis function structure
    ## FORMAT [xBF] = spm_get_bf(xBF)
    ##
    ## xBF.dt      - time bin length {seconds}
    ## xBF.name    - description of basis functions specified
    ## xBF.length  - window length (seconds)
    ## xBF.order   - order
    ## xBF.bf      - Matrix of basis functions
    ##
    ## xBF.name  'hrf'
    ##       'hrf (with time derivative)'
    ##       'hrf (with time and dispersion derivatives)'
    ##       'Fourier set'
    ##       'Fourier set (Hanning)'
    ##       'Gamma functions'
    ##       'Finite Impulse Response'}
    ##
    ## (any other specification will default to hrf)
    ##__________________________________________________________________________
    ##
    ## spm_get_bf prompts for basis functions to model event or epoch-related
    ## responses.  The basis functions returned are unitary and orthonormal
    ## when defined as a function of peri-stimulus time in time-bins.
    ## It is at this point that the distinction between event and epoch-related
    ## responses enters.
    ##_______________________________________________________________________
    ## Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

    ## Karl Friston
    ## $Id: spm_get_bf.m 3934 2010-06-17 14:58:25Z guillaume $

    ## length of time bin
    ##--------------------------------------------------------------------------
    if (missing(xBF)) {
        xBF <- list()
        while (is.null(xBF$dt) || is.na(xBF$dt)) {
            xBF$dt <- readline('time bin for basis functions in seconds (.0625): ')
            if (nchar(xBF$dt) == 0L) {
                xBF$dt <- .0625  #default
            } else {
                xBF$dt <- as.numeric(xBF$dt)
            }            
        }
    }
    dt   = xBF$dt

    ## assemble basis functions
    ##==========================================================================

    ## model event-related responses
    ##--------------------------------------------------------------------------
    if (!exists('name', where=xBF)) {
        mnames <- c("hrf", "hrf (with time derivative)", "hrf (with time and dispersion derivatives)",
                    "Fourier set", "Fourier set (Hanning)", "Gamma functions", "Finite Impulse Response")
                    
        while (is.null(xBF$name) || is.na(xBF$name) || nchar(xBF$name) == 0L) {
            cat("Hemodynamic Basis functions...\n")
            cat(paste0("  (1) ", mnames[1], "\n  (2) ", mnames[2], "\n  (3) ", mnames[3], 
                       "\n  (4) ", mnames[4], "\n  (5) ", mnames[5], "\n  (6) ", mnames[6], "\n  (7) ", mnames[7], "\n"))
            msel <- as.numeric(readline("Select basis set (1-7): "))
            if (is.na(msel) || msel < 1 || msel > 7) {
                xBF$name <- NULL
            } else {
                xBF$name <- mnames[msel]
            }
        }
    }

    ## get order and length parameters
    ##--------------------------------------------------------------------------
    if (xBF$name %in% c("Fourier set", "Fourier set (Hanning)", "Gamma functions", "Finite Impulse Response")) {
        while (is.null(xBF$length) || is.na(xBF$length) || xBF$length < 1) {
            xBF$length <- readline("window length in seconds (32): ")
            if (nchar(xBF$length) == 0L) { xBF$length <- 32 } else { xBF$length <- as.numeric(xBF$length) }  #default to 32s
        }

        while (is.null(xBF$order) || is.na(xBF$order) || xBF$order < 1) {            
            xBF$order <- readline("order (4): ")
            if (nchar(xBF$order) == 0L) { xBF$order <- 4 } else { xBF$order <- as.numeric(xBF$order) }  #default to fourth-order
        }
    }
    h   = xBF$order
    l   = xBF$length
    
    ## create basis functions
    ##--------------------------------------------------------------------------
    if (xBF$name %in% c("Fourier set", "Fourier set (Hanning)")) {
        pst   = seq(from=0, to=l, by=xBF$dt)
        pst   = pst/max(pst)

        ## hanning window
        ##----------------------------------------------------------------------
        if (xBF$name == "Fourier set (Hanning)") {
            g  = (1 - cos(2*pi*pst))/2
        } else {
            g  = rep(1, length(pst))
        }

        ## zeroth and higher Fourier terms
        ##----------------------------------------------------------------------
        bf    = g
        for (i in 1:h) {
            bf = cbind(bf, g * sin(i*2*pi*pst))
            bf = cbind(bf, g * cos(i*2*pi*pst))
        }
    } else if (xBF$name == "Gamma functions") {
        pst   = seq(from=0, to=l, by=xBF$dt)
        bf    = spm_gamma_bf(pst,h)
    } else if (xBF$name == "Finite Impulse Response") {
        bin   = l/h
        bf    = kronecker(diag(h), rep(1, round(bin/dt)))
    } else if (xBF$name == "NONE") {
        bf = 1
    } else {
        ## canonical hemodynamic response function
        ##----------------------------------------------------------------------
        hrf            = spm_hrf(dt)
        bf             = hrf$hrf
        p              = hrf$p

        ## add time derivative
        ##----------------------------------------------------------------------
        if (grepl("time", xBF$name, ignore.case = TRUE)) {
            dp     = 1
            p[6]   = p[6] + dp
            D      = (bf[,1] - spm_hrf(dt,p)$hrf)/dp
            bf     = cbind(bf, D)
            p[6]   = p[6] - dp

            ## add dispersion derivative
            ##------------------------------------------------------------------
            if (grepl("dispersion", xBF$name, ignore.case = TRUE)) {
                dp    = 0.01
                p[3]  = p[3] + dp
                D     = (bf[,1] - spm_hrf(dt,p)$hrf)/dp
                bf    = cbind(bf, D)
            }
        }

        ## length and order
        ##----------------------------------------------------------------------
        xBF$length = dim(bf)[1]*dt
        xBF$order  = dim(bf)[2]         
    }

    xBF$bf <- spm_orth(bf)
    return(xBF)

}


## compute Gamma functions
##--------------------------------------------------------------------------
spm_gamma_bf <- function(u,h) {
    ## returns basis functions used for Volterra expansion
    ## FORMAT bf = spm_gamma_bf(u,h)
    ## u   - times {seconds}
    ## h   - order
    ## bf  - basis functions (mixture of Gammas)
    ##__________________________________________________________________________
    
    u     = as.vector(u)
    bf    = c()
    for (i in 2:(1 + h)) {
        m   = 2^i
        s   = sqrt(m)
        bf  = cbind(bf, dgamma(u, shape=(m/s)^2, rate=m/s^2)) #replace spm_Gpdf in favor of R built-in dgamma
    }
   
    return(bf)
}

spm_orth <- function(X, OPT) {
    ## Recursive Gram-Schmidt orthogonalisation of basis functions
    ## FORMAT X = spm_orth(X,OPT)
    ##
    ## X   - matrix
    ## OPT - 'norm' for Euclidean normalisation
    ##     - 'pad'  for zero padding of null space [default]
    ##
    ## Serial orthogonalisation starting with the first column
    ##
    ## Reference:
    ## Golub, Gene H. & Van Loan, Charles F. (1996), Matrix Computations (3rd
    ## ed.), Johns Hopkins, ISBN 978-0-8018-5414-9.
    ##__________________________________________________________________________
    ## Copyright (C) 2002-2012 Wellcome Trust Centre for Neuroimaging

    ## Karl Friston
    ## $Id: spm_orth.m 4708 2012-04-02 11:28:14Z guillaume $

    if (missing(OPT)) { OPT = "pad" }

    ##-Recursive Gram-Schmidt orthogonalisation
    ##--------------------------------------------------------------------------
    ## sw = warning('off','all'); ##suppress warnings?
    d     = dim(X)
    n     = d[1]; m <- d[2]
    X     = X[, apply(X, 2, function(col) { any(col != 0.0) }), drop=FALSE] ##Only retain non-zero columns
    rankX = qr(X)$rank

    x = j = NULL
    
    tryCatch({
        x     <<- X[,1, drop=FALSE] #assign in global environment
        j     <<- 1

        if (ncol(X) > 1L) {
            for (i in 2:ncol(X)) {
                D = X[,i, drop=FALSE]
                D = D - x%*%(MASS::ginv(x) %*% D)
                if (norm(D,"1") > exp(-32)) {
                    x                <<-  cbind(x, D)
                    j[length(j) + 1] <<- i                    
                }
                if (length(j) == rankX) { break }
            }
        }
    }, error=function(e) {
        print(e)
        x     <<- matrix(0, nrow=n, ncol=0)
        j     <<- c()        
    })

    ## and normalisation, if requested
    ##--------------------------------------------------------------------------
    if (OPT == "pad") {
        X      = matrix(0, nrow=n, ncol=m)
        X[,j]  = x
    } else {
        ## Euclidean norm each column in X        
        X <- apply(x, 2, function(col) { if (any(x != 0.0)) { col/sqrt(sum(col^2)) } else { col } })
    }

    return(X)
}