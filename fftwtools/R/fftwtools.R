##     The fftwtools R package
##     fftw2d and general mvfftw tools in R
##     Copyright (C) 2013 Karim Rahim 

##     This file is part of the fftwtools package for R.

##     The fftwtools package is free software: you can redistribute it and
##     /or modify
##     it under the terms of the GNU General Public License as published by
##     the Free Software Foundation, either version 2 of the License, or
##     any later version.

##     The fftwtools package is distributed in the hope that it will be 
##     useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
##     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##     GNU General Public License for more details.

##     You should have received a copy of the GNU General Public License
##     along with Foobar.  If not, see <http://www.gnu.org/licenses/>.

##     If you wish to report bugs please contact the author. 
##     karim.rahim@gmail.com
##     Jeffery Hall, Queen's University, Kingston Ontario
##     Canada, K7L 3N6


##generic function to call fftw3
## we recommend the FFTW library for setting the plan
fftw <- function(data, inverse=0, HermConj=1, n=NULL) {
    res <- NULL
    if(inverse==0) {
        if(!is.complex(data)) {
            res <- fftw_r2c(data, HermConj=HermConj)
        } else {
            res <- fftw_c2c(data, inverse=0)
        }
    } else {
        if(HermConj==0 && is.complex(data) && is.numeric(n)) {
            res <- fftw_c2r(data, HermConj=HermConj, n=n)
        } else{
            res <- fftw_c2c(data, inverse=1)
        }
    }
    return(res)
}

##generic 2d fft
fftw2d <- function(data, inverse=0, HermConj=1) {
    res <- NULL
    if(inverse==0) {
        if(!is.complex(data)) {
            res <- fftw_r2c_2d(data, HermConj=HermConj)
        } else {
            res <- fftw_c2c_2d(data, inverse=0)
        }
    } else {
        res <- fftw_c2c_2d(data, inverse=1)
    }
    return(res)
}

##generic function to call multicol fftw3
mvfftw <- function(data, inverse=0, HermConj=1, n=NULL, fftplanopt=0) {
    res <- NULL
    if(inverse==0) {
        if(!is.complex(data)) {
            res <- mvfftw_r2c(data, HermConj=HermConj, fftplanopt=fftplanopt)
        } else {
            res <- mvfftw_c2c(data, inverse=0, fftplanopt=fftplanopt)
        }
    } else {
        if(HermConj==0 && is.complex(data) && is.numeric(n)) {
            res <- mvfftw_c2r(data, HermConj=HermConj, n=n,
                              fftplanopt=fftplanopt)
        } else {
            res <- mvfftw_c2c(data, inverse=1, fftplanopt=fftplanopt)
        }
    }
    return(res)
}


fftw_r2c <- function(data, HermConj=1) {

    n <- length(data)

    nc <- NULL
    if(HermConj ==1) {
        nc <-  n
    } else {
        nc <- as.integer(n/2) +1
    }    

    out <- .C("fft_r2c", as.integer(n), as.double(data),
              res=complex(nc),
              as.integer(HermConj))
    return(out$res)
}


fftw_c2r <- function(data, HermConj=1, n=NULL) {

    ##when HermConj == 0, we need the length of the original
    ## data set.

    len <- length(data)
    nc <- NULL
    if(HermConj==0) {
        stopifnot(is.numeric(n))
        nc <- as.integer(n/2) + 1
        stopifnot(nc == len)
    } else {
        n <- len
        nc <- as.integer(n/2) + 1
    }
      
    ## only pass what is needed to the C function
    out <- .C("fft_c2r", as.integer(n), as.complex(data[1:nc]),
              res=double(n))
    
    return(out$res)
}


fftw_c2c <- function(data, inverse=0) {

    n <- length(data)
    
    out <- .C("fft_c2c", as.integer(n), as.complex(data),
              res=complex(n), as.integer(inverse))
    
    return(out$res)
}

mvfftw_r2c <- function(data, HermConj=1, fftplanopt=0) {

    data <- as.matrix(data)
    n <- dim(data)[1]
    m <- dim(data)[2]  ## ncol

    nc <- as.integer(n/2) +1
       
    
    out <- .C("mvfft_r2c", as.integer(n), as.integer(m),
              as.double(data),
              res=matrix(as.complex(0), nc, m),
              as.integer(fftplanopt))

    res <- as.matrix(out$res)

    ## moving this routine to C requires making space at the end each column,
    ## but the data is stored in a contiguous format when returned from
    ## fftw, and it is easier to manipulate this in R.
    if(HermConj==1) {
        ## if( n %% 2 == 0) {
        ##     res <- rbind(res, Conj(res[(nc-1):2,]))
        ## } else {
        ##     res <- rbind(res, Conj(res[nc:2,]))
        ## }
        idx <- n-nc+1
        res <- rbind(res, Conj(res[idx:2,]))
    }

    return(res)
}


mvfftw_c2r <- function(data, HermConj=1, n=NULL, fftplanopt=0) {

    ##HermConj is used to determine the true length, n,  of the
    ##input, fftw does not use the the Hermitian conjugate
    ## in this case the Hermatian conjugate must be stripped from the data.

    data <- as.matrix(data)
    ##n <- length(data[,1])
    len <- dim(data)[1]
    m <- dim(data)[2]
    nc <- NULL

      
    if(HermConj==0) {
        stopifnot(is.numeric(n))
        nc <- as.integer(n/2) + 1
        stopifnot(nc == len)
    } else{
        n <- len
        nc <- as.integer(n/2) + 1
    }
        
    out <- .C("mvfft_c2r", as.integer(n), as.integer(m),
              as.complex(data[1:nc,]),
              res=matrix(as.double(0), n, m), as.integer(fftplanopt))

    return(out$res)
}

mvfftw_c2c <- function(data, inverse=0, fftplanopt=0) {

    data <- as.matrix(data)
    n <- dim(data)[1]
    m <- dim(data)[2]
    
    out <- .C("mvfft_c2c", as.integer(n), as.integer(m),
              as.complex(data),
              res=matrix(as.complex(0), n, m),
              as.integer(inverse), as.integer(fftplanopt))
    
    return(out$res)
}

fftw_r2c_2d <- function(data, HermConj=1) {

    ## can be done with two mvfft's t(mvfft(t(mvfft(a))))
    ## == mvfft(t(mvfft(t(a))))

    data <- as.matrix(data)

    nR <- dim(data)[1]
    nC <- dim(data)[2]
    nRc <- floor(nR/2) +1
    
    isEven <- 1 - (nR %% 2)
    idxRowAppend <- NULL
    
    if(isEven) {
        idxRowAppend <- (nRc -1):2
    } else {
        idxRowAppend <- nRc:2
    }
    
    ##correct for the fact the c call is column-major

    out <- .C("fft_r2c_2d", as.integer(nC), as.integer(nR),
              as.double(data), res=matrix(as.complex(0), nRc , nC))

    res <- as.matrix(out$res)
    if(HermConj==1 && nR > 2) {
        if(length(idxRowAppend) == 1) {
            ##If the number of rows to append is one, R
            ##will create a column vector for cbind unless
            ##we transpose.

            res <- rbind(res, Conj(cbind(res[2,1],
                                         t(res[2,nC:2]))))
         } else {
             ##With the exception of the first row, the
             ## resulting matrix is Hermatian --conjugate symmetric         
             res <- rbind(res, Conj(cbind(res[idxRowAppend,1],
                                          res[idxRowAppend,nC:2])))
         }
    }
    
    return(res)
}

fftw_c2c_2d <- function(data, inverse=0) {

    ## can be done with two mvfft's t(mvfft(t(mvfft(a))))
    ## == mvfft(t(mvfft(t(a))))

    data <- as.matrix(data)

    nR <- dim(data)[1]
    nC <- dim(data)[2]

    ##we correct for the fact the c call is column-major

    out <- .C("fft_c2c_2d", as.integer(nC), as.integer(nR),
              as.complex(data),
              res=matrix(as.complex(0), nR , nC),
              as.integer(inverse))

    return(out$res)
}
