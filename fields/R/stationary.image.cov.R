# fields  is a package for analysis of spatial data written for
# the R software environment .
# Copyright (C) 2016
# University Corporation for Atmospheric Research (UCAR)
# Contact: Douglas Nychka, nychka@ucar.edu,
# National Center for Atmospheric Research, PO Box 3000, Boulder, CO 80307-3000
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with the R software environment if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
# or see http://www.r-project.org/Licenses/GPL-2    
stationary.image.cov <- function(ind1, ind2, Y, cov.obj = NULL, 
    setup = FALSE, grid, M = NULL, N = NULL, cov.function="stationary.cov",delta=NULL, cov.args=NULL, ...) {
    #
    # if cov object is missing then create
    # basically need to enlarge domain and find the FFT of the
    # covariance
    #
    cov.args<-c( cov.args, list(...))
    if (is.null(cov.obj)) {
        dx <- grid$x[2] - grid$x[1]
        dy <- grid$y[2] - grid$y[1]
        m <- length(grid$x)
        n <- length(grid$y)
        #
        # determine size of padding
        # default is twice domain and will then yeild exact results
        # delta indicates that covariance is zero beyond a distance delta
        # so using a smaller grid than twice domain will stil give exact results.
        if(!is.null(delta)){
           M<- ceiling(m + 2*delta/dx)
           N<- ceiling(n + 2*delta/dy)
         }
        if (is.null(M)) 
            M <- (2 * m)
        if (is.null(N)) 
            N <- (2 * n)
        xg <- make.surface.grid(list((1:M) * dx, (1:N) * dy))
        center <- matrix(c((dx * M)/2, (dy * N)/2), nrow = 1, 
            ncol = 2)
        #
        # here is where the actual covariance form is used
        # note passed arguments from call for parameters etc.
        #
        out<- do.call(cov.function, c(cov.args, list(x1 = xg, x2 = center)))  
        # check if this is a sparse result and if so expand to full size
        if( class( out)=="spam"){
           out <- spam2full(out)
         }
        # coerce to a matrix (image)
           out<- matrix( c(out), nrow = M, ncol = N)
        print( dim( out))
        temp <- matrix(0, nrow = M, ncol = N)
        #
        # a simple way to normalize. This could be avoided by
        # translating image from the center ...
        #
        temp[M/2, N/2] <- 1
        wght <- fft(out)/(fft(temp) * M * N)
        #
        # wght is the discrete FFT for the covariance suitable for fast
        # multiplication by convolution.
        #
        cov.obj <- list(m = m, n = n, grid = grid, N = N, M = M, 
            wght = wght, call = match.call())
        if (setup) {
            return(cov.obj)
        }
    }
    temp <- matrix(0, nrow = cov.obj$M, ncol = cov.obj$N)
    if (missing(ind1)) {
        temp[1:cov.obj$m, 1:cov.obj$n] <- Y
        Re(fft(fft(temp) * cov.obj$wght, inverse = TRUE)[1:cov.obj$m, 
            1:cov.obj$n])
    }
    else {
        if (missing(ind2)) {
            temp[ind1] <- Y
        }
        else {
            temp[ind2] <- Y
        }
        #
        # as promised this is a single clean step
        #
        Re(fft(fft(temp) * cov.obj$wght, inverse = TRUE)[ind1])
    }
}
