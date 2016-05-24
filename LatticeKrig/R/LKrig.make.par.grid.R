# LatticeKrig  is a package for analysis of spatial data written for
# the R software environment .
# Copyright (C) 2012
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

LKrig.make.par.grid <- function(par.grid = NULL, LKinfo = NULL) {    
    # if par.grid is missing find all the information for the LKinfo list
    if (is.null(par.grid)) {
        par.grid <- list()
        par.grid$a.wght <- list(LKinfo$a.wght)
        par.grid$alpha <- list(LKinfo$alpha)
        par.grid$llambda <- ifelse(is.na(LKinfo$lambda), 0, log(LKinfo$lambda))
        return(par.grid)
    }
    if (is.null(par.grid$llambda)) {
        par.grid$llambda <- NA
    }
   # if needed create alpha from nu
    if( !is.null(par.grid$nu)){
        nlevel<- LKinfo$nlevel
        M<- length( par.grid$nu)
        alpha<- matrix( NA,nrow=M, ncol=LKinfo$nlevel)
        for ( k in 1:M){
          alphaTemp <- exp(-2 * (1:nlevel) * par.grid$nu[k])
          alpha[k,] <- alphaTemp/sum(alphaTemp)
        }
        par.grid$alpha<- alpha
     } 
    # if needed create alpha from gamma
    if (!is.null(par.grid$gamma)) {
        if (!is.matrix(par.grid$gamma)) {
            par.grid$gamma <- cbind(par.grid$gamma)
        }
        par.grid$alpha <- cbind(1, exp(par.grid$gamma))
        for (j in 1:nrow(par.grid$gamma)) {
            par.grid$alpha[j, ] <- par.grid$alpha[j, ]/sum(par.grid$alpha[j, 
                ])
        }
    }
    #convert alpha to a list format 
    if (is.matrix(par.grid$alpha)) {
        M <- nrow(par.grid$alpha)
        temp.list <- list()
        for (k in (1:M)) {
            temp.list <- c(temp.list, list(par.grid$alpha[k, 
                ]))
        }
        par.grid$alpha <- temp.list
    }
    #convert a.wght to list format
    if (is.matrix(par.grid$a.wght)) {
        M <- nrow(par.grid$a.wght)
        temp.list <- list()
        for (k in (1:M)) {
            temp.list <- c(temp.list, list(par.grid$a.wght[k, 
                ]))
        }
        par.grid$a.wght <- temp.list
    }
    # some checks
    NG <- length(par.grid$alpha)
    if (length(par.grid$llambda) != NG) {
        stop("llambda values not same length as alpha")
    }
    if (length(par.grid$a.wght) != NG) {
        stop("a.wght values not same length as alpha")
    }
    return(par.grid)
}
