###############################################################################
## Time Series Model Specification for Air Pollution and Health
## Copyright (C) 2004-2010, Roger D. Peng <rpeng@jhsph.edu>
##     
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2 of the License, or
## (at your option) any later version.
## 
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
###############################################################################

## Copied/modified from 'split<-.data.frame'
"split<-.matrix" <- function (x, f, value) {
        ix <- split(seq(NROW(x)), f)
        n <- length(value)
        j <- 0
        for (i in ix) {
                j <- j%%n + 1
                x[i, ] <- value[[j]]
        }
        x
}

## Return a length(v) x (length(k) * (length(k)-1))/2 matrix of lag
## interactions
xLag <- function(v, k, group = NULL) {
        if(length(k) < 2)
                stop("'k' should be > 1")
        xlag.f <- function(x) {
                lagmat <- Lag(x, k, group = NULL)
                M <- matrix(nrow = length(x), ncol = length(k) * (length(k) - 1) / 2)
                nc <- NCOL(lagmat)
                
                Mlist <- lapply(seq(to = nc - 1), function(i) { 
                        m <- lagmat[, i] * lagmat[, -seq(to = i), drop = FALSE]
                        colnames(m) <- paste(i - 1, seq(i, nc - 1), sep = ":")
                        m
                })
                do.call("cbind", Mlist)
        }
        D <- if(!is.null(group)) {
                groupLag <- tapply(v, group, xlag.f)

                ## Need to get the column names
                cn <- colnames(groupLag[[1]])

                lagmat <- matrix(nrow = length(v), ncol = length(k) * (length(k) - 1) / 2)
                split(lagmat, group) <- groupLag
                colnames(lagmat) <- cn
                lagmat
        }
        else
                xlag.f(v)
        ## colnames(D) <- as.character(k)
        drop(D)
}

Lag <- function(v, k, group = NULL) {
        stopifnot(length(k) > 0)
        v <- as.numeric(v)
        
        if(max(abs(k)) >= length(v))
                stop("largest lag in 'k' must be less than 'length(v)'")

        lag.f <- function(x) { 
                lagmat <- matrix(nrow = length(x), ncol = length(k))
                n <- length(x)
                
                for(i in seq(along = k)) {
                        lag <- k[i]
                        
                        if(lag > 0)            
                                lagmat[, i] <- c(rep(NA, lag), x[1:(n - lag)])
                        else if(lag < 0)
                                lagmat[, i] <- c(x[(-lag + 1):n], rep(NA, -lag))
                        else
                                lagmat[, i] <- x
                }
                lagmat
        }
        lagmat <- if(!is.null(group)) {
                groupLag <- tapply(v, group, lag.f)
                lagmat <- matrix(nrow = length(v), ncol = length(k))
                split(lagmat, group) <- groupLag
                lagmat
        }
        else
                lag.f(v)
        colnames(lagmat) <- as.character(k)
        drop(lagmat)
}

dLag <- function(v, k, group = NULL) {
        stopifnot(length(k) > 0)
        v <- as.numeric(v)

        dlag.f <- function(x) {
                lagmat <- Lag(x, k, group = NULL)

                if(length(k) == 1)
                        return(matrix(lagmat))
                use <- complete.cases(lagmat)
                xqr <- qr(lagmat[use, ], LAPACK = FALSE)
                
                if(xqr$rank < ncol(lagmat))  
                        stop("problem with rank of lag matrix")
                d <- qr.R(xqr)[1, 1]
                X <- qr.Q(xqr, Dvec = rep(d, length(k)))
                
                if(any(!use)) {
                        M <- matrix(NA, nrow = length(use), ncol = length(k))
                        M[use, ] <- X
                }
                else
                        M <- X
                M
        }
        D <- if(!is.null(group)) {
                groupdLag <- tapply(v, group, dlag.f)
                lagmat <- matrix(nrow = length(v), ncol = length(k))
                split(lagmat, group) <- groupdLag
                lagmat
        }
        else
                dlag.f(v)
        colnames(D) <- as.character(k)
        drop(D)
}

runMean <- function(v, lags = 0, group = NULL, filter = NULL) {
        if(is.null(filter)) {
                first.lag <- min(lags)
                lags <- sort(lags) - first.lag + 1
                filter <- rep(0, max(lags) - min(lags))
                filter[lags] <- 1
                filter <-  filter / sum(filter)
        }
        rm.f <- function(x) { 
                x <- as.vector(Lag(x, k = first.lag, group = NULL))
                f <- filter(x, filter = filter, sides = 1, circular = FALSE)
                as.numeric(f)
        }
        if(!is.null(group)) {
                groupFilter <- tapply(v, group, rm.f)
                ## unname(unlist(groupFilter))
                unname(unsplit(groupFilter, group))
        }
        else
                rm.f(v)
}

## This function does not allow arbitrary filters, just running averages
runMean2 <- function(v, lags = 0, group = NULL, na.rm = FALSE) {
        lagMat <- Lag(v, lags, group)
        rowMeans(lagMat, na.rm = na.rm)
}

harmonic <- function(x, nfreq, period, intercept = FALSE) {
        stopifnot(nfreq > 0)
        pi <- base::pi  ## Just in case someone has redefined pi!
        x <- as.numeric(x)

        N <- seq(0, nfreq - 1)
        k <- 2^N * 2 * pi / period
        M <- outer(x, k)
        sinM <- apply(M, 2, sin)
        cosM <- apply(M, 2, cos)
        if(!intercept) 
                cbind(sinM, cosM)
        else
                cbind(1, sinM, cosM)
}





######################################################################

## Not exported yet
## constrDL <- function(x, f, by = NULL) {
##         stopifnot(is.list(f))
##         x <- as.numeric(x)
##         np <- length(f)
##         lagmat <- matrix(nrow = length(x), ncol = np)
## 
##         for(i in seq(along = f)) {
##                 mark <- f[[i]]
##                 if(length(mark) > 1)
##                         lagmat[, i] <- runMean(x, len = mark[2]-mark[1]+1, lag = mark[1])
##                 else if(length(mark) == 1)
##                         lagmat[, i] <- Lag(x, mark)
##                 else
##                         stop("Something wrong")
##         }
##         lagmat
## }


## adjustTimeDF <- function(object, dfseq, timeVar = "time", smoothType = "ns") {
##         stopifnot(inherits(object, "tsModel"))
##         stopifnot(length(dfseq) > 0)
##         
##         confounderFormula <- as.formula(object$call$confounders)
## 
##         stopifnot(timeVar %in% all.vars(confounderFormula))
## 
##         origTimeVar <- grep(timeVar, attr(terms(confounderFormula), "term.labels"),
##                             value = TRUE)
##         stopifnot(length(origTimeVar) == 1)
##         
##         timeVec <- paste("ns(", timeVar, ", ", dfseq, ")", sep = "")
##         results <- vector("list", length = length(dfseq))
## 
##         for(i in seq(along = dfseq)) {
##                 newTimeVar <- timeVec[i]
##                 newFormula <- as.formula(paste("~ . -", origTimeVar, "+", newTimeVar))
##                 results[[i]] <- update(object, confounders = newFormula)
##         }
##         if(length(results) == 1)
##                 results <- results[[1]]
##         structure(results, class = "adjustTimeDF", dfseq = dfseq)
## }



