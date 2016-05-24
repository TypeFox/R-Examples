## Copyright (C) 2013 Marius Hofert
##
## This program is free software; you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation; either version 3 of the License, or (at your option) any later
## version.
##
## This program is distributed in the hope that it will be useful, but WITHOUT
## ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
## FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
## details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.


##' @title Iteratively oppositely order a column of a matrix with respect to
##'        the sum of all other colums
##' @param x matrix
##' @return matrix as x after one iteration over all columns
##' @author Marius Hofert
##' @note Iteration here means that the jth column already involves the changed
##'       values in the (j-1)th column
oppositely.order <- function(x)
{
    stopifnot(is.matrix(x), (d <- ncol(x)) >= 2)
    for(j in 1:d) x[,j] <- sort(x[,j], decreasing=TRUE)[order( order(rowSums(x[,-j])) )]
    x
}

##' @title Computing lower and upper bounds for the (smallest or largest) VaR
##' @param alpha confidence level
##' @param N tail discretization parameter
##' @param qmargins list of marginal quantile functions
##' @param bound character string indicating the VaR bound to
##'        be approximated (largest (default) or smallest)
##' @param verbose logical indicating whether progress information is displayed
##' @return lower and upper bounds for the (lower or upper) VaR bound
##' @author Marius Hofert
##' @note The 'steps' refer to the rearrangement algorithms on p. 9 in
##'       Embrechts, Puccetti, Rueschendorf
VaRbound <- function(alpha, N, qmargins, bound = c("upper", "lower"), verbose = FALSE)
{
    ## checks
    stopifnot(0 < alpha, alpha < 1, N == as.integer(N), N >= 1,
              (d <- length(qmargins)) >= 2, is.list(qmargins))
    bound <- match.arg(bound)

    ## step 2
    q <- switch(bound, # length N+1
                "upper" = alpha+(1-alpha)*0:N/N,
                "lower" = alpha*0:N/N,
                stop("wrong argument 'bound'"))
    Z <- sapply(1:d, function(j) qmargins[[j]](q)) # (N+1, d) matrix
    X <- Z[1:N,] # (N, d) matrix; lower bound for the VaR bound
    Y <- Z[2:(N+1),] # (N, d) matrix; upper bound for the VaR bound

    ## step 3
    X. <- apply(X, 2, sample)
    Y. <- apply(Y, 2, sample)

    ## steps 4 and 5
    while (TRUE) {
        X.. <- oppositely.order(X.)
        if(verbose)
            switch(bound,
                   "upper" = cat("Min of row-wise sums of X: ", min(rowSums(X..)),"\n"),
                   "lower" = cat("Max of row-wise sums of X: ", max(rowSums(X..)),"\n"),
                   stop("wrong argument 'bound'"))
        if(all(X.. == X.)) break else X. <- X..
    }

    ## step 6
    while (TRUE) {
        Y.. <- oppositely.order(Y.)
        if(verbose)
            switch(bound,
                   "upper" = cat("Min of row-wise sums of Y: ", min(rowSums(Y..)),"\n"),
                   "lower" = cat("Max of row-wise sums of Y: ", max(rowSums(Y..)),"\n"),
                   stop("wrong argument 'bound'"))
        if(all(Y.. == Y.)) break else Y. <- Y..
    }

    ## step 7
    switch(bound,
           "lower" = c(lower=max(rowSums(X..)), upper=max(rowSums(Y..))), # lower and upper bounds for the smallest VaR
           "upper" = c(lower=min(rowSums(X..)), upper=min(rowSums(Y..))), # lower and upper bounds for the largest VaR
           stop("wrong argument 'bound'"))
}
