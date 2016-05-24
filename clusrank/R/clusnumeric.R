################################################################################
##
##   R package clusrank by Mei-Ling Ting Lee, Jun Yan, and Yujing Jiang
##   Copyright (C) 2015
##
##   This file is part of the R package clusrank.
##
##   The R package clusrank is free software: you can redistribute it and/or
##   modify it under the terms of the GNU General Public License as published
##   by the Free Software Foundation, either version 3 of the License, or
##   (at your option) any later version.
##
##   The R package clusrank is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##
##   You should have received a copy of the GNU General Public License
##   along with the R package reda. If not, see <http://www.gnu.org/licenses/>.
##
################################################################################

#' The Wilcoxon Signed Rank Test for Clustered Data
#'
#' Performs one-sample Wilcoxon test on vectors of data using
#' large sample.
#'
#' @note This function is able to deal with data with
#' clusterentitical or variable cluster size. When the data
#' is unbalanced, adjusted signed rank statistic is used.
#' Ties are dropped in the test.
#' @examples
#' data(crsd)
#' cluswilcox.test(z, cluster = id, data = crsd)
#' data(crsdUnb)
#' cluswilcox.test(z, cluster = id, data = crsdUnb)
#' @author Yujing Jiang
#' @references
#' Bernard Rosner, Robert J. Glynn, Mei-Ling Ting Lee(2006)
#' \emph{The Wilcoxon Signed Rank Test for Paired Comparisons of
#'  Clustered Data.} Biometrics, \bold{62}, 185-192.
#' @describeIn cluswilcox.test numeric interface for signed rank test.
#' @importFrom  stats complete.cases
#' @export

cluswilcox.test.numeric <- function(x, y = NULL,
                                    cluster = NULL,
                                    data = parent.frame(),
                                    alternative = c("two.sided", "less", "greater"),
                                    mu = 0, permutation = FALSE,
                                    n.rep = 500, ...) {
  ## Process the input arguments before feeding them to
  ## signed rank test . Assign a class (better to
  ## be S4) to the processed arguments for them to be
  ## sent to the corresponding functions.
  ##
  ## Inputs:
  ##   The same as cluswilcox.test.
  ##   x: numeric vector of data values. Non-finite
  ##     (e.g., infinite or missing) values will be omitted.
  ##
  ##
  ##   y: an optional numeric vector of data values:
  ##     as with x non-finite values will be omitted.
  ##
  ##
  ##   cluster:  an integer vector. Cluster cluster.If not provclustered,
  ##     assume there is no cluster.
  ##
  ##   data: an optional matrix or data frame
  ##     (or similar: see model.frame) containing the variables.
  ##     By default the variables are taken from environment(formula).
  ##
  ##
  ##   alternative: a character string specifying the
  ##     alternative hypothesis, must be one of
  ##    "two.sclustered" (default), "greater" or "less".
  ##     You can specify just the initial letter.
  ##
  ##   mu: a number specifying an optional parameter
  ##       used to form the null hypothesis.
  ##
  ##   paired: a logical indicating whether you want a paired test.
  ##
  ## permuation:

  METHOD <- "Wilcoxon signed rank test for clutered data"

  pars <- as.list(match.call()[-1])

  ## If data name existed, take out the x (and y) observations,
  ## group cluster, cluster cluster, stratum cluster, otherwise, no need to
  ## take values from a data frame.

  if(!is.null(pars$data)) {
    x <- data[, as.character(pars$x)]
    DNAME <- (pars$x)

    if(!is.null(pars$y)) {
      y <- data[, as.character(pars$y)]
      DNAME <- paste(DNAME, "and", pars$y)
    } else {
      y <- NULL
    }


    if(!is.null(pars$cluster)) {
      cluster <- data[, as.character(pars$cluster)]
      DNAME <- paste0(DNAME, ", cluster: ", pars$cluster)
    } else {
      cluster <- NULL
    }


    DNAME <- paste0(DNAME, " from ", pars$data)


  } else {

    DNAME <- deparse(substitute(x))

    if(!is.null(y)) {
      DNAME <- paste(DNAME, "and", deparse(substitute(y)))
    }


    if(!is.null(cluster)) {
      DNAME <- paste0(DNAME, ", cluster id: ", deparse(substitute(cluster)))
    }

  }


  ## Check and initialize cluster if not given,
  ## transform it to numeric if given as characters.

  l.x <- length(x)

  if( is.null(cluster)) {
    cluster <- c(1 : l.x)
  } else {
    if(!is.numeric(cluster)) {
      if(!is.character(cluster)) {
        stop("'cluster' has to be numeric or characters")
      }
      if(length(cluster) != l.x) {
        stop("'cluster' and 'x' must have the same lengths")
      }
      uniq.cluster <- unique(cluster)
      l.uniq.cluster <- length(uniq.cluster)
      cluster <- as.numeric(recoderFunc(cluster, uniq.cluster, c(1 : l.uniq.cluster)))
    }
  }



  ## Check x.
  if ( !is.numeric(x))
    stop("'x' must be numeric")


  ## Check data for paired test, paired test
  ## do not deal with stratified data.


  if( !is.null(y)) {
    if (!is.numeric(y))
      stop("'y' must be numeric")
    l.y <- length(y)
    if( l.y != l.x) {
      stop("'x' and 'y' must have the same
           lengths for signed rank test.")
    }

    OK <- complete.cases(x, y, cluster)
    x <- x[OK] - y[OK] - mu
    cluster <- cluster[OK]
    finite.x <- is.finite(x)
    x <- x[finite.x]
    cluster <- cluster[finite.x]

    if(length(x) < 1L) {
      stop("not enough (finite) 'x' observations")
    }

    } else {

      ## If only x is given, it is the difference score.

      OK <- complete.cases(x, cluster)
      x <- x[OK]
      cluster <- cluster[OK]
      finite.x <- is.finite(x)
      x <- x[finite.x] - mu
      cluster <- cluster[finite.x]
      if(length(x) < 1L) {
        stop("not enough (finite) 'x' observations")
      }
    }





  alternative <- match.arg(alternative)
  if(!missing(mu) && ((length(mu) > 1L) || !is.finite(mu))) {
    stop("'mu' must be a single number")
  }

  if(permutation == FALSE) {
    return(cluswilcox.test.signedrank(x, cluster, alternative, mu, DNAME, METHOD))
  } else {
    METHOD <- paste(METHOD, "using permutation")
    return(cluswilcox.test.signedrank.permutation(x, cluster, alternative, mu,
                                                  n.rep,
                                                  DNAME, METHOD))
  }
}
