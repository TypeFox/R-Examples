#' Calibrate sample weights
#'
#' Calibrate sample weights according to known marginal population totals.
#' Based on initial sample weights, the so-called \emph{g}-weights are computed
#' by generalized raking procedures.
#'
#' The methods return a list containing both the \emph{g}-weights (slot
#' \code{g_weights}) as well as the final weights (slot \code{final_weights})
#' (initial sampling weights adjusted by the \emph{g}-weights.
#'
#' @name calibSample
#' @aliases calibSample,df_or_dataObj_or_simPopObj,dataFrame_or_Table-method
#' @docType methods
#' @note This is a faster implementation of parts of
#' \code{\link[sampling]{calib}} from package \code{sampling}. Note that the
#' default calibration method is raking and that the truncated linear method is
#' not yet implemented.
#' @section Methods: \describe{The function provides methods with the following
#' signatures.  \item{list("signature(inp=\"df_or_dataObj_or_simPopObj\",
#' totals=\"dataFrame_or_Table\",...)")}{ Argument 'inp' must be an object of
#' class \code{data.frame}, \code{\linkS4class{dataObj}} or
#' \code{\linkS4class{simPopObj}} and the totals must be specified in either
#' objects of class \code{table} or \code{data.frame}. If argument 'totals' is
#' a data.frame it must be provided in a way that in the first columns
#' n-columns the combinations of variables are listed. In the last column, the
#' frequency counts must be specified. Furthermore, variable names of all but
#' the last column must be available also from the sample data specified in
#' argument 'inp'. If argument 'total' is a table (e.g. created with function
#' \code{\link{tableWt}}, it must be made sure that the dimnames match the
#' variable names (and levels) of the specified input data set.  } }
#' @author Andreas Alfons and Bernhard Meindl
#' @references Deville, J.-C. and Saerndal, C.-E. (1992)
#' Calibration estimators in survey sampling. \emph{Journal of the American
#' Statistical Association}, \bold{87}(418), 376--382.  Deville, J.-C.,
#' Saerndal, C.-E. and Sautory, O. (1993) Generalized raking
#' procedures in survey sampling. \emph{Journal of the American Statistical
#' Association}, \bold{88}(423), 1013--1020.
#' @keywords survey methods
#' @export calibSample
#' @examples
#' data(eusilcS)
#' eusilcS$agecut <- cut(eusilcS$age, 7)
#' inp <- specifyInput(data=eusilcS, hhid="db030", hhsize="hsize", strata="db040", weight="db090")
#'
#' ## for simplicity, we are using population data directly from the sample, but you get the idea
#' totals1 <- tableWt(eusilcS[, c("agecut","rb090")], weights=eusilcS$rb050)
#' totals2 <- tableWt(eusilcS[, c("rb090","agecut")], weights=eusilcS$rb050)
#' totals3 <- tableWt(eusilcS[, c("rb090","agecut","db040")], weights=eusilcS$rb050)
#' totals4 <- tableWt(eusilcS[, c("agecut","db040","rb090")], weights=eusilcS$rb050)
#'
#' weights1 <- calibSample(inp, totals1)
#' totals1.df <- as.data.frame(totals1)
#' weights1.df <- calibSample(inp, totals1.df)
#' identical(weights1, weights1.df)
#'
#' # we can also use a data.frame and an optional weight vector as input
#' df <- as.data.frame(inp@@data)
#' w <- inp@@data[[inp@@weight]]
#' weights1.x <- calibSample(df, totals1.df, w=inp@@data[[inp@@weight]])
#' identical(weights1, weights1.x)
#'
#' weights2 <- calibSample(inp, totals2)
#' totals2.df <- as.data.frame(totals2)
#' weights2.df <- calibSample(inp, totals2.df)
#' identical(weights2, weights2.df)
#'
#' \dontrun{
#' ## approx 10 seconds computation time ...
#' weights3 <- calibSample(inp, totals3)
#' totals3.df <- as.data.frame(totals3)
#' weights3.df <- calibSample(inp, totals3.df)
#' identical(weights3, weights3.df)
#'
#' ## approx 10 seconds computation time ...
#' weights4 <- calibSample(inp, totals4)
#' totals4.df <- as.data.frame(totals4)
#' weights4.df <- calibSample(inp, totals4.df)
#' identical(weights4, weights4.df)
#' }
NULL

setClassUnion("df_or_dataObj_or_simPopObj", c("data.frame", "dataObj", "simPopObj"))
setClassUnion("dataFrame_or_Table", c("data.frame", "table"))
setGeneric("calibSample", function(inp, totals, ...) {
  standardGeneric("calibSample")
})

#' @export
setMethod("calibSample", c(inp="df_or_dataObj_or_simPopObj", totals="dataFrame_or_Table"), function(inp, totals, ...) {
  if ( class(inp) == "data.frame" ) {
    samp <- data.table(inp)
  }
  if ( class(inp) == "dataObj" ) {
    samp <- inp@data
  }
  if ( class(inp) == "simPopObj" ) {
    samp <- inp@sample@data
  }

  if ( class(totals) == "table" ) {
    totals <- as.data.frame(totals)
  }

  if ( ncol(totals) < 2 ) {
    stop("we need at least one dimension variable and one column for frequencies!\n Check your input!\n")
  }

  freqs <- totals[,ncol(totals)]
  totals <- totals[,-ncol(totals), drop=F]
  vnames <- colnames(totals)

  #  check if vars exist
  if ( !all(vnames %in% colnames(samp)) ) {
    print(vnames)
    stop("not all slotnames of argument 'totals' exist in the sample slot of argument 'inp'!\n")
  }

  # handle additional arguments
  args <- list(...)
  if ( is.null(args$method) ) {
    method <- "raking"
  }
  if ( is.null(args$bounds) ) {
    bounds <- c(0,10)
  }
  if ( is.null(args$maxit) ) {
    maxit <- 500
  }
  if ( is.null(args$tol) ) {
    tol <- 1e-06
  }
  if ( is.null(args$q) ) {
    q <- NULL
  }
  if ( is.null(args$eps) ) {
    eps <- .Machine$double.eps
  }

  # everything ok so far
  tmp <- samp[,vnames,with=FALSE]
  fac <- apply(tmp, 1, function(x) { paste(x, collapse="-") })

  # check if all combinations are available in the dataset
  grid <- expand.grid(lapply(tmp, function(x) {
    unique(x)
  }))

  # check consistency
  if ( nrow(grid) > nrow(totals) ) {
    stop("some combinations of characteristics are missing from input argument 'totals'!\n")
  }
  if ( nrow(grid) < nrow(totals) ) {
    stop("in input argument 'totals' some combinations are listed that are not available from the sample!\n")
  }

  grid_fac <- apply(grid, 1, function(x) {
    paste(x, collapse="-")
  })
  totals_fac <- apply(totals, 1, function(x) {
    paste(x, collapse="-")
  })
  ii <- match(grid_fac, totals_fac)
  if ( any(is.na(ii)) ) {
    stop("some characteristings in argument 'totals' differ from those in the actual data!")
  }

  # create binary factors
  X <- calibVars(fac)

  # order totals
  ii <- match(colnames(X), totals_fac)

  grid <- grid[ii,]
  freqs <- freqs[ii]

  # initial sample weights
  if ( class(inp) == "dataObj" ) {
    w <- samp[[inp@weight]]
  }
  if ( class(inp) == "simPopObj" ) {
    w <- samp[[inp@sample@weight]]
  }
  if ( class(inp) == "data.frame" ) {
    if ( !is.null(args$w) ) {
      w <- args$w
      if ( length(w) != nrow(samp) ) {
        stop("if argument 'w' was provided, then its dimension must match the number of rows from argument 'inp'!\n")
      }
    } else {
      w <- rep(1, nrow(samp))
    }
  }

  totals <- freqs

  # g-weights
  g_weights <- calibSample_work(X, d=w, totals=totals, q=q, method=method, bounds=bounds, maxit=maxit, tol=tol, eps=eps)

  # final-weights
  final_weights <- g_weights*w

  invisible(list(g_weights=g_weights, final_weights=final_weights))
})

calibSample_work <- function(X, d, totals, q=NULL,
  method=c("raking", "linear", "logit"),
  bounds=c(0, 10), maxit=500, tol=1e-06,
  eps=.Machine$double.eps) {

  ## initializations and error handling
  X <- as.matrix(X)
  d <- as.numeric(d)
  totals <- as.numeric(totals)
  haveNA <- c(any(is.na(X)), any(is.na(d)),
    any(is.na(totals)), !is.null(q) && any(is.na(q)))
  if ( any(haveNA) ) {
    argsNA <- c("'X'", "'d'", "'totals'", "'q'")[haveNA]
    stop("missing values in the following arguments", paste(argsNA, collapse=", "),"\n")
  }
  n <- nrow(X)  # number of rows
  if ( length(d) != n ) {
    stop("length of 'd' not equal to number of rows in 'X'!\n")
  }
  p <- ncol(X)  # number of columns
  if ( length(totals) != p ) {
    stop("length of 'totals' not equal to number of columns in 'X'!\n")
  }
  if ( is.null(q) ) {
    q <- rep.int(1, n)
  } else {
    q <- as.numeric(q)
    if ( length(q) != n ) {
      stop("length of 'q' not equal to number of rows in 'X'!\n")
    }
    if ( any(is.infinite(q)) ) {
      stop("infinite values in 'q'")
    }
  }
  method <- match.arg(method)

  ## computation of g-weights
  if ( method == "linear" ) {
    ## linear method (no iteration!)
    lambda <- ginv(t(X * d * q) %*% X, tol=eps) %*% (totals - as.vector(t(d) %*% X))
    g <- 1 + q * as.vector(X %*% lambda)  # g-weights
  } else {
    ## multiplicative method (raking) or logit method
    lambda <- matrix(0, nrow=p)  # initial values
    # function to determine whether teh desired accuracy has
    # not yet been reached (to be used in the 'while' loop)
    tolNotReached <- function(X, w, totals, tol) {
      max(abs(crossprod(X, w) - totals)/totals) >= tol
    }
    if ( method == "raking" ) {
      ## multiplicative method (raking)
      # some initial values
      g <- rep.int(1, n)  # g-weights
      w <- d  # sample weights
      ## iterations
      i <- 1
      while ( !any(is.na(g)) && tolNotReached(X, w, totals, tol) && i <= maxit ) {
        # here 'phi' describes more than the phi function in Deville,
        # Saerndal and Sautory (1993); it is the whole last term of
        # equation (11.1)
        phi <- t(X) %*% w - totals
        T <- t(X * w)
        dphi <- T %*% X  # derivative of phi function (to be inverted)
        lambda <- lambda - ginv(dphi, tol=eps) %*% phi  # update 'lambda'
        g <- exp(as.vector(X %*% lambda) * q)  # update g-weights
        w <- g * d  # update sample weights
        i <- i + 1  # increase iterator
      }
      ## check wether procedure converged
      if( any(is.na(g)) || i > maxit ) {
        warning("no convergence!\n")
        g <- NULL
      }
    } else {
      ## logit (L, U) method
      ## error handling for bounds
      if ( length(bounds) < 2 ) {
        stop("'bounds' must be a vector of length 2")
      } else {
        bounds <- bounds[1:2]
      }
      if ( bounds[1] >= 1 ) {
        stop("the lower bound must be smaller than 1")
      }
      if ( bounds[2] <= 1 ) {
        stop("the lower bound must be larger than 1")
      }
      ## some preparations
      A <- diff(bounds)/((1 - bounds[1]) * (bounds[2] - 1))
      # function to bound g-weights
      getG <- function(u, bounds) {
        (bounds[1] * (bounds[2]-1) + bounds[2] * (1-bounds[1]) * u) /
        (bounds[2]-1 + (1-bounds[1]) * u)
      }
      ## some initial values
      g <- getG(rep.int(1, n), bounds)  # g-weights
      # in the procedure, g-weights outside the bounds are moved to the
      # bounds and only the g-weights within the bounds are adjusted.
      # these duplicates are needed since in general they are changed in
      # each iteration while the original values are also needed
      X1 <- X
      d1 <- d
      totals1 <- totals
      q1 <- q
      g1 <- g
      indices <- 1:n
      # function to determine which g-weights are outside the bounds
      anyOutOfBounds <- function(g, bounds) {
        any(g < bounds[1]) || any(g > bounds[2])
      }
      ## iterations
      i <- 1
      while ( !any(is.na(g)) && (tolNotReached(X, g*d, totals, tol ) || anyOutOfBounds(g, bounds)) && i <= maxit) {
        # if some of the g-weights are outside the bounds, these values
        # are moved to the bounds and only the g-weights within the
        # bounds are adjusted
        if ( anyOutOfBounds(g, bounds) ) {
          g[g < bounds[1]] <- bounds[1]
          g[g > bounds[2]] <- bounds[2]
          # values within the bounds
          tmp <- which(g > bounds[1] & g < bounds[2])
          if ( length(tmp) > 0 ) {
            indices <- tmp
            X1 <- X[indices,]
            d1 <- d[indices]
            if ( length(indices) < n ) {
              totals1 <- totals - as.vector(t(g[-indices] * d[-indices]) %*% X[-indices, , drop=FALSE])
            }
            q1 <- q[indices]
            g1 <- g[indices]
          }
        }
        w1 <- g1 * d1  # current sample weights
        # here 'phi' describes more than the phi function in Deville,
        # Saerndal and Sautory (1993); it is the whole last term of
        # equation (11.1)
        phi <- t(X1) %*% w1 - totals1
        T <- t(X1 * w1)
        dphi <- T %*% X1  # derivative of phi function (to be inverted)
        lambda <- lambda - ginv(dphi, tol=eps) %*% phi  # update 'lambda'
        # update g-weights
        u <- exp(A * as.vector(X1 %*% lambda) * q1)
        g1 <- getG(u, bounds)
        g[indices] <- g1
        i <- i+1  # increase iterator
      }
      ## check wether procedure converged
      if ( any(is.na(g)) || i > maxit ) {
        warning("no convergence!\n")
        g <- NULL
      }
    }
  }
  ## return g-weights
  invisible(g)
}
