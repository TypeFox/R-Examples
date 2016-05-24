# $Id: fast.prcomp.R 1025 2006-11-28 22:38:11Z warnes $

# The current implementation of the function svd() in S-Plus and R is
# much slower when operating on a matrix with a large number of
# columns than on the transpose of this matrix, which has a large
# number of rows. R's La.svd does not suffer from this problem.
#
# For R, the simple solution is to use La.svd instead of svd.  For
# S-Plus the solution is to check if the matrix is wider than tall,
# and to SVD on the transpose when this is the case.

if(exists("is.R") && is.R()==TRUE)
  {
    # This fast.prcomp() function is a slight modification of the
    # standard R prcomp function from package mva.  It uses La.svd instead
    # of the standard svd.  Consequently, it can be used with matrices
    # with many rows without a performance penalty.

    fast.prcomp <- function (x, retx = TRUE, center = TRUE, scale. = FALSE,
                        tol = NULL)
    {
        x <- as.matrix(x)
        x <- scale(x, center = center, scale = scale.)
        s <- La.svd(x, nu = 0)
        if (!is.null(tol)) {
            rank <- sum(s$d > (s$d[1] * tol))
            if (rank < ncol(x))
                s$vt <- s$vt[, 1:rank, drop = FALSE]
        }
        s$d <- s$d/sqrt(max(1, nrow(x) - 1))

        dimnames(s$vt) <- list(paste("PC", seq(len = nrow(s$vt)), sep = ""),
                               colnames(x) )
        r <- list(sdev = s$d, rotation = t(s$vt) )
        if (retx)
            r$x <- x %*% t(s$vt)
        class(r) <- "prcomp"
        r
    }

    fast.svd <- function( x, nu = min(n, p), nv = min(n, p), ...)
      {
        x <- as.matrix(x)
        dx <- dim(x)
        n <- dx[1]
        p <- dx[2]

        retval <- La.svd(x, nu=nu, nv=nv,  ... )
        retval$v <- t(retval$vt)
        retval$vt <- NULL
        retval
      }

  } else
  {
      # The fast.svd() function checks if the number of columns is
      # larger than the number of rows.  When this is the case, it
      # transposes the matrix, calles svd, and then flips the returned
      # u and v matrixes. Otherwise it just calls svd.
      #
      # This permits an SVD to be computed efficiently regardless of whether
      # n >>p or vice versa.

      # The fast.prcomp() function is simply a copy of the standard R
      # prcomp function from package mva which calls fast.svd instead of the
      # standard svd.  Consequently, it can be used with matrices with many
      # rows without a performance penalty.


    fast.prcomp <- function (x, retx = TRUE, center = TRUE, scale. = FALSE,
                        tol = NULL)
    {
        x <- as.matrix(x)
        x <- scale(x, center = center, scale = scale.)
        s <- fast.svd(x, nu = 0)
        if (!is.null(tol)) {
            rank <- sum(s$d > (s$d[1] * tol))
            if (rank < ncol(x))
                s$v <- s$v[, 1:rank, drop = FALSE]
        }
        s$d <- s$d/sqrt(max(1, nrow(x) - 1))
        dimnames(s$v) <- list(colnames(x), paste("PC", seq(len = ncol(s$v)),
            sep = ""))
        r <- list(sdev = s$d, rotation = s$v)
        if (retx)
            r$x <- x %*% s$v
        class(r) <- "prcomp"
        r
    }


    fast.svd <- function( x, nu = min(n, p), nv = min(n, p), ...)
      {
        x <- as.matrix(x)
        dx <- dim(x)
        n <- dx[1]
        p <- dx[2]

        if(  p <= n )
          return( svd( x, nu, nv, ... ) )
        else
          {
          s <- svd( t(x), nu=nv, nv=nu, ...)
          retval <- list()
          retval$d <- s$d
          retval$u <- s$v
          retval$v <- s$u
          return(retval)
        }
      }

    NULL

   }
