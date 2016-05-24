bpca.prcomp <- function(x,
                        d=1:2, ...)
{
  stopifnot(class(x) == 'prcomp')

  if (!length(x$x))
    stop(gettextf("object '%s' has no objects coordinates!",
                  deparse(substitute(x))), domain=NA)

  if (is.complex(x$x))
    stop("bpca is not defined for complex PCA!")

  li <- d[1]
  le <- d[length(d)]
  n.lambda <- (le - li + 1)

  if(n.lambda < 2 || n.lambda > 3)
    stop('Please, check the parameter d:\n',
         'The (d[1] - d[length(d)] + 1) must equal to 2 (for bpca.2d) or 3 (for bpca.3d).\n\n')

  # Due to necessity of different type of factorization it will go back
  # from prcom, i.e, it will regenerate the data x already scaled!
  xreg <- x$x %*%
  (solve(t(x$rotation) %*%
         x$rotation) %*%
  t(x$rotation))
  # xreg <- x$x %*%
  # ginv(x$rotation) # another option (require MASS)
  bpca.default(xreg,
               d,
               center=ifelse(x$center[1] == FALSE, 0, 2),
               scale=ifelse(x$scale[1] == FALSE, 0, 1), ...)
}
