##  Multivariate normal density
#o  .dmnorm <- function (X,mu,sigma) ((2 * pi)^(-length(mu) / 2)) * (det(sigma)^(-1/ 2)) * exp (-0.5 * mahalanobis (X, mu, sigma))
.dmnorm <- function (X,mu,sigma)
{
  ((2 * pi)^(-length(mu) / 2)) *
  (det(sigma)^(-1/ 2)) *
  exp (-0.5 * .evmaha (X, mu, sigma))
}

.evmaha <- function (X, mu, sigma)    ##  calculate mahalanobis distances 
                                      ##  using the eigenvalues and eigenvectors. 
                                      ##  thus no singularity problems arise.
{                                     ##  Mahalanobis distance == Inf is possible.
    v <- eigen (sigma)
    Zt <- t (v$vectors) %*% (t (X) - mu)
    colSums ((Zt * v$values^(-.5))^2)
}

#.dmnorm <-
#function(X,mu,sigma)
#{
#  ((2*pi)^(-length(mu)/2))*(det(sigma)^(-1/2))*exp(-0.5*mahalanobis(X,mu,sigma))
#}

.doEllipses <-
function (eigval, eigvec, eigen, center, cov, n = 100, size = 1, ...)
{
  if (!missing (cov))
    eigen <- base::eigen (cov)
  if (!missing (eigen))
  {
    eigval <- eigen$values
    eigvec <- eigen$vectors
  }

  eigval[eigval < 0] <- 0

                      ##  check dimensionality of eigenvalues
  if (!is.numeric (eigval) || !length (eigval) == 2)
    stop ("argument eigval has to be a numeric vector of length 2.")

                      ##  check dimensionality of center
  if (!is.numeric (center) || !length (center) == 2)
    stop ("argument center has to be a numeric vector of length 2.")

                      ##  check dimensionality of eigenvectors
  if (!is.matrix (eigvec) || any (dim (eigvec) != 2))
    stop ("argument eigvec has to be a numeric mamtrix of dimension 2x2.")

  r <- seq (0, 2 * pi, length.out = n)    ##  create rad rep. of circle

  uc <- rbind (sin(r), cos (r))      ##  create a unit circle
                           ##  "stretch" circle corresponding to ev & sizefact
  uc = t(uc * sqrt(eigval) * size)
                           ##  rotate resulting ellipses from PC into XY coords
  XY = uc %*% t(eigvec)          

  XY = t( t(XY) + center)  ##  move ellipses to the specified center

  lines (XY[,1], XY[,2], ...) ##  draw ellipses
}

.getsubmatrix <-
function (x, idx)  matrix (x[drop = FALSE,,,idx], nrow = dim (x)[1])

.stretch <-
function (x, fact)
{
  range = diff (sort (x))
  c (x + range * fact * c(-1,1))
}

.vline <-
function (x, yfact = 2, col = 1, lty = 1, lwd = 1, ...)
{
   ylim = par ("usr")[3:4]
   ylim <- ylim + diff (ylim) * (1 - 1 / yfact) / 2 * c(1,-1)
   
   col = rep (col, length (x))
   lty = rep (lty, length (x))
   lwd = rep (lwd, length (x))
   
   for (i in 1:length (x))
     lines (rep (x[i], 2), ylim, col = col[i], lty = lty [i], lwd = lwd [i]
             , ...)
}
# 
# .setargs <- function (ARGS, FORMALS, DEFAULT = FALSE, ...)
# {
#   args <- list (...)
#   names.args <- names (args)
#   names.call <- names (ARGS)
#   
#   idx.call <- pmatch (names.call, FORMALS)
#   idx.arg <- pmatch (names.args, FORMALS)
# 
#   idx <- match (idx.arg, idx.call)
# 
#   idx.NA <- is.na (idx)
#   if (!DEFAULT)
#     ARGS[names.call [!idx.NA]] <- args[!idx.NA]
#     
#   ARGS[names.args[idx.NA]] <- args[idx.NA]
#   ARGS
# }
# 
# .getarg <- function (ARGS, FORMALS, arg)
# {
#   idx.arg <- which (FORMALS == arg)
#   stopifnot (length (idx.arg) == 1)
# 
#   idx.call <- which (pmatch (names (ARGS), FORMALS) == idx.arg)
#   stopifnot (length (idx.call) == 1)
#   
#   ARGS[[idx.call]]
# }
# 
# .formalsnames <- function (x) names (formals (x))
# 
# .getformals <- function (f)
# {
#   args <- unlist (sapply (f, .formalsnames, simplify = TRUE))
#   args [args != "..."]
# }



.Conv2Matrix <- function (x, sx = substitute (x))
{
    if(is.matrix(x))
        return (x)
    if(is.data.frame(x))
        return (data.matrix(x))
    return (matrix(x, nrow = length(x), ncol = 1, dimnames = list(names(x), deparse(sx))))
}
