
tclust <-
function (x, k = 3, alpha = 0.05, nstart = 50, iter.max = 20,
          restr = c ("eigen", "deter", "sigma"), restr.fact = 12,
          equal.weights = FALSE, center = 0, scale = 1, store.x = TRUE,
          drop.empty.clust = TRUE, trace = 0, warnings = 3, zero.tol = 1e-16
         )
{
#  Disabled arguments:
#  restr = c ("eigen", "deter", "sigma", "dir.eigen", "dir.deter", "prop"),
#  iter.tune,
  ovv <- 0      ##  optimization parameter for "optVectors"
  fuzzy <- FALSE
  m <- 2
  iter.tune <- 10
  x.s <- substitute(x)

  par <- list (x = x, x.s = x.s, k = k, alpha = alpha, nstart = nstart,
               iter.max = iter.max, restr = restr[1], restr.fact = restr.fact,
               equal.weights = equal.weights, center = center, scale = scale,
			   store.x = store.x, drop.empty.clust = drop.empty.clust,
			   trace = trace, warnings = warnings, zero.tol = zero.tol,
			   ovv = ovv, fuzzy = fuzzy, m = m, iter.tune = iter.tune)

  par <- .tclust.preproc (par)

  ret.C <- .C ( "tclust", PACKAGE = "tclust", DUP = TRUE
        , as.integer (c (dim (par$x), par$k, par$fuzzy, par$nstart, par$iter.max,
		            par$equal.weights, par$restr.C, par$deter.C, par$usetrace,
					par$iter.tune, par$ovv))
        , parN = integer (5)
        , as.double (c (par$alpha, par$restr.fact, par$m, par$zero.tol))
        , parD = double (2)
        , as.double (par$x)
        , center = double (par$p * par$k)
        , cov = double (par$p * par$p * par$k)
        , cluster = integer (nrow (par$x))
        , size = double (par$k)
        , weights = double (par$k)
        , z = double (par$z.size)
        , er.obj = double (par$nstart)
        , er.conv = integer (par$nstart)
  )

  par <- .tclust.postproc (par, ret.C)

  return (par$ret)
}

.tclust.int <-
function (x, k = 3, alpha = 0.05, nstart = 50, iter.max = 20,
          restr = c ("eigen", "deter", "sigma"), restr.fact = 12,
          equal.weights = FALSE, center = 0, scale = 1, store.x = TRUE,
          drop.empty.clust = TRUE, trace = 0, warnings = 3, zero.tol = 1e-16,
		  fuzzy = FALSE, m = 2, iter.tune = 10, ovv = 0
         )
{
#  Disabled arguments:
#  restr = c ("eigen", "deter", "sigma", "dir.eigen", "dir.deter", "prop"),
#  iter.tune,
#  ovv <- 0      ##  optimization parameter for "optVectors"
#  fuzzy <- FALSE
#  m <- 2
#  iter.tune <- 10
  x.s <- substitute(x)

  par <- list (x = x, x.s = x.s, k = k, alpha = alpha, nstart = nstart,
               iter.max = iter.max, restr = restr[1], restr.fact = restr.fact,
               equal.weights = equal.weights, center = center, scale = scale,
			   store.x = store.x, drop.empty.clust = drop.empty.clust,
			   trace = trace, warnings = warnings, zero.tol = zero.tol,
			   ovv = ovv, fuzzy = fuzzy, m = m, iter.tune = iter.tune)

  par <- .tclust.preproc (par)

  ret.C <- .C ( "tclust", PACKAGE = "tclust", DUP = TRUE
        , as.integer (c (dim (par$x), par$k, par$fuzzy, par$nstart, par$iter.max,
		            par$equal.weights, par$restr.C, par$deter.C, par$usetrace,
					par$iter.tune, par$ovv))
        , parN = integer (5)
        , as.double (c (par$alpha, par$restr.fact, par$m, par$zero.tol))
        , parD = double (2)
        , as.double (par$x)
        , center = double (par$p * par$k)
        , cov = double (par$p * par$p * par$k)
        , cluster = integer (nrow (par$x))
        , size = double (par$k)
        , weights = double (par$k)
        , z = double (par$z.size)
        , er.obj = double (par$nstart)
        , er.conv = integer (par$nstart)
  )

  par <- .tclust.postproc (par, ret.C)

  return (par$ret)
}



.tkmeans.preproc <- function (O)
{

  O$x <- .Conv2Matrix (O$x, O$x.s)
  if (!is.numeric (O$x))
    stop ("parameter x: numeric matrix/vector expected")

  O$usetrace <- .tr (O$trace, 0)

  O$p <- ncol (O$x)
  O$n <- nrow (O$x)

  scaled <- .ScaleAdv (O$x, O$center, O$scale)
  O$x <- scaled$x
  O$scale <- scaled$scale
  O$center <- scaled$center
  O
}

.tclust.preproc <- function (O)
{
  O <- .tkmeans.preproc (O)

  rd <- .match.restr (O$restr)
  O$restr.C <- rd$restr
  O$deter.C <- rd$deter

  O$iter.tune <- rep (O$iter.tune, len = 3)

  O$z.size <- ifelse (O$fuzzy, O$n * O$k, 0)
  O
}

.tkmeans.postproc <- function (O, ret.C)
{

  if (ret.C$parN[4])	##	error excecution
    stop ()

  parlist <- list (k = O$k, alpha = O$alpha, nstart = O$nstart,
                   iter.max = O$iter.max, equal.weights = O$equal.weights)

  if (O$store.x)
    parlist$x <- O$x

  if (O$drop.empty.clust)
    idxuse <- which (ret.C$size > 0)
  else
    idxuse <- 1:O$k

  idxuse <- idxuse [order (ret.C$size[idxuse], decreasing = TRUE)]

  ClusterIDs <- rep (0, O$k + 1)
  ClusterIDs[idxuse + 1] <- 1:length (idxuse)

  k.real <- length (idxuse)

  int <- list (
      iter.converged = ret.C$parN[1]
    , iter.successful = ret.C$parN[2]
    , dim = dim (O$x)
    , code = ret.C$parN[3]
    , er.obj = ret.C$er.obj
    , er.conv = ret.C$er.conv
  )

  ret <- list (
    centers = array (ret.C$center, c(O$p, O$k)) [,idxuse, drop = FALSE]
    , cluster = ClusterIDs [ret.C$cluster + 2]
    , par = parlist
    , k = sum (ret.C$size > 0)
    , obj = ret.C$parD[1]
    , size = ret.C$size[idxuse]
    , weights = ret.C$weights[idxuse]
    , int = int
  )

  dn.x <- dimnames (O$x)

  rownames (ret$centers) <-
    if (is.null (colnames (O$x)))
      paste ("X", 1:ncol (O$x))
    else
      colnames (O$x)
#    if (is.null (dn.x[[2]])) paste ("X", 1:ncol (O$x)) else dn.x[[2]]
  colnames (ret$centers) <- paste ("C", 1:k.real)

  ret <- .tkmeans.warnings (O, ret)
  .tkmeans.warn (O, ret)

  for (i in 1:ret$k)
    ret$centers[,i] <- ret$centers[,i] * O$scale + O$center

  class (ret) <- "tkmeans"

  O$ret <- ret
  O$idxuse <- idxuse
  O$ClusterIDs <- ClusterIDs

  return (O)
}

.get.Mm.eigen <- function (x, O)
{
	idx <- which (x == O$ret$cluster)
	if (length (idx) <= 1)
		return (rep (NA, ncol (O$x)))
	eigen  (cov (O$x[idx, , drop = FALSE]))$value
}

.get.Mm.det <- function (x, O)
{
	idx <- which (x == O$ret$cluster)
	if (length (idx) <= 1)
		return (NA)
	det (cov (O$x[idx, , drop = FALSE]))
}

.tclust.postproc <- function (O, ret.C)
{
  O <- .tkmeans.postproc (O, ret.C)

  O$ret$par$restr.fact <- O$restr.fact
  O$ret$par$restr      <- O$restr
  O$ret$par$restr.C    <- O$restr.C
  O$ret$par$deter.C    <- O$deter.C

  O$ret$cov <- array (ret.C$cov, c (O$p, O$p, O$k))[, , O$idxuse, drop = FALSE]

  dimnames (O$ret$cov) <- dimnames (O$ret$centers)[c (1, 1, 2)]

	## calculate the "unrestr.fact"
#  if (O$deter.C)
#    get.Mm <- function (x, O) det (cov (O$x[x == O$ret$cluster, ]))
#  else
#    get.Mm <- function (x, O) eigen (cov (O$x[x == O$ret$cluster, ]))$values

  get.Mm <- if (O$deter.C) .get.Mm.det else .get.Mm.eigen

  EV <- sapply (1:O$ret$k, get.Mm, O = O)

  if (all (is.na (EV)))
	O$ret$unrestr.fact <- 1
  else
	  O$ret$unrestr.fact <- ceiling (max (EV, na.rm = TRUE) / min (EV, na.rm = TRUE))

  O$ret <- .tclust.warnings (O, O$ret)
  .tclust.warn (O, O$ret)



  if (O$fuzzy)
    O$ret$z <- matrix (ret.C$z, O$n, O$k)[,O$idxuse]

  cmm <- O$scale %*% t (O$scale)
       #  matrix (O$scale, nrow = O$p, ncol = 1) %*%
       #  matrix (O$scale, nrow = 1, ncol = O$p)

  O$ret$mah <- array (dim = O$n)
  O$ret$mah[!O$ret$cluster] <- NA

  for (i in 1:O$ret$k)
  {
    idx <- O$ret$cluster == i
    O$ret$mah[idx] <- mahalanobis (O$x[idx, , drop = FALSE], center = O$ret$centers[, i], cov = O$ret$cov[,, i])
    O$ret$cov[,,i] <- O$ret$cov[,,i] * cmm
  }

  class (O$ret) <- c (class (O$ret), "tclust")

  return (O)
}

#.tkmeans.err.exec <- function (ret.C) ret.C$parN[4]

.tkmeans.warnings <- function (O, ret)
{
  ret$warnings <- list (
                  singular = ret$int$code == 2
                , iter = ret$int$iter.successful && ret$int$iter.converged /
                         ret$int$iter.successful < 0.5
                , drop = O$k > ret$k
                , size = any (ret$size < O$n / 50)
                , sizep = min (ret$size) <= O$p
                , smallobj = ret$obj < (-1e+20)
                )
  return (ret)
}

.tclust.warnings <- function (O, ret)
{

  ret$warnings$restr.lo = ret$unrestr.fact > O$restr.fact
  ret$warnings$restr.hi = ret$unrestr.fact * 2 < O$restr.fact

  if (ret$warnings$sizep || ret$warnings$size)
    ret$warnings$restr.lo <- FALSE

  return (ret)
}


.tkmeans.warn <- function (O, ret)
{
  if (O$warnings >= 1)
  {
    if (ret$warnings$iter)
      warning (paste ("Less than 50% of the iterations (",
        round (ret$int$iter.converged / ret$int$iter.successful * 100, 1),
               "%) converged - please increase iter.max.", sep = ""))
    if (ret$warnings$smallobj)
      warning ("Due to a very small objective function's value, the final solution does not seem to be reliable.\n  More iterations are probably needed (-> increase \"nstart\").")
  }

  if (O$warnings >= 2)
  {
    if (ret$warnings$singular)    ## not all iterations could be executed
      if (ret$par$deter.C)
        warning ("All observations are concentrated in k subspaces after trimming.")
      else
        warning ("All observations are concentrated in k points after trimming.")
#        warning ("Points in the data set are concentrated in k subspaces after trimming.")
#      else
#        warning ("points in the data set are concentrated in k points after trimming")

    if (ret$warnings$drop)
    {
      n.drop <- ret$par$k - ret$k
      if (n.drop == 1)
          warning (paste (n.drop, "empty cluster has been detected - try reducing k."))
      else
          warning (paste (n.drop, "empty clusters have been detected - try reducing k."))
    }
    else if (ret$warnings$size)
      warning ("Clusters with size < n * 0.02 found - try reducing k.")
    else if (ret$warnings$size)
     warning ("Clusters with size <= p found - try reducing k.")

  }

}

.tclust.warn <- function (O, ret)
{
  if (O$warnings >= 3)
  {
    if (ret$par$restr != "sigma")
    {
      if (ret$warnings$restr.lo)
         warning (paste ("The result is artificially constrained due to ",
         "restr.fact = ",ret$par$restr.fact,".", sep = ""))
#       warning (paste ("The chosen restriction factor (", ret$par$restr.fact,
#       ") artificially restricts the solution.\n",
#       "  This solution implies a restriction factor of ",
#       ceiling (ret$unrestr.fact), ".", sep = ""))

#      if (ret$warnings$restr.hi)    ##  warning currently disabled..
#        warning (paste ("The restriction factor (", ret$par$restr.fact,
#        ") has been chosen too large for this solution.\n",
#        "  This solution implies a restriction factor of ",
#        ceiling (ret$unrestr.fact), ".", sep = ""))
    }
  }
}

.match.restr <- function (restr)
{
  stopifnot (is.character (restr))

  restr <- restr[1]
  restr <- match.arg (restr, c ("eigen", "deter", "sigma", "dir.eigen",
                               "dir.deter", "prop", "none"))
  deter <- (restr == "deter" || restr == "dir.deter")

  if (restr == "eigen" || restr == "deter")
    restr <- 0
  else if (restr == "dir.eigen" || restr == "dir.deter")
    restr <- 1
  else if (restr == "sigma")
    restr <- 2
  else if (restr == "prop")
    restr <- 3
  else if (restr == "none")
    restr <- 4
  else
    stop ("unknown value for parameter restr")

  list (restr = restr, deter = deter)
}

#convplot <- function (x)
#{  ##  function for evaluating the convergence of a run of tclust
#   ord <- order (x$int$er.obj)
#   plot (x$int$er.obj[ord], col = 2-x$er.conv[ord],
#   xlab = "Iteration (ordered)",
#   ylab = "Objective Function", main = "Convergence Plot")
#}
