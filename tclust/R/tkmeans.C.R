
tkmeans <-
function (x, k = 3, alpha = 0.05, nstart = 50, iter.max = 20,
          equal.weights = FALSE, center = 0, scale = 1, store.x = TRUE,
		  drop.empty.clust = TRUE, trace = 0, warnings = 2, zero.tol = 1e-16)
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
               iter.max = iter.max, equal.weights = equal.weights,
			   center = center, scale = scale, store.x = store.x,
			   drop.empty.clust = drop.empty.clust, trace = trace,
			   warnings = warnings, zero.tol = zero.tol)

  par <- .tkmeans.preproc (par)

  ret.C <- .C ( "tkmeans", PACKAGE = "tclust", DUP = TRUE
        , as.integer (c (dim (par$x), par$k, par$equal.weights, par$usetrace, par$nstart, par$iter.max))
        , parN = integer (5)
        , as.double (c (par$alpha, par$zero.tol))
        , parD = double (2)
        , as.double (par$x)
        , center = double (par$p * par$k)
        , cluster = integer (nrow (par$x))
        , size = double (par$k)
        , weights = double (par$k)
        , er.obj = double (par$nstart)
        , er.conv = integer (par$nstart)
  )

  par <- .tkmeans.postproc (par, ret.C)

  return (par$ret)
}
