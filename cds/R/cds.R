#' Constrained Dual Scaling for Successive Categories with Groups
#' 
#' Uses an alternating nonnegative least squares algorithm combined with a
#' k-means-type algorithm to optimize the constrained group dual scaling
#' criterion outlined in the reference. Parallel computations for random starts
#' of the grouping matrix is supported via package \pkg{parallel}.
#' 
#' See the reference for more details.
#' 
#' @param x an object of class \code{"dsdata"} (see \code{cds.sim()}),
#' or a matrix (or object coercible to a matrix) containing the data for n
#' individuals on m objects. The data does not yet contain any additional
#' columns for the rating scale.
#' @param K The number of response style groups to look for. If a vector of 
#' length greater than one is given, the algorithm is run for each element
#' and a list of class \code{cdslist} is returned.
#' @param q The maximum rating (the scale is assumed to be \code{1:q}).
#' @param eps.ALS Numerical convergence criterion for the alternating least
#' squares part of the algorithm (updates for row and column scores).
#' @param eps.G Numerical convergence criterion for the k-means part of the
#' algorithm.
#' @param nr.starts.G Number of random starts for the grouping matrix.
#' @param nr.starts.a Number of random starts for the row scores.
#' @param maxit.ALS Maximum number of iterations for the ALS part of the
#' algorithm. A warning is given if this maximum is reached. Often it is not a
#' concern if this maximum is reached.
#' @param maxit Maximum number of iterations for the k-means part of the
#' algorithm.
#' @param Gstarts Facility to supply a list of explicit starting values for the
#' grouping matrix G. Each start consists of a two element list: \code{i} giving
#' and integer number the start, and \code{G} giving the starting configuration
#' as an indicator matrix.
#' @param astarts Supply explicit starts for the a vectors, as a list.
#' @param parallel logical. Should parallelization over starts for the grouping
#' matrix be used?
#' @param random.G logical. Should the k-means part consider the individuals in
#' a random order?
#' @param times.a.multistart The number of times that random starts for the row
#' scores are used. If == 1, then random starts are only used once for each
#' start of the grouping matrix.
#' @param info.level Verbosity of the output. Options are 1, 2, 3 and 4.
#' @param mc.preschedule Argument to mclapply under Unix.
#' @param seed Random seed for random number generators. Only partially
#' implemented.
#' @param LB logical. Load-balancing used in parallelization or not? Windows only.
#' @param reorder.grps logical. Use the Hungarian algorithm to reorder group
#' names so that the trace of the confusion matrix is maximized.
#' @param rescale.a logical. Rescale row score to length sqrt(2n) if TRUE
#' (after the algorithm has converged).
#' @param tol tolerance \code{tol} passed to \code{\link{lsei}} of the
#' \pkg{limSolve} package. Defaults to \code{sqrt(.Machine$double.eps)}
#' @param update.G Logical indicating whether or not to update the G matrix
#' from its starting configuration. Useful when clustering is known apriori or
#' not desired.
#' 
#' @return Object of class \code{ds} with elements: \item{G}{Grouping indicator
#' matrix.} \item{K}{Number of groups K.} \item{opt.crit}{Optimum value of the
#' criterion.} \item{a}{The 2n-vector of row scores.} \item{bstar}{The m-vector
#' of object scores.} \item{bkmat}{The matrix of group-specific boundary scores
#' for the ratings.} \item{alphamat}{The estimated spline coefficients for each
#' group.} \item{iter}{The number of iterations used for the optimal random
#' start wrt the grouping matrix.} \item{time.G.start}{The number of seconds it
#' took for the algorithm to converge for this optimal random start.}
#' \item{grp}{The grouping of the individuals as obtained by the algorithm.}
#' \item{kloss}{Loss value from G update (not equivalent to that of ALS
#' updates).} \item{hitrate, confusion}{Confusion and hitrates of original data
#' object contained a grouping vector.} \item{loss.G}{Optimality criterion
#' values for the random starts of G.} \item{q}{The number of ratings in the
#' Likert scale \code{1:q}} \item{time.total}{Total time taken for the
#' algorithm over all random starts} \item{call}{The function call.}
#' \item{data}{The input data object.}
#' 
#' @author Pieter C. Schoonees
#' @references Schoonees, P.C., Velden, M. van de & Groenen, P.J.F. (2013).
#' Constrained Dual Scaling for Detecting Response Styles in Categorical Data.
#' (EI report series EI 2013-10). Rotterdam: Econometric Institute.
#' @importFrom parallel makePSOCKcluster detectCores clusterSetRNGStream parLapplyLB parLapply 
#' clusterExport stopCluster clusterExport
#' @importFrom graphics abline barplot legend matplot par plot points rug symbols
#' @importFrom methods is
#' @importFrom stats cor pnorm qnorm quantile rnorm runif
#' @keywords multivariate
#' @examples
#' 
#' set.seed(1234)
#' dat <- cds.sim()
#' out <- cds(dat)
#' 
#' @export
cds <- function (x, K = 4, q = NULL, eps.ALS = 1e-3, eps.G = 1e-7, 
                 nr.starts.G = 20, nr.starts.a = 5, maxit.ALS = 20, 
                 maxit = 50, Gstarts = NULL, astarts = NULL, parallel = FALSE, 
                 random.G = FALSE, times.a.multistart = 1, info.level = 1, 
                 mc.preschedule = TRUE, seed = NULL, LB = FALSE, 
                 reorder.grps = TRUE, rescale.a = TRUE, 
                 tol = sqrt(.Machine$double.eps), update.G = TRUE)
{
  time1 <- proc.time()[3]
  cll <- match.call()
  
  ## Recurse if K is a vector
  if (length(K) > 1) {
    res <- lapply(K, cds, x = x, q = q, eps.ALS = eps.ALS, eps.G = eps.G, 
                  nr.starts.G = nr.starts.G, nr.starts.a = nr.starts.a, maxit.ALS = maxit.ALS, 
                  maxit = maxit, Gstarts = Gstarts, astarts = astarts, parallel = parallel, 
                  random.G = random.G, times.a.multistart = times.a.multistart, info.level = info.level, 
                  mc.preschedule = mc.preschedule, seed = seed, LB = LB, 
                  reorder.grps = reorder.grps, rescale.a = rescale.a, 
                  tol = tol, update.G = update.G)
    class(res) <- "cdslist"
    return(res)
  }
  
  ## Some sanity checks
  stopifnot(inherits(x, "cdsdata"))
  if (K == 1) {
    nr.starts.G <- 1
    message("Setting parallel = FALSE and nr.starts.G = 1 since there is only one group")
  }
  if (nr.starts.G == 1) parallel <- FALSE
  
  ## Some preliminaries
  n <- nrow(x$postrs)
  m <- x$m
  scales <- x$scales
  q <- max(scales)
  
  ## Calculate spline basis matrix
  x.bounds <- scales[-1] - 0.5
  tvec <- c(min(scales) + 0.5, mean(x.bounds), max(scales - 0.5))
  Mmat <- ispline(x.bounds, tvec = tvec)
  
  ## Extract Fr and calculate constant term in loss
  Fr.cent.rs <- x$Fr.cent.rs
  const <-  sum(Fr.cent.rs*Fr.cent.rs)
  
  ## Create G starts if not supplied
  if (!is.null(Gstarts))	{
  	  nr.starts.G <- length(Gstarts)
  	  Glst <- Gstarts
  	  ## Checks 
  	  if (any(sapply(Gstarts, function(x) ncol(x$G)) != K)) 
  	    stop("Argument 'K' does not match the number of classes implied by 'Gstarts'.")
  } else {
    ## Generate random starts
    Glst <- lapply(1:nr.starts.G, function(x) 
      list(i = x, G = diag(K)[sample(x = 1:K, size = n, replace = TRUE), , drop = FALSE]))
  }
  
  ## Evaluate starts in parallel if required  
  if (parallel) {
    RNGkind("L'Ecuyer-CMRG")
    if(!is.null(seed)) set.seed(seed)
    if(.Platform$OS.type != "unix") {
        cl <- parallel::makePSOCKcluster(parallel::detectCores())
        parallel::clusterExport(cl, c("Glst","group.ALS","Lfun", "Lfun.G.upd", "updateG", 
                            "G.start", "nr.starts.a", "maxit", "maxit.ALS", "n", 
                            "m", "q", "Fr.cent.rs", "Mmat", "info.level", "eps.ALS", 
                            "eps.G", "const", "times.a.multistart", "K", "random.G", "update.G"),
                            envir = environment())
#         clusterEvalQ(cl, {require(limSolve)
#                             require(clue)})
        if (!is.null(seed)) parallel::clusterSetRNGStream(cl, seed)
        if (LB) out.lst <- parallel::parLapplyLB(cl, Glst, G.start, nr.starts.a = nr.starts.a, 
                                      astarts = astarts, maxit = maxit, 
                                      maxit.ALS = maxit.ALS, n = n, m = m, q = q, 
                                      Fr.cent = Fr.cent.rs, Mmat = Mmat, 
                                      info.level = info.level, eps.ALS = eps.ALS, 
                                      eps.G = eps.G, const = const, 
                                      times.a.multistart = times.a.multistart, 
                                      K = K, random.G = random.G, tol = tol, update.G = update.G)
        else out.lst <- parallel::parLapply(cl, Glst, G.start, nr.starts.a = nr.starts.a, 
                                  astarts = astarts, maxit = maxit, 
                                  maxit.ALS = maxit.ALS, n = n, m = m, q = q, 
                                  Fr.cent = Fr.cent.rs, Mmat = Mmat, 
                                  info.level = info.level, eps.ALS = eps.ALS, 
                                  eps.G = eps.G, const = const, 
                                  times.a.multistart = times.a.multistart, K = K, 
                                  random.G = random.G, tol = tol, update.G = update.G)
        parallel::stopCluster(cl)
        } else {
      if (!is.null(seed)) parallel::mc.reset.stream()  
      out.lst <- parallel::mclapply(Glst, G.start, nr.starts.a = nr.starts.a, 
                          astarts = astarts, maxit = maxit, maxit.ALS = maxit.ALS,
                          n = n, m = m, q = q, Fr.cent = Fr.cent.rs, 
                          Mmat = Mmat, info.level = info.level, eps.ALS = eps.ALS, 
                          eps.G = eps.G, const = const, 
                          times.a.multistart = times.a.multistart, K = K, 
                          random.G = random.G, tol = tol, update.G = update.G, 
                          mc.cores = parallel::detectCores(), 
                          mc.preschedule = mc.preschedule)
      }
    } else {
    if (!is.null(seed)) set.seed(seed)
    out.lst <- lapply(Glst, G.start, nr.starts.a = nr.starts.a, astarts = astarts, maxit = maxit, 
                         maxit.ALS = maxit.ALS, n = n, m = m, q = q, Fr.cent = Fr.cent.rs, 
                         Mmat = Mmat, info.level = info.level, eps.ALS = eps.ALS, eps.G = eps.G,
                         const = const, times.a.multistart = times.a.multistart, K = K, 
                         random.G = random.G, tol = tol, update.G = update.G)
    }
  
  ## Determine best loss value and retain that one
  crit.G <- sapply(out.lst, '[[', "minloss")
  which.G <- which.min(crit.G)
  out <- out.lst[[which.G]]

  # Add info to object so that it can be returned
  K <- out$K

  # Determine response style clustering vector
  out$grp <- apply(out$G, 1, which.max)

  # Compare grp to true grouping if available
  if(!is.null(x$grp.rs)) {
    confmat <- table(x$grp.rs, out$grp)
    if(nrow(confmat) == K && reorder.grps && K > 1) {
        ord <- clue::solve_LSAP(confmat, maximum = TRUE)
        confmat <- confmat[,ord]
        out$G <- out$G[,ord]
        out$alphamat <- out$alphamat[ord,]
        out$bkmat <- out$bkmat[,ord]
        out$grp <- apply(out$G, 1, which.max)
        confmat <- table(x$grp.rs, out$grp)
        out$hitrate <- sum(diag(confmat))/n
        }
    out$confusion <- confmat
    }  
  out$loss.G <- crit.G
  
  # Rescale a to length 2*n and apply inverse scaling to b
  if(rescale.a){
    a.len <- as.numeric(sqrt(crossprod(out$a)/(2*n)))
    out$a <- out$a/a.len
    out$bkmat <- out$bkmat*a.len
    out$alphamat <- out$alphamat*a.len
    out$bstar <- out$bstar*a.len
  }
  
  # Set dimnames for alphamat
  rownames(out$alphamat) <- colnames(out$bkmat) <- paste0("g", 1:out$K)
  colnames(out$alphamat) <- c("mu", "a1","a2","a3")
  
  out$q <- q
  out$postrs <- x$postrs
  time2 <- proc.time()[3]
  out$time.total <- time2 - time1
  out$call <- cll
  class(out) <- c("cds","list")
  out
}
