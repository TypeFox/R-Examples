#' Constrained Dual Scaling for a Single Random G Start
#' 
#' Run algorithm for a single G matrix.
#' 
#' @param X List of two elements, namely \code{i} giving the number of the start
#' and \code{G} given the starting configuration
#' @param nr.starts.a The number or random starts for \code{a} to use in the ALS.
#' @param astarts Explicit starts for a, if applicable.
#' @param maxit The maximum number of iterations with respect to G. 
#' @param n The number of respondents.
#' @param m The number of items.
#' @param q The maximum rating category such that the rating scale is \code{1:q}.
#' @param Fr.cent The centred Fr matrix.
#' @param maxit.ALS The maximum number of ALS iterations.
#' @param Mmat The basis matrix for the quadratic monotone splines.
#' @param eps.G The absolute error tolerance for the G updates.
#' @param info.level Integer controlling the amount of information printed.
#' @param times.a.multistart The number of times random starts for \code{a} is used.
#' @param eps.ALS The absolute error tolerance for the ALS.
#' @param const The constant part of the loss function.
#' @param K The number of groups.
#' @param random.G The \code{random} argument passed to \code{\link{updateG}}.
#' @param tol tolerance \code{tol} passed to \code{\link{lsei}} of the
#' \pkg{limSolve} package)
#' @param update.G Logical indicating whether or not to update the starting configuration
#' \code{G} in \code{X}
#' @keywords multivariate
#' @export G.start
G.start <- function(X, nr.starts.a, astarts, maxit, n, m, q, Fr.cent, maxit.ALS, Mmat, eps.G,
                      info.level, times.a.multistart, eps.ALS, const, K, random.G, tol, update.G) {
  time.start <- proc.time()[3]
  G <- X$G
  i <- X$i
  Gmoves <- -1 ## For convergence check when update.G is FALSE
  iter <- 1
  overall.crit <- rep(NA, maxit)
  
  # Allow for start(s) for a to be supplied explicitly
  if(!is.null(astarts)){
    a.lst <- as.list(astarts)
    nr.starts.a <- length(a.lst)
    times.a.multistart <- 1L
    cat(nr.starts.a, "explicit start(s) for a supplied\n")
  }
  
  ### Iterate between ALS and k-means 
  
  while(iter <= maxit)
  {
    ### Update a and X given G 
    
    # Do random starts for a, only times.a.multistart times
    # Keep the best start from the previous iteration for all iter > 1
    
    if(iter == 1 & is.null(astarts)) a.lst <- lapply(rep(2*n, nr.starts.a), rnorm)
    if(iter > 1 && iter <= times.a.multistart) {
      a.lst <- lapply(rep(2*n, nr.starts.a), rnorm)
      a.lst[[1]] <- a.cur
    }
    if(iter > times.a.multistart) a.lst <- list(a.cur)
    a.out.lst <- lapply(a.lst, group.ALS, m = m, q = q, G = G, Fr.cent = Fr.cent, 
                        maxit = maxit.ALS, eps = eps.ALS, Mmat = Mmat, 
                        info.level = info.level, const= const, K = K, n = n, tol = tol)
    
    crit.a <- sapply(a.out.lst, '[[', "crit.opt")
    which.a <- which.min(crit.a)
    out <- a.out.lst[[which.a]]
    a.cur <- out$a
    alphamat <- out$alpha
    bkmat <- out$bk
    awts2 <- out$awts2
    bwts2 <- out$bwts2
    Fr.bk <- out$Fr.bk
    
    ### Update G given a and X only if required
    if (update.G){
      out.G <- updateG(G = G, a = a.cur, bwts2 = bwts2, Fr.bk = Fr.bk, 
                       const = const, n = n, m = m, q = q, random = random.G, 
                       info.level = info.level)
      G <- out.G$G
      K <- out.G$K
      bkmat <- bkmat[, out.G$classes.left]
      alphamat <- alphamat[, out.G$classes.left]
      if(K == 1) {
        bkmat <- matrix(bkmat, ncol = 1)
        alphamat <- matrix(alphamat, ncol = 1)
      }
      bwts2 <- bwts2[out.G$classes.left]
      Gmoves <- out.G$moves
    }
    
    ## Get criterion value
    overall.crit[iter] <- Lfun(a.cur = a.cur, bkmat = bkmat, G = G, Fr.cent = Fr.cent,
                               n = n, m = m, q = q, const = const, K = K)
    
    if(iter > 1 && abs(overall.crit[iter - 1] - overall.crit[iter]) < eps.G || Gmoves == 0) break
    if(iter == maxit) {
        # warning("k-means: Maximum number of iterations reached.")
        break
        }
    iter <- 1 + iter
  }
  crit.G <- overall.crit[iter]
  time.now <- proc.time()[3]
  if(info.level > 1) cat("= = = =\n")    
  if(info.level > 0) cat("Start", sprintf("%3d", i), 
                         " | loss = ", sprintf("%.10f", crit.G), " |", 
                         sprintf("%5.2f", round(time.now - time.start, 2)), "s | i =", 
                         sprintf(paste0("%", nchar(as.character(maxit)), "d"), iter), "\n") #, ifelse(iter == maxit, "***\n", "\n") )
  if(info.level >  1) cat("= = = =\n\n")
  
  list(G = G, K = K, minloss = crit.G, a = a.cur, bstar = out$b1, bkmat = bkmat, alphamat = t(alphamat), 
       time.G.start = time.now - time.start, iter.G = iter, loss.G.update = ifelse(update.G, out.G$kloss, Inf))
}
