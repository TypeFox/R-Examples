#' Fit a mixture weighted distance-based model
#' 
#' This function computes fitting of mixture weighted distance-based 
#' model for the given data set of complete rankings.
#' 
#' @param dset data set of complete rankings
#' @param G number of clusters
#' @param dset.agg whether the data set is in the aggregated form
#'                 (default as FALSE)
#' @param dtype type of the weighted distance measure 
#' Kendall or K(default) : "Weighted Kendall's tau", SqrtSpearman 
#' or SS : "Square root of weighted Spearman", Spearman or S : 
#' "Weighted Spearman", Footrule or F : "Weighted Spearman's 
#' footrule" 
#' @param noise whether a noise cluster is contained (default as FALSE)
#' @param iter number of iterations of the EM algorithm (default as 100)
#' @return a list of the fitting result, containing the following objects:
#'         $clusterNum number of clusters (excluding the noise)
#'         $dtype type of the distance measure
#'         $noise whether a noise cluster is contained
#'         $iterNum actual number of iterations of the EM algorithm
#'         $convergence whether the complete-data loglikelihood converges
#'         $clusterProb probability of each cluster
#'         $modalRank modal rankings
#'         $weight weight vectors for clusters
#'         $trueLoglik the true loglikelihood by the fitted model
#'         $squaredPearsonStat the sum of squares of Pearson residuals
#' @export
#' @author Yumin Zhang <zymneo@@gmail.com>
#' @examples
#' data(Croon)
#' # Time comu
#' # mwdbm(Croon, 3)
mwdbm <- function(dset, G, dset.agg=TRUE, 
                  dtype="Kendall", noise=FALSE, iter=100)
{
  # Fit a mixture weighted distance-based model
  # 
  # This function computes fitting of mixture weighted distance-based 
  # model for the given data set of complete rankings.
  # Author 	: Yumin Zhang
  # Email 	: zymneo@gmail.com
  # Created	: Oct 16, 2013
  #
  # 1. Preparation
  #
  # The checking for (1) single-observation data set, (2) the matrix type of 
  # dset, and (3) the too-many #clusters is omitted here for simplicity.
  if (!dset.agg) {
    dset.original <- dset
    dset <- rankagg(dset.original)
  }
  rownames(dset) <- rownames(dset, do.NULL=FALSE, prefix="Obs.")
  colnames(dset) <- colnames(dset, do.NULL=FALSE, prefix="pi")
  colnames(dset)[ncol(dset)] <- "n"
  nitem <- ncol(dset) - 1
  dset.unique <- dset[, 1:nitem]
  nobs.unique <- nrow(dset)
  freq <- dset[, ncol(dset)]
  nobs <- sum(freq)
  if (!noise) {
    nclst <- G
  } else {
    nclst <- G + 1
  }
  seq.cplt <- generate.perms(nitem)
  colnames(seq.cplt) <- colnames(seq.cplt, do.NULL=FALSE, prefix="pi")
  # Containers for calculation
  cd.loglikelihood <- rep(0, iter)
  convergence <- FALSE
  z <- matrix(0, nrow=nclst, ncol=nobs.unique)
  rownames(z) <- rownames(z, do.NULL=FALSE, prefix="G")
  colnames(z) <- colnames(z, do.NULL=FALSE, prefix="Obs.")
  dmat <- matrix(0, nrow=nclst, ncol=nobs.unique)
  rownames(dmat) <- rownames(dmat, do.NULL=FALSE, prefix="G")
  colnames(dmat) <- colnames(dmat, do.NULL=FALSE, prefix="Obs.")
  dmat.pc <- matrix(0, nrow=nclst, ncol=factorial(nitem))
  rownames(dmat.pc) <- rownames(dmat.pc, do.NULL=FALSE, prefix="G")
  colnames(dmat.pc) <- colnames(dmat.pc, do.NULL=FALSE, prefix="Seq.")
  pc <- rep(0, nclst)
  p.dist <- matrix(0, nrow=nclst, ncol=nobs.unique) # may not be used
  rownames(p.dist) <- rownames(p.dist, do.NULL=FALSE, prefix="G")
  colnames(p.dist) <- colnames(p.dist, do.NULL=FALSE, prefix="Obs.")
  p.clst <- rep(1 / nclst, nclst)
  w <- matrix(runif(nclst * nitem, min=0, max=1), nrow=nclst, ncol=nitem)
  rownames(w) <- rownames(w, do.NULL=FALSE, prefix="G")
  colnames(w) <- colnames(w, do.NULL=FALSE, prefix="w")
  if (noise) {
    w[nclst, ] <- rep(0, nitem)
    pc[nclst] <- factorial(nitem)
  }
  sorted.meanrank <- 1:nitem
  meanrank <- destat(dset)$mean.rank
  for (i in 1:nitem) {
    sorted.meanrank[which(meanrank == sort(meanrank)[i])] <- i
  }
  Pi <- matrix(sorted.meanrank, nrow=G, ncol=nitem, byrow=T)
  rownames(Pi) <- rownames(Pi, do.NULL=FALSE, prefix="G")
  colnames(Pi) <- colnames(Pi, do.NULL=FALSE, prefix="Item")
  #
  # 2. Function definition
  #
  cd.loglik.clst <- function(weight, modal, Gindex) { 
    # Note that for the use of optimization functions in R, we define the 
    # complete-data loglikelihood computed by this function to be negative.
    clst.dmat <- wdmat(dset1=dset.unique, dset2=modal, 
                       dset1.agg=FALSE, dset2.agg=FALSE,
                       dtype=dtype, weight=weight, modal=modal)$mat
    clst.dmat.cplt <- wdmat(dset1=seq.cplt, dset2=modal, 
                            dset1.agg=FALSE, dset2.agg=FALSE,
                            dtype=dtype, weight=weight, modal=modal)$mat
    clst.pc <- sum(exp(-clst.dmat.cplt))
	penalty <- ifelse(any(weight < 0), 100000000, 0)
    cd.loglik <- -(sum(freq * z[Gindex, ]) * log(p.clst[Gindex] / clst.pc)
                   - sum(freq * z[Gindex, ] * clst.dmat)) + penalty
    return(cd.loglik)
  }
  #
  # 3. Iterations
  #
  for (g in 1:G) {
    dmat[g, ] <- wdmat(dset1=dset.unique, dset2=Pi[g, ],
                       dset1.agg=FALSE, dset2.agg=FALSE,
                       dtype=dtype, weight=w[g, ], modal=Pi[g, ])$mat
    dmat.pc[g, ] <- wdmat(dset1=seq.cplt, dset2=Pi[g, ],
                          dset1.agg=FALSE, dset2.agg=FALSE,
                          dtype=dtype, weight=w[g, ], modal=Pi[g, ])$mat
  }
  pc <- rowSums(exp(-dmat.pc))
  i <- 1
  while (i <= iter) {
    exp.dmat <- exp(-dmat)
    p.dist <- apply(exp.dmat, 2, function(arg) arg / pc)
    z <- apply(p.dist, 2, function(arg) (arg * p.clst) / sum(arg * p.clst))
	p.clst <- apply(z, 1, function(arg) sum(arg * freq)) / nobs
    for (g in 1:G) {
      seq.cplt.cdlc <- apply(seq.cplt, 1, cd.loglik.clst, 
                             weight=w[g, ], Gindex=g)
      Pi[g, ] <- seq.cplt[which(seq.cplt.cdlc == min(seq.cplt.cdlc)), ]
      w.update <- optim(w[g, ], cd.loglik.clst, gr=NULL, modal=Pi[g, ],
	               Gindex=g, method="Nelder-Mead", hessian=FALSE) 
      w[g, ] <- w.update$par
      seq.cplt.cdlc <- apply(seq.cplt, 1, cd.loglik.clst, 
                             weight=w[g, ], Gindex=g)
      Pi[g, ] <- seq.cplt[which(seq.cplt.cdlc == min(seq.cplt.cdlc)), ]
      dmat[g, ] <- wdmat(dset1=dset.unique, dset2=Pi[g, ],
                         dset1.agg=FALSE, dset2.agg=FALSE,
                         dtype=dtype, weight=w[g, ], modal=Pi[g, ])$mat
      dmat.pc[g, ] <- wdmat(dset1=seq.cplt, dset2=Pi[g, ],
                            dset1.agg=FALSE, dset2.agg=FALSE,
                            dtype=dtype, weight=w[g, ], modal=Pi[g, ])$mat
    }
    pc <- rowSums(exp(-dmat.pc))
    cd.loglikelihood[i] <- -sum(sapply(1:nclst, function(arg) 
                                      cd.loglik.clst(w[arg, ], Pi[arg, ], arg)))
    cd.loglikelihood.final <- cd.loglikelihood[i]
	niter <- i
    if (i > 1) {
      if (abs((cd.loglikelihood[i] - cd.loglikelihood[i - 1]) 
               / cd.loglikelihood[i - 1]) < (1e-6)) {
        convergence <- TRUE
        i <- iter
      }
    }
    i <- i + 1
  }
  true.loglikelihood <- sum(log(apply(exp(-dmat), 2, function(arg) 
                                      sum((arg / pc) * p.clst))) * freq)
  # compute the sum of squares Pearson residuals
  expect.cplt <- apply(exp(-dmat.pc), 2, 
                       function(arg) sum((arg / pc) * p.clst)) * nobs
  freq.cplt <- rep(0, factorial(nitem))
  for (a in 1:factorial(nitem)) {
    for (b in 1:nobs.unique) {
      if (all(seq.cplt[a, ] == dset.unique[b, ])) {
        freq.cplt[a] <- freq[b]
      }
    }
  }
  ss <- sum(((freq.cplt - expect.cplt)^2) / expect.cplt)
  #
  # 4. Output
  #
  message("Maximum likelihood estimation of the mixture weighted distance-based model")
  dtype.full <- switch(dtype, 
                       SqrtSpearman = "Square root of weighted Spearman", 
                       SS           = "Square root of weighted Spearman", 
                       Spearman     = "Weighted Spearman", 
                       S            = "Weighted Spearman", 
                       Footrule     = "Weighted Spearman's footrule", 
                       F            = "Weighted Spearman's footrule", 
                       Kendall      = "Weighted Kendall's tau", 
                       K            = "Weighted Kendall's tau")
  message("Distance type: ", dtype.full, ".")
  message("Number of clusters: ", G) 
  noise.YN <- ifelse(noise, "Yes", "No")
  message("With additional noice cluster: ", noise.YN)
  cov.YN <- ifelse(convergence, "Yes", "No")
  message("Convergent or not: ", cov.YN)
  message("Number of iterations: ", niter)
  message("Cluster probability: ")
  print(p.clst)
  message("\nModal ranking: ")
  print(Pi)
  message("\nWeights: ")
  print(w)
  message("\nComplete-data loglikelihood: ")
  print(cd.loglikelihood.final)
  message("True loglikelihood: ")
  print(true.loglikelihood)
  message("\nThe sum of squared Pearson residuals: ")
  print(ss)
  mwdbm.fit <- list(clusterNum=G, distanceType=dtype.full, noise=noise.YN,
                    iterNum=niter, convergence=cov.YN, clusterProb=p.clst, 
					modalRank=Pi, weight=w, trueLoglik=true.loglikelihood,
                    squaredPearsonStat=ss)
  return(mwdbm.fit)
}
