#' Fits a Multi-Modal Mallows' model to ranking data.
#' 
#' Fits the Multi-Modal Mallows' model to partial or full ranking data, using
#' Kendall's metric and an EM algorithm.  This is essentially metric sequence
#' clustering.
#' 
#' @param datas Matrix of partial or fully-ranked data.
#' @param G Number of modes, 2 or greater.
#' @param iter Maximum number of iterations.
#' @param hyp Hypothesis sequence vector, to initialize one of the cluster
#' centers at.
#' @param plot.like Should the likelihood be printed at each iteration?
#' @return See output of FormatOut
#' @author Erik Gregory
#' @references "Mixtures of distance-based models for ranking data". Thomas 
#' Brendan Murphy & Donal Martin. 1 April 2002. Computational Statistics & 
#' Data Analysis 41 (2003) 645-655.
#' @keywords cluster Mallow
 
Mallows <-
function(datas, G, iter = 10,hyp = NULL, plot.like = FALSE) {
  top <- 20
  # Number of subjects.
  N <- nrow(datas)
  # Number of items being ranked.
  abils <- ncol(datas)
  
  dists.table <- DistanceDistribution(abils)
  # Initialize the p-value of membership in each cluster.
  p <- rep(1/G, G)
  # Initialize the modal sequences.
  R <- Rgen(G, hyp, abils)
  # Initialize the lambda values.
  lambda <- runif(G)
  cat("Solving...\n")
  likelihood <- 0*(1:iter)
  best.like <- 0
  all.dists.data <- NULL
  if (plot.like == TRUE) {
  	x11()
  }
  i <- 1
  infos <- KendallInfo(datas)
  while (i <= iter) {
    # Calculate the normalizing coefficients for 
    # lambda
    C.lam <- unlist(lapply(lambda, 
                           function(i) C_lam(i, dists.table = dists.table)))
    # E Step
    z <- EStep(R, datas, p, lambda, G, N, C.lam, all.dists.data)
    # M Step
    #R <- UpdateR(R, seqs, datas, z[[i]], G, all.dists.data)
    R <- UpdateR(datas, z, infos)

    p <- UpdateP(z)

    all.dists.data <- AllKendall(datas, do.call("rbind", R), infos)
    lambda <- UpdateLambda(datas, R, z, G, all.dists.data,
                           dists.table = dists.table)
    likelihood[i] <- Likelihood(z, p, C.lam, lambda, all.dists.data)
    if (plot.like == TRUE) {
	    if (i %%20 == 0) {
	      top <- top + 20
	    }
	    plot(1:i, likelihood[1:i], main = "Likelihood", xlab = "Iteration", ylab = "Log-Likelihood",
		 xlim = c(0, top), "l", cex = 0.5, col = "red")
    }
    if (i > 2) {
      if (likelihood[i] == likelihood[i -1]) {
        print("Algorithm converged")
        i <- iter
      }
    }
    i <- i + 1
  }
  cat("Formatting output")
  out <- FormatOut(R, p, lambda, z, datas, likelihood)
  if (plot.like == TRUE) {
     dev.off()
  }
  return(out)
}
