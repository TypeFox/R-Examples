#----------------------------------------------------------------------#
# Generate a random network where both the network structure and the   #
# partial correlation coefficients are random. The data matrices are   #
# generated from multivariate normal distribution with the covariance  #
# matrix corresponding to the network.                                 #
#----------------------------------------------------------------------#
# Inputs :                                                             #
#  G    The number of variables (vertices)                             #
#  etaA The proportion of non-null edges among all the G(G-1)/2 edges. #
#  n    The sample size.                                               #
#  r    The number of replicated G by n data matrices.                 #
#  dist A function which indicates the distribution of sample.         #
#       "mvnorm" is multivariate normal distribution and               #
#       "mvt" is multivariate t distribution with df=2.                #
#       The default is set by "mvnorm".                                #
# Outputs :                                                            #
#  data            A list of length rep, each element a simulated      #
#                  n x G data matrix.                                  #
#  true.partialcor The partial correlation matrix which the datasets   #
#                  are generated from.                                 #
#  truecor.scaled  The covariance matrix calculted from the partial    #
#                  correlation matrix.                                 #
#  sig.node        The indices of nonzero upper triangle elements of   #
#                  partial correlation matrix.                         #
#----------------------------------------------------------------------#
simulateData <- function(G, etaA, n, r, dist = "mvnorm") {

  expression <- list()
  dist <- tolower(dist)

  partialcov <- matrix(data = 0.0, nrow = G, ncol = G)

  tri.w <- which(upper.tri(partialcov))

  no.sig <-  ceiling(etaA * length(tri.w))

  sig.node <- sample(tri.w, no.sig)

  partialcov[sig.node] <- runif(n = no.sig, min = -1.0, max = 1.0)

  partialcov <- t(partialcov)  + partialcov

  diag(partialcov) <- colSums(abs(partialcov)) + 0.0001

  dpc <- diag(partialcov)
  partialcor <- partialcov / sqrt(dpc %o% dpc)
  invcor <- - partialcor
  diag(invcor) <- 1.0
  truecor <- try(solve(invcor), silent = TRUE)
  if( is(truecor, "try-error") ) {
    stop("Unable to invert correlation matrix.", call. = FALSE)
  }

  dtc <- diag(truecor)
  truecor.scaled <- truecor/ sqrt(dtc %o% dtc)

  if( dist == "mvnorm" ) {
    for( i in 1L:r ) {
      expression[[i]] <- rmvnorm(n = n, mean = rep(0,G), truecor)
    }
  } else if (dist == "mvt") {
    for( i in 1L:r ) {
      expression[[i]] <- rmvt(n = n, sigma = truecor, df = 2L)
    }
  }

  return(list("data" = expression,
              "true.partialcor" = partialcor,
              "truecor.scaled" = truecor.scaled,
              "sig.node" = sig.node))
}

