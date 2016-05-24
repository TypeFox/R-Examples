################################################################################
# Pure Rank Metrics
################################################################################

#' Calculates the Kendall Tau correlation between two ranks
#' 
#' @param rank1 two rankings. Order does not matter
#' @param rank2 two rankings. Order does not matter
#' @return The Kendall Tau correlation
#' @export
#' @examples
#' rank1 <- scramble(1:10)
#' rank2 <- scramble(1:10)
#' Evaluation.KendallTau(rank1, rank2)
Evaluation.KendallTau <- function(rank1, rank2) {
  cor(rank1, rank2, method="kendall")
}

#' Calculates the location of the True winner in the estimated ranking
#' 
#' @param EstimatedRank estimated ranking
#' @param TrueRank true ranking
#' @return The location of the true best in the estimated rank
#' @export
#' @examples
#' rank1 <- scramble(1:10)
#' rank2 <- scramble(1:10)
#' Evaluation.LocationofWinner(rank1, rank2)
Evaluation.LocationofWinner <-  function(EstimatedRank, TrueRank) {
  location.of.true.winner <- which(TrueRank == min(TrueRank))
  EstimatedRank[location.of.true.winner]
}

################################################################################
# Meta-Search Metrics
################################################################################
#From Paper 1-A Flexible Generative model for Preference Aggregation

#' Calculates the Normalized Discounted Cumluative Gain
#' 
#' @param EstimatedRank estimated ranking
#' @param RelevanceLevel score for the document
#' @return The NDCG for this estimation and relevance level
#' @export
#' @examples
#' EstimatedRank <- scramble(1:10)
#' RelevanceLevel <- runif(10)
#' Evaluation.NDCG(EstimatedRank, RelevanceLevel)
Evaluation.NDCG = function(EstimatedRank, RelevanceLevel) {
  Ordered.RelevanceLevel <- RelevanceLevel[order(-EstimatedRank)]
  Actual.RelevanceLevel <- RelevanceLevel[order(-RelevanceLevel)]
  dcg <- 0
  normalization.constant <- 0
  for(i in 1:length(Ordered.RelevanceLevel)) dcg <- dcg + (2^Ordered.RelevanceLevel[i] - 1) / log2(i+1)
  for(i in 1:length(Actual.RelevanceLevel)) normalization.constant <- normalization.constant + (2^Actual.RelevanceLevel[i] - 1) / log2(i+1)
  dcg / normalization.constant
}

#' Calculates the Average Precision
#' 
#' @param EstimatedRank estimated ranking
#' @param RelevanceLevel score for the document
#' @return The AP for this estimation and relevance level
#' @export
#' @examples
#' EstimatedRank <- scramble(1:10)
#' RelevanceLevel <- runif(10)
#' Evaluation.AveragePrecision(EstimatedRank, RelevanceLevel)
Evaluation.AveragePrecision <- function(EstimatedRank, RelevanceLevel) {
  Ordered.RelevanceLevel <- RelevanceLevel[order(-EstimatedRank)]
  numerator <- 0
  for(k in 1:length(Ordered.RelevanceLevel)) numerator <- numerator + Evaluation.Precision.at.k(EstimatedRank, RelevanceLevel, k) * Ordered.RelevanceLevel[k]
  numerator / sum(Ordered.RelevanceLevel)  
}

#' Calculates the Average Precision at k
#' 
#' @param EstimatedRank estimated ranking
#' @param RelevanceLevel score for the document
#' @param k positive that we want to run this algorithm for
#' @return The AP at k for this estimation and relevance level
#' @export
#' @examples
#' EstimatedRank <- scramble(1:10)
#' RelevanceLevel <- runif(10)
#' Evaluation.Precision.at.k(EstimatedRank, RelevanceLevel, 5)
Evaluation.Precision.at.k <- function(EstimatedRank, RelevanceLevel, k) {
  Ordered.RelevanceLevel <- RelevanceLevel[order(-EstimatedRank)]
  mean(Ordered.RelevanceLevel[1:k])
}

################################################################################
# Predictive Modeling Metrics
################################################################################

#' Calculates KL divergence between empirical pairwise preferences
#' and modeled pairwise preferences
#' 
#' @param Data.pairs data broken up into pairs using Breaking function
#' @param m number of alternatives
#' @param Estimate estimation object from an Estimate function
#' @param pairwise.prob Function that given two alternatives from 
#' the the Parameters argument, returns back a model probability that one is larger than the other 
#' @param prior prior weight to put in pairwise frequency matrix
#' @param nonparametric indicator that model is nonparametric (default FALSE)
#' @param ... additional arguments passed to generateC.model
#' @return the KL divergence between modeled and empirical pairwise preferences, 
#' thinking of the probabilities as a probability distribution over the (n choose 2) pairs
#' @export
#' @examples
#' data(Data.Test)
#' Data.Test.pairs <- Breaking(Data.Test, "full")
#' m <- 5
#' Estimate <- Estimation.PL.GMM(Data.Test.pairs, m)
#' Evaluation.KL(Data.Test.pairs, m, Estimate, PL.Pairwise.Prob)
Evaluation.KL <- function(Data.pairs, m, Estimate, pairwise.prob = NA, prior = 0, nonparametric = FALSE, ...) { 
  C.empirical <- generateC(Data.pairs, m, prior = prior)
  C.model     <- generateC.model(Estimate, pairwise.prob, nonparametric, ...)
  KL(C.empirical, C.model)
}

#' Calculates KL Divergence between non-diagonal entries of two matrices
#' 
#' @param A first matrix, this is the "true" distribution
#' @param B second matrix, this is the "estimated" distribution
#' @return KL divergence
#' @export
#' @examples
#' KL(matrix(runif(25), nrow=5), matrix(runif(25), nrow=5))
KL <- function(A, B) {
  non.diag    <- col(A) != row(A)
  P <- array(A[non.diag]) / sum(A[non.diag])
  Q <- array(B[non.diag]) / sum(B[non.diag])
  sum(log(P[P != 0]/Q[P != 0]) * P[P != 0])
}

#' Calculates MSE between non-diagonal entries of two matrices
#' if the diagonal elements are 0s
#' 
#' @param A first matrix
#' @param B second matrix
#' @return MSE divergence
#' @export
#' @examples
#' MSE(matrix(runif(25), nrow=5), matrix(runif(25), nrow=5))
MSE <- function(A, B) {
  m <- dim(A)[1]
  sum((A - B)^2) / (m^2 - m)
}

#' Calculates MSE between empirical pairwise preferences
#' and modeled pairwise preferences
#' 
#' @param Data.pairs data broken up into pairs using Breaking function
#' @param m number of alternatives
#' @param Estimate estimation object from an Estimate function
#' @param pairwise.prob Function that given two alternatives from 
#' @param prior prior weight to put in pairwise frequency matrix
#' @param nonparametric indicator that model is nonparametric (default FALSE)
#' the the Parameters argument, returns back a model probability that one is larger than the other 
#' @param ... additioanal parameters passed into generateC.model
#' @return the KL divergence between modeled and empirical pairwise preferences, 
#' thinking of the probabilities as a probability distribution over the (n choose 2) pairs
#' @export
#' @examples
#' data(Data.Test)
#' Data.Test.pairs <- Breaking(Data.Test, "full")
#' m <- 5
#' Estimate <- Estimation.PL.GMM(Data.Test.pairs, m)
#' Evaluation.MSE(Data.Test.pairs, m, Estimate, PL.Pairwise.Prob)
Evaluation.MSE <- function(Data.pairs, m, Estimate, pairwise.prob = NA, prior = 0, nonparametric = FALSE, ...) { 
  C.empirical <- generateC(Data.pairs, m, prior = prior)
  C.model     <- generateC.model(Estimate, pairwise.prob, nonparametric, ...)
  MSE(C.empirical, C.model)
}

#' Calculates TVD between empirical pairwise preferences
#' and modeled pairwise preferences
#' 
#' @param Data.pairs data broken up into pairs using Breaking function
#' @param m number of alternatives
#' @param Estimate estimation object from an Estimate function
#' @param pairwise.prob Function that given two alternatives from 
#' @param prior prior weight to put in pairwise frequency matrix
#' @param nonparametric indicator that model is nonparametric (default FALSE)
#' the the Parameters argument, returns back a model probability that one is larger than the other 
#' @param ... additional arguments passed to generateC.model
#' @return the TVD between modeled and empirical pairwise preferences, 
#' thinking of the probabilities as a probability distribution over the (n choose 2) pairs
#' @export
#' @examples
#' data(Data.Test)
#' Data.Test.pairs <- Breaking(Data.Test, "full")
#' m <- 5
#' Estimate <- Estimation.PL.GMM(Data.Test.pairs, m)
#' Evaluation.TVD(Data.Test.pairs, m, Estimate, PL.Pairwise.Prob)
Evaluation.TVD <- function(Data.pairs, m, Estimate, pairwise.prob = NA, prior = 0, nonparametric = FALSE, ...) { 
  C.empirical <- generateC(Data.pairs, m, prior = prior)
  C.model     <- generateC.model(Estimate, pairwise.prob, nonparametric, ...)
  TVD(C.empirical, C.model)
}

#' Calculates TVD between two matrices
#' 
#' @param A first matrix
#' @param B second matrix
#' @return Total variation distance
#' @export
#' @examples
#' TVD(matrix(runif(25), nrow=5), matrix(runif(25), nrow=5))
TVD <- function(A, B) {
  sum(abs(A - B))/2
}