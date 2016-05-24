################################################################################
## Data Generation
################################################################################

#' Generates possible scores for a Zemel model
#' 
#' @param m Number of alternatives
#' @return a set of scores, all whose logs sum to 1
#' @export
#' @examples
#' Generate.Zemel.Parameters(10)
Generate.Zemel.Parameters <- function(m) {
  scores <- rnorm(m)
  list(Score = scores - log(sum(exp(scores))))
}

#' Generates pairwise ranks from a Zemel model given a set of scores
#' 
#' @param scores a vector of scores
#' @param m Number of alternatives
#' @param n Number of pairwise alternatives to generate
#' @return simulated pairwise comparison data
#' @export
#' @examples
#' scores <- Generate.Zemel.Parameters(10)$Score
#' Generate.Zemel.Ranks.Pairs(scores, 10, 10)

Generate.Zemel.Ranks.Pairs <- function(scores, m, n)
{
  Generate.Zemel.Pair <- function() {
    pair <- sample(m, 2)
    a <- pair[1]
    b <- pair[2]
    p <- exp(scores[a] - scores[b]) / (exp(scores[a] - scores[b]) + exp(scores[b] - scores[a]))
    if(rbinom(n = 1, prob = p, size = 1) == 0) pair <- c(b, a)
    pair
  }
  t(replicate(n, Generate.Zemel.Pair()))
}

################################################################################
## Likelihood
################################################################################

#' Gives Zemel pairwise Log-likelihood with data and scores
#' 
#' Calculates the log-likelihood in the pairwise Zemel model
#' 
#' @param Data.pairs data broken up into pairwise comparisons
#' @param Estimate Inference object from Estimate function
#' @return a log-likelihood of the data under the Zemel model
#' @export
#' @examples
#' Estimate <- Generate.Zemel.Parameters(10)
#' pairs <- Generate.Zemel.Ranks.Pairs(Estimate$Score, 10, 10)
#' Likelihood.Zemel(pairs, Estimate)
Likelihood.Zemel <- function(Data.pairs, Estimate) {
  m <- length(Estimate$Score)
  Tee <- nrow(Data.pairs)
  
  # in the case of the likelihood, I do not factor in the rank distances in the Count matrix
  Count.matrix <- generateC(Data.pairs, m, weighted = FALSE, normalized = FALSE)
  
  # calculate the denominator in the fraction
  denominator <- 0
  for(k in 1:m){
    for(l in 1:m){
      if(k != l)
        denominator <- denominator + exp(Estimate$Score[k] - Estimate$Score[l])
    }
  }
  
  # calculate the part that involves the scores
  scorepart <- 0
  for(i in 1:m){
    for(j in 1:m){
      if(i != j)
        scorepart <- scorepart + Count.matrix[i, j] * (log(exp(Estimate$Score[i] - Estimate$Score[j])) - log(denominator))
    }
  }
  
  lfactorial(Tee) - sum(lfactorial(Count.matrix)) + scorepart
}


################################################################################
## Parameter Estimation
################################################################################


#' Estimates Zemel Parameters via Gradient Descent
#' 
#' This function takes in data broken into pairs, and estimates
#' the parameters of the Zemel mode via Gradient Descent
#' 
#' @param Data.pairs data broken up into pairwise comparisons
#' @param m how many alternatives
#' @param threshold turning parameter for gradient descent
#' @param learning.rate turning parameter for gradient descent
#' @export
#' @return a set of scores for the alternatives, normalized such that the sum 
#' of the log scores is 0
#' scores <- Generate.Zemel.Parameters(10)$Score
#' pairs <- Generate.Zemel.Ranks.Pairs(scores, 10, 10)
#' Estimation.Zemel.MLE(pairs, 10, threshold = .1)
Estimation.Zemel.MLE = function(Data.pairs, m, threshold = 0.0001, learning.rate = 1/30000)
{
  
  # Tee is the total number of pairwise data points
  # spelled Tee instead of T since T is reserved in R
  Tee <- nrow(Data.pairs)
  
  t0 <- proc.time()
  
  # issues with using the rank-reweighting DOES NOT CONVERGE
  C <- generateC(Data.pairs, m, weighted = FALSE, normalized = FALSE)
  
  # random initialization of scores
  scores <- rnorm(m)
  oldscores <- rep(0, m)
  
  gradient <- function(i) {
    denominator <- 0
    numerator <- 0
    for(k in 1:m) for(l in 1:m) if(k != l) denominator <- denominator + exp(scores[k] - scores[l])
    for(j in 1:m) if(j != i) numerator <- numerator - exp(scores[j] - scores[i]) + exp(scores[i] - scores[j])
    sum(C[i,]) - sum(C[,i]) - Tee * numerator / denominator
  }
  
  # learning rate for gradient descent
  
  while(max(abs(scores - oldscores)) > threshold) {
    oldscores <- scores
    for(i in 1:m) scores[i] <- scores[i] + learning.rate * gradient(i)
    scores <- scores - log(sum(exp(scores)))
  }
  
  # Zemel scores are invariant to addition so I subtract a number from them that will make
  new.scores <- scores - log(sum(exp(scores)))
  list(m = m, order = order(-new.scores), Score = new.scores, Time = proc.time() - t0, Parameters = convert.vector.to.list(new.scores, name = "Score"))
  
}