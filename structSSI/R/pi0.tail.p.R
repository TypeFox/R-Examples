  # This function estimates the value of pi0 given a set of p-values
  # using storey's tail p-value method. The idea is to set
  # pi0 = (# p vales > lambda)/(1 - lambda). See the final
  # report for an intuitive reason for this procedure.
  #
  # Input: 1) lambda - [single numeric between 0 and 1] - A parameter
  # that we will compare p-values to in order to find the proportion of true. null hypotheses. The
  # convention is to use lambda = 0.5. The only requirement is
  # that lambda lies between 0 and 1.
  # 2) p.values - [vector of numeric] - The unadjusted p-values
  # that resulted from a
  # multiple testing procedure.
  #
  # Output: 1) pi0 - An estimate of the proportion of null hypotheses
  # within the original p-values.


pi0.tail.p <- function(lambda, p.values){
    num <- length(which(p.values >= lambda))
    denom <- length(p.values)*(1 - lambda)
    return(min(num/denom, 1))
}
