# @title Log-likelihood of a CUB model with shelter effect
# @aliases loglikcubshe
# @description Compute the log-likelihood of a CUB model with a shelter effect
#  for the given absolute frequency distribution.
# @usage loglikcubshe(m, freq, pai1, pai2, csi, shelter)
# @param m Number of ordinal categories
# @param freq Vector of the absolute frequency distribution
# @param pai1 Mixing coefficient for the shifted Binomial component of the mixture distribution
# @param pai2 Mixing coefficient for the discrete Uniform component of the mixture distribution
# @param csi Feeling parameter
# @param shelter Category corresponding to the shelter choice
#' @keywords internal



loglikcubshe <-
function(m,freq,pai1,pai2,csi,shelter)
{t(freq)%*%log(probcubshe1(m,pai1,pai2,csi,shelter))}
