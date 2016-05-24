#' RADanalysis: A package for normalization of abundance tables to desired
#' number of ranks using MaxRank Normalization method.
#'
#' @description RADanalysis package has tools for normalizing rank abundance
#'     distributions (RAD) to a desired number of ranks using MaxRank
#'     Normalization method.
#'     RADs are commonly used in biology/ecology and mathematically equivalent
#'     to complementary cumulative distributions (CCDFs) which are used in
#'     physics, linguistics and sociology and more generally in data science.
#'
#' @section Rank Abundance Distributions (RAD):
#' Rank Abundance Distributions (RADs) are a way to capture the distribution
#' of biological species in communities, where we use the term "species" for
#' all types of distinct biological entities, e.g. microbial species in a
#' microbiome, viral strains in a quasi-species, the diverse variants B cells
#' in a person, etc. A RAD can be thought of as a plot with the two axes rank
#' (x-axis) and abundance (y-axis). For the most abundant species we draw a
#' point at the (x,y) coordinates  (1,a1) , with  a1  the abundance of this
#' most abundant species. For the second most abundant species we draw a point
#' at  (2,a2).
#'
#' @section MaxRank Normalization:
#' MaxRank normalization is the method to normalize RADs. MaxRank normalization
#' maps all rank abundance vectors to the same rank range from 1 to a common
#' maximum rank R. First we chose the maximum rank or "MaxRank" or "R". Second generated
#' for each sample s a pool of N_s of all individuals in s. From this pool we drew
#' individuals at random with
#' uniform probability and without replacement as long as the number of sampled
#' ranks of the original RAD did not exceed R. In this way we generated a new,
#' reduced abundance vector of R ranks, with a reduced number of individuals.
#' Division of these reduced abundances by sum of reduced abundances transforms the
#' reduced abundance
#' vector to a probability distribution for the R ranks with rank probabilities
#' summing up to 1. If R < total number of ranks in the original sample , the random
#' drawing of individuals from the pool
#' in general introduces a sampling error in the abundances. To control this error,
#' one should repeat the procedure several times (typically 10-100 times) and
#' averaged over all sampled abundance distributions.
#'
#'
#' @source Saeedghalati et al. 2016 "Quantitative comparison of abundance structures of genetic communities", submitted
#'
#' @docType package
#' @name RADanalysis
NULL

