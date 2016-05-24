#' averisk: Calculation of Average Population Attributable Fractions and Confidence Intervals.
#'
#' Average population attributable fractions are calculated for a set of risk factors (either binary or ordinal valued) for both prospective and case-control designs.  Confidence intervals are found by Monte Carlo simulation. The method can be applied to either prospective or case control designs, provided an estimate of disease prevalence is provided.  In addition to an exact calculation of AF, an approximate calculation, based on randomly sampling permutations has been implemented to ensure the calculation is computationally tractable when the number of risk factors is large.   
#' 
#' @section averisk functions:
#' getAF
#'
#' @docType package
#' @name averisk
NULL
#> NULL