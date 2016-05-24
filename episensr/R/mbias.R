#' Sensitivity analysis to correct for selection bias caused by M bias.
#'
#' Simple sensitivity analysis to correct for selection bias caused by M bias using
#' estimates of the odds ratios relating the variables.
#'
#' @param or Vector defining the input bias parameters, in the following order:
#' \enumerate{
#' \item Odds ratio between A and the exposure E,
#' \item Odds ratio between A and the collider C,
#' \item Odds ratio between B and the collider C,
#' \item Odds ratio between B and the outcome D,
#' \item Odds ratio observed between the exposure E and the outcome D.
#' }
#' @param var Vector defining variable names, in the following order:
#' \enumerate{
#' \item Outcome,
#' \item Exposure,
#' \item A,
#' \item B,
#' \item Collider.
#' }
#' 
#' @return A list with elements:
#' \item{mbias.parms}{Maximum bias parameters.}
#' \item{adj.measures}{Selection bias corrected measures.}
#' \item{bias.parms}{Input bias parameters.}
#'
#' @references Greenland S. Quantifying biases in causal models: classical
#' confounding vs. collider-stratification bias. Epidemiology 2003;14:300-6. 
#' @examples
#' mbias(or = c(2, 5.4, 2.5, 1.5, 1),
#' var = c("HIV", "Circumcision", "Muslim", "Low CD4", "Participation"))
#' @export
mbias <- function(or,
                  var) {
    if(is.null(or))
        stop('Missing input bias parameters.')
    else or <- or
    if(!is.vector(or))
        stop('The argument or should be a vector of length 5')
    if(length(or) != 5)
        stop('The argument or should be made of 5 components in the following order: (1) Odds ratio between A and the exposure E, (2) Odds ratio between A and the collider C, (3) Odds ratio between B and the collider C, (4) Odds ratio between B and the outcome D, and (5) Odds ratio observed between the exposure E and the outcome D.')
    if(!all(or >= 0))
        stop('Odds ratios should be greater than 0.')

    or.ae <- or[1]
    or.ac <- or[2]
    or.bc <- or[3]
    or.bd <- or[4]
    or.ed <- or[5]

    mbias.ce <- (or.ae * or.ac * (1 / (1 + sqrt(or.ae * or.ac))) +
                     1 - (1 / (1 + sqrt(or.ae * or.ac)))) /
        ((or.ae * (1 / (1 + sqrt(or.ae * or.ac))) +
              1 - (1 / (1 + sqrt(or.ae * or.ac)))) *
             (or.ac * (1 / (1 + sqrt(or.ae * or.ac))) +
                  1 - (1 / (1 + sqrt(or.ae * or.ac)))))
    mbias.cd <- (or.bc * or.bd * (1 / (1 + sqrt(or.bc * or.bd))) +
                     1 - (1/(1 + sqrt(or.bc * or.bd)))) /
        ((or.bc * (1 / (1 + sqrt(or.bc * or.bd))) +
              1 -(1 / (1 + sqrt(or.bc * or.bd)))) *
             (or.bd * (1 / (1 + sqrt(or.bc * or.bd))) +
                  1 - (1 / (1 + sqrt(or.bc * or.bd)))))
    mbias.ed <- (mbias.ce * mbias.cd * (1 / (1 + sqrt(mbias.ce * mbias.cd))) +
                     1 - (1 / (1 + sqrt(mbias.ce * mbias.cd)))) /
        ((mbias.ce * (1 / (1 + sqrt(mbias.ce * mbias.cd))) +
              1 - (1 / (1 + sqrt(mbias.ce * mbias.cd)))) *
             (mbias.cd * (1 / (1 + sqrt(mbias.ce * mbias.cd))) +
                  1 - (1 / (1 + sqrt(mbias.ce * mbias.cd)))))

    or.corr <- or.ed / mbias.ed
   
    res <- list(mbias.parms = c(mbias.ce, mbias.cd, mbias.ed),
                adj.measures = or.corr,
                bias.parms = or,
                labels = var)
    class(res) <- c("mbias", "list")
    res
}
