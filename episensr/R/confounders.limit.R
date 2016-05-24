#' Bounding the bias limits of unmeasured confounding.
#'
#' Function to elicit the limits on measures of effect corrected for an unmeasured
#' confounder when only some of the bias parameters are known.
#'
#' @param p Proportion with the confounder among the unexposed group.
#' @param RR Relative risk between the confounder and the outcome.
#' @param OR Odds ratio between the confounder and the outcome.
#' @param crude.RR Crude relative risk between the exposure and the outcome.
#' @param dec Number of decimals in the printout.
#' @param print A logical scalar. Should the results be printed?
#' 
#' @return A list with elements:
#' \item{conf.limits}{Limits on confounding.}
#' \item{bias.parms}{Input bias parameters p, RR, OR, and crude RR.}
#'
#' @references
#' Lash, T.L., Fox, M.P, Fink, A.K., 2009 \emph{Applying Quantitative
#' Bias Analysis to Epidemiologic Data}, pp.59--78, Springer.
#'
#' Flanders, W. Dana, Khoury, Muin J., 1990. Indirect Assessment of
#' Confounding: Graphic Description and Limits on Effect of Adjusting for
#' Covariates. \emph{Epidemiology} 1(3): 239--246.
#' 
#' @examples
#' confounders.limit(OR = 1.65, crude.RR = 1.5)
#' @export
confounders.limit <- function(p = NA,
                              RR = NA,
                              OR = NA,
                              crude.RR = NULL,
                              dec = 4,
                              print = TRUE){
    if(is.null(crude.RR))
        stop('Please provide crude relative risk.')
    if(is.null(p) & is.null(RR) & is.null(OR))
        stop('Not enough information.')

    q <- ifelse(is.null(p), NA, 1 - p)
    lower.bound <- crude.RR /
        min(RR, OR, 1/p, RR/(q+RR*p), OR/(q+OR*p), na.rm = TRUE)
    upper.bound <- crude.RR

    if (print)
        cat("\nLimits on adjusted RR:", round(lower.bound, dec),
            "<= RRadj <=", round(upper.bound, dec), "\n")
    if (print)
        cat("\nInput Bias Parameters:",
            "\n----------------\n\n")
    if (print)
        cat("  p(Confounder+|Exposure-):", p,
            "\n    RR(Confounder-Disease):", RR,
            "\n   OR(Confounder-Exposure):", OR,
            "\nCrude RR(Exposure-Disease):", crude.RR, "\n")
    invisible(list(conf.limits = c(lower.bound, upper.bound),
                   bias.parms = c(p, RR, OR, crude.RR)))
}
