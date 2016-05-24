#' Subject level data from SOLVD trial
#'
#' Dataset for use in \pkg{beanz} examples and vignettes.
#'
#' Subject level data from SOLVD trial. SOLVD is a randomized controlled
#' trial of the effect of an Angiotensin-converting-enzyme inhibitor (ACE
#' inhibitor) called enalapril on survival in patients with reduced left
#' ventricular ejection fraction and congestive heart failure (CHF).
#'
#' @name solvd.sub
#' @aliases solvd.sub
#'
#' @format A dataframe with 6 variables:
#' \describe{
#'   \item{trt}{treatment assignment}
#'   \item{y}{time to death or first hospitalization}
#'   \item{censor}{censoring status}
#'   \item{sodium}{level of sodium}
#'   \item{lvef}{level of lvef}
#'   \item{any.vasodilator.use}{level of use of vasodilator}
#' }
#'
#' @references
#'
#' Solvd Investigators and others, Effect of enalapril on survival in patients
#'     with reduced left ventricular ejection fraction and congestive heart
#'     failure. N Engl J Med. 1991, 325:293-302
#'
NULL
