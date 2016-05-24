

############# documentation ###############

#' EffectLiteR
#'
#' @name EffectLiteR
#' @docType package
NULL


#' Dataset nonortho.
#' 
#' A simulated dataset. The variables are:
#' 
#' \itemize{
#'   \item y. Continuous dependent variable depression.
#'   \item x. Treatment variable with values 0 (control), 1 (treat1), and 2 (treat2).
#'   \item z. Categorical covariate with values 0 (low neediness), 1 (medium neediness) and 2 (high neediness).
#' }
#' 
#' @docType data
#' @keywords datasets
#' @format A data frame with 500 rows and 3 variables
#' @name nonortho
NULL



#' Dataset example01.
#' 
#' A simulated dataset. The variables are:
#' 
#' \itemize{
#'   \item x. Treatment variable with values control, treat1, and treat2.
#'   \item k1. Categorical covariate with values male and female.
#'   \item kateg2. Categorical covariate with values 1 and 2.
#'   \item z1-z3. Continuous covariates.
#'   \item dv. Coninuous dependent variable.
#' }
#' 
#' @docType data
#' @keywords datasets
#' @format A data frame with 2000 rows and 7 variables.
#' @name example01
NULL



#' Dataset example02lv.
#' 
#' A simulated dataset with latent variables. The variables are:
#' 
#' \itemize{
#'   \item CPM11. First indicator of latent covariate.
#'   \item CPM21. Second indicator of latent covariate.
#'   \item CPM12. First indicator of latent outcome.
#'   \item CPM22. Second indicator of latent outcome.
#'   \item x. Dichotomous treatment variable with values 0 (control), and 1 (treatment).
#'   \item k. Categorical covariate with values first, second, and third.
#' }
#' 
#' @docType data
#' @keywords datasets
#' @format A data frame with 300 rows and 6 variables.
#' @name example02lv
NULL


#' Dataset example_multilevel.
#' 
#' A simulated dataset with a cluster ID and sampling weights to test multilevel options. The variables are:
#' 
#' \itemize{
#'   \item y. Coninuous dependent variable.
#'   \item x. Treatment variable with values 0, 1.
#'   \item z. Continuous covariate.
#'   \item xz. Product of x and z.
#'   \item cid. Cluster ID.
#'   \item weights. Sampling weights.
#'   \item iptw. Classic inverse probability of treatment weights based on a logistic regression of x on z. Use with care (only for average effects).
#' }
#' 
#' @docType data
#' @keywords datasets
#' @format A data frame with 800 rows and 7 variables.
#' @name example_multilevel
NULL




############## namespace ###########

#' @importFrom methods new is
NULL

#' @importMethodsFrom methods show 
NULL

#' @importFrom stats as.formula ftable model.frame model.matrix pnorm relevel var
NULL

#' @importFrom utils capture.output read.csv read.table combn
NULL



