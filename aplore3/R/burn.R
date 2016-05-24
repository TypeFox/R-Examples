#' BURN1000 data
#' 
#' burn1000 dataset.
#' 
#' @format A data.frame with 1000 rows and 9 variables:
#' \describe{
#' \item{id}{Identification code (1 - 1000)}
#' \item{facility}{Burn facility (1 - 40)}
#' \item{death}{Hospital discharge status (1: Alive, 2: Dead)}
#' \item{age}{Age at admission (Years)}
#' \item{gender}{Gender (1: Female, 2: Male)}
#' \item{race}{Race (1: Non-White, 2: White)}
#' \item{tbsa}{Total burn surface area (0 - 100\%)}
#' \item{inh_inj}{Burn involved inhalation injury (1: No, 2: Yes)}
#' \item{flame}{Flame involved in burn injury (1: No, 2: Yes)}
#' }
#' @source Hosmer, D.W., Lemeshow, S. and Sturdivant, R.X. (2013) Applied
#' Logistic Regression, 3rd ed., New York: Wiley
#' @examples
#' head(burn1000, n = 10)
#' summary(burn1000)
#'
#' ## Table 3.15 p. 80
#' summary(mod3.15 <- glm(death ~ tbsa + inh_inj + age + gender + flame + race,
#'                        family = binomial, data = burn1000 ))
"burn1000"


#' BURN_EVAL_1 data
#' 
#' burn_eval_1 dataset.
#' 
#' @format A data.frame with 500 rows and 9 variables: the covariate are
#' the same as those from \code{\link{burn1000}}.
#' @source Hosmer, D.W., Lemeshow, S. and Sturdivant, R.X. (2013) Applied
#' Logistic Regression, 3rd ed., New York: Wiley
#' @examples
#' head(burn_eval_1, n = 10)
#' summary(burn_eval_1)
"burn_eval_1"

#' BURN_EVAL_2 data
#' 
#' burn_eval_2 dataset.
#' 
#' @format A data.frame with 500 rows and 9 variables: the covariate are
#' the same as those from \code{\link{burn1000}}.
#' @source Hosmer, D.W., Lemeshow, S. and Sturdivant, R.X. (2013) Applied
#' Logistic Regression, 3rd ed., New York: Wiley
#' @examples
#' head(burn_eval_2, n = 10)
#' summary(burn_eval_2)
"burn_eval_2"

#' BURN13M data
#' 
#' burn13m dataset.
#' 
#' @format A data.frame with 388 rows and 11 variables: the covariate are
#' the same as those from \code{\link{burn1000}} with the addition of
#' \describe{
#' \item{pair}{Pair Identification Code (1-119)}
#' \item{pairid}{Subject Identification Code within pair (1-4)}
#' }
#' @source Hosmer, D.W., Lemeshow, S. and Sturdivant, R.X. (2013) Applied
#' Logistic Regression, 3rd ed., New York: Wiley
#' @examples
#' head(burn13m, n = 10)
#' summary(burn13m)
"burn13m"
