#' An hypothetical dataset used to demonstrate functions.
#'
#' A hypothetical dataset with continuous and binary outcome.
#'
#'
#' @format Data frame with 270 rows and 5 variables:
#' \itemize{
#'   \item continuous: continuous outcome
#'   \item binary: binary outcome
#'   \item educ: education time(years)
#'   \item female: 0=male, 1=female
#'   \item treat: 0=control group, 1=treatment group
#' }
"tippingdata"



#'  Imputation results under different methods
#'
#' Imputation results based on Missing At Random(MAR) and Missing Completely At Random(MCAR) assumption for treatment and control group.
#' @format Data frame with 500 rows and 8 variables:
#' \itemize{
#'   \item MAR_T1: Average value of nonrespondents for continuous outcome in treatment group under MAR assumption.
#'   \item MAR_C1: Average value of nonrespondents for continuous outcome in control group under MAR assumption.
#'   \item MAR_T2: Number of success of nonrespondents for binary outcome in treatment group under MAR assumption.
#'   \item MAR_C2: Number of success of nonrespondents for binary outcome in control group under MAR assumption.
#'   \item MCAR_T1: Average value of nonrespondents for continuous outcome in treatment group under MCAR assumption.
#'   \item MCAR_C1: Average value of nonrespondents for continuous outcome in control group under MCAR assumption.
#'   \item MCAR_T2: Number of success of nonrespondents for binary outcome in treatment group under MCAR assumption.
#'   \item MCAR_C2:  Number of success of nonrespondents for binary outcome in control group under MCAR assumption.
#' }
#'
"imputedata"
