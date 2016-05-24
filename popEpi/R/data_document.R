


# sire - simulated survival data ------------------------------------------

#' sire - a simulated cohort of Finnish female rectal cancer patients
#'
#' \code{sire} is a simulated cohort pertaining female Finnish rectal cancer patients
#' diagnosed between 1993-2012. Instead of actual original dates, the dates are masked
#' via modest randomization within several time windows. 
#'
#' The closing date for the pertinent data was 2012-12-31, meaning status information was
#' available only up to that point --- hence the maximum possible \code{ex_date} is \code{2012-12-31}.
#'
#' @source The Finnish Cancer Registry
#' @format data.table with columns
#' \itemize{
#'  \item sex - gender of the patient (1 = female)
#'  \item bi_date - date of birth
#'  \item dg_date - date of cancer diagnosis
#'  \item ex_date - date of exit from follow-up (death or censoring)
#'  \item status  - status of the person at exit; 0 alive; 1 dead due to pertinent cancer; 2 dead due to other causes
#'  \item dg_age  - age at diagnosis expressed as fractional years
#' }
#' @author Karri Seppa
#' @name sire
NULL


# sibr - simulated survival data ------------------------------------------

#' sibr - a simulated cohort of Finnish female breast cancer patients
#'
#' \code{sibr} is a simulated cohort pertaining female Finnish breast cancer patients
#' diagnosed between 1993-2012. Instead of actual original dates, the dates are masked
#' via modest randomization within several time windows. The dataset is additionally
#' a random sample of 10 000 cases from the pertaining time window. 
#'
#' The closing date for the pertinent data was 2012-12-31, meaning status information was
#' available only up to that point --- hence the maximum possible \code{ex_date} is \code{2012-12-31}.
#'
#' @source The Finnish Cancer Registry
#' @format data.table with columns
#' \itemize{
#'  \item sex - gender of the patient (1 = female)
#'  \item bi_date - date of birth
#'  \item dg_date - date of cancer diagnosis
#'  \item ex_date - date of exit from follow-up (death or censoring)
#'  \item status  - status of the person at exit; 0 alive; 1 dead due to pertinent cancer; 2 dead due to other causes
#'  \item dg_age  - age at diagnosis expressed as fractional years
#' }
#' @author Karri Seppa
#' @name sibr
NULL




# International standard weights ------------------------------------------
#' Age standardisation weights from the ICSS scheme.
#'
#' Contains three sets age-standardisation weights for age-standardized survival (net, relative or observed).
#' 
#'
#' @source 
#' \href{http://seer.cancer.gov/stdpopulations/survival.html}{ICSS weights (US National Cancer Institute website)}
#' 
#' Corazziari, Isabella, Mike Quinn, and Riccardo Capocaccia. "Standard cancer patient population for age standardising survival ratios." European Journal of Cancer 40.15 (2004): 2307-2316.
#' @format data.table with columns
#' \itemize{
#'  \item age - lower bound of the age group
#'  \item ICSS1 - first set of weights, sums to 100 000
#'  \item ICSS2 - second set of weights, sums to 100 000
#'  \item ICSS3 - third set of weights, sums to 100 000
#' }
#' @name ICSS
#' @examples 
#' ## aggregate weights to a subset of age groups
#' data(ICSS)
#' cut <- c(0, 30, 50, 70, Inf)
#' agegr <- cut(ICSS$age, cut, right = FALSE)
#' aggregate(ICSS1~agegr, data = ICSS, FUN = sum)
NULL


# popmort -------------------------------------------------------------------
#' Population mortality rates in Finland 1951 - 2013 in 101 age groups and
#' by gender
#'
#'
#' @source Statistics Finland
#' @format \code{data.table} with columns
#' \itemize{
#'  \item sex gender coded as male, female (0, 1)
#'  \item year calendar year
#'  \item agegroup - coded 0 to 100; one-year age groups
#'  \item haz the average population mortality rate per person year
#' }
#' @name popmort
NULL




# stdpop18 ------------------------------------------------------------------

#' Standard populations from 2000: world, europe and nordic.
#'
#' World, europe, nordic standard populations by 18 age categories. Sums to 100000.
#'
#' @source Nordcan, 2000
#' @format data.table with columns
#' \itemize{
#'  \item agegroup (character)
#'  \item world std (numeric)
#'  \item europe std (numeric)
#'  \item nordic std (numeric)
#' }
#' @name stdpop18
NULL


# stdpop101 -----------------------------------------------------------------

#' World standard population by 1 year agegroups from 1 to 101. Sums to 100 000.
#'
#'
#' @source Standard population is from:
#' \href{http://seer.cancer.gov/stdpopulations/stdpop.singleages.html}{world standard pop "101of1"}
#' 
#' @format data.table with columns
#' \itemize{
#'  \item \code{world_std} weight in 100000 years (numeric)
#'  \item \code{agegroup} age group from 1 to 101 (numeric)
#' }
#' @name stdpop101
NULL


