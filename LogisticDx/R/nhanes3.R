#' @name nhanes3
#' @docType data
#' @title NHANES III data
#'
#' @format A \code{data.frame} with
#' \eqn{17030} observations (rows)
#' and \eqn{16} variables (columns).
#'
#' @details
#' A subset of data from the National Health and Nutrition
#' Examination Study (NHANES) III. Subjects age >=20 are included.
#' \cr
#' A sample of 39,695 subjects was selected, representing more
#' than 250 million people living in the USA. Data was collected
#' 1988-1994.
#' \cr \cr
#' 49 pseudo strata were created with 2 pseudo-PSU's in each
#' stratum (primary sampling units).
#' \cr \cr
#' This is a subset of the original dataset.
#' \cr \cr
#' Columns are:
#' \describe{
#'  \item{SEQN}{Respondent sequence number.}
#'  \item{SDPPSU6}{Pseudo-PSU (primary sampling unit).}
#'  \item{SDPSTRA6}{Pseudo stratum.}
#'  \item{WTPFHX6}{Statistical weight.
#'                 Range \eqn{225.93} to \eqn{139744.9}.}
#'  \item{HSAGEIR}{Age (years).}
#'  \item{HSSEX}{Gender (\code{factor}):
#'   \describe{
#'    \item{0}{female}
#'    \item{1}{male}}}
#'  \item{DMARACER}{Race (\code{factor}):
#'   \describe{
#'    \item{1}{white}
#'    \item{2}{black}
#'    \item{3}{other}}}
#'  \item{BMPWTLBS}{Body weight (lbs).}
#'  \item{BMPHTIN}{Standing height (inches).}
#'  \item{PEPMNK1R}{Average Systolic BP.}
#'  \item{PEPMNK5R}{Average Diastolic BP.}
#'  \item{HAR1}{Has respondent smoked >100 cigarettes
#'              in life (\code{factor}):
#'   \describe{
#'    \item{1}{yes}
#'    \item{2}{no}}}
#'  \item{HAR3}{Does respondent smoke cigarettes now?
#'              (\code{factor}):
#'   \describe{
#'    \item{1}{yes}
#'    \item{2}{no}}}
#'  \item{SMOKE}{Smoking (\code{factor}):
#'   \describe{
#'    \item{1}{never (HAR1 = 2)}
#'    \item{2}{>100 cigs (HAR1 = 1 & HAR3 = 2)}
#'    \item{3}{current (HAR1 =1 & HAR3 = 1)}}}
#'  \item{TCP}{Serum cholesterol (mg/100ml).}a
#'  \item{HBP}{High blood pressure? (\code{factor}):
#'   \describe{
#'    \item{1}{yes (PEPMNK1R > 140)}
#'    \item{2}{no (PEPMNK1R <= 140)}}}
#' }
#'
#' @note
#' Taken from:
#' \cr
#' ANALYTIC AND REPORTING GUIDELINES:
#' The Third National Health and Nutrition Examination Survey,
#' NHANES III (1988-94).
#' \cr \cr
#' In the NHANES III, 89 survey locations were randomly divided
#' into 2 sets or phases, the first consisting of
#' 44 and the other, 45
#' locations. One set of primary sampling units (PSUs) was allocated
#' to the first 3-year survey period (1988-91) and the other set to
#' the second 3-year period (1991-94).
#' \cr
#' Therefore, unbiased national
#' estimates of health and nutrition characteristics can be
#' independently produced for each phase as well
#' as for both phases
#' combined. Computation of national estimates from both phases
#' combined (i.e. total NHANES III) is the preferred option;
#' individual phase estimates may be highly variable. In addition,
#' individual phase estimates are not statistically independent.
#' \cr \cr
#' It is also difficult to evaluate whether
#' differences in individual
#' phase estimates are real or due to methodological
#' differences. That is, differences may be due to changes
#' in sampling methods or data
#' collection methodology over time.
#' At this time, there is no valid
#' statistical test for examining differences between phase 1 and
#' phase 2.
#' \cr \cr
#' NHANES III is based on a complex multistage probability sample
#' design. Several aspects of the NHANES design must be taken into
#' account in data analysis, including the sampling weights and the
#' complex survey design. Appropriate sampling weights are needed to
#' estimate prevalence, means, medians, and other statistics.
#' Sampling weights are used to produce correct population estimates
#' because each sample person does not have an equal probability of
#' selection. The sampling weights incorporate the differential
#' 3 probabilities of selection and include
#' adjustments for noncoverage
#' and nonresponse.
#' \cr \cr
#' With the large
#' oversampling of young children, older persons, black persons, and
#' Mexican Americans in NHANES III, it is essential that
#' the sampling
#' weights be used in all analyses. Otherwise, misinterpretation of
#' results is highly likely.
#' \cr \cr
#' Other aspects of the design that must be
#' taken into account in data analyses are the strata
#' and PSU pairings
#' from the sample design. These pairings should be used to estimate
#' variances and test for statistical significance.
#' \cr \cr
#' For weighted analyses, analysts can use special computer
#' software packages that
#' use an appropriate method for estimating variances for complex
#' samples such as SUDAAN (Shah 1995) and WesVarPC (Westat 1996).
#' \cr \cr
#' Although initial exploratory analyses may be performed on
#' unweighted data with standard
#' statistical packages assuming simple
#' random sampling, final analyses should be done on weighted data
#' using appropriate sampling weights.
#'
#' @keywords datasets
#'
#' @examples
#' ## use simpler column names
#' data("nhanes3", package="LogisticDx")
#' n1 <- c("ID", "pStrat", "pPSU", "sWt", "age", "sex",
#'         "race", "bWt", "h", "sysBP", "diasBP", "sm100",
#'         "smCurr", "smok", "chol", "htn")
#' names(nhanes3) <- n1
#'
#' @source
#' \href{ftp://ftp.wiley.com/public/sci_tech_med/logistic}{Wiley FTP}
#' @references
#' \bold{H&L 2nd ed.} Page 215. Table 6.3.
#'
#' National Center for Health Statistics (US) and others 1996.
#' NHANES III reference manuals and reports.
#' \emph{National Center for Health Statistics}.
#' \href{http://www.cdc.gov/nchs/nhanes/nh3rrm.htm}{CDC (free)}
NULL

