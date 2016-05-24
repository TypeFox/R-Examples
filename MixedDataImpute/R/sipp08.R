#' Housholder earnings from the SIPP
#'
#' A dataset extracted from the first wave of the 2008 Survey of Income and Program Participation,
#' constructed  by retaining
#' all heads of household with positive earnings from employment and complete
#' cases. This is useful as a realistic population for studying imputation routines.
#'
#' @format A data frame with 30507 rows and 13 variables:
#' \describe{
#'   \item{total_earnings}{Total monthly earnings from employment (USD)}
#'   \item{age}{Age (years)}
#'   \item{sex}{Sex}
#'   \item{race}{Single race/ethnicity}
#'   \item{marital_status}{Marital Status}
#'   \item{born_us}{Born in or outside the United States}
#'   \item{own_kid}{Number of respondent's own children living in the household}
#'   \item{edu_level}{Level of education}
#'   \item{occ}{Occupation code for primary job (recoded from variable TJBOCC1 in the SIPP)}
#'   \item{worker_class}{Private/Non-profit/Government worker}
#'   \item{union}{Union member}
#'   \item{hourly}{Primary job is paid hourly}
#'   \item{hrs}{Usual hours worked per week}
#' }
#' @source Extracted from the SIPP public use files using Anthony Damico's scripts \url{http://www.asdfree.com/}
"sipp08"