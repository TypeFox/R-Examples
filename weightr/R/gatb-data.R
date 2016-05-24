#' Studies of the Predictive Validity of the General Ability Subscale of the
#' General Aptitude Test Battery (GATB)
#'
#' Results from 755 studies on the General Aptitude Test Battery's predictive validity
#' of job perfomance (General Ability subscale).
#'
#' @docType data
#'
#' @usage dat.gatb
#'
#' @format A data frame containing the following columns:
#' \describe{
#'  \item{\code{z}}{Fisher's z-transformed correlation coefficients}
#'  \item{\code{v}}{corresponding sampling variance}
#'  }
#'
#' @keywords datasets
#'
#' @details The General Aptitude Test Battery (GATB) is designed to measure nine
#' cognitive, perceptual, and psychomotor skills thought relevant to the prediction
#' of job performance. From 1947 to 1993, a total of 755 studies were completed in
#' order to assess the validity of the GATB and its nine scales, and the GATB has
#' been found to be a moderately valid predictor of job performance. This dataset
#' consists of validity coefficients for the General Ability scale of the GATB.
#'
#' @references Vevea, J. L., Clements, N. C., & Hedges, L. V. (1993). Assessing the
#' effects of selection bias on validity data for the General Aptitude Test Battery.
#' Journal of Applied Psychology, 78(6), 981-987.
#'
#' U.S. Department of Labor, Division of Counseling and Test Development, Employment
#' and Training Administration. (1983a). The dimensionality of the General Aptitude
#' Test Battery (GATB) and the dominance of general factors over specific factors in
#' the prediction of job performance for the U.S. Employment Service (U.S. Employment
#' Service Test Research Rep. No. 44). Washington, DC.
#'
#' U.S. Department of Labor, Division of Counseling and Test Development, Employment
#' and Training Administration. (1983b). Test validity for 12,000 jobs: An application
#' of job classification and validity generalization analysis to the General Aptitude
#' Test Battery (U.S. Employment Service Test Research Rep. No. 45). Washington, DC.
#'
#' @source U.S. Department of Labor, Division of Counseling and Test Development,
#'  Employment and Training Administration. (1983a). The dimensionality of the
#'  General Aptitude Test Battery (GATB) and the dominance of general factors over
#'  specific factors in the prediction of job performance for the U.S. Employment
#'  Service (U.S. Employment Service Test Research Rep. No. 44). Washington, DC.
#'
#'  U.S. Department of Labor, Division of Counseling and Test Development, Employment
#'  and Training Administration. (1983b). Test validity for 12,000 jobs: An application
#'  of job classification and validity generalization analysis to the General Aptitude
#'  Test Battery (U.S. Employment Service Test Research Rep. No. 45). Washington, DC.
#'
#' @examples
#' \dontrun{
#' dat.gatb
#' effect <- dat.gatb$z
#' v <- dat.gatb$v
#' weightfunct(effect, v)
#' }
"dat.gatb"
