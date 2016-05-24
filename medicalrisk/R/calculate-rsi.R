# Code for calculating Sessler RSI

#' Returns the covariate coefficient for a particular diagnosis
#' or procedure code, along with the actual code that was found in
#' the internal database of coefficients. If a child code is supplied
#' but its parent is in the database, the coefficient for the parent
#' will be returned, along with that parent code. For example, if
#' D1231 is supplied but only D123 is available, D123 will be used.
#' This is apparently how the SPSS sample code works.
#' 
#' @param code A single ICD-9-CM code
#' @param betalist One of the rsi_beta_* datasets (supplied with this package)
#' @return Covariate coefficient. You must sum all of these for a given patient and then
#'     subtract the appropriate population beta (e.g. rsi_beta_1yrpod$popbeta)
#' @importFrom hash has.key
#' @examples
#' # get coefficient for hypercholesterolemia
#' sessler_get_single_beta('D2720', rsi_beta_inhosp)
#' # Also works with extra 0 on the end
#' sessler_get_single_beta('D27200', rsi_beta_inhosp)
#' @export
sessler_get_single_beta <- function(code, betalist) {
  code <- as.character(code)
  if(nchar(code) < 4) return(list(coeff=0, code=NULL)) # Sanity check
  
  # Actually find the entry. For some reason R does not have
  # real hashes built in, so we have to resort to the 'hash'
  # add-on package available in CRAN.
  if(has.key(code, betalist)) {
    return(list(coeff=betalist[[code]], code=code))
  } else {
    # Didn't find a match? Perhaps only this code's parent is in the
    # betalist. Chop off a character and try again. We do this down
    # to 4 characters.
    if(nchar(code) > 4) {
      code <- substr(code, 1, nchar(code)-1)
      return(sessler_get_single_beta(code, betalist))
    }
    # Couldn't find it, so give up and return zero
    return(list(coeff=0, code=NULL))
  }
}

#' Returns composite Sessler risk stratification index, given a list of ICD-9-CM codes.
#'
#' ICD-9-CM codes must have periods removed.  Diagnostic codes are prefixed with 
#' 'D' while procedure codes are prefixed with 'P'. So, diagnostic code
#' \code{404.03} should be \code{"D40403"}.
#' 
#' Note: A subsequent publication (Sigakis, 2013) found the following:
#' "Calibration "in-the-large" for RSI in-hospitalmortality illustrated a discrepancy 
#' between actual (1.5%) and predicted (51.7%) in-hospital mortality. The authors 
#' identified a regression constant (-2.198) in the published RSI 
#' "all-covariates.xls" file that was not used in the published SPSS implementation 
#' macro."
#' 
#' @author Tom Joseph <thomas.joseph@@mountsinai.org>, 
#' Patrick McCormick <patrick.mccormick@@mountsinai.org>
#' @param icd9 a unique character vector of ICD-9-CM codes
#' @return The risk stratification index score
#' @references 1. Sessler DI, Sigl JC, Manberg PJ, Kelley SD, Schubert A, Chamoun NG.
#'     Broadly applicable risk stratification system for predicting duration of hospitalization and mortality.
#'     Anesthesiology. 2010 Nov;113(5):1026-37. doi: 10.1097/ALN.0b013e3181f79a8d.   
#' 
#' 2. Sigakis MJG, Bittner EA, Wanderer JP: Validation of a risk stratification 
#'     index and risk quantification index for predicting patient outcomes: 
#'     in-hospital mortality, 30-day mortality, 1-year mortality, and length-of-stay. 
#'     Anesthesiology 2013; 119:525-40
#' @examples
#' # Calculate RSI for each patient ("id") in dataframe
#' cases <- data.frame(id=c(1,1,1,2,2,2),
#'   icd9cm=c("D20206","D24220","D4439","D5064","DE8788","D40403"))
#' library(plyr)
#' ddply(cases, .(id), function(x) { icd9cm_sessler_rsi(x$icd9cm) } )
#' @export
icd9cm_sessler_rsi <- function(icd9) {
  icdlist <- unique(icd9)    
  data.frame(
    rsi_1yrpod = sum(sapply(sapply(icdlist, sessler_get_single_beta, medicalrisk::rsi_beta_1yrpod)['coeff',], function(x) x[[1]])) - medicalrisk::rsi_beta_1yrpod$popbeta,
    rsi_30dlos = sum(sapply(sapply(icdlist, sessler_get_single_beta, medicalrisk::rsi_beta_30dlos)['coeff',], function(x) x[[1]])) - medicalrisk::rsi_beta_30dlos$popbeta,
    rsi_30dpod = sum(sapply(sapply(icdlist, sessler_get_single_beta, medicalrisk::rsi_beta_30dpod)['coeff',], function(x) x[[1]])) - medicalrisk::rsi_beta_30dpod$popbeta,
    rsi_inhosp = sum(sapply(sapply(icdlist, sessler_get_single_beta, medicalrisk::rsi_beta_inhosp)['coeff',], function(x) x[[1]])) - medicalrisk::rsi_beta_inhosp$popbeta
  ) 
}

#' Validates this Sessler RSI implementation against reference data
#' 
#' Requires that "sample data rev2.csv" and "sample results rev2.csv" be available 
#' in datasrc directory
#'
#' @author Patrick McCormick <patrick.mccormick@@mountsinai.org>
#' @return Table of patients with scores >0.001 difference.
#' @seealso rsi_sample_data, rsi_sample_results
#' @references 1. Sessler DI, Sigl JC, Manberg PJ, Kelley SD, Schubert A, Chamoun NG.
#'     Broadly applicable risk stratification system for predicting duration of hospitalization and mortality.
#'     Anesthesiology. 2010 Nov;113(5):1026-37. doi: 10.1097/ALN.0b013e3181f79a8d.
#' @importFrom plyr ddply
#' @importFrom reshape2 melt
#' @examples
#' \dontrun{
#' verify_sessler_rsi()
#' }
#' @export
verify_sessler_rsi <- function() {
  sessler_sample_melted <- melt(medicalrisk::rsi_sample_data, c("patient_id"), value.name="icd9cm")
  sessler_medicalrisk_results <- ddply(sessler_sample_melted, "patient_id", function(x) { icd9cm_sessler_rsi(x$icd9cm) })
  test_grid <- merge(sessler_medicalrisk_results,
                     medicalrisk::rsi_sample_results[,c('patient_id','RSI_LOS','RSI_INHOSP','RSI_30Day','RSI_1YR')])
  test_grid$diff_1yrpod <- test_grid$rsi_1yrpod - test_grid$RSI_1YR
  test_grid$diff_30dlos <- test_grid$rsi_30dlos - test_grid$RSI_LOS
  test_grid$diff_30dpod <- test_grid$rsi_30dpod - test_grid$RSI_30Day
  test_grid$diff_inhosp <- test_grid$rsi_inhosp - test_grid$RSI_INHOSP
  test_grid$diff_sum <- test_grid$diff_1yrpod + test_grid$diff_30dlos + test_grid$diff_30dpod + test_grid$diff_inhosp
  test_grid$diff_problem <- test_grid$diff_sum > 0.0001
  
  test_grid[test_grid$diff_problem,]
}
