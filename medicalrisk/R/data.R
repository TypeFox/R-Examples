#' List of Charlson comorbidities
#' 
#' @docType data
#' @name charlson_list
#' @usage charlson_list
#' @format A list, with one column for each comorbidity; value is a textual description 
#' @keywords datasets
#' @examples
#' # List the strings used to refer to Charlson comorbidities
#' names(charlson_list)
#' 
#' # List descriptions of comorbidities
#' charlson_list
NULL

#' List of Elixhauser comorbidities
#' 
#' @docType data
#' @name elixhauser_list
#' @usage elixhauser_list
#' @format A list, with one column for each comorbidity; value is a textual description 
#' @keywords datasets
#' @examples
#' # List the strings used to refer to Elixhauser comorbidities
#' names(elixhauser_list)
#' 
#' # List descriptions of comorbidities
#' elixhauser_list
NULL

#' List of ICD-9-CM diagnostic and procedural codes
#' 
#' ICD-9-CM codes have the periods removed.  Diagnostic codes are prefixed with 
#' 'D' while procedure codes are prefixed with 'P'. So, diagnostic code
#' \code{404.03} appears as \code{"D40403"}.
#' 
#' Obsolete codes not active in 2012 are not present, and may cause this dataset
#' to miss certain classifications when applied to older datasets.  For example, 
#' codes 043 and 044 (both obsolete AIDS codes) are not included.
#' 
#' @docType data
#' @name icd9cm_list
#' @usage icd9cm_list
#' @format A string vector
#' @references 1. \url{https://www.cms.gov/Medicare/Coding/ICD9ProviderDiagnosticCodes/codes.html}
#'   
#' @keywords datasets
#' @examples
#' # Count procedural codes
#' length(icd9cm_list[grep('^P',icd9cm_list)])
NULL

#' Map of Charlson comorbidity categories to weights
#' 
#' List that links the Charlson comorbidity categories to the original 
#' weights (specified in the original Charlson paper, Table 3)
#' 
#' Original Weights:
#' 
#' 1 = MI, CHF, PVD, CVD, Dementia, Chronic pulm dz, Connective tissue dz,
#' Ulcer, Mild liver dz, Diabetes
#' 
#' 2 = Hemiplegia, Mod or severe renal dz, Diabetes with end organ damage, Any
#' tumor, Leukemia, Lymphoma
#' 
#' 3 = Moderate or severe liver dz
#' 
#' 6 = Metastatic solid tumor, AIDS
#' 
#' @docType data
#' @name charlson_weights_orig
#' @usage charlson_weights_orig
#' @format A list, with Charlson comorbidities as names and weight as value
#' @references 1. Charlson ME, Pompei P, Ales KL, MacKenzie CR: A new method of
#'   classifying prognostic comorbidity in longitudinal studies: development and
#'   validation. Journal of chronic diseases 1987; 40:373-83 
#'   \url{http://www.ncbi.nlm.nih.gov/pubmed/3558716} 
#' @seealso \code{\link{charlson_weights}},
#' \code{\link{icd9cm_charlson_deyo}}, 
#' \code{\link{icd9cm_charlson_romano}},
#' \code{\link{icd9cm_charlson_quan}},
#' \code{\link{melt_icd9list}}
#' @keywords datasets
#' @examples
#' charlson_weights_orig["aids"]
NULL

#' Map of Charlson comorbidity categories to revised weights
#' 
#' List that links the Charlson comorbidity categories to revised 
#' weights as calculated by Schneeweiss in Table 4 of his paper.
#' 
#' Revised Schneeweiss weights:
#' 
#' 0 = Connective tissue dz, Ulcer
#' 
#' 1 = MI, PVD, CVD, Diabetes, Hemiplegia
#' 
#' 2 = CHF, Chronic pulm dz, Mild liver dz, Diabetes with end organ damage, Any 
#' tumor, Leukemia, Lymphoma
#' 
#' 3 = Dementia, Mod or severe renal dz
#' 
#' 4 = Moderate or severe liver dz, AIDS
#' 
#' 6 = Metastatic solid tumor
#' 
#' @docType data
#' @name charlson_weights
#' @usage charlson_weights
#' @format A list, with Charlson comorbidities as names and weight as value
#' @references 1. Schneeweiss S, Wang PS, Avorn J, Glynn RJ: Improved comorbidity adjustment 
#'   for predicting mortality
#'   in Medicare populations. Health services research 2003; 38:1103
#'   \url{http://www.ncbi.nlm.nih.gov/pubmed/12968819}
#' @seealso \code{\link{charlson_weights_orig}},
#' \code{\link{icd9cm_charlson_deyo}}, 
#' \code{\link{icd9cm_charlson_romano}},
#' \code{\link{icd9cm_charlson_quan}},
#' \code{\link{melt_icd9list}}
#' @keywords datasets
#' @examples
#' charlson_weights["dementia"]
NULL

#' First 100 patients and their ICD-9-CM codes from the Vermont Uniform Hospital
#' Discharge Data Set for 2011, Inpatient.
#' 
#' Diagnostic ICD-9 codes are prefixed with 'D', while procedural ICD-9 codes are prefixed with 'P'.
#' 
#' @docType data
#' @name vt_inp_sample
#' @usage vt_inp_sample
#' @format A data frame, with column "id" (numeric), some descriptive columns,
#'   "dx" (factor), and "icd9cm" (factor)
#' @references 
#' \url{http://healthvermont.gov/research/hospital-utilization/RECENT_PU_FILES.aspx}
#' @seealso \code{\link{icd9cm_charlson_deyo}}, 
#' \code{\link{icd9cm_charlson_quan}}, 
#' \code{\link{icd9cm_charlson_romano}}, 
#' \code{\link{icd9cm_elixhauser_quan}}, 
#' \code{\link{icd9cm_elixhauser_ahrq37}}
#' @keywords datasets
#' @examples
#' max(vt_inp_sample$scu_days)
NULL

#' Values for calculating RSI for in-hospital mortality.
#' @docType data
#' @name rsi_beta_inhosp
#' @usage rsi_beta_inhosp
#' @format A hash (see package "hash"), where key is icd9cm code, and value is beta.
#' Special key "popbeta" has the population beta for the entire table.
#' @references 1. Sessler DI, Sigl JC, Manberg PJ, Kelley SD, Schubert A, Chamoun NG.
#'     Broadly applicable risk stratification system for predicting duration of hospitalization and mortality.
#'     Anesthesiology. 2010 Nov;113(5):1026-37. doi: 10.1097/ALN.0b013e3181f79a8d.
#' \url{http://my.clevelandclinic.org/anesthesiology/outcomes-research/risk-stratification-index.aspx}
NULL

#' Values for calculating RSI for 30-day mortality
#' @docType data
#' @name rsi_beta_30dpod
#' @usage rsi_beta_30dpod
#' @format A hash (see package "hash"), where key is icd9cm code, and value is beta.
#' Special key "popbeta" has the population beta for the entire table.
#' @references 1. Sessler DI, Sigl JC, Manberg PJ, Kelley SD, Schubert A, Chamoun NG.
#'     Broadly applicable risk stratification system for predicting duration of hospitalization and mortality.
#'     Anesthesiology. 2010 Nov;113(5):1026-37. doi: 10.1097/ALN.0b013e3181f79a8d.
#' \url{http://my.clevelandclinic.org/anesthesiology/outcomes-research/risk-stratification-index.aspx}
NULL

#' Values for calculating RSI for 1 year mortality
#' @docType data
#' @name rsi_beta_1yrpod
#' @usage rsi_beta_1yrpod
#' @format A hash (see package "hash"), where key is icd9cm code, and value is beta.
#' Special key "popbeta" has the population beta for the entire table.
#' @references 1. Sessler DI, Sigl JC, Manberg PJ, Kelley SD, Schubert A, Chamoun NG.
#'     Broadly applicable risk stratification system for predicting duration of hospitalization and mortality.
#'     Anesthesiology. 2010 Nov;113(5):1026-37. doi: 10.1097/ALN.0b013e3181f79a8d.
#' \url{http://my.clevelandclinic.org/anesthesiology/outcomes-research/risk-stratification-index.aspx}
NULL

#' Values for calculating RSI for 30-day length of stay
#' @docType data
#' @name rsi_beta_30dlos
#' @usage rsi_beta_30dlos
#' @format A hash (see package "hash"), where key is icd9cm code, and value is beta.
#' Special key "popbeta" has the population beta for the entire table.
#' @references 1. Sessler DI, Sigl JC, Manberg PJ, Kelley SD, Schubert A, Chamoun NG.
#'     Broadly applicable risk stratification system for predicting duration of hospitalization and mortality.
#'     Anesthesiology. 2010 Nov;113(5):1026-37. doi: 10.1097/ALN.0b013e3181f79a8d.
#' \url{http://my.clevelandclinic.org/anesthesiology/outcomes-research/risk-stratification-index.aspx}
NULL

#' Sample data for validating RSI
#' @docType data
#' @name rsi_sample_data
#' @usage rsi_sample_data
#' @seealso \code{\link{rsi_sample_results}}, \code{\link{verify_sessler_rsi}}
#' @format A data table with a patient ID and several columns with ICD-9-CM codes.
#' @references 1. Sessler DI, Sigl JC, Manberg PJ, Kelley SD, Schubert A, Chamoun NG.
#'     Broadly applicable risk stratification system for predicting duration of hospitalization and mortality.
#'     Anesthesiology. 2010 Nov;113(5):1026-37. doi: 10.1097/ALN.0b013e3181f79a8d.
#' \url{http://my.clevelandclinic.org/anesthesiology/outcomes-research/risk-stratification-index.aspx}
NULL

#' Sample results for validating RSI
#' @docType data
#' @name rsi_sample_results
#' @usage rsi_sample_results
#' @seealso \code{\link{rsi_sample_data}}, \code{\link{verify_sessler_rsi}}
#' @format A data table with a patient ID, principal diagnosis, principal procedure, and RSI result columns.
#' @references 1. Sessler DI, Sigl JC, Manberg PJ, Kelley SD, Schubert A, Chamoun NG.
#'     Broadly applicable risk stratification system for predicting duration of hospitalization and mortality.
#'     Anesthesiology. 2010 Nov;113(5):1026-37. doi: 10.1097/ALN.0b013e3181f79a8d.
#' \url{http://my.clevelandclinic.org/anesthesiology/outcomes-research/risk-stratification-index.aspx}
NULL
