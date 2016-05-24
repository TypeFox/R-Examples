#' ICU data
#' 
#' icu dataset.
#' 
#' @format A data.frame with 200 rows and 21 variables:
#' \describe{
#' \item{id}{Identification code (ID Number)}
#' \item{sta}{Vital Status at hospital discharge (1: Died, 2: Lived)}
#' \item{age}{Age (Years)}
#' \item{gender}{Gender (1: Male, 2: Female)}
#' \item{race}{Race (1: White, 2: Black, 3: Other)}
#' \item{ser}{Service at ICU admission (1: Medical, 2: Surgical)}
#' \item{can}{Cancer part of present problem (1: No, 2: Yes)}
#' \item{crn}{History of chronic renal failure (1: No, 2: Yes)}
#' \item{inf}{Infection probable at ICU admission (1: No, 2: Yes)}
#' \item{cpr}{CPR prior to ICU admission (1: No, 2: Yes)}
#' \item{sys}{Systolic blood pressure at ICU admission (mm Hg)}
#' \item{hra}{Heart rate at ICU admission (Beats/min)}
#' \item{pre}{Previous admission to an ICU within 6 months (1: No, 2: Yes)}
#' \item{type}{Type of admission (1: Elective, 2: Emergency)}
#' \item{fra}{Long bone, multiple, neck, single area, or hip fracture (1:
#' No, 2: Yes)}
#' \item{po2}{PO2 from initial blood gases (1: > 60, 2: <= 60)}
#' \item{ph}{PH from initial blood gases (1: >= 7.25, 2: < 7.25)}
#' \item{pco}{PCO2 from initial blood gases (1: <= 45, 2: > 45)}
#' \item{bic}{Bicarbonate from initial blood gases (1: >= 18, 2: < 18)}
#' \item{cre}{Creatinine from initial blood gases (1: <= 2.0, 2: > 2.0)}
#' \item{loc}{Level of consciousness at ICU admission (1: No coma or deep
#' stupor, 2: Deep stupor, 3: Coma)}
#' }
#' @source Hosmer, D.W., Lemeshow, S. and Sturdivant, R.X. (2013) Applied
#' Logistic Regression, 3rd ed., New York: Wiley
#' @examples
#' head(icu, n = 10)
#' summary(icu)
"icu"
