#' Synthetic Hospital Admissions
#'
#' A dataset containing synthetic hospital admissions in the classic longitudinal format.
#' The dataset counts imaginary 10 patients who undergo different (re)admission into a hospital.
#' Some demographic and clinical variables are also included.
#'
#' @format A \code{data.table} with 53 rows and 12 variables:
#' \describe{
#'   \item{subj}{Subject ID (integer)}
#'   \item{adm_number}{Hospital admissions counter (integer)}
#'   \item{gender}{Gender of patient (factor with 2 levels: "F" = females, "M" = males)}
#'   \item{age}{Age of patient in years at the given observation (integer)}
#'   \item{rehab}{Rehabilitation flag: if the admission has been in rehabilitation,
#'   then rehab = 1, else = 0 (integer)}
#'   \item{it}{Intensive Therapy flag: if the admission has been in intensive therapy,
#'   then it = 1, else = 0 (integer)}
#'   \item{rehab_it}{String which in one place marks the hospital admission types based on
#'   rehab and it. The standard admission is coded as "df" (default). If admission was in
#'   rehabilitation or in intensive therapy, rehab_it = "rehab" or "it", respectively (character)}
#'   \item{label_2}{Subject status at the end of the study. It takes 2 values: "alive" and "dead"
#'   (character)}
#'   \item{label_3}{Subject status at the end of the study. It takes 3 values: "alive" and "dead_in"
#'   and "dead_out" (character)}
#'   \item{dateIN}{Exact admission date (date)}
#'   \item{dateOUT}{Exact discharge date (date)}
#'   \item{dateCENS}{Either censoring time or exact death time (date)}
#' }
"hosp"
