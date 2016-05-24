#' Metabolic profiles in tuberculosis.
#'
#' Relative abundances of metabolites from serum samples of three groups of
#' individuals
#'
#' A data frame with 136 observations on 425 metabolic variables.
#' 
#' Serum samples from three groups of individuals were compared: tuberculin
#' skin test negative (NEG), positive (POS) and clinical tuberculosis (TB).
#' 
#' @references
#'    Weiner J 3rd, Parida SK, Maertzdorf J, Black GF, Repsilber D, et al.
#'    (2012) Biomarkers of Inflammation, Immunosuppression and Stress 
#'    Are Revealed by Metabolomic Profiling of Tuberculosis
#'    Patients. PLoS ONE 7(7): e40221. doi:10.1371/journal.pone.0040221
#' 
#' @examples
#'  data(metabo)
#'  # maybe str(metabo) ; plot(metabo) ...
#'  pca <- prcomp( metabo[,-1] )
#' @keywords datasets
#' @name metabo
NULL
