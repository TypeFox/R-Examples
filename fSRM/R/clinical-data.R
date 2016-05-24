#' Data set on negativity in problematic and nonproblematic families

#' A simulated dataset, mimicking the study performed by Eichelsheim et al. (2011) who investigated whether there are
#' differences in patterns of negativity between families with and without an adolescent with externalizing problem
#' behavior. In the study of Eichelsheim et al. (2011), four members of the same family (mother, father, target 
#' adolescent and sibling) reported on the amount of negativity they experienced in relation to each other. 
#' In total, these authors studied 120 Dutch four-person families with a target adolescent scoring above the 
#' externalizing behavior clinical norm scores on either the Child Behavior Check List (N = 47; CBCL; Achenbach, 1991) 
#' or the Youth Self Report (N = 73; YSR; Achenbach, 1991). Because of confidentiality reasons, not the original data 
#' but mimicked data are used here, so results of social relations analyses will deviate from the original paper. 
#' 
#' In sum, this dataset contains a measures of negativity for each of the 12 relationships.
#' Four roles are present: Mothers "M", fathers "F", the asolescent with externalizing problem behavior "T", and the adolescent sibling without problem behavior "S".
#' A wide version of the same data set is in \code{clinical.wide}.
#'
#' The variables are as follows:
#'
#' \itemize{	
#' \item family.id An indicator for the family.
#' \item actor.id An indicator for the perceiver, either "M", "F", "T", or "S".
#' \item partner.id An indicator for the target, either "M", "F", "T", or "S".
#' \item neg Negativity measure.
#' }
#'

#' @docType data
#' @keywords datasets
#' @name clinical
#' @aliases clinical.wide
#' @usage data(clinical)
NULL
