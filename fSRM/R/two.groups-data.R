#' Data set on negative interactions

#' A simulated dataset, mimicking the study performed by Eichelsheim et al. (2011) who investigated whether there are
#' differences in patterns of negativity between families with and without an adolescent with externalizing problem
#' behavior. The problematic and nonproblematic group consist of 120 and 153 four-person families, respectively.
#' This dataset contains a measures of negativity for each of the 12 relationships.
#' Four roles are present: Mothers "M", fathers "F", the asolescent with externalizing problem behavior "T", and the adolescent sibling without problem behavior "S".
#' A wide version of the same data set is in \code{two.groups.wide}.
#'
#' The variables are as follows:
#'
#' \itemize{	
#' \item family.id An indicator for the family.
#' \item actor.id An indicator for the perceiver, either "M", "F", "T", or "S".
#' \item partner.id An indicator for the target, either "M", "F", "T", or "S".
#' \item group An indicator for the group, 1 represents the problematic families and 2 the nonproblematic families.
#' \item neg Negativity measure.
#' }
#'
#' @references Eichelsheim, V. I., Buist, K. L., Dekovic, M., Cook, W. L., Manders, W., Branje, S. J. T., et al. (2011). Negativity in problematic and nonproblematic families: A multigroup social relations model analysis with structured means. \emph{Journal of Family Psychology}, 25, 152-6. DOI: 10.1037/a0022450.

#' @format A data frame with 3276 rows and 6 variables (273 families with 4 members each, round-robin design)
#' @docType data
#' @keywords datasets
#' @name two.groups
#' @aliases two.groups.wide
#' @usage data(two.groups)
NULL
