############################################################################################################################

# Function: CSE
# Argument: ....
# Description: This function applies the metrics specified in the evaluation model to the test results (p-values) and
# summaries to the statistic results.
#' @export
CSE = function(data, analysis, evaluation, simulation) {
  UseMethod("CSE")
}