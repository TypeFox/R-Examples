######################################################################################################################

# Function: GenerateReport.
# Argument: Results returned by the CSE function and presentation model and Word-document title and Word-template.
# Description: This function is used to create a summary table with all results
#' @export
GenerateReport = function(presentation.model = NULL, cse.results, report.filename, report.template = NULL){
  UseMethod("GenerateReport")
}