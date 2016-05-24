#'@title
#'Extract results of an experiment
#'@description
#'Query Optimizely API to extract top-level results experiment. Results returned are computed by Optimizely Stats Engine.
#'Metrics  for each experiment are listed for every combination of variations and goals defined for that experiment.
#' 
#'@export GetExperimentResults GetResults
#'
#'@param  experiment.list list of experiment identifier. 
#'
#'
#'@examples
#'\dontrun{
#'# Extract results 
#'# Assign token before getting results 
#' set_token('abcdefghihjklmnopqrs:54321')
#' results.df<-GetExperimentResults(c('123123','1234567'))
#'}
#'
#'@return data frame with experiments results. A data frame representing every combination of variations and goals that have been defined for each experiment in list. For example, if there are three variations and two goals defined for an experiment, the response will contain six rows representing each variation_id and goal_id combination.
#'
#'
#'@seealso \code{\link{GetResults}}

GetExperimentResults <- function(experiment.list) {
    
    # Validate status code and extract data
    if (length(experiment.list) != 0) {
        experiment.results.df <- do.call(rbind, lapply(experiment.list, GetResults))
    }
    return(experiment.results.df)
} 
