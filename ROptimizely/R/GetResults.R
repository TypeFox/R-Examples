#'@title 
#'Optimizely API to extract results of individual experiment
#'
#'@description
#' \code{GetResults} extracts results for an individual experiment.
#'Metrics  for each experiment are listed for every combination of variations and goals defined for that experiment. All experiments started on or after January 21 2015 have statistics computed by Optimizely Stats Engine.
#'@export 
#'
#'
#'@param  experiment.id experiment identifier. 
#'
#'@param result.type Type of results to be extracted (valid options 'stats' or 'results')
#' Setting this argument to results  will extract top-level results.This option will be consistent with Optimizely Results page for experiments started before Jan 21, 2015. Default is set to 'stats'
#' Default is set to stats to extract results from  Optimizely Stats engine.
#'@param audience.id Audience identifier. Setting this parameter can be used to filter experiment results by particular audience. 
#'@param dimension.id Dimension identifier. Setting this parameter together with dimension.value to filter results to visitors with custom dimension value. Default to NA
#'@param dimension.value value of a dimension. Use this in combination with dimension.id
#'@examples
#'\dontrun{
#' # Extract results of a single experiment 
#' # Assign token before getting results 
#' # set_token('abcdefghihjklmnopqrs:54321')
#' # Default results 
#' exp.df<-GetResults('123456')
#' 
#' # Change result type to gather results
#' exp.df<-GetResults('123456',result.type='results)
#'
#' # Filter results by audience
#' exp.df<-GetResults('123456',audience.id='123141')
#'
#' # filter results by dimension
#' 
#' exp.df<-GetResults('123456',dimension.id='456', dimension.value='state')
#'
#'   
#'}
#'
#'@return data frame with experiment results. A data frame representing every combination of variations and goals that have been defined for the experiment. For example, if there are three variations and two goals defined for an experiment, the response will contain six rows representing each variation_id and goal_id combination
#'
#'
#'@seealso \code{\link{GetExperimentResults}}

GetResults <- function(experiment.id, result.type = "stats", audience.id = NA, dimension.id = NA, dimension.value = NA) {
    token.id <- get_token()
    base.args <- list()
    if (!is.na(audience.id)) {
        base.args <- c(base.args, audience_id = audience.id)
    } else if (!is.na(dimension.id) & is.na(dimension.value)) {
        base.args <- c(base.args, dimension_id = dimension.id, dimension_value = dimension.value)
    }
    
    base.url <- "https://www.optimizelyapis.com/experiment/v1/experiments"
    # Currently Depreciated experiment.result.url<-paste(base.url,experiment.id,'results',sep='/')
    if (result.type == "results") {
        experiment.result.url <- paste(base.url, experiment.id, "results", sep = "/")
    } else {
        experiment.result.url <- paste(base.url, experiment.id, "stats", sep = "/")
    }
    if (length(base.args) > 0) {
        experiment.result.response <- GET(experiment.result.url, add_headers(Token = token.id), query = base.args)
    } else {
        experiment.result.response <- GET(experiment.result.url, add_headers(Token = token.id))
    }
    
    if (experiment.result.response$status_code == 200) {
        experiment.results.df <- data.frame(fromJSON(content(experiment.result.response, as = "text")))
    } else {
        stop(paste("Status Code :", experiment.result.response$status_code, " Error Message :", content(experiment.result.response), 
            "\n Validate your token at https://app.optimizely.com/tokens or check permissions for the entity"))
    }
    return(experiment.results.df)
} 
