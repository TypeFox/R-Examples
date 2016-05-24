#'@title Experiment information
#'@description
#'Query Optimizely  to read metadata information for a single experiment
#'  
#'
#'@export 
#'
#'
#'@param  experiment.id valid experiment identifier
#'
#'
#'@examples
#'\dontrun{
#' experiment.info <- GetExperimentInfo('1234545')
#' #use this to extract information on a list of experiments
#'do.call(rbind, lapply(c(1231231,14335333'), GetExperimentInfo))
#'}
#'
#'@return data frame of metadata information for an experiment
#'
#'
#'@seealso \code{\link{GetExperimentList}}

GetExperimentInfo <- function(experiment.id) {
    # Obtain token information
    token.id <- get_token()
    
    # Construct request url
    base.url <- paste("https://www.optimizelyapis.com/experiment/v1/experiments", experiment.id, sep = "/")
    
    # Get response in json forma
    experiment.response <- GET(base.url, add_headers(Token = token.id))
    
    # Validate status code and extract data
    if ((experiment.response$status_code == 200)) {
        experiment.list <- fromJSON(content(experiment.response, as = "text"))
        experiment.list$display_goal_order_lst <- paste(experiment.list$display_goal_order_lst, collapse = ",")
        experiment.df <- data.frame(t(experiment.list))
    } else {
        stop(paste("Status Code :", experiment.response$status_code, " Error Message :", content(experiment.response), 
            "\n Validate your token at https://app.optimizely.com/tokens"))
    }
    
    return(experiment.df)
} 
