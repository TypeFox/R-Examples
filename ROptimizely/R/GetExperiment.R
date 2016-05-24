#'@export GetExperiment
GetExperiment <- function(project.id) {
    # Obtain token information
    token.id <- get_token()
    
    # Construct request url
    base.url <- paste("https://www.optimizelyapis.com/experiment/v1/projects", project.id, "experiments/", 
        sep = "/")
    
    # Get response in json format
    experiment.response <- GET(base.url, add_headers(Token = token.id))
    
    # Validate status code and extract data
    if ((experiment.response$status_code == 200)) {
        experiment.df <- data.frame(fromJSON(content(experiment.response, as = "text")))
    } else {
        stop(paste("Status Code :", experiment.response$status_code, " Error Message :", content(experiment.response), 
            "\n Validate your token at https://app.optimizely.com/tokens or check permission for the id you are trying to view"))
    }
    
    return(experiment.df)
} 
