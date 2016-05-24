#'@title Get project information
#'@description
#'Query Optimizely API to read a project. A project is a collection of experiments, goals and audiences. A javascript file associated with each project is included. 
#'
#'@export 
#'
#'
#'@param  project.id valid project identifier
#'
#'
#'@examples
#'\dontrun{
#' project.info <- GetProjectInfo('1234545')
#'}
#'
#'@return data frame of metadata for a single project
#'
#'
#'@seealso \code{\link{GetProjectList}}
#'
#' 
GetProjectInfo <- function(project.id) {
    # Obtain token information
    token.id <- get("token", envir = roptimizely_cache)
    # Construct request url
    base.url <- paste("https://www.optimizelyapis.com/experiment/v1/projects", project.id, sep = "/")
    
    # Get response in json forma
    project.response <- GET(base.url, add_headers(Token = token.id))
    
    
    # Validate status code and extract data
    if ((project.response$status_code == 200)) {
        project.list <- lapply(fromJSON(content(project.response, as = "text")), function(x) ifelse(is.null(x), 
            NA, x))
        project.df <- data.frame(project.list)
    } else {
        stop(paste("Status Code :", project.response$status_code, " Error Message :", content(project.response), 
            "\n Validate your token at https://app.optimizely.com/tokens"))
    }
    
    return(project.df)
} 
