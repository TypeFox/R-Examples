#'@title List Optimizely entities
#'@description
#'Query Optimizely API to get project entities listed in a project 
#'  \itemize{
#' \item Goals
#' \item Audiences
#' \item Targeting list
#' \item Dimensions
#'}
#'@param  project.id identifier for project.
#'@param  param  Identifies the entity to list
#'  \itemize{
#' \item G - Goals
#' \item A - Audiences
#' \item T - Targeting list (as on June 2015 user lists is part of paid feature. Confirm if you have access to user lists)
#' \item D - Dimensions
#'}
#'@export
#'@examples
#'\dontrun{
#' # Information for an dimension
#' # Use set_token to assign account token 
#' # set_token('123456789101112:zzzzz')
#' # Lists all audiences in a project 
#' final.df<-ListProjectParam('1234567',param='A')
#' 
#' # Lists goals in a project (default)
#' final.df<-ListProjectParam('1234567')
#'}
#'@return  data frame - of goals/audiences/targeting lists/dimensions in a project
#'@seealso \code{\link{GetProjectList}}

ListProjectParam <- function(project.id, param = "G") {
    # Obtain token information
    token.id <- get_token()
    
    # Construct request url from param definition
    base.url <- "https://www.optimizelyapis.com/experiment/v1/projects"
    
    if (toupper(param) == "G") {
        final.url <- paste(base.url, project.id, "goals/", sep = "/")
    } else if (toupper(param) == "A") {
        final.url <- paste(base.url, project.id, "audiences/", sep = "/")
    } else if (toupper(param) == "T") {
        final.url <- paste(base.url, project.id, "targeting_lists/", sep = "/")
    } else if (toupper(param) == "D") {
        final.url <- paste(base.url, project.id, "dimensions/", sep = "/")
    }
    
    
    
    
    # Get response in json format
    response <- GET(final.url, add_headers(Token = token.id))
    
    # Validate status code and extract data
    if ((response$status_code == 200)) {
        final.df <- data.frame(fromJSON(content(response, as = "text")))
    } else {
        stop(paste("Status Code :", response$status_code, " Error Message :", content(response), "\n Verify your identifier.\n Validate your token at https://app.optimizely.com/tokens"))
    }
    return(final.df)
} 
