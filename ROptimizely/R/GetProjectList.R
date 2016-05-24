#'@title
#'Query Optimizely API to list projects in an account 
#'@description
#'\code{GetProjectList} Returns a data frame of projects listed in account authorized by the token
#' Obtain  API token from \url{https://app.optimizely.com/tokens}
#'@export 
#'
#'
#'
#'@examples
#'\dontrun{
#'  # use set_token(Token) to assign your API token
#' project.df <- GetProjectList()
#'}
#'
#'@return  data frame containing the response from Optimizely API
#'
#'@seealso \code{\link{GetProjectInfo}}
#'
GetProjectList <- function() {
    # Obtain token information
    token.id <- get_token()
    
    # Get response in json forma
    project.response <- GET("https://www.optimizelyapis.com/experiment/v1/projects", add_headers(Token = token.id))
    
    # Validate status code and extract data
    if ((project.response$status_code == 200)) {
        project.df <- data.frame(fromJSON(content(project.response, as = "text")))
    } else {
        stop(paste("Status Code :", project.response$status_code, " Error Message :", content(project.response), 
            "\n Validate your token at https://app.optimizely.com/tokens"))
    }
    
    return(project.df)
} 
