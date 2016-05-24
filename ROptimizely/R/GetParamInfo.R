
#'@title Information on  Optimizely entities
#'@description
#'Query Optimizely to obtain information on project entities
#' \describe{
#'  \item{ Goal}{A user list is a set of user identifiers that you have uploaded to Optimizely. Membership in a targeting list can be used to define an Optimizely audience, so you can target experiments only to a particular set of users.}
#' \item{ Audience }{ An Audience is a group of visitors that match set conditions. You can target an experiment to one or more audiences, or you can segment experiment results to see how different audiences performed.}
#' \item{ Targeting(User) list}{ A user list is a set of user identifiers that you have uploaded to Optimizely. Membership in a targeting list can be used to define an Optimizely audience}
#' \item{Dimension}{ Dimensions are attributes of visitors to your website or mobile app, such as demographic data, behavioral characteristics, or any other information particular to a visitor. Dimensions can be used to construct audiences and segment experiment results.}
#' \item{Variation}{Every experiment contains a set of variations that each change the visitor's experience in a different way. Variations define the code that should be applied on a page to change the experience, and the percentage of visitors who should see that code. A standard 'A/B' test has two variations (including the original), and Optimizely supports adding many more variations.}
#' }
#' \url{https://help.optimizely.com}
#'@param  id Identifier for project entity i.e. goal id, audience id, user list id, dimension id
#'@param  param   alphabet character  to identify the parameter list 
#' \itemize{
#' \item G - Goals
#' \item T - Targeting (user) lists
#' \item V - Variations
#' \item A - Audience 
#' \item D - Dimensions
#'}
#'@return  data frame - of goals/audiences/targeting lists/dimensions in a project
#'@examples 
#'\dontrun{
#' # Information for an dimension
#' # Use set_token to assign account token 
#' set_token('123456789101112:zzzzz')
#' # Lists all audiences in a project 
#' final.df<-GetParamInfo('1234567',param='A')
#' # Lists goals in a project (default)
#' final.df<-GetParamInfo('1234567')
#'}
#'
#'@seealso \code{\link{GetProjectList}}
#'@export GetParamInfo
#'@export empty_list_check

GetParamInfo <- function(id, param) {
    # Obtain token information
    token.id <- get_token()
    
    # Construct request url from param definition
    base.url <- "https://www.optimizelyapis.com/experiment/v1"
    
    if (toupper(param) == "G") {
        final.url <- paste(base.url, "goals", id, sep = "/")
    } else if (toupper(param) == "A") {
        final.url <- paste(base.url, "audiences", id, sep = "/")
    } else if (toupper(param) == "T") {
        final.url <- paste(base.url, "targeting_lists", id, sep = "/")
    } else if (toupper(param) == "D") {
        final.url <- paste(base.url, "dimensions", id, sep = "/")
    } else if (toupper(param) == "V") {
        final.url <- paste(base.url, "variations", id, sep = "/")
    } else {
        stop("Enter a valid option for entity param")
    }
    # Get response in json format
    response <- GET(final.url, add_headers(Token = token.id))
    
    # Validate status code and extract data
    if ((response$status_code == 200)) {
        if (toupper(param) == "G") {
            json_response = fromJSON(content(response, as = "text"))
            final.df <- data.frame(lapply(json_response, empty_list_check))
        } else {
            final.df <- data.frame(fromJSON(content(response, as = "text")))
        }
    } else {
        stop(paste("Status Code :", response$status_code, " Error Message :", content(response), "\n Verify your identifier.\n Validate your token at https://app.optimizely.com/tokens"))
    }
    return(final.df)
}

 
