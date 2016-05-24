#'@title API token Authorizations
#'@description 
#'All API request needs to be authenticated. To authenticate, use token generated from 
#'\url{https://app.optimizely.com/tokens}
#'@param token copy token string from 'Manage Tokens'
#'@examples 
#' set_token('abcdefghijklmnopqrstuvwxyz:123456')
#'@export set_token
#'@import jsonlite
#'@import httr
set_token <- function(token) {
    assign("Token", token, envir = roptimizely_cache)
}
#'@export get_token
get_token <- function() {
    if (!exists("Token", envir = roptimizely_cache)) {
        stop("Token Not assigned for this session")
    }
    
    return(get("Token", envir = roptimizely_cache))
    
} 
