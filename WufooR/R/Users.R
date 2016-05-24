#' @title Return information about the user
#' 
#' @description Shows information about the user. 
#' 
#' @references wufoo_api - User's API Key. Default: \code{\link{auth_key}}.
#' @param wufoo_name - User's Name. Default: \code{\link{auth_name}}
#' @param domain - for thatever reason domain's ccTLD may change. With this param you can change 
#' default \code{wufoo.com} to another one such as \code{wufoo.eu}
#' @param showRequestURL - use only for DEBUG purposes
#' @param debugConnection - samilar as above. Either "0L" (false; default) or "1L"
#' 
#' @return \url{http://help.wufoo.com/articles/en_US/SurveyMonkeyArticleType/The-Users-API}
#' 
#' @examples 
#' user_info(debugConnection = 1L, domain = "wufoo.eu")
#' 
#' @export
user_info <- function(wufoo_name = auth_name(NULL), domain = "wufoo.com", showRequestURL = FALSE, debugConnection = 0L) {
  
  user_url <- paste0("https://", wufoo_name, ".", domain, "/api/v3/users.json")
  
  executedUserGetRst <- doRequest(user_url, showURL = showRequestURL, debugConnection = debugConnection)
  df_user <- executedUserGetRst$Users
  
  return(df_user)
}




