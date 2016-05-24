#' Setup connection with LiveChat's REST API using email and an api key
#'
#' @param email Email address used to login from \href{https://my.livechatinc.com/agents/api-key}{Live Chat API Docs}
#' @param api_key Your API Key from \href{https://my.livechatinc.com/agents/api-key}{Live Chat API Docs}
#'
#' @return A list of 2, the first item is called email and the second item is called api_key
#' @export
#'
#' @importFrom magrittr '%>%'
#' @importFrom dplyr tbl_df select
#' @importFrom jsonlite fromJSON
#' @importFrom data.table rbindlist setkey
#' @importFrom purrr map flatten
#'
#' @examples
#' # account <- livechatCreateAccount("email_here", "api_key_here")
#' # account$email
#' # account$api_key
livechatCreateAccount <- function(
  email,
  api_key
  ) {
  if(!is.null(email) & !is.null(api_key)) {
    obj = list(email = email, api_key = api_key)
    obj
  } else {
    warning("Email and API Key both must be provided in order to access your organization's LiveChat data.")
  }
}


