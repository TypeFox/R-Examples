#' Create a list of greetings with LiveChat credentials as parameters
#'
#' @param account defined by the livechatCreateAccount function
#' @param group defaults to 0, see \href{https://developers.livechatinc.com/rest-api/#!get-greetings}{LiveChat REST API}
#'
#' @return a data frame of greetings
#' @export
#'
#' @examples
#' # account <- livechatCreateAccount("email_here", "api_key_here")
#' # livechatGetGreetings(account, group = 0)
livechatGetGreetings <- function(
  account,
  group = 0
) {

  email <- account$email
  api_key <- account$api_key
  email_api_key <- paste0(email, ":", api_key)

  query <- paste0("curl 'https://api.livechatinc.com/greetings?group=",
                  group,"'",
                  " -u ",
                  email_api_key,
                  " -H X-API-Version:2")
  raw_json <- system(query, intern = TRUE)
  raw_json_char <- paste(raw_json, collapse = "")
  df <- fromJSON(raw_json_char)
  df

}

