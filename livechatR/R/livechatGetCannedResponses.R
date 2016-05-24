#' Create a list of canned responses with LiveChat credentials as parameters
#'
#' @param account defined by the livechatCreateAccount function
#' @param group defaults to 0, see \href{https://developers.livechatinc.com/rest-api/#!canned-responses}{LiveChat REST API}
#'
#' @return a data frame of canned responses
#' @export
#'
#' @examples
#' # account <- livechatCreateAccount("email_here", "api_key_here")
#' # livechatGetCannedResponses(account, group = 0)
livechatGetCannedResponses <- function(
  account,
  group = 0
) {

  email <- account$email
  api_key <- account$api_key
  email_api_key <- paste0(email, ":", api_key)

  query <- paste0("curl 'https://api.livechatinc.com/canned_responses?group=",
                  group,"'",
                  " -u ",
                  email_api_key,
                  " -H X-API-Version:2")
  cr_raw_json <- system(query, intern = TRUE)
  raw_json_char <- paste(cr_raw_json, collapse = "")
  cr <- fromJSON(raw_json_char) %>% tbl_df()
  cr

}

