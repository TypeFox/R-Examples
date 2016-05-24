#' Create a list of agents with LiveChat credentials as parameters
#'
#' @param account defined by the livechatCreateAccount function
#'
#' @return A date frame of agents
#' @export
#'
#' @examples
#' # account <- livechatCreateAccount("email_here", "api_key_here")
#' # livechatGetAgents(account)
livechatGetAgents <- function(
  account
) {

  email <- account$email
  api_key <- account$api_key
  email_api_key <- paste0(email, ":", api_key)

  query <- paste0("curl 'https://api.livechatinc.com/agents'",
                  " -u ",
                  email_api_key,
                  " -H X-API-Version:2")
  agents_raw_json <- system(query, intern = TRUE)
  raw_json_char <- paste(agents_raw_json, collapse = "")
  agents <- fromJSON(raw_json_char) %>% tbl_df()
  agents

}

