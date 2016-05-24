#' Create a list of goals with LiveChat credentials as parameters
#'
#' @param account defined by the livechatCreateAccount function
#'
#' @return a data frame of goals
#' @export
#'
#' @examples
#' # account <- livechatCreateAccount("email_here", "api_key_here")
#' # livechatGetGoals(account)
livechatGetGoals <- function(
  account
) {

  email <- account$email
  api_key <- account$api_key
  email_api_key <- paste0(email, ":", api_key)

  query <- paste0("curl 'https://api.livechatinc.com/goals'",
                  " -u ",
                  email_api_key,
                  " -H X-API-Version:2")
  raw_json <- system(query, intern = TRUE)
  raw_json_char <- paste(raw_json, collapse = "")

  if(raw_json_char == "[]") { # no goals set
    message("There are currently no goals set in your LiveChat account.")
  } else {
    df <- fromJSON(raw_json_char) %>% tbl_df()
    df
  }

}

