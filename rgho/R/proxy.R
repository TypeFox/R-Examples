#' Find Local Proxy Settings
#'
#' Tries to find the local proxy settings to access a given url.
#'
#' @param url The url to access, as a character string.
#'
#' @return A list of proxy, or \code{NULL} if not proxy are found.
#'
get_proxy_list <- function(url) {
  proxy_string <- curl::ie_get_proxy_for_url(target_url = url)

  if (is.null(proxy_string)) {
    return(NULL)
  } else {
    return(parse_proxy_string(proxy_string))
  }
}

#' Parse Proxy IP From a Character String
#'
#' Given character string containig a list of proxy,
#' returns the list in a more computer readable form.
#'
#' @param proxy_string A list of proxy IP as a concatenated string.
#'
#' @return A list of proxy settings.
#'
parse_proxy_string <- function(proxy_string) {
  proxy_string %>%
    strsplit(";") %>%
    extract2(1) %>%
    strsplit(":") %>%
    lapply(
      function(x) {
        httr::use_proxy(x[1], as.integer(x[2]))
      }
    )
}
