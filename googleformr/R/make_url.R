#' Make Either Get or Post Url for Google Forms
#'
#' Make Either Get or Post Url for Google Forms
#'
#' @param form Can be either the form_url or form_id
#' @param do Can be either "get" or "post"
#' @export
#' @include get_form_id.R
#' @examples
#' \dontrun{
#'   url %>% make_url("get") -> get_url
#'   url %>% make_url("post") -> post_url
#' }
make_url <- function (form, do) {
    if (grepl("^https://docs.google.com", form)) {
        form <- get_form_id(form)
    }
    switch(tolower(do),
           "get"  = httr::modify_url("https://docs.google.com/",
                                    path = paste0("forms/d/", form, "/viewform") ),
           "post" = httr::modify_url("https://docs.google.com/",
                                     path = paste0("forms/d/", form, "/formResponse") ) )
}

