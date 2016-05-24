#' Search
#'
#' @export
#' @param query (character) Query terms. Required.
#' @param order (character) One of views, latest, likes
#' @param status (character) One of open, closed, archived, noreplies, or single_user
#' @param category (character) Category to search for
#' @param username (character) User name
#' @param group (character) Groupo name
#' @param badge (character) Badge name
#' @param in_ (character) One of likes, posted, watching, tracking, private, bookmarks, first
#' @param posts_count (integer) Number of posts per topic
#' @param min_age (integer) Minimum age
#' @param max_age (integer) Maximum age
#' @template args
#' @examples \dontrun{
#' dg_search(query = "poo")
#' dg_search(posts_count = 1)
#' dg_search(in_ = "posted")
#' dg_search(status = "open")
#' }
dg_search <- function(query = NULL, order = NULL, status = NULL, category = NULL,
  username = NULL, group = NULL, badge = NULL, in_ = NULL, posts_count = NULL,
  min_age = NULL, max_age = NULL, url = NULL, key = NULL, user = NULL, ...) {

  args <- dc(list(api_key = check_key(key), api_username = check_user(user),
                  term = query))
  other_args <- dc(list(order = order, status = status, category = category,
                        user = username, group = group, badge = badge, `in` = in_,
                        posts_count = posts_count, min_age = min_age, max_age))
  oa <- paste0(paste(names(other_args), unlist(unname(other_args)), sep = ":"),
               collapse = "+")
  args$term <- paste0(args$term, "+", oa)
  disc_GET(check_url(url), "search/query", args, ...)
}
