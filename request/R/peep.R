#' Peek at a query
#'
#' @export
#' @param .data (list) input, using higher level interface
#' @examples \dontrun{
#' api('https://api.github.com/') %>% peep
#' api('https://api.github.com/') %>%
#'   api_path(repos, ropensci, rgbif, issues) %>%
#'   peep
#'
#' repo_info <- list(username = 'craigcitro', repo = 'r-travis')
#' api('https://api.github.com/') %>%
#'   api_template(template = 'repos/{{username}}/{{repo}}/issues', data = repo_info) %>%
#'   peep
#'
#' api("http://api.plos.org/search") %>%
#'   api_query(q = ecology, wt = json, fl = id, fl = journal) %>%
#'   peep
#' }
peep <- function(.data) {
  pipe_autoexec(toggle = FALSE)
  structure(.data, class = "req")
}
