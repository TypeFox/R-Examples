#' Find out about http status codes
#'
#' @name http
#' @param code (character) An http status code, or a regex search for HTTP
#' status codes.
#' @param text (character) A text string to search the messages or descriptions
#' of HTTP status codes.
#' @examples
#' # search by code
#' http_code(100)
#' http_code(400)
#' http_code(503)
#'
#' # fuzzy code search
#' http_code('1xx')
#' http_code('3xx')
#' http_code('30[12]')
#' http_code('30[34]')
#' http_code('30[34]')
#'
#' # search by text message
#' http_search("request")
#' http_search("forbidden")
#' http_search("too")
#'
#' @examples \dontrun{
#' http_search("birds")
#' http_code(999)
#' }

#' @export
#' @rdname http
http_code <- function(code = NULL){
  code <- as.character(code)
  if (is.null(code)) {
    print_codes()
  } else {
    if (is_three_digit_code(code)) {
      print_code(code)
    } else {
      print_filtered_codes(code)
    }
  }
}

#' @export
#' @rdname http
http_search <- function(text = NULL) print_search(text)


# helpers -------------------------------
print_filtered_codes <- function(code){
  code2 <- paste0(gsub("x", "\\\\d", code), "$")
  found_codes <- nozero(sapply(names(status_codes), function(x) grep(code2, x, value = TRUE)))
  if (length(found_codes) == 0) stopcode('No code found corresponding to', code)
  print_codes(found_codes)
}

print_codes <- function(codes=names(status_codes)) lapply(sort(codes), print_code)

print_code <- function(code){
  twocodes <- status_codes[code]
  if (length(twocodes[[1]]) != 2) stopcode('No description found for code', code)
  structure(msg_list(code, twocodes[[1]][[1]], twocodes[[1]][[2]]), class = "http_code")
}

#' @export
print.http_code <- function(x, ...){
  cat(sprintf("<Status code: %s>", x$status_code), sep = "\n")
  cat(sprintf("  Message: %s", x$message), sep = "\n")
  cat(sprintf("  Explanation: %s", x$explanation), sep = "\n")
}

print_search <- function(text){
  found_codes <- vapply(status_codes, function(x){
    any(vapply(x, function(y) grepl(text, y, ignore.case = TRUE), logical(1)))
  }, logical(1))
  if (any(found_codes)) {
    print_codes(names(status_codes[found_codes]))
  } else {
    stopcode('No status code found for search: ', text)
  }
}
