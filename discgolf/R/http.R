disc_GET <- function(url, endpt, args = list(), ...){
  res <- GET(as.url(url, endpt), query = args, dg_head(), ...)
  check_res(res)
  parse_json(res)
}

disc_POST <- function(url, endpt, args = list(), ...){
  res <- POST(as.url(url, endpt), query = args, dg_head(), ...)
  check_res(res)
  parse_json(res)
}

disc_PUT <- function(url, endpt, args = list(), ...){
  res <- PUT(as.url(url, endpt), query = args, dg_head(), ...)
  check_res(res)
  parse_log(res)
}

disc_DELETE <- function(url, endpt, args = list(), ...){
  res <- DELETE(as.url(url, endpt), query = args, dg_head(), ...)
  check_res(res)
  parse_log(res)
}

check_res <- function(x) {
  if (x$status_code > 201) {
    err <- err_handle(x)
    stop(err$status, " - ",
         err$mssg,
         call. = FALSE)
  }
}

parse_json <- function(x) {
  stopifnot(x$headers$`content-type` == "application/json; charset=utf-8")
  res <- content(x, as = "text", encoding = "UTF-8")
  jsonlite::fromJSON(res)
}

parse_log <- function(x) {
  stopifnot(
    x$headers$`content-type` %in%
      c("application/json; charset=utf-8",
        "text/plain; charset=utf-8",
        "text/html; charset=utf-8"))
  x$status_code == 200
}

dg_head <- function() {
  add_headers(Accept = 'application/json', user_agent = "discgolf R client")
}

as.url <- function(x, y) file.path(x, y)

err_handle <- function(y) {
  z <- httr::content(y, as = "text", encoding = "UTF-8")
  if (nchar(z) == 0) {
    list(status = y$status_code,
         mssg = http_status(y)$reason)
  } else {
    if (grepl("html", y$headers$`content-type`)) {
      html <- xml2::read_html(z)
      list(status = y$status_code,
           mssg = xml2::xml_text(xml2::xml_find_one(html, "//h1"))
      )
    } else {
      bb <- jsonlite::fromJSON(z)
      list(status = y$status_code, mssg = unlist(bb$errors))
    }
  }
}
