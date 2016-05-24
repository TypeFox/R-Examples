ec <- function(l) Filter(Negate(is.null), l)

error_handler <- function(x){
  res_info <- content(x)$info
  if (x$status_code %in% c(400, 500)) {
    stop(sprintf("%s : %s", res_info$message, gsub('\"', "'", res_info$additional)), call. = FALSE)
  }
  stopifnot(x$headers$`content-type` == 'application/json; charset=utf-8')
  dat <- content(x, as = "text", encoding = 'utf-8')
  jsonlite::fromJSON(dat, FALSE)
}

enigma_GET <- function(url, args, ...){
  if (length(args) == 0) args <- NULL
  res <- GET(url, query = args, ...)
  error_handler(res)
}

check_dataset <- function(dataset){
  if (is.null(dataset)) stop("You must provide a dataset") else dataset
}

check_key <- function(x){
  tmp <- if (is.null(x)) Sys.getenv("ENIGMA_KEY", "") else x
  if (tmp == "") getOption("enigmaKey", stop("need an API key for the Enigma API")) else tmp
}

en_base <- function() 'https://api.enigma.io/v2'
