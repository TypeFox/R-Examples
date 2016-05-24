pluck <- function(x, name, type) {
  if (missing(type)) {
    lapply(x, "[[", name)
  } else {
    vapply(x, "[[", name, FUN.VALUE = type)
  }
}

is_url <- function(x){
  grepl("https?://", x, ignore.case = TRUE) || grepl("localhost:[0-9]{4}", x, ignore.case = TRUE)
}

is_port <- function(x) {
  # strip other characters
  x <- strextract(x, "[[:digit:]]+")
  if (length(x) == 0) {
    FALSE
  } else {
    grepl("[[:digit:]]", x) && nchar(x) == 4
  }
}

add_http <- function(x) {
  if (!grepl("https?://", x, ignore.case = TRUE)) {
    paste0("http://", x)
  } else {
    x
  }
}

comp <- function(l) {
  Filter(Negate(is.null), l)
}

empty <- function(l) {
  is_length_zero <- function(z) {
    length(z) == 0
  }
  tmp <- Filter(Negate(is_length_zero), l)
  if (length(tmp) == 1 && is(tmp, "list")) {
    tmp[[1]]
  } else {
    tmp
  }
}

strextract <- function(str, pattern) {
  regmatches(str, regexpr(pattern, str))
}

strtrim <- function(str) {
  gsub("^\\s+|\\s+$", "", str)
}

trimslash <- function(str) {
  gsub("\\/+$", "", str)
}

combconfig <- function(x) {
  if (is.null(x)) {
    NULL
  } else {
    req <- do.call("c", x[vapply(x, class, "") == "request"])
    c(req, x[vapply(x, class, "") != "request"])
  }
}

gather_paths <- function(x) {
  x$url <- trimslash(x$url)
  if (!is.null(x$paths) && !is.null(x$template)) {
    stop("Cannot pass use both api_template and api_path", call. = FALSE)
  }
  if (!is.null(x$paths)) {
    file.path(x$url, paste(unlist(x$paths), collapse = "/"))
  } else if (!is.null(x$template)) {
    file.path(x$url, x$template)
  } else {
    x$url
  }
}

make_ua <- function() {
  versions <- c(curl = curl::curl_version()$version,
                curl = as.character(packageVersion("curl")),
                httr = as.character(packageVersion("httr")),
                request = as.character(packageVersion("request")))
  paste0(names(versions), "/", versions, collapse = " ")
}

get_links <- function(w) {
  lk <- w$link
  if (is.null(lk)) {
    NULL
  } else {
    if (is(lk, "character")) {
      links <- strtrim(strsplit(lk, ",")[[1]])
      lapply(links, each_link)
    } else {
      nms <- sapply(w, "[[", "name")
      tmp <- unlist(w[nms %in% "next"])
      grep("http", tmp, value = TRUE)
    }
  }
}

each_link <- function(z) {
  tmp <- strtrim(strsplit(z, ";")[[1]])
  nm <- gsub("\"|(rel)|=", "", tmp[2])
  url <- gsub("^<|>$", "", tmp[1])
  list(name = nm, url = url)
}
