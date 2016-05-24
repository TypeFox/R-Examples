# as request functions and print method for request S3 class
as.req <- function(x) {
  UseMethod("as.req")
}

as.req.req <- function(x) {
  x
}

as.req.endpoint <- function(x){
  req(x$url)
}

as.req.character <- function(x){
  if (is_url(tryCatch(as.url(x), error = function(e) e))) {
    req(x)
  } else {
    stop("error ...")
  }
}

req <- function(x){
  structure(list(url = as.url(x)), class = "req")
}

#' @export
print.req <- function(x, ...){
  cat("<http request> ", sep = "\n")
  cat(paste0("  url: ", x$url), sep = "\n")
  cat(paste0("  paths: ",
             paste(unlist(x$paths), collapse = "/")), sep = "\n")
  cat(paste0("  query: ",
             paste(names(x$query), unname(unlist(x$query)), sep = "=", collapse = ", ")), sep = "\n")
  cat(paste0("  body: ",
             print_body(x$body)), sep = "\n")
  cat(paste0("  paging: ",
             paste(names(x$paging), unname(unlist(x$paging)), sep = "=", collapse = ", ")), sep = "\n")
  cat(paste0("  headers: ",
             print_heads(x$headers)), sep = "\n")
  cat(paste0("  rate limit: ",
             print_rate(x$rate_limit)), sep = "\n")
  cat(paste0("  retry (n/delay (s)): ",
             paste0(x$retry$n, "/", x$retry$time)), sep = "\n")
  cat(paste0("  error handler: ",
             names(x$error)), sep = "\n")
  cat("  config: ", sep = "\n")
  if (!is.null(x$config)) print(x$config, sep = "\n")
}

print_rate <- function(z) {
  if (!is.null(z)) {
    z2 <- unname(unlist(z))
    sprintf("%s @ %s - on_limit: %s", z2[1], z2[2], z2[3])
  }
}

print_heads <- function(x) {
  if (is.logical(x) || is.null(x) || is.character(x)) {
    return(x)
  } else {
    print_lazy(as.list(x$headers))
  }
}

print_body <- function(x) {
  if ("body_value" %in% names(x) && length(x) == 1) x <- unlist(unname(x))
  if (is.logical(x) || is.null(x) || is.character(x)) {
    return(x)
  } else if (any(grepl("upload_file", x[[1]]))) {
    "  File Upload"
  } else {
    print_lazy(x)
  }
}

print_lazy <- function(x) {
  out <- list()
  for (i in seq_along(x)) {
    val <- if (is(x[[i]], "name")) {
      deparse(x[[i]])
    } else {
      x[[i]]
    }
    out[[i]] <- sprintf("    %s: %s", names(x)[i], val)
  }
  return(paste0("\n", paste0(out, collapse = "\n")))
}
