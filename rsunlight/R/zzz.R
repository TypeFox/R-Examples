sc <- function(x){
  Filter(Negate(is.null), x)
}

return_obj <- function(x, y){
  y <- err_hand(y)
  x <- match.arg(x, c('response', 'list', 'table', 'data.frame'))
  if (x == 'response') {
    y
  } else {
    if (x == 'list') {
      jsonlite::fromJSON(y, simplifyVector = FALSE, flatten = TRUE)
    } else {
      jsonlite::fromJSON(y, flatten = TRUE)
    }
  }
}

# check if stupid single left bracket returned
err_hand <- function(z) {
  tmp <- httr::content(z, "text")
  if (identical(tmp, "[")) {
    q <- httr::parse_url(z$request$opts$url)$query
    q <- paste0("\n - ", paste(names(q), q, sep = "="), collapse = "")
    stop("The following query had no results:\n", q, call. = FALSE)
  } else {
    tmp
  }
}

give_noiter <- function(as, url, endpt, args, ...) {
  tmp <- return_obj(as, query(paste0(url, endpt), args, ...))
  switch(as,
         table = structure(tmp, class = c("sunlight", "data.frame")),
         list = tmp,
         response = tmp)
}

give <- function(as, url, endpt, args, ...) {
  iter <- get_iter(args)
  if (length(iter) == 0) {
    tmp <- return_obj(as, query(paste0(url, endpt), args, ...))
  } else {
    tmp <- lapply(iter[[1]], function(w) {
      args[[ names(iter) ]] <- w
      return_obj(as, query(paste0(url, endpt), args, ...))
    })
    if (as == "table") {
      tmp <- tmp[vapply(tmp, length, numeric(1)) != 0]
      tmp <- plyr::rbind.fill(tmp)
    }
  }
  switch(as,
         table = structure(tmp, class = c("sunlight", "data.frame")),
         list = tmp,
         response = tmp)
}

give_cg <- function(as, url, endpt, args, ...) {
  iter <- get_iter(args)
  if (length(iter) == 0) {
    tmp <- return_obj(as, query(paste0(url, endpt), args, ...))
  } else {
    tmp <- lapply(iter[[1]], function(w) {
      args[[ names(iter) ]] <- w
      return_obj(as, query(paste0(url, endpt), args, ...))
    })
    if (as == "table") {
      res <- lapply(tmp, "[[", "results")
      res <- Filter(function(x) !(is.null(x) || length(x) == 0), res)
      for (i in seq_along(res)) {
        res[[i]][names(iter)] <- iter[[1]][i]
      }
      res <- res[vapply(res, length, numeric(1)) != 0]
      tmp <- rbind.fill(res)
    }
  }
  switch(as,
         table = structure(if (length(iter) == 0) tmp$results else tmp,
                           class = c("sunlight", "data.frame")),
         list = tmp,
         response = tmp)
}

one_vec <- function(x) {
  lens <- x[vapply(x, length, 1) > 1]
  if (length(lens) > 1) {
    stop("Only one parameter can be vectorized per function call", call. = FALSE)
  }
}

get_iter <- function(z) {
  z[vapply(z, length, 1) > 1]
}

#' @export
print.sunlight <- function(x, ..., n = 10){
  cat("<Sunlight data>", sep = "\n")
  cat(sprintf("   Dimensions:   [%s X %s]\n", NROW(x), NCOL(x)), sep = "\n")
  trunc_mat(x, n = n)
}

check_key <- function(x){
  tmp <- if (is.null(x)) {
    Sys.getenv("SUNLIGHT_LABS_KEY", "")
  } else {
    x
  }

  if (tmp == "") {
    getOption("SunlightLabsKey", stop("need an API key for Sunlight Labs"))
  } else {
    tmp
  }
}

cgurl <- function() 'https://congress.api.sunlightfoundation.com'
cwurl <- function() 'http://capitolwords.org/api'
ieurl <- function() 'http://transparencydata.com/api/1.0'
osurl <- function() 'http://openstates.org/api/v1'
rtieurl <- function() "http://realtime.influenceexplorer.com/api/"

