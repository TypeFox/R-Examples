RequestIterator <- R6::R6Class("RequestIterator",
  public = list(
  result = list(),
  limit = NA,
  limit_max = NA,
  links = list(),
  initialize = function(result, limit, limit_max, links) {
    if (!missing(result)) self$result <- result
    if (!missing(limit)) self$limit <- limit
    if (!missing(limit_max)) self$limit_max <- limit_max
    if (!missing(links)) self$links <- links
  },
  GET = function(.data, ...) {
    if (length(self$links) == 0) {
      .data <- as.req(.data)
      .data$config <- c(httr::user_agent(make_ua()), .data$config, .data$headers)
      .data$url <- gather_paths(.data)
      .data$query <- if (is.null(.data$query)) NULL else as.list(.data$query)
      res <- suppressWarnings(httr::GET(.data$url[1], .data$config, .data$write,
                                        query = .data$query, ...))
    } else {
      .data <- as.req(self$links[[1]]$url)
      res <- suppressWarnings(httr::GET(.data$url[1], .data$config, .data$write, ...))
    }
    # error catching
    self$handle_errors(.data, res)
    # cache links
    self$links <- get_links(res$headers)
    # give back result
    self$result <- empty(list(self$result, res))
  },
  POST = function(.data, ...) {
    if (length(self$links) == 0) {
      .data <- as.req(.data)
      .data$config <- c(httr::user_agent(make_ua()), .data$config)
      .data$url <- gather_paths(.data)
      res <- suppressWarnings(httr::POST(.data$url[1], .data$config, body = .data$body, ...))
    } else {
      .data <- as.req(self$links[[1]]$url)
      res <- suppressWarnings(httr::POST(.data$url[1], .data$config, ...))
    }
    # error catching
    self$handle_errors(.data, res)
    # cache links
    self$links <- get_links(res$headers)
    # give back result
    self$result <- empty(list(self$result, res))
  },
  body = function() {
    self$result
  },
  status = function() {
    if (is(self$result, "response")) {
      self$result$status_code
    } else {
      lapply(self$result, function(z) {
        if (is(z, "response")) {
          z$status_code
        } else {
          NULL
        }
      })
    }
  },
  parse = function(parse = TRUE) {
    x <- self$result
    if (is(x, "response")) {
      httr_parse(x)
    } else {
      lapply(x, httr_parse)
    }
  },
  count = function() {
    if (is(self$result, "response")) {
      length(httr::content(self$result))
    } else {
      sum(sapply(self$result, function(x) length(httr::content(x))))
    }
  },
  handle_errors = function(.data, x) {
    if (!is.null(.data$retry)) {
      i <- 0
      while (x$status_code > 201 && i < .data$retry$n) {
        i <- i + 1
        message("Retrying request\n")
        x <- self$GET(.data)
        Sys.sleep(.data$retry$time)
      }
      return(x)
    }
    if (is.null(.data$error)) {
      # httr::stop_for_status(x)
      try_error(x)
    } else {
      .data$error[[1]](x)
    }
  }
))

try_error <- function(x) {
  if (x$status_code > 201) {
    one <- tryCatch(content(x), error = function(e) e)
    if (!is(one, "error")) {
      two <- tryCatch(one$error, error = function(e) e)
      if (!is(two, "error")) {
        msg <- sprintf("%s - %s", x$status_code, two)
      } else {
        msg <- http_status(x)$message
      }
    } else {
      msg <- http_status(x)$message
    }
    stop(msg, call. = FALSE)
  }
}

httr_parse <- function(x) {
  if (grepl("json", x$headers$`content-type`)) {
    if (!is.null(x$request$output$path)) {
      return(x$request$output$path)
    } else {
      txt <- httr::content(x, "text")
      jsonlite::fromJSON(txt, parse, flatten = TRUE)
    }
  } else {
    content(x)
  }
}

# dd <- api('https://api.github.com/') %>%
#   api_path(repos, ropensci, rgbif, commits) %>%
#   api_paging(limit = 220, limit_max = 100)
# # dd <- api('https://api.github.com/') %>%
# #   api_path(repos, ropensci, rplos, commits) %>%
# #   api_paging(limit = 220, limit_max = 100)
# rr <- GetIter$new(limit = dd$paging$limit, limit_max = dd$paging$limit_max)
# rr$GET(dd)
# rr$count()
# rr$result
# rr$links
# rr$parse()
