
#' Client size load balancing
#'
#' @docType package
#' @name spareserver-package
NULL

spare_services <- new.env()

clean_services <- function() {
  rm(list = ls(spare_services), envir = spare_services)
}

.onLoad <- function(libname, pkgname) {
  ub <- unlockBinding
  ub("spare_services", asNamespace(pkgname))
  invisible()
}

.onAttach <- function(libname, pkgname) {
  ub <- unlockBinding
  ub("spare_services", asNamespace(pkgname))
  invisible()
}

#' Create a spare server object
#'
#' @param base_url Base URL of the server.
#' @param priority Non-negative initial odds of using this server, it has
#'   no units. Note that a server with zero odds can be used if all other
#'   servers are down.
#' @param timeout Timeout in seconds. After this time, we give up on
#'   the server.
#' @return A spare server object
#'
#' @export

server <- function(base_url, priority = 1, timeout = 1) {
  structure(
    list(
      base_url = base_url,
      priority = priority,
      timeout = timeout,
      state = "unknown",
      timestamp = now()
    ),
    class = "spare_server"
  )
}

#' Create a redundant service
#'
#' @param name Name of the service.
#' @param ... Server objects, zero or more. Servers can also be added
#'   later, using \code{add_server}. See \code{server}.
#' @return The name of the service, invisibly.
#'
#' @export
#' @importFrom assertthat assert_that is.string

add_service <- function(name, ...) {
  assert_that(is.string(name))
  if (name %in% ls(spare_services)) {
    stop("Service already exists", call. = FALSE)
  }
  if (any(! sapply(list(...), is, class2 = "spare_server"))) {
    stop("Not a server object", call. = FALSE)
  }

  assign(name, list(servers = list(...)), envir = spare_services)

  invisible(name)
}

#' List redundant services
#'
#' @return List of services, and all servers for them.
#'
#' @export

services <- function() {
  as.list(spare_services)
}

#' Remove a redundant service
#'
#' @param name The name of the service to remove.
#' @return Nothing.
#'
#' @export
#' @importFrom assertthat assert_that is.string

remove_service <- function(name) {
  assert_that(is.string(name))
  if (! name %in% ls(spare_services)) {
    stop("Unknown service", call. = FALSE)
  }

  rm(list = name, envir = spare_services)

  invisible()
}

#' Add a spare server to a service
#'
#' @param service Name of the service.
#' @param ... Server objects to add, or modify. To remove a server,
#'   set it to \code{NULL} here.
#' @return Name of the service, invisibly.
#'
#' @export
#' @importFrom assertthat assert_that is.string
#' @importFrom methods is
#' @importFrom utils modifyList

add_server <- function(service, ...) {
  assert_that(is.string(service))
  if (! service %in% ls(spare_services)) {
    stop("Unknown service", call. = FALSE)
  }
  if (any(! (sapply(..., is, class2 = "spare_server")) |
            (sapply(..., is.null)))) {
    stop("Not a server object or NULL", call. = FALSE)
  }

  new_servers <- modifyList(get(service, envir = spare_services),
                            list(...))
  assign(service, new_servers, envir = spare_services)

  invisible(service)
}

#' Make a web query that makes use of redundant servers
#'
#' @param service Name of the service.
#' @param url URL, within the base URL of the server(s). This will be
#'   appended to the base URLs.
#' @param fun The function to call. It should have an argument called
#'   \code{url}.
#' @param ... Additional arguments to pass to \code{fun}.
#' @return The return value of \code{fun}.
#'
#' @export

spare_q <- function(service, url, fun, ...) {
  assert_that(is.string(service))
  if (! service %in% ls(spare_services)) {
    stop("Unknown service", call. = FALSE)
  }
  servers <- get(service, envir = spare_services)$servers
  odds <- sapply(servers, "[[", "priority")
  robust_q(service, url, fun, list(...), servers, odds)
}

## Make a robust query. We can check the internet connection in
## general, but we cannot rely on that completely, because the
## specified server can be still available. E.g. if it is on
## localhost.
##
## Goals:
## * minimal latency for the user, obviously.
## * get out of the way, sparing should go unnoticed (if possible)
## * deterministic (?)
##
## Algorithm
## * Each server has a state, with a time label. The state is simply 'on',
##   'off' or 'unknown'.
## * States expire, relatively quickly, e.g. in 5 or 10 minutes.
##   Then they become 'unknown', effectively.
## * Sort the servers according to their priorities.
## * Find the first server with an 'on' state. If needed, ping
##   servers in an 'unknown' state, to see if they are up.
## * Try the server with the 'on' state. If it works, good, update
##   its time stamp.
## * Otherwise set its state to 'off' with the current timestamp,
##   and continue with the next server.
##
## If all servers are down, then we start over, with fives times bigger
## timeout. We only do two rounds currently.

robust_q <- function(service, url, fun, args, servers, odds,
                     timeout_multiplier = 1.0, no_rounds = 2) {
  order <- order(odds, decreasing = TRUE)
  servers <- servers[order]
  for (i in seq_along(servers)) {
    ## Unknown or expired?
    if (servers[[i]]$state == "unknown" ||
        now() - servers[[i]]$timestamp > as.difftime(3, units = "mins")) {
      servers[[i]] <- ping_server(servers[[i]], mult = timeout_multiplier)
    }
    if (servers[[i]]$state == "on") {
      res <- try_server(servers[[i]], url, fun, args,
                        mult = timeout_multiplier)
      if (!inherits(res, "try-error")) {
        servers[[i]]$state <- "on"
        servers[[i]]$timestamp <- now()
        obj <- get(service, spare_services)
        obj$servers <- servers
        assign(service, obj, envir = spare_services)
        return(res)
      } else {
        servers[[i]]$state <- "off"
        servers[[i]]$timestamp <- now()
      }
    }
  }
  obj <- get(service, spare_services)
  obj$servers <- servers
  assign(service, obj, envir = spare_services)
  if (no_rounds == 1) {
    stop("Cannot do query '", url, "'.")
  } else {
    for (i in seq_along(servers)) servers[[i]]$state <- "unknown"
    robust_q(service, url, fun, args, servers, odds,
             timeout_multiplier = timeout_multiplier * 50,
             no_rounds = no_rounds - 1)
  }
}

## Check if it is up

#' @importFrom pingr ping_port
#' @importFrom httr parse_url

ping_server <- function(server, mult) {
  parsed_url <- parse_url(server$base_url)
  protocol <- parsed_url$scheme
  host <- parsed_url$hostname
  port <- parsed_url$port

  if (protocol == "") { protocol <- "http" }
  if (is.null(port) && protocol == "http") { port <- 80 }
  if (is.null(port) && protocol == "https") { port <- 443 }

  resp_time <- ping_port(host, port = port, count = 1,
                         timeout = server$timeout * mult)
  server$state <- if (is.na(resp_time)) "off" else "on"
  server$timestamp <- now()
  server
}

## Try to make a query

#' @importFrom httr timeout

try_server <- function(server, url, fun, args, mult) {
  full_url <- paste0(server$base_url, url)
  all_args <- c(list(full_url, httr::timeout(server$timeout * mult)), args)
  try(silent = TRUE, do.call(fun, all_args))
}

regexp_to_df <- function(text, match) {
  positive <- match != -1
  g_text <- text[positive]
  g_start <- attr(match, "capture.start")[positive, , drop = FALSE]
  g_length <- attr(match, "capture.length")[positive, , drop = FALSE]

  lapply(seq_len(sum(positive)), function(i) {
    data.frame(row.names = attr(match, "capture.names"),
               start = g_start[i,],
               length = g_length[i,],
               match = substring(text[i], g_start[i,],
                 g_start[i,] + g_length[i,] - 1),
               stringsAsFactors = FALSE)
  })
}

now <- function() {
  as.POSIXct(Sys.time())
}
