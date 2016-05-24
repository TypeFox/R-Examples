#' Ping a url to time the request
#'
#' @export
#'
#' @param .request A httr response object
#' @param count integer, Number of requests to do.
#' @param delay integer, Seconds to delay successive calls by. Default: 0.5 seconds.
#' @param flood logical; If TRUE, no delay between requests. If FALSE, delay by 0.5
#' second.
#' @param verbose logical; If TRUE, print progress.
#' @param ... Further args passed on to functions in \code{httr}
#' @examples \dontrun{
#' GET("http://httpbin.org/get") %>% time()
#' GET("https://api.github.com") %>% time()
#' GET("http://google.com") %>% time()
#' }

time <- function(.request, count=10, delay = 0.5, flood = FALSE, verbose=TRUE, ...) {
  stopifnot(is(.request, "response"))
  if (flood) delay <- 0 else stopifnot(is(as.numeric(delay), "numeric"))
  if (verbose) cat(sprintf("%s kb - %s code:%s time(ms):%s", get_kb(.request), .request$url, .request$status_code, .request$times[["total"]]*1000), sep = "\n")
  if (count < 2) stop("count parameter must be greater than 1", call. = FALSE)
  if (count > 1) count_ <- count - 1
  reps <- replicate(count_, rerequest_(.request, delay, verbose), simplify = FALSE)
  all <- do.call(c, list(list(.request), reps))
  times <- pluck(all, "times")
  nmz <- names(times[[1]])
  avgs <- setNames(lapply(nmz, function(z){
    tt <- pluck(times, z, double(1))
    setNames(sapply(list(min, max, mean), function(x) calc(tt, x)), c("min","max","mean"))
  }), nmz)
  structure(list(times = times, averages = avgs, request = .request),
            count = count,
            delay = delay,
            class = "http_time")
}

calc <- function(x, fxn) format(fxn(x*1000, na.rm = TRUE), scientific = FALSE)

rerequest_ <- function(x, delay, verbose=TRUE){
  Sys.sleep(delay)
  tmp <- rerequest(x)
  if (verbose) cat(sprintf("%s kb - %s code:%s time(ms):%s", get_kb(tmp), tmp$url, tmp$status_code, tmp$times[["total"]]*1000), sep = "\n")
  return(tmp)
}

get_kb <- function(x) (pryr::object_size(x)/1000)[[1]]

#' @export
print.http_time <- function(x, ...){
  cat("<http time>", sep = "\n")
  cat(paste0("  Avg. min (ms):  ", x$averages[['total']][['min']]), sep = "\n")
  cat(paste0("  Avg. max (ms):  ", x$averages[['total']][['max']]), sep = "\n")
  cat(paste0("  Avg. mean (ms): ", x$averages[['total']][['mean']]), sep = "\n")
}

#' @export
summary.http_time <- function(object, ...){
  cat("<http time, averages (min max mean)>", sep = "\n")
  cat(paste0("  Total (s):           ", p0(object$averages[['total']])), sep = "\n")
  cat(paste0("  Tedirect (s):        ", p0(object$averages[['redirect']])), sep = "\n")
  cat(paste0("  Namelookup time (s): ", p0(object$averages[['namelookup']])), sep = "\n")
  cat(paste0("  Connect (s):         ", p0(object$averages[['connect']])), sep = "\n")
  cat(paste0("  Pretransfer (s):     ", p0(object$averages[['pretransfer']])), sep = "\n")
  cat(paste0("  Starttransfer (s):   ", p0(object$averages[['starttransfer']])), sep = "\n")
}

p0 <- function(x) paste0(x, collapse = " ")
