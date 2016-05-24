#' Parse Rprof output.
#'
#' Parses the output of \code{\link{Rprof}} into an alternative format
#' described in \code{\link{profr}}. This produces a flat data frame, which is
#' somewhat easier to summarise and visualise.
#'
#' @param path path to \code{\link{Rprof}} output
#' @param interval real-time interval between samples (in seconds)
#' @keywords debugging
#' @return \code{\link{data.frame}} of class \code{profr}
#' @seealso \code{\link{profr}} for profiling and parsing
#' @import stringr plyr
#' @export
#' @examples
#' nesting_ex <- system.file("samples", "nesting.rprof", package="profr")
#' nesting <- parse_rprof(nesting_ex)
#'
#' reshape_ex <- system.file("samples", "reshape.rprof", package="profr")
#' diamonds <- parse_rprof(reshape_ex)
parse_rprof <- function(path, interval=0.02) {
  lines <- readLines(path)[-1]

  calls <- str_split(lines, " ")
  calls <- lapply(calls, function(x) rev(str_replace_all(x, "\"", ""))[-1])

  df <- .simplify(calls)

  times <- c("time", "start", "end")
  df[times] <- df[times] * interval

  df
}

group_id <- function(x, y) {
  n <- length(x)
  cumsum(c(TRUE, x[-1] != x[-n]))
}


.simplify <- function(calls) {
  df <- ldply(seq_along(calls), function(i) {
    call <- calls[[i]]
    call_info(call, i - 1)
  })
  df$hist <- id(list(df$hist))

  # A group consists of all calls with the same history, in a
  # consecutive block of time
  levels <- ddply(df, "level", function(df) {
    mutate(df,
      g_id = group_id(hist),
      t_id = cumsum(c(TRUE, diff(start) != 1))
    )
  })

  collapsed <- ddply(levels, c("level", "g_id", "t_id"), summarise,
    f = f[1],
    start = min(start),
    end = max(end),
    n = length(f),
    leaf = leaf[1]
  )
  collapsed <- mutate(collapsed,
    time = end - start,
    source = function_source(f)
  )
  # subset(collapsed, time != n)

  structure(collapsed, class = c("profr", "data.frame"))
}

call_info <- function(call, i) {
  n <- length(call)
  history <- unlist(lapply(seq_along(call), function(i) {
    paste(call[seq_len(i)], collapse = "")
  }))

  quickdf(list(
    f = call,
    level = seq_along(call),
    start = rep(i, n),
    end = rep(i + 1, n),
    leaf = c(rep(FALSE, n - 1), TRUE),
    hist = history
  ))
}

function_source <- function(f) {
  pkgs <- search()
  names(pkgs) <- pkgs
  all_objs <- ldply(pkgs, as.data.frame(ls))
  names(all_objs) <- c("package", "f")
  all_objs$package <- str_replace_all(all_objs$package, "package:", "")

  all_objs$package[match(f, all_objs$f)]
}
