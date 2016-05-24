
#' Countdown before deadline.
#' @param deadline     POSIXct. deadline
#' @param complete_str string. displayed when deadline is passed
#' @export
countdown = function(deadline, complete_str = intToUtf8(10004)) {
  days_left = difftime(deadline, Sys.time(), units="days")
  if (days_left>1)
    return(paste("J -", floor(days_left)))
  else if (days_left>0)
    return(paste("H -", floor(days_left*24)))
  else
    return(complete_str)
}

#' Formatted last update date before deadline.
#' @param deadline POSIXct. deadline
#' @param format   string. see \code{\link{format.POSIXct}}
#' @export
last_update <- function(deadline, format = "%d %b %Y %H:%M") {
  last = Sys.time()
  if (last>deadline+1)
    last = deadline+1
  return(format(last, format))
}

#' Print read errors.
#' @param read_err list of read errors returned by \code{\link{store_new_submissions}}
#' @export
#' @return \code{NULL}
print_readerr <- function(read_err = list()) {
  if (length(read_err)==0)
    cat("No read error.\n")
  for (i in seq(along=read_err)) {
    cat("Team", names(read_err)[i], ":\n")
    for (j in seq(along=read_err[[i]])) {
      cat("   ", paste(basename(names(read_err[[i]])[j]), ": "))
      cat(read_err[[i]][[j]]$message, "\n")
    }
  }
  invisible(NULL)
}

#' String displayed for the rank.
#' 
#' Concatenates the rank number with symbols indicating the progress since the last change.
#' 
#' @param r     integer. rank
#' @param r_d   integer. rank difference
#' @param symb  named list of characters. symbols used for the progress in ranking:
#'   no change (\code{const}), ascent (\code{up}) and descent (\code{down})
#' @keywords internal
str_rank <- function(r, r_d, symb = list(const = html_img(glyphicon("right_arrow"), "10px"),
                                         up = html_img(glyphicon("up_arrow"), "10px"),
                                         down = html_img(glyphicon("down_arrow"), "10px"))) {
  paste0(r, '. ', ifelse(r_d==0,
                         symb["const"], 
                         ifelse(r_d<0,
                                paste(rep(symb["up"], -r_d), collapse=""), 
                                paste(rep(symb["down"], r_d), collapse=""))))
}


#' Format the leaderboard in Markdown.
#' 
#' @param best    list of the best submissions per team and per metric as returned
#'   by \code{\link{get_best}}.
#' @param metric  string. name of the metric considered
#' @param test_name string. name of the test set used: \code{"quiz"} or \code{"test"}
#' @param ... further parameters to pass to \code{\link[knitr]{kable}}er
#' 
#' @return \code{print_leaderboard} returns a character vector of the table source code
#'   to be used in a Markdown document.
#'   
#' @note Chunk option \code{results='asis'} has to be used
#' 
#' @export
#' @seealso \code{\link[knitr]{kable}}
#' @importFrom knitr kable
print_leaderboard <- function(best, metric, test_name = "quiz", ...) {
  metric_column = paste(metric, test_name, sep=".")
  leaderboard = data.frame(Rank = mapply(FUN = str_rank, best[[metric]]$rank, best[[metric]]$rank_diff),
                           Team = best[[metric]]$team,
                           Submissions = paste(best[[metric]]$n_submissions),
                           Date = format(best[[metric]]$date, format="%d/%m/%y %H:%M"),
                           Score = format(best[[metric]][[metric_column]], digits=3))
  knitr::kable(leaderboard, ...)
}

#' Path to glyphicon image file.
#' @param name string. name of the glyphicon.
#' @param path string. folder of search.
#' @return the path to the file.
#' @export
glyphicon <- function(name, path = system.file('glyphicons', package = 'rchallenge')) {
  file <- list.files(path, pattern = paste("glyphicons_[0-9]+_", name, ".png", sep=""))
  if (length(file)==0)
    file <- list.files(path, pattern = paste("glyphicons_social_[0-9]+_", name, ".png", sep=""))
  if (length(file)==0)
    stop("glyphicon ", name, " not found")
  return(file.path(path, file))
}

#' html code for an image.
#' @param file string. image file.
#' @param width string. width of display.
#' @export
html_img <- function(file, width = "10px") {
  paste('<img src="', file, '" style="width: ', width, ';"/>', sep="")
}
