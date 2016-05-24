#' Convert Eprime Frames into data-frames
#'
#' @details Individual EprimeFrames are converted to a data-frame using
#'   \code{as.data.frame}. (Strings are not converted to factors.)
#'
#'   Each of the individual data-frames are then \code{rbind}ed together, with
#'   missing columns being filled with NA.
#'
#' @param x an EprimeFrame object, or a FrameList object (a list of
#'   EprimeFrames)
#' @return all of the EprimeFrames combined into a single data frame.
#' @export
#' @seealso \link{rbind.fill}
to_data_frame <- function(x) {
  UseMethod("to_data_frame")
}

#' @export
to_data_frame.default <- function(x) {
  data_frames <- lapply(x, as.data.frame.list, stringsAsFactors = FALSE)
  rbind.fill(data_frames)
}

#' @export
to_data_frame.FrameList <- function(x) {
  data_frames <- lapply(x, to_data_frame)
  rbind.fill(data_frames)
}

#' @export
to_data_frame.EprimeFrame <- function(x) {
  as.data.frame.list(x, stringsAsFactors = FALSE)
}

