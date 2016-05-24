#' Make a filtering predicate
#' @keywords internal
make_filter <- function(key, values) {
  assert_that(length(key) == 1)
  function(frame_list) {
    # Convert NULLs to NA
    plucked <- pluck_apply(key, frame_list)
    plucked[sapply(plucked, is.null)] <- NA
    is.element(unlist(plucked), values)
  }
}



#' Filter levels in or out of a FrameList based on attribute values
#'
#' @param frame_list a list of \code{EprimeFrame} objects
#' @param key the name of the attribute to filter in or out
#' @param values the whitelisted or blacklisted values of the attribute
#' @return for \code{filter_in}, only log-frames where \code{key} is one of the
#'   \code{values} are kept. for \code{filter_out}, log-frames where \code{key}
#'   is one of the \code{values} are omitteed.
#' @export
filter_in <- function(frame_list, key, values) {
  has_key_value <- make_filter(key, values)
  as.FrameList(frame_list[has_key_value(frame_list)])
}

#' @rdname filter_in
#' @export
filter_out <- function(frame_list, key, values) {
  lacks_key_value <- Negate(make_filter(key, values))
  as.FrameList(frame_list[lacks_key_value(frame_list)])
}

#' Filter levels in or out of a FrameList based on Eprime.Level values
#'
#' These functions are shortcuts for calls to \code{filter_in} or
#' \code{filter_out}.
#'
#' Note that the meaning of Eprime.Level value in a log-frame ultimately is
#' equal to one plus the number of tabs before each line in the log-frame.
#'
#' @inheritParams filter_in
#' @param level_numbers the whitelisted or blacklisted values for Eprime.Level
#' @return for \code{keep_levels}, only log-frames where the level matches one
#'   of the \code{level_numbers} are kept. for \code{keep_levels}, log-frames
#'   where the level matches one of the \code{level_numbers} are omitted.
#' @export
keep_levels <- function(frame_list, level_numbers) {
  filter_in(frame_list, rprime_cols$level, as.character(level_numbers))
}

#' @rdname keep_levels
#' @export
drop_levels <- function(frame_list, level_numbers) {
  filter_out(frame_list, rprime_cols$level, as.character(level_numbers))
}
