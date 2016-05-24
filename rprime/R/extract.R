
#' Extract log-frames from an Eprime log file
#'
#' Almost all of the information in an Eprime file comes in chunks of text
#' bracketed by the lines \code{*** LogFrame Start ***} and \code{*** LogFrame
#' End ***}. The exception is the header information which is bracketed by
#' \code{*** Header Start ***} and \code{*** Header End ***}.
#'
#' \code{extract_chunks} extracts the bracketed text, storing each log-frame of
#' text in a list. The lists also include some additional lines of text as
#' metadata: \code{Eprime.FrameNumber} and \code{Eprime.Basename} (the name of
#' the source file). The header log-frame also gets dummy lines:
#' \code{Procedure: Header} and \code{Running: Header}.
#'
#' These chunks of colon-separated lines are converted into lists by
#' \code{FrameList(...)}.
#'
#' @param eprime_log a character vector containing the lines of text from Eprime
#'   txt file
#' @return a list of character vectors, where each vector contains the lines of
#'   a log-frame
#' @export
extract_chunks <- function(eprime_log) {
  parsed <- parse_chunks(eprime_log)
  basename <- ifelse(!has_attr(eprime_log, "basename"), NA,
                     attr(eprime_log, "basename"))
  parsed <- parsed %>%
    update_header %>%
    insert_frame_numbers %>%
    insert_basename(basename)
  lapply(parsed, as.EprimeChunk)
}

as.EprimeChunk <- function(x) {
  class(x) <- c("EprimeChunk", class(x))
  x
}

# Add "Running", "Procedure" lines to the header log-frame
update_header <- function(chunked) {
  if (has_header(chunked)) {
    header_position <- Position(is_header, chunked)
    header <- chunked[[header_position]]
    row_run <- new_line("Running", "Header")
    row_prc <- new_line("Procedure", "Header")
    header <- insert_line(header, c(row_run, row_prc))
    chunked[[header_position]] <- header
  }
  chunked
}

# Add "Eprime.FrameNumber" lines to every frame
insert_frame_numbers <- function(chunked) {
  rows <- new_line(rprime_cols$frame, seq_along(chunked))
  Map(insert_line, chunked, rows)
}

# Add "Eprime.Basename" lines to every log-frame
insert_basename <- function(chunked, basename) {
  rows <- new_line(rprime_cols$basename, basename)
  Map(insert_line, chunked, rows)
}

# Insert a line in the second-to-last position in a log-frame
insert_line <- function(xs, ys) {
  c(but_last(xs), ys, last(xs))
}

#' Extract the lines of text from each log-frame
#'
#' @param eprime_log a character vector of lines from an Eprime log file
#' @return a list of character vectors, where each vector contains the lines of
#'   a log-frame
#' @keywords internal
parse_chunks <- function(eprime_log) {
  # Find all the texts between LogFrame boundaries
  starts <- str_which(eprime_log, patterns$bracket_start)
  ends <- str_which(eprime_log, patterns$bracket_end)
  ranges <- make_ranges(starts, ends, eprime_log)

  pull_lines <- function(lines) eprime_log[lines]
  frames <- lapply(ranges, pull_lines)
  frames
}


#' Check Eprime log-frame line ranges
#'
#' @param starts the line numbers of the log-frame start lines
#' @param ends the line numbers of the log-frame end lines
#' @param eprime_log a character vector of lines from an Eprime log file
#' @return a list of sequences. Each list contains the line numbers contained by
#'   an Eprime log-frame. If there is a log frame without an end-line, that
#'   partial frame is excluded and its contents are previed in a warning
#'   message.
#' @keywords internal
make_ranges <- function(starts, ends, eprime_log) {
  # There should be the same number of starts and ends
  min_chunks <- min(length(starts), length(ends))
  old_starts <- starts
  starts <- starts[seq_len(min_chunks)]
  ends <- ends[seq_len(min_chunks)]

  # Warn if there is an incomplete frame (more old_starts than ends)
  bad_lines <- setdiff(old_starts, starts)
  if (!length_zero(bad_lines)) {
    last_bad_line <- max(bad_lines)

    # Give a preview of the incomplete chunk in the warning
    last_preview_line <- min(last_bad_line + 10, length(eprime_log))
    sample_range <- seq(last_bad_line, last_preview_line)
    warning_header <- paste0("Incomplete Log Frame found on line ", bad_lines)
    lines <- paste0(c(warning_header, eprime_log[sample_range]), collapse = "\n")
    warning(lines)
  }

  # Each start should come before its corresponding end
  well_ordered_pairs <- starts < ends
  starts <- starts[well_ordered_pairs]
  ends <- ends[well_ordered_pairs]

  Map(seq, starts, ends)
}

