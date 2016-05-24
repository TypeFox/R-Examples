#' Convert lines from an Eprime file into EprimeFrame objects
#'
#' Convert character vectors of implicit key-value pairs (e.g., \code{c("key1:
#' value1", "key2: value2")}), into  lists of explicit key-value pairs,
#' \code{list(key1 = "value1", key2 = "value2")}.
#'
#' @details During the conversion, if \code{Running: x}, then the
#'   \code{x.Sample} and \code{x.Cycle} lines are simplified into \code{Sample}
#'   and \code{Cycle} lines. The \code{x: value} line is recoded as
#'   \code{Eprime.LevelName: x_value}. The purpose of this tidying is to force
#'   the same set of key names (eventually, column names) onto frames with
#'   different values for "Running".
#'
#' @param x a character vector with lines of the form \code{"key: value"}, or a
#'   list of vectors of colon-separated text
#' @return When passed a list of character vectors of \code{"key: value"} lines,
#'   a FrameList object (a list of EprimeFrames) is returned. when passed a
#'   single vector vector of \code{"key: value"} lines, a single EprimeFrame
#'   object is returned inside of a FrameList object.
#' @export
#' @examples
#' lines <- c("\t*** LogFrame Start ***",
#'            "\tProcedure: FamTask",
#'            "\titem1: bear",
#'            "\titem2: chair",
#'            "\tCorrectResponse: bear",
#'            "\tImageSide: Left",
#'            "\tDuration: 885",
#'            "\tFamiliarization: 1",
#'            "\tFamInforcer: 1",
#'            "\tReinforcerImage: Bicycle1",
#'            "\tFamiliarization.Cycle: 1",
#'            "\tFamiliarization.Sample: 1",
#'            "\tRunning: Familiarization",
#'            "\tFamTarget.RESP: ",
#'            "\tCorrect: True",
#'            "\t*** LogFrame End ***")
#' # List of 1
#' # $ :List of 17
#' # ..$ Eprime.Level      : num 2
#' # ..$ Eprime.LevelName  : chr "Familiarization_1"
#' # ..$ Eprime.Basename   : chr "NA"
#' # ..$ Eprime.FrameNumber: chr "1"
#' # ..$ Procedure         : chr "FamTask"
#' # ..$ Running           : chr "Familiarization"
#' # ..$ item1             : chr "bear"
#' # ..$ item2             : chr "chair"
#' # ..$ CorrectResponse   : chr "bear"
#' # ..$ ImageSide         : chr "Left"
#' # ..$ Duration          : chr "885"
#' # ..$ FamInforcer       : chr "1"
#' # ..$ ReinforcerImage   : chr "Bicycle1"
#' # ..$ Cycle             : chr "1"
#' # ..$ Sample            : chr "1"
#' # ..$ FamTarget.RESP    : chr ""
#' # ..$ Correct           : chr "True"
#' # ..- attr(*, "class")= chr [1:2] "EprimeFrame" "list"
#' # - attr(*, "class")= chr [1:2] "list" "FrameList"
FrameList <- function(x) UseMethod("FrameList")

#' @export
FrameList.character <- function(x) {
  FrameList.list(extract_chunks(x))
}

#' @export
FrameList.list <- function(x) {
  assert_that(is_list_of(x, "EprimeChunk") | is_list_of(x, "character"))
  as.FrameList(lapply(x, EprimeFrame))
}



#' Create an EprimeFrame object
#'
#' This constructor function converts a character vector into an
#' \code{EprimeFrame} object, which is just a list with some special metadata
#' values. Strings with the format \code{"key: value"} are parsed into \code{key
#' = value} list items (via \code{listify}).
#'
#' @param keys_values a character vector of containing some \code{"key: value"}
#'   strings.
#' @return a list with the class \code{EprimeFrame} and with special
#'   \code{Eprime.} metadata, \code{Running} and \code{Procedure} values, all
#'   set to NA by default.
#' @export
#' @examples
#' # Default metadata values
#' lines <- c(
#'   "key: value",
#'   "question: answer",
#'   "garbage text")
#'
#' EprimeFrame(lines)
#' # List of 8
#' # $ Eprime.Level      : num 1
#' # $ Eprime.LevelName  : logi NA
#' # $ Eprime.Basename   : logi NA
#' # $ Eprime.FrameNumber: logi NA
#' # $ Procedure         : logi NA
#' # $ Running           : logi NA
#' # $ key               : chr "value"
#' # $ question          : chr "answer"
#'
#' # Normalize [Running] related lines
#' keys_values <- c(
#'   "Running: Demo",
#'   "Demo: ExampleCode",
#'   "Demo.Cycle: 1",
#'   "Demo.Sample: 1",
#'   "Key: Value")
#'
#' EprimeFrame(keys_values)
#' # List of 9
#' # $ Eprime.Level      : num 1
#' # $ Eprime.LevelName  : chr "Demo_ExampleCode"
#' # $ Eprime.Basename   : logi NA
#' # $ Eprime.FrameNumber: logi NA
#' # $ Procedure         : logi NA
#' # $ Running           : chr "Demo"
#' # $ Cycle             : chr "1"
#' # $ Sample            : chr "1"
#' # $ Key               : chr "Value"
EprimeFrame <- function(keys_values) UseMethod("EprimeFrame")

#' @export
EprimeFrame.EprimeChunk <- function(keys_values = character(0)) {
  EprimeFrame.character(keys_values)
}

#' @export
EprimeFrame.character <- function(keys_values = character(0)) {
  keys_values <- merge_lists(listify(keys_values), count_tabs(keys_values))
  frame <- as.EprimeFrame(keys_values)
  tidy(frame)
}



#' Convert a list of EprimeFrames into a FrameList object
#' @param xs a list of EprimeFrames
#' @return the original list as a \code{FrameList} object
#' @export
as.FrameList <- function(xs) {
  assert_that(is_list_of(xs, "list"), is_list_of(xs, "EprimeFrame"))
  class(xs) <- unique(c(class(xs), "FrameList", "list"))
  xs
}

#' Convert a list into an EprimeFrame object
#' @param xs a list
#' @return the original list as an \code{EprimeFrame} object (along with dummy
#'   Eprime metadata fields)
#' @export
as.EprimeFrame <- function(xs) {
  assert_that(is.list(xs))
  with_defaults <- merge_lists(default_metadata, xs)
  structure(with_defaults, class = c("EprimeFrame", "list"))
}



#' @export
print.FrameList <- function(...) str(...)

#' @export
print.EprimeFrame <- function(...) str(...)



#' Clean up `Running`-related attributes
#'
#' A fresh EprimeFrame object has the structure:
#'
#' \code{
#'   Eprime.LevelName: NA
#'   Eprime.Level: [Level]
#'   Running: [Key]
#'   [Key]: [Value]
#'   [Key].Cycle: [Cycle]
#'   [Key].Sample: [Sample]
#' }
#'
#' These \code{[Key]} values make it harder to merge together data-frames later
#' on, since each unique \code{[Key]} gets its own column name. Therefore, we
#' normalize these field names early on. The end result has the structure:
#'
#' \code{
#'   Eprime.LevelName: [Key]_[Value]
#'   Eprime.Level: [Level]
#'   Running: [Key]
#'   Cycle: [Cycle]
#'   Sample: [Sample]
#' }
#'
#' @keywords internal
tidy <- function(x) {
  eprime_frame <- x
  level_depth <- eprime_frame[[rprime_cols$level]]
  level_label <- eprime_frame[[rprime_cols$running]]
  level_name  <- eprime_frame[[rprime_cols$level_name]]

  # Tidy if Running field is used and level name is still NA
  has_level_label <- !is.na(level_label)
  needs_level_name <- is.na(level_name)

  if (needs_level_name & has_level_label) {
    # Store [Key]_[Value] as "Eprime.LevelName"
    level_index <- eprime_frame[[level_label]]
    level_name  <- paste0(level_label, "_", level_index)
    new_list <- structure(as.list(level_name), names = rprime_cols$level_name)

    # Remove "[Key]: [Value]" item and then remove "[Key]." from names
    eprime_frame[level_label] <- NULL
    key_dot <- paste0(level_label, "\\.")
    names(eprime_frame) <- str_replace(names(eprime_frame), key_dot, "")
    eprime_frame <- merge_lists(eprime_frame, new_list)
  }

  class(eprime_frame) <- class(x)
  eprime_frame
}


#' Convert a vector of colon-separated text lines into a list of named elements
#'
#' @details Some minor cleaning of the input is performed:
#'   \itemize{
#'     \item Lines without a colon-space separator \code{": "} are filtered out.
#'     \item Once the strings are split at the separator, white-space on the
#'           left and right sides of each half-string is omitted.}
#' @param colon_sep_xs a character vector with lines of the form
#'   \code{"key: value"}
#' @return a named list of the values in the colon-separated lines.
#'   \code{"key: value"} yields \code{list(key = "value")}
#' @export
listify <- function(colon_sep_xs) {
  colon_sep_xs <- Filter(is_row, colon_sep_xs)
  splits <- str_split_fixed(colon_sep_xs, pattern = ": ", 2)
  # Trim after splitting so "X: " lines are correctly parsed
  splits <- apply(splits, 2, str_trim)
  # apply reduces single row matrix into a vector
  splits <- if (!is.matrix(splits)) matrix(splits, ncol = 2) else splits
  structure(as.list(splits[, 2]), names = splits[, 1])
}


#' Infer level of nesting (for a log-frame) by counting tabs
#'
#' The number of tabs before the "key: value" information in a log-frame tells
#' where the frame is nested in the experiment's structure.
#' @keywords internal
count_tabs <- function(colon_sep_xs) {
  colon_sep_xs <- ifelse(length_zero(colon_sep_xs), "", colon_sep_xs)
  # Add one because there is no level 0
  level <- str_count(first(colon_sep_xs), "\\t") + 1
  structure(list(level), names = rprime_cols$level)
}

