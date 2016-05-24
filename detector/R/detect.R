#' Detect if a data object contains PII.
#'
#' @param .x A data object.
#' @return A logical value indicating if that data object contains PII.
#' @examples
#' # atomic vectors
#' detect(letters)
#' detect(1:10)
#' detect(as.Date("2014-01-01"))
#'
#' # data.frames
#' detect(mtcars)
#' @export
detect <- function(.x) {
  UseMethod("detect")
}

#' @describeIn detect Method for default vectors.
#' @export
detect.default <- function(.x) {
  return(detect(as.character(.x)))
}

#' @describeIn detect Method for character vectors.
#' @export
detect.character <- function(.x) {
  temp <- data.frame()
  temp[1, "has_email_addresses"] <- has_email_addresses(.x)
  temp[, "has_phone_numbers"] <- has_phone_numbers(.x)
  temp[, "has_national_identification_numbers"] <- has_national_identification_numbers(.x)
  return(temp)
}

detect_column <- function(.nm, .df) {
  return(cbind(data.frame(column_name = .nm, stringsAsFactors = FALSE), detect(.df[, .nm])))
}

#' @describeIn detect Method for data.frames.
#' @export
detect.data.frame <- function(.x) {
  temp <- lapply(colnames(.x), detect_column, .df = .x)
  return(do.call(rbind, temp))

}
