#
# Print benchmark result
#
#' @export
print.benchmark <- function(x, digits = 2, ...) {
  # Validate arguments
  new_args <- .validate_print_benchmark_args(x, digits, ...)

  # Print
  print(new_args$x$tab, digits = new_args$digits)
}

#
# Validate arguments and return updated arguments
#
.validate_print_benchmark_args <- function(x, digits, ...) {

  if (!methods::is(x, "benchmark")) {
    stop("Ivalid object type", call. = FALSE)
  }

  assertthat::assert_that(assertthat::is.number(digits))

  list(x = x, digits = digits)
}

#
# Print curve evaluation result
#
#' @export
print.evalcurve <- function(x, data_type = "summary", ...) {
  # Validate arguments
  new_args <- .validate_print_evalcurve_args(x, data_type, ...)

  # Print
  if (new_args$data_type == "summary") {
    newdf <- new_args$x$testsum[, c("testset", "toolset", "toolname", "label")]
    names(newdf) <- c("testset", "toolset", "toolname", "score")
    print(newdf)
  } else if (new_args$data_type == "all") {
    print(new_args$x$testscores)
  } else if (new_args$data_type == "basepoints") {
    print(new_args$x$basepoints)
  } else if (new_args$data_type == "predictions") {
    print(new_args$x$predictions)
  } else if (new_args$data_type == "rawsummary") {
    print(new_args$x$testsum)
  }
}

#
# Validate arguments and return updated arguments
#
.validate_print_evalcurve_args <- function(x, data_type, ...) {

  if (!methods::is(x, "evalcurve")) {
    stop("Ivalid object type", call. = FALSE)
  }

  assertthat::assert_that(assertthat::is.string(data_type))
  idx <- pmatch(data_type, c("summary", "all", "basepoints", "predictions",
                             "rawsummary"))
  if (is.na(idx)){
    stop("Invalid data_type", call. = FALSE)
  }

  if (idx == 1) {
    data_type <- "summary"
  } else if (idx == 2) {
    data_type <- "all"
  } else if (idx == 3) {
    data_type <- "basepoints"
  } else if (idx == 4) {
    data_type <- "predictions"
  } else if (idx == 5) {
    data_type <- "rawsummary"
  }

  list(x = x, data_type = data_type)
}
