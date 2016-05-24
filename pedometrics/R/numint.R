#' Tests for data types
#' 
#' Evaluate the data type contained in an object.
#' 
#' @param x Object to be tested.
#' 
#' @return
#' \code{TRUE} or \code{FALSE} depending on whether \code{x} contains a given
#' data type.
#' 
#' @seealso \code{\link[base]{is.numeric}}, \code{\link[base]{is.integer}},
#' \code{\link[base]{is.factor}}.
#' 
#' @author Alessandro Samuel-Rosa \email{alessandrosamuelrosa@@gmail.com}
#' @aliases isNumint allNumint anyNumint allInteger anyInteger allFactor
#' anyFactor allNumeric anyNumeric uniqueClass
#' @examples
#' # Vector of integers
#' x <- 1:10
#' isNumint(x) # FALSE
#' 
#' # Vector of numeric integers
#' x <- as.numeric(x)
#' isNumint(x) # TRUE
#' 
#' # Vector of numeric values
#' x <- c(1.1, 1, 1, 1, 2)
#' isNumint(x) # FALSE
#' allNumint(x) # FALSE
#' anyNumint(x) # TRUE
#' whichNumint(x)
#' 
#' # Single numeric integer
#' isNumint(1) # TRUE
#' 
#' # Single numeric value
#' isNumint(1.1) # FALSE
# FUNCTION - NUMERIC INTEGERS ##################################################
#' @rdname numint
#' @export
isNumint <-
  function (x) {
    if (is.integer(x) || is.factor(x) || is.character(x)) return (FALSE)
    if (is.numeric(x) && length(x) > 1) {
      res <- ifelse(round(x, digits = 0) == x, TRUE, FALSE)
      res <- ifelse(length(unique(res)) == 1, TRUE, FALSE)
    } else {
      res <- ifelse(round(x, digits = 0) == x, TRUE, FALSE)
    }
    return (res)
  }
#' @rdname numint
#' @export
allNumint <-
  function (x) {
    res <- sapply(x, pedometrics::isNumint)
    res <- all(res == TRUE)
    return (res)
  }
#' @rdname numint
#' @export
anyNumint <-
  function (x) {
    res <- sapply(x, pedometrics::isNumint)
    res <- any(res == TRUE)
    return (res)
  }
#' @rdname numint
#' @export
whichNumint <-
  function (x) {
    res <- sapply(x, pedometrics::isNumint)
    res <- which(res == TRUE)
    return (res)
  }
# FUNCTION - INTEGERS ##########################################################
#' @rdname numint
#' @export
allInteger <-
  function (x) {
    res <- sapply(x, is.integer)
    res <- all(res == TRUE)
    return (res)
  }
#' @rdname numint
#' @export
anyInteger <-
  function (x) {
    res <- sapply(x, is.integer)
    res <- any(res == TRUE)
    return (res)
  }
#' @rdname numint
#' @export
whichInteger <-
  function (x) {
    res <- sapply(x, is.integer)
    res <- which(res == TRUE)
    return (res)
  }
# FUNCTION - FACTORS ###########################################################
#' @rdname numint
#' @export
allFactor <-
  function (x) {
    res <- sapply(x, is.factor)
    res <- all(res == TRUE)
    return (res)
  }
#' @rdname numint
#' @export
anyFactor <-
  function (x) {
    res <- sapply(x, is.factor)
    res <- any(res == TRUE)
    return (res)
  }
#' @rdname numint
#' @export
whichFactor <-
  function (x) {
    res <- sapply(x, is.factor)
    res <- which(res == TRUE)
    return (res)
  }
# FUNCTION - NUMERIC ###########################################################
#' @rdname numint
#' @export
allNumeric <-
  function (x) {
    res <- sapply(x, is.numeric)
    res <- all(res == TRUE)
    return (res)
  }
#' @rdname numint
#' @export
anyNumeric <-
  function (x) {
    res <- sapply(x, is.numeric)
    res <- any(res == TRUE)
    return (res)
  }
#' @rdname numint
#' @export
whichNumeric <-
  function (x) {
    res <- sapply(x, is.numeric)
    res <- which(res == TRUE)
    return (res)
  }
# FUNCTION - IS ONE TYPE #######################################################
#' @rdname numint
#' @export
uniqueClass <-
  function (x) {
    res <- sapply(x, class)
    res <- length(unique(res))
    res <- ifelse(res == 1, TRUE, FALSE)
    return (res)
  }
