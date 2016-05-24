#' @title Is single
#' 
#' @description Tests if an object is single (i.e. of length 1)
#' 
#' @param x an R object
#' @seealso \code{\link{is_single_number}}, \code{\link{is_single_string}},
#' \code{\link{is_single_logical}}
#' @export
#' @examples
#' is_single("hoskdflksfd")  # TRUE
#' is_single("1.0")  # TRUE
#' is_single(1:5)  # FALSE
#' is_single(matrix(runif(4), 2, 2))  # FALSE
is_single <- function(x) {
  (length(x) == 1)
}


#' @title Is single string
#' 
#' @description Tests if an object is a single string
#' 
#' @param x an R object
#' @seealso \code{\link{is_single}}
#' @export
#' @examples
#' is_single_string(1.0)  # FALSE
#' is_single_string("hoskdflksfd")  # TRUE
#' is_single_string(c("1.0", "sd"))  # FALSE
is_single_string <- function(x) {
  if (is_single(x)) {
    is_string(x)
  } else FALSE
}


#' @title Is single number
#' 
#' @description Tests if an object is a single number
#' 
#' @param x an R object
#' @seealso \code{\link{is_single}}
#' @export
#' @examples
#' is_single_number(1.0)  # TRUE
#' is_single_number("hoskdflksfd")  # FALSE
#' is_single_number("1.0")  # FALSE
#' is_single_number(1:5)  # FALSE
is_single_number <- function(x) {
  if (is_single(x)) {
    is.numeric(x)
  } else FALSE
}


#' @title Is single positive number
#' 
#' @description Tests if an object is a single positive number
#' 
#' @param x an R object
#' @seealso \code{\link{is_single}}, \code{\link{is_single_negative}}
#' @export
#' @examples
#' is_single_positive(1.0)  # TRUE
#' is_single_positive(c(1.0,2))  # FALSE
#' is_single_positive(-1.0)  # FALSE
#' is_single_positive(0)  # FALSE
#' is_single_positive(NA)  # FALSE
is_single_positive <- function(x) {
  if (is_single(x)) {
    is_positive(x)
  } else FALSE
}


#' @title Is single negative number
#' 
#' @description Tests if an object is a single negative number
#' 
#' @param x an R object
#' @seealso \code{\link{is_single}}, \code{\link{is_single_positive}}
#' @export
#' @examples
#' is_single_negative(1.0)  # FALSE
#' is_single_negative(-1.0)  # TRUE
#' is_single_negative(c(-1.0,-2))  # FALSE
#' is_single_negative(0)  # FALSE
#' is_single_negative(NA)  # FALSE
is_single_negative <- function(x) {
  if (is_single(x)) {
    is_negative(x)
  } else FALSE
}


#' @title Is single decimal
#' 
#' @description Tests if an object is a single decimal number
#' 
#' @param x an R object
#' @seealso \code{\link{is_single}}
#' @export
#' @examples
#' is_single_decimal(0.01)  # TRUE
#' is_single_decimal(-3/4)  # TRUE
#' is_single_decimal("hoskdflksfd")  # FALSE
#' is_single_decimal("1.0")  # FALSE
#' is_single_decimal(1:5)  # FALSE
is_single_decimal <- function(x) {
  if (is_single(x)) {
    is_decimal(x)
  } else FALSE
}


#' @title Is single positive decimal
#' 
#' @description Tests if an object is a single positive decimal
#' 
#' @param x an R object
#' @seealso \code{\link{is_single}}, \code{\link{is_single_positive}},
#' \code{\link{is_single_negative_decimal}}
#' @export
#' @examples
#' is_single_positive_decimal(0.01)  # TRUE
#' is_single_positive_decimal(-3/4)  # FALSE
#' is_single_positive_decimal("hoskdflksfd")  # FALSE
#' is_single_positive_decimal("1.0")  # FALSE
#' is_single_positive_decimal(1:5)  # FALSE
is_single_positive_decimal <- function(x) {
  if (is_single(x)) {
    is_positive_decimal(x)
  } else FALSE
}


#' @title Is single negative decimal
#' 
#' @description Tests if an object is a single positive decimal
#' 
#' @param x an R object
#' @seealso \code{\link{is_single}}, \code{\link{is_single_negative}},
#' \code{\link{is_single_positive_decimal}}
#' @export
#' @examples
#' is_single_negative_decimal(-3/4)  # TRUE
#' is_single_negative_decimal(0.01)  # FALSE
#' is_single_negative_decimal("hoskdflksfd")  # FALSE
#' is_single_negative_decimal("1.0")  # FALSE
#' is_single_negative_decimal(1:5)  # FALSE
is_single_negative_decimal <- function(x) {
  if (is_single(x)) {
    is_negative_decimal(x)
  } else FALSE
}


#' @title Is single positive integer
#' 
#' @description Tests if an object is a single positive integer
#' 
#' @param x an R object
#' @seealso \code{\link{is_single}}, \code{\link{is_single_negative_integer}}
#' @export
#' @examples
#' is_single_positive_integer(1.0)  # TRUE
#' is_single_positive_integer(c(1.0,2))  # FALSE
#' is_single_positive_integer(-1.0)  # FALSE
#' is_single_positive_integer(0)  # FALSE
#' is_single_positive_integer(NA)  # FALSE
is_single_positive_integer <- function(x) {
  if (is_single(x)) {
    is_positive_integer(x)
  } else FALSE
}


#' @title Is single negative integer
#' 
#' @description Tests if an object is a single negative integer
#' 
#' @param x an R object
#' @seealso \code{\link{is_single}}, \code{\link{is_single_positive_integer}}
#' @export
#' @examples
#' is_single_negative_integer(-1.0)  # TRUE
#' is_single_negative_integer(1.0)  # FALSE
#' is_single_negative_integer(c(1.0,2))  # FALSE
#' is_single_negative_integer(0)  # FALSE
#' is_single_negative_integer(NA)  # FALSE
is_single_negative_integer <- function(x) {
  if (is_single(x)) {
    is_negative_integer(x)
  } else FALSE
}


#' @title Is single odd
#' 
#' @description Tests if an object is a single odd number
#' 
#' @param x an R object
#' @seealso \code{\link{is_single}}, \code{\link{is_single_even}}
#' @export
#' @examples
#' is_single_odd(1.0)  # TRUE
#' is_single_odd(2)  # FALSE
#' is_single_odd(c(1.0,2))  # FALSE
#' is_single_odd(2)  # FALSE
#' is_single_odd(0)  # FALSE
#' is_single_odd(NA)  # FALSE
is_single_odd <- function(x) {
  if (is_single(x)) {
    is_odd(x)
  } else FALSE
}


#' @title Is single even
#' 
#' @description Tests if an object is a single even number
#' 
#' @param x an R object
#' @seealso \code{\link{is_single}}, \code{\link{is_single_odd}}
#' @export
#' @examples
#' is_single_even(2)  # TRUE
#' is_single_even(5)  # FALSE
#' is_single_even(c(1.0,2))  # FALSE
#' is_single_even(-1.0)  # FALSE
#' is_single_even(0)  # TRUE
#' is_single_even(NA)  # FALSE
is_single_even <- function(x) {
  if (is_single(x)) {
    is_even(x)
  } else FALSE
}


#' @title Is single logical
#' 
#' @description Tests if an object is a single logical
#' 
#' @param x an R object
#' @seealso \code{\link{is_single}}, \code{\link{is_single_true}},
#' \code{\link{is_single_false}}
#' @export
#' @examples
#' is_single_logical(TRUE)  # TRUE
#' is_single_logical(FALSE)  # TRUE
#' is_single_logical(c(TRUE, FALSE))  # FALSE
#' is_single_logical(-1.0)  # FALSE
#' is_single_logical(0)  # FALSE
#' is_single_logical(NA)  # FALSE
is_single_logical <- function(x) {
  if (is_single(x)) {
    if (is.na(x)) {
      FALSE
    } else {
      is.logical(x)      
    }
  } else FALSE
}


#' @title Is single true
#' 
#' @description Tests if an object is a single TRUE
#' 
#' @param x an R object
#' @seealso \code{\link{is_single}}, \code{\link{is_single_false}}
#' @export
#' @examples
#' is_single_true(TRUE)  # TRUE
#' is_single_true(FALSE)  # FALSE
#' is_single_true(c(TRUE, FALSE))  # FALSE
#' is_single_true(-1.0)  # FALSE
#' is_single_true(0)  # FALSE
#' is_single_true(NA)  # FALSE
is_single_true <- function(x) {
  if (is_single(x)) {
    is_TRUE(x)
  } else FALSE
}


#' @title Is single false
#' 
#' @description Tests if an object is a single FALSE
#' 
#' @param x an R object
#' @seealso \code{\link{is_single}}, \code{\link{is_single_true}}
#' @export
#' @examples
#' is_single_false(FALSE)  # TRUE
#' is_single_false(TRUE)  # FALSE
#' is_single_false(c(TRUE, FALSE))  # FALSE
#' is_single_false(-1.0)  # FALSE
#' is_single_false(0)  # FALSE
#' is_single_false(NA)  # FALSE
is_single_false <- function(x) {
  if (is_single(x)) {
    is_FALSE(x)
  } else FALSE
}
