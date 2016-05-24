#' Get the generic length of an object 
#' 
#' This function gets the generic length of an object.
#'
#' @section Usage: 
#' anylength(data)
#'
#' @section Details:
#' This function consolidates size dimensions for one and two dimensional 
#' data structures. The idea is that many operations require knowing 
#' either how long a vector is or how many rows are in a matrix.
#' So rather than switching between \code{length} and \code{nrow},
#' \code{anylength} provides the appropriate polymorphism to return the 
#' proper value. 
#'
#' When working with libraries, it is easy to forget the return type of a
#' function, particularly when there are a lot of switches between vectors,
#' matrices, and other data structures. This function along with its
#' \code{\link{anynames}} counterpart provides a single interface for
#' accessing this information across objects
#' 
#' The core assumption is that in most cases length is semantically synonymous 
#' with \code{nrow} such that the number of columns in two-dimensional 
#' structures is
#' less consequential than the number of rows. This is particularly true of
#' time-based objects, such as zoo or xts where the number of observations is
#' equal to the number of rows in the structure.
#'
#' When working with functions that are polymorphic, \code{lambda.r} function
#' clauses with guard conditions on the length of the input data structure
#' can use \code{anylength} instead of using \code{length} or \code{nrow},
#' which preserves polymorphism and reduces the number of function clauses 
#' necessary. For example, instead of one clause to check \code{length}
#' and another to check \code{nrow}, \code{anylength} can test for both
#' situations in a single clause.
#'
#' \preformatted{
#' slice(x, expression) \%::\% a : logical : list
#'
#' slice(x, expression) \%when\% \{ length(expression) == length(x) \}
#'
#' slice(x, expression) \%::\% a : logical : 
#'
#' slice(x, expression) \%when\% \{ length(expression) == nrow(x) \}}
#'
#' These two clauses can be replaced with
#'
#' \preformatted{
#' slice(x, expression) \%::\% a : logical : .
#'
#' slice(x, expression) \%when\% \{ length(expression) == anylength(x) \}}
#'
#' Another use of \code{anylength} is when working with \code{sapply}.
#' The output value is governed by the result of the higher-order
#' function, so it is difficult to know a priori whether the result
#' will be a vector or a matrix. With \code{anylength} it doesn't 
#' matter since the same function is used in either case.
#'
#' @name anylength 
#' @param data Any indexable data structure
#' @return The conceptual 'length' of a data structure.
#'
#' @author Brian Lee Yung Rowe
#'
#' @examples
#' # Get the rows of the matrix
#' anylength(matrix(c(1,2,3,4,5,6), ncol=2))
#'
#' # Get the length of the vector
#' anylength(c(1,2,3,4,5))
#
anylength(data) %when% { ! is.null(nrow(data)) } %as% nrow(data)
anylength(data) %as% length(data)

#' Get the names of a data structure. This attempts to 
#' create some polymorphism around lists, vectors, and data.frames
#' 
#' This function gets the useful names of a data structure. This attempts to 
#' create some polymorphism around lists, vectors, and data.frames.
#'
#' @section Usage:
#' anynames(data)
#'
#' @section Details:
#' Depending on the type of structure utilized in code, one needs to 
#' call either names or \code{colnames} to get information related to 
#' the data sets within the structure. The use of two separate functions 
#' can cause errors and slows development time as data structures 
#' passed from intermediate functions may change over time,
#' resulting in a broken interface.
#'
#' By providing a thin layer over underlying accessors, 
#' this function attempts to expedite development and add a 
#' bit of polymorphism to the semantics of names.
#' The explicit assumption is that data sets in 
#' two dimensional structures are organized by column, 
#' as this is compatible with time-series objects such as
#' zoo and xts. 
#'
#' @name anynames
#' @aliases anynames<-
#' @param data Any indexable data structure
#' @return Returns the names for a data structure.
#'
#' @author Brian Lee Yung Rowe
#' @keywords attribute
#'
#' @examples
#' m <- matrix(c(1,2,3,4,5,6), ncol=2)
#' anynames(m) <- c('d','e')
#' anynames(m)
#'
#' v <- c(a=1,b=2,c=3,d=4,e=5)
#' anynames(v)
#'
#' l <- list(a=1,b=2,c=3,d=4,e=5)
#' anynames(l)
#' 
#' df <- data.frame(a=1:10, b=1:10,c=1:10,d=1:10,e=1:10)
#' anynames(df)
anynames(data) %when% { ! is.null(names(data)) } %as% names(data)
anynames(data) %when% { ! is.null(colnames(data)) } %as% colnames(data)
anynames(data) %as% NULL

"anynames<-" <- function(data, value)
{
  if (is.null(names(data))) colnames(data) <- value
  else names(data) <- value
  invisible(data)
}

#' Show the types of a list or data.frame
#'
#' This function shows the types of the columns in a data.frame or the 
#' elements of a list.
#'
#' @section Usage:
#' anytypes(data, fn) %::% list : Function : character
#' anytypes(data, fn=class) 
#' 
#' anytypes(data, fn) %::% data.frame : Function : character
#' anytypes(data, fn=class) 
#'
#' @section Details:
#' This is a convenience function to see the types associated with the
#' elements of a list or the columns of a data.frame.
#'
#' @name anytypes
#' @param data A data.frame
#' @param fn The function used to get the types. Defaults to class, although
#'    type or mode, etc. could be used
#' @return A vector containing the types of the columns of a data structure
#'
#' @author Brian Lee Yung Rowe
#' @keywords attribute
#'
#' @examples
#' x <- data.frame(ints=1:3, chars=c('a','b','c'), nums=c(.1,.2,.3))
#' anytypes(x)
#'
#' x <- list(ints=1:4, chars=c('a','b','c'), nums=c(.1,.2,.3))
#' anytypes(x)
anytypes(data, fn) %::% list : Function : character
anytypes(data, fn=class) %as% {
  ts <- sapply(data, fn)
  names(ts) <- anynames(data)
  return(ts)
}

anytypes(data, fn) %::% data.frame : Function : character
anytypes(data, fn=class) %as% {
  ts <- apply(matrix(anynames(data), ncol=1), 1, function(x) fn(data[,x]))
  names(ts) <- anynames(data)
  return(ts)
}


#' Check whether data is bad or empty
#'
#' These functions quickly test whether data within an object has 
#' bad values or if the object is defined (i.e. not null) but has no data.
#'
#' @section Usage:
#' is.empty(x) %::% a : logical
#' is.empty(x)
#'
#' is.bad(x)
#'
#' @section Details:
#' Depending on the type of an object, knowing whether an object contains a
#' valid value or not is different. These functions unify the interfaces across
#' different data types quickly indicating whether an object contains bad
#' values and also whether an object has a value set.
#'
#' For example, a data.frame may be initialized with no data. This results in
#' an object that is non-null but also unusable. Instead of checking whether
#' something is both non-null and has positive length, just check is.bad().
#'
#' If you know that an object is non-null, then you can call is.empty() which
#' is a shortcut for checking the length of an object.
#' 
#' @name is.empty
#' @aliases is.bad
#' @param x An object containing the data to test
#' @return Logical values that indicate whether the test was successful 
#' or not. For matrices and data.frames, a matrix of logical values 
#' will be returned.
#'
#' @author Brian Lee Yung Rowe
#' @keywords attribute
#'
#' @examples
#' a <- data.frame(a=NULL, b=NULL)
#' is.bad(a)
#' 
#' b <- list(a=1:3, b=NULL, c=NA, d='foo')
#' is.bad(b)
#' 
#' c <- list()
#' is.empty(c)
is.empty(x) %::% a : logical
is.empty(x) %as% { length(x) < 1 }


is.bad(x) %::% a : logical
is.bad(NULL) %as% TRUE
is.bad(EMPTY) %as% TRUE

is.bad(x) %::% list : list
is.bad(x) %as% { lapply(x, is.bad) }

is.bad(x) %::% data.frame : matrix
is.bad(x) %as% { sapply(x, is.bad) }

is.bad(x) %::% matrix : matrix
is.bad(x) %as% { apply(x,1, is.bad) }

is.bad(x) %::% a : logical
is.bad(x) %as% { is.na(x) }
