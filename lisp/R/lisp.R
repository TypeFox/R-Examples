##' The empty list
##' @export
nil <- list()

##' Whether a list is empty.
##' @param list the list to test
##' @return Whether the list is empty
##' @export
is.nil <- function(list)
  length(list) == 0 || is.null(car(list))

##' First element of a list
##' @param list the list to first
##' @return The first element
##' @export
car <- function(list)
  list[[1]]

##' Return elements after the first of a list.
##' @param list the list from which to extract
##' @return The elements after the first, or \code{nil}
##' if only one
##' @export
cdr <- function(list) {
  if (is.nil(list))
    stop('CDRing a null list')
  length <- length(list)
  if (length == 1)
    nil
  else
    list[2:length]
}

##' Composite \code{car}/\code{cdr}
##' @param list the list from which to extract
##' @return The extracted elements
##' @export
caar <- function(list) {
  car(car(list))
}

##' Composite \code{car}/\code{cdr}
##' @param list the list from which to extract
##' @return The extracted elements
##' @export
cadr <- function(list) {
  car(cdr(list))
}

##' Composite \code{car}/\code{cdr}
##' @param list the list from which to extract
##' @return The extracted elements
##' @export
cddr <- function(list) {
  cdr(cdr(list))
}

##' Composite \code{car}/\code{cdr}
##' @param list the list from which to extract
##' @return The extracted elements
##' @export
caddr <- function(list) {
  car(cddr(list))
}

##' Composite \code{car}/\code{cdr}
##' @param list the list from which to extract
##' @return The extracted elements
##' @export
cadar <- function(list) {
  cadr(car(list))
}

##' Composite \code{car}/\code{cdr}
##' @param list the list from which to extract
##' @return The extracted elements
##' @export
cdddr <- function(list) {
  cddr(cdr(list))
}

##' Is a number even?
##' @param a the number to test
##' @return Whether the number is even
##' @export
is.even <- function(a) {
  a %% 2 == 0
}

##' Is a number odd?
##' @param a the number to test
##' @return Whether the number is odd
##' @export
is.odd <- function(a) {
  Negate(is.even)(a)
}

##' Zip \emph{n} lists together into tuplets of
##' length \emph{n}.
##' @param zipper the zipping function
##' @param \dots the lists to be zipped
##' @return A list of tuplets
##' @export
zip <- function(zipper, ...) {
  if (length(list(...)) == 1)
    Map(zipper, c(...))
  else {
    m <- mapply(zipper, ...)
    split(m, col(m))
  }
}

##' Zip using \code{c}.
##' @param \dots the lists to be zipped
##' @return A list of tuplets
##' @seealso \code{\link{zip}}
##' @export
zip.c <- function(...) {
  zip(c, ...)
}

##' Zip using \code{list}.
##' @param \dots the lists to be zipped
##' @return A list of tuplets
##' @seealso \code{\link{zip}}
##' @export
zip.list <- function(...) {
  zip(list, ...)
}

##' Do a less efficient zip whilst preserving names.
##' @param ... lists to be zipped whilst preserving names
##' @export
zip.with.names <- function(...) {
  if (length(list(...)) == 1)
    Map(list, c(...))
  else {
    lists <- list(...)
    iter <- function(zipped, lists, names) {
      if (Reduce(`||`, Map(is.nil, lists), FALSE))
        zipped
      else {
        this.list <- Map(car, lists)
        these.names <- Map(car, names)
        names(this.list) <- these.names
        iter(append(zipped, list(this.list)),
             Map(cdr, lists),
             Map(cdr, names))
      }
    }
    structure(iter(NULL, lists, Map(names, lists)),
              names=names(lists))
  }
}

##' Combine a list into pairwise elements; lists should
##' be of the same length. In case of odd numbers of members,
##' the last will be removed.
##' @param list the list to be pairwise decomposed
##' @return A list of pairwise elements
##' @export
pairwise <- function(list) {
  length <- length(list)
  if (length < 2)
    return(nil)
  length <- ifelse(is.odd(length),
                   length - 1,
                   length)
  odds <- seq(1, length, 2)
  evens <- seq(2, length, 2)
  zip.c(list[odds], list[evens])
}

##' Apply \code{f} to the successive elements of \code{...}.
##' @param f the function to apply, whose arity should match the
##' cardinality of \code{...}
##' @param ... lists upon which to apply \code{f} successively
##' @return NULL
##' @export
for.each <- function(f, ...) {
  args <- list(...)
  while (!do.call(any, Map(is.nil, args))) {
    applicanda <- Map(car, args)
    do.call(f, applicanda)
    args <- Map(cdr, args)
  }
}

##' Try to get the \code{cdrs}; otherwise, return \code{nil}.
##' @param ... lists to \code{cdr}
##' @return the \code{cdr} of the lists
cdrs <- function(...)
  tryCatch(Map(cdr, list(...)),
           error=function(e) nil)

##' Last element in a list.
##' @param list The list to last
##' @export
last <- function(list) car(tail(list, 1))

##' pair-fold-right from SRFI-1.
##' @param f function to apply over the list-tails
##' @param nil the default value
##' @param ... the lists whose tails fold over
##' @TODO one-list fast-path
##' @export
pair.fold.right <- function(f, nil, ...) {
  lists <- list(...)
  iter <- function(lists) {
    cdrs <- do.call(cdrs, lists)
    if (is.nil(cdrs))
      nil
    else
      do.call(f, append(lists, list(iter(cdrs))))
  }
  iter(lists)
}
