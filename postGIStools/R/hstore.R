#' Create a empty hstore
#'
#' This function creates an empty list of lists, which can be appended to a
#' data frame as a hstore column.
#'
#' @param nr Number of records.
#' @return A empty hstore (list of lists) of length \code{nr}.
#'
#' @examples
#' contacts <- data.frame(name = c("Anne", "Bert", "Chris"))
#' contacts$phone <- new_hstore(3)
#' contacts$phone %->% "home" <- c("555-123-4567", "555-923-9134", "555-276-1123")
#' contacts$phone[2] %->% "cell" <- "555-889-9134"
#' str(contacts)
#'
#' @seealso \code{\link{\%->\%}} to read or edit the resulting data structure.
#' @export
new_hstore <- function(nr) {
    if (!is.numeric(nr) || length(nr) > 1) stop("nr must be a single number")
    rep(list(list()), nr)
}

#' Extract or replace hstore values by key
#'
#' Operator to get or set values corresponding to a given key for all records
#' in a hstore.
#'
#' Based on the hstore "\code{->}" operator in PostgreSQL, the \code{\%->\%}
#' operator returns values associated with a given key for all records in the
#' hstore. The assignment version of the operator i.e.
#' \code{hstore \%->\% key <- value} either creates a new key-value pair or,
#' if they key exists, update the associated value. It can also delete a key
#' by assigning its value to NULL.
#'
#' Note that to subset the records in \code{hstore} to which the operator applies,
#' you must use single brackets so that the result remains a list of lists. See
#' below for usage examples.
#'
#' @param hstore A hstore (i.e. list of lists).
#' @param key Character string corresponding to a key in \code{hstore}.
#' @return For the extract version, a vector of the same length as \code{hstore},
#'  containing the value correponding to \code{key} for each record
#'  (or \code{NA} if none). For the replace version, the modified hstore.
#'
#' @examples
#' contacts <- data.frame(name = c("Anne", "Bert", "Chris"))
#' contacts$phone <- new_hstore(3)
#' contacts$phone %->% "home" <- c("555-123-4567", "555-923-9134", "555-276-1123")
#' contacts$phone[2:3] %->% "home"
#'
#' contacts$phone[2] %->% "home" <- NULL
#' contacts$phone %->% "home"
#' contacts$phone[2] %->% "cell" <- "555-889-9134"
#' contacts$phone %->% "cell"
#'
#' @seealso \code{\link{new_hstore}} to create a empty hstore.
#' @export
#' @rdname extract-hstore
`%->%` <- function(hstore, key) {
    if (!all(vapply(hstore, is.list, FALSE))) {
        stop(paste(deparse(substitute(hstore)), "is not a valid hstore"))
    }
    test_single_str(key)
    lst <- lapply(hstore, `[[`, key)
    lst[vapply(lst, is.null, FALSE)] <- NA
    unlist(lst)
}

#' @param value Vector of values of the same length as \code{hstore}.
#' @usage hstore \%->\% key <- value
#' @export
#' @rdname extract-hstore
`%->%<-` <- function(hstore, key, value) {
    if (!all(vapply(hstore, is.list, FALSE))) {
        stop(paste(deparse(substitute(hstore)), "is not a valid hstore"))
    }
    test_single_str(key)
    if (!is.null(value) & length(hstore) != length(value)) {
        stop("length of value does not match number of records in hstore")
    }
    for (i in 1:length(hstore)) hstore[[i]][[key]] <- value[i]
    hstore
}



