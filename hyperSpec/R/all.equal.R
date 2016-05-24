.all.equal <- function (target, current, ..., check.attributes = FALSE, check.names = FALSE, 
                        check.column.order = FALSE, check.label = FALSE){
  validObject (target)
  validObject (current)
  
  result <- character (0)
  
  cmp <- all.equal (target = target@wavelength, current = current@wavelength, ...,
                    check.attributes = check.attributes, check.names = check.names)
  if (! isTRUE (cmp)) result <- c("@wavelength:", cmp)
  
  if (check.column.order)
    cmp <- all.equal (target = target@data, current = current@data, ...,
                      check.attributes = check.attributes)
  else
    cmp <- all.equal (target  = target@data  [order (colnames ( target@data))],
                      current = current@data [order (colnames (current@data))],
                      ...,
                      check.attributes = check.attributes, check.names = check.names)
  if (! isTRUE (cmp)) result <- c (result, "@data:", cmp)
  
  if (check.label){
    cmp <- all.equal (target  = target@label  [order (names (target@label))],
                      current = current@label [order (names (current@label))],
                      ...,
                      check.attributes = check.attributes, check.names = check.names)
    if (! isTRUE (cmp)) result <- c (result, "@label:", cmp)
  }
  
  if (length (result) == 0)
    TRUE
  else
    result
}


.test (.all.equal) <- function () {
  checkTrue (all.equal (flu, --flu))
  checkTrue (all.equal (flu, --flu, check.attributes = TRUE)) 
  checkTrue (all.equal (flu, --flu, check.names = TRUE))
  
  ## labels
  checkTrue (all.equal (flu, --flu, check.label = TRUE))
  ## order of labels does *not* matter
  tmp <- flu
  tmp@label <- rev (tmp@label)
  checkTrue (all.equal (flu, tmp, check.label = TRUE))
  ## whether character or expression does matter:
  tmp@label <- lapply (tmp@label, as.expression)
  checkTrue (! isTRUE (all.equal (flu, tmp, check.label = TRUE)))
  
  ## column order
  checkTrue ( all.equal (flu, --flu, check.column.order = TRUE))
  checkTrue (! isTRUE (all.equal (flu, flu [, rev (colnames (flu))], check.column.order = TRUE)))
  checkTrue (          all.equal (flu, flu [, rev (colnames (flu))], check.column.order = FALSE))
}


##' @aliases all.equal  all.equal,hyperSpec,hyperSpec-method
##' @rdname Comparison
##' @param target,current two \code{hyperSpec} objects that are tested for
##'   equality
##' @param ... handed to \code{\link[base]{all.equal}} when testing the slots of the
##'   \code{hyperSpec} objects
##' @param check.column.order If two objects have the same data, but the order
##'   of the columns (determined by the names) differs, should they be regarded
##'   as different?
##' @param check.label Should the slot \code{label} be checked? \cr If the
##'   labels differ only in the order of their entries, they are conidered
##'   equal.
##' @param check.attributes,check.names see \code{\link[base]{all.equal}}
##' @return \code{all.equal} returns either \code{TRUE}, or a character vector describing the
##' differences. In conditions, the result must therefore be tested with
##' \code{\link[base]{isTRUE}}.
##' @seealso \code{\link[base]{all.equal}} and \code{\link[base]{isTRUE}}
##' @export
setMethod ("all.equal", signature (target = "hyperSpec", current = "hyperSpec"), .all.equal)
