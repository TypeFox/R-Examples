# Copyright 2009-2015 Meik Michalke <meik.michalke@hhu.de>
#
# This file is part of the R package klausuR.
#
# klausuR is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# klausuR is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with klausuR.  If not, see <http://www.gnu.org/licenses/>.


#' Sort method for S4 objects of class klausuR
#'
#' Returns the given object, with global results and anonymized results sorted by the given variable.
#'
#' @param x An object of class \code{klausuR}
#' @param decreasing Logical, whether sorting should be sone increasing or decreasing.
#' @param ... Additional arguments.
#' @author m.eik michalke \email{meik.michalke@@uni-duesseldorf.de}
#' @seealso \code{\link[klausuR:klausur]{klausur}}
#' @keywords methods
#' @examples
#' \dontrun{
#' klsr.obj <- klausur(data.obj)
#' sort(klsr.obj, sort.by="Points")
#' }
#' @export
#' @docType methods
#' @rdname sort-methods
setGeneric("sort", function(x, decreasing=FALSE, ...) standardGeneric("sort"))

#' @rdname sort-methods
#' @docType methods
#' @aliases sort,-methods sort,klausuR-method sort.klausuR,klausuR-method
#' @param sort.by An optional character string naming a variable to sort the results by. Defaults to \code{c()}, i.e. no re-ordering.
setMethod("sort", signature(x="klausuR"), function(x, decreasing=FALSE, sort.by=c()){
  if(length(x@results) == 0){
    return(invisible(NULL))
  } else {}

  global.results <- x@results
  anon.results <- x@anon
  # sort results
  if(length(sort.by) > 0){
    if(!sort.by %in% names(global.results)){
      stop(simpleError(paste("Can't sort by '",sort.by,"', there's no such variable!", sep="")))
    } else {}
    new.order <- order(global.results[[sort.by]], decreasing=decreasing)
    global.results <- global.results[new.order,]
    anon.results <- anon.results[new.order,]
    dimnames(global.results)[[1]] <- 1:nrow(global.results)
    dimnames(anon.results)[[1]] <- 1:nrow(anon.results)
  } else {}

  x@results <- global.results
  x@anon <- anon.results
  return(x)
})
