#' Is the input a data.table?
#'
#' Checks to see if the input is a data.table.
#'
#' @param x Input to check.
#' @param .xname Not intended to be used directly.
#' @param severity How severe should the consequences of the assertion be?  
#' Either \code{"stop"}, \code{"warning"}, \code{"message"}, or \code{"none"}.
#' @return \code{is_data.table} wraps \code{is.data.table}, providing more 
#' information on failure. The \code{assert_*} functions return nothing but
#' throws an error if the corresponding \code{is_*} function returns 
#' \code{FALSE}.
#' @seealso \code{\link[data.table]{is.data.table}}.
#' @examples
#' if(requireNamespace("data.table"))
#' {
#'   assert_is_data.table(data.table::data.table(x = 1:5))
#'   #These examples should fail.
#'   assertive.base::dont_stop(assert_is_data.table(list(1,2,3)))
#' } else
#' {
#'   message("This example requires the data.table package to be installed.")
#' }
#' @importFrom assertive.base is2
#' @export
is_data.table <- function(x, .xname = get_name_in_parent(x))
{
  if(!(ok <- is_data.frame(x, .xname)))
  {
    return(ok)
  }
  is2(x, "data.table", .xname)
}
