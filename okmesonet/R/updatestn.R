#' Update list of Oklahoma Mesonet stations
#'
#' The list of Oklahoma Mesonet stations was updated November 28, 2014. 
#' Use \code{updatestn} to manually update the list.
#'
#' @export
#' @name updatestn
#' @seealso \code{\link{okstations}}
#' 
#' @examples
#' \dontrun{
#' ## Update Oklahoma Mesonet station list
#' okstations <- updatestn()
#' }

updatestn <- function() {
  ## Update okstations by calling downloadstn()
  ##
  ## Arguments: none
  ## Returns: updated okstations objects
  downloadstn()
}