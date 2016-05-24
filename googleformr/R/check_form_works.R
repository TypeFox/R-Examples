#' Check If The Form Is Sending Correctly
#'
#' Check the HTTP status message is below 400, meaning that the transmission worked.
#'
#' @param r response from HTTP POST request
#' @export
#' @examples
#' \dontrun{
#' why_R_u_opensource("R is Cool!") %>% check_form_works()
#' # "All good."
#' }

check_form_works <- function(r) {
   if (httr::status_code(r) < 400) {
     message("All good.")
   } else {
     message("There seems to be an error message: ",httr::message_for_status(r$status))
   }
}