#' Placement and Searches for Carcasses at Altamont
#' 
#' An event-level dataset containing information pertinent to the ACME model, 
#' from Warren-Hicks et al. (2012), "Improving Methods for Estimating 
#' Fatality of Birds and Bats at Wind Energy Facilities."
#' @docType data
#' @format A data frame with 3984 observations and 6 variables:
#' \describe{
#'  \item{Date}{\code{character} - Data of the event}
#'  \item{Time}{\code{character} - Time of the event}
#'  \item{ID}{\code{character} - ID of the carcass}
#'  \item{Species}{\code{character} - 4-letter AOU species code}
#'  \item{Event}{\code{character} - Type of event - Check, Place, or Search}
#'  \item{Found}{\code{logical} - TRUE if carcass was discovered}
#' }
#' 
#' @source Warren-Hicks, W., Newman, J., Wolpert, R. L., Karas, B., and Tran, L. 
#' (2012), Improving Methods for Estimating Fatality of Birds and Bats at Wind 
#' Energy Facilities, California Wind Energy Association publication 
#' CEC-500-2012-086.
"altamont"
