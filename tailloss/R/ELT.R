#' Event Loss Table
#'
#' Function to create an ELT object
#'
#' @param X Data frame containing at least two numeric columns. The column \code{Loss} contains the expected losses from each single occurrence of event. The column \code{Rate} contains the arrival rates of a single occurrence of event.
#' @param Rate Positive numeric vector of arrival rates
#' @param Loss Positive numeric vector of losses
#' @param ID Vector event ID.

#'  @return An object ELT, a data frame with 3 columns. The column \code{ID} contains the ID of each event.  The column \code{Rate} contains the arrival rates of a single occurrence of event. The column \code{Loss} contains the expected losses from each single occurrence of event.

#' @seealso \code{\link{data.frame}}
#' @export
#' @examples
#' rate <- c(.1, .02, .05)
#' loss <- c(2, 5, 7)
#'
#' ELT(Rate = rate, Loss = loss)
#' # Same as
#' rl <- data.frame(Rate = rate, Loss = loss)
#' ELT(rl)

ELT <- function(X = NULL, Rate = NULL, Loss = NULL, ID = NULL) {
	if (is.data.frame(X)) {
		stopifnot(c("Rate", "Loss") %in% names(X), X$Rate > 0, X$Loss > 0)
		if (is.null(X$ID))
			X$ID <- 1L:length(X$Rate)
		X <- data.frame(ID = X$ID, Rate = X$Rate, Loss = X$Loss)
		} else {
	stopifnot(is.null(X), Rate > 0, Loss > 0, length(Rate) == length(Loss))
	if (is.null(ID))
		ID <- 1L:length(Rate)
	else stopifnot(length(ID) == length(Loss))
	X <- data.frame(ID = ID, Rate = Rate, Loss = Loss)
	}
	class(X) <- c("ELT", "data.frame")
	X
}
