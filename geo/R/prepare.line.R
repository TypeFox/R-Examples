#' Prepare a line ?
#' 
#' Prepare a line for some action ?
#' 
#' 
#' @param x Line vector ?
#' @return Returns a list with components: \item{lx1}{lx1 ?} \item{lx2}{lx2 ?}
#' \item{nlx}{nlx ?}
#' @note Needs elaboration.
#' @seealso Called by \code{\link{cut_multipoly}} and \code{\link{findline}}.
#' @keywords manip
#' @export prepare.line
prepare.line <-
function(x)
{
	n <- length(x)
	x1 <- x[2:n]
	x2 <- x[1:(n - 1)]
	ind <- c(1:(n - 1))
	ind1 <- ind[!is.na(x1) & is.na(x2)]
	ind2 <- ind[is.na(x1) & !is.na(x2)]
	if(length(ind1) > 0)
		lx1 <- ind1 + 1
	if(length(ind2) > 0)
		lx2 <- ind2
	if(length(ind1) == 0)
		lx1 <- 1
	if(length(ind2) == 0)
		lx2 <- n
	if(!is.na(x[1]))
		lx1 <- unique(c(1, lx1))
	if(!is.na(x[n]))
		lx2 <- unique(c(lx2, n))
	nlx <- length(lx1)
	return(list(lx1 = lx1, lx2 = lx2, nlx = nlx))
}

