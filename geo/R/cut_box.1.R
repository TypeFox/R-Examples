#' Cut box ?
#' 
#' Cut box in connection with overlaying a Lambert projected plot with a grid.
#' 
#' 
#' @param x,y Longitude and latitude of gridlines
#' @param xb,yb Longitude and latitude limits
#' @return List with compoents: \item{x}{Longitudes??} \item{y}{Latitudes??}
#' \item{ind}{Indices??}
#' @note Internal in \code{gridaxes.Lambert}, needs elaboration.
#' @seealso \code{\link{gridaxes.Lambert}}
#' @keywords manip
#' @export cut_box.1
cut_box.1 <-
function(x, y, xb, yb)
{
	ind <- c(1:length(x))
	ind <- ind[is.na(x)]
	x <- matrix(x,  , ind[1], byrow = T)
	y <- matrix(y,  , ind[1], byrow = T)
	n <- ind[1] - 1
	t1 <- (yb[1] - y[, 1])/(y[, n] - y[, 1])
	t2 <- (yb[2] - y[, 1])/(y[, n] - y[, 1])
	x1 <- y1 <- matrix(NA, nrow(x), 3)
	x1[, 1] <- x[, 1] + t1 * (x[, n] - x[, 1])
	x1[, 2] <- x[, 1] + t2 * (x[, n] - x[, 1])
	y1[, 1] <- y[, 1] + t1 * (y[, n] - y[, 1])
	y1[, 2] <- y[, 1] + t2 * (y[, n] - y[, 1])
	ind2 <- cut(x1[, 1], xb,labels=FALSE) # labels=FALSE R ver.
	ind <- c(1:length(ind2))
	ind2 <- ind[!is.na(ind2)]
	ind <- cut(x1[, 1], c(-9999999, xb),labels=FALSE) # labels=FALSE R ver.
	ind1 <- c(1:length(ind))
	ind1 <- ind1[!is.na(ind) & ind == 1]
	t1 <- (xb[1] - x[ind1, 1])/(x[ind1, n] - x[ind1, 1])
	x1[ind1, 1] <- x[ind1, 1] + t1 * (x[ind1, n] - x[ind1, 1])
	y1[ind1, 1] <- y[ind1, 1] + t1 * (y[ind1, n] - y[ind1, 1])
	ind1 <- c(1:length(ind))
	ind1 <- ind1[is.na(ind)]
	t1 <- (xb[2] - x[ind1, 1])/(x[ind1, n] - x[ind1, 1])
	x1[ind1, 1] <- x[ind1, 1] + t1 * (x[ind1, n] - x[ind1, 1])
	y1[ind1, 1] <- y[ind1, 1] + t1 * (y[ind1, n] - y[ind1, 1])
	ind <- cut(x1[, 2], c(-9999999, xb),labels=FALSE )# labels=FALSE R ver.
	ind1 <- c(1:length(ind))
	ind1 <- ind1[ind == 1 | is.na(ind)]
	x1[ind1,  ] <- NA
	y1[ind1,  ] <- NA
	return(list(x = t(x1), y = t(y1), ind = ind2))
}

