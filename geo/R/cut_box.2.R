#' Cut box type 2???
#' 
#' Cut box type 2???.
#' 
#' 
#' @param x Longtude ?
#' @param y Latitude ?
#' @param xb Limits of longitude?
#' @param yb Limits of latitude?
#' @return List with components: \item{x, y}{Longitude and latitude of ???}
#' \item{x1, y1}{Some other Longitude and latitude of ???}
#' @note Internal to \code{gridaxes.Lambert}, needs elaboration.
#' @seealso \code{\link{gridaxes.Lambert}}
#' @keywords manip
#' @export cut_box.2
cut_box.2 <-
function(x, y, xb, yb)
{
	ind <- c(1:length(x))
	inds <- ind[is.na(x)]
	ind <- inds
	xx <- matrix(x,  , ind[1], byrow = T)
	yy <- matrix(y,  , ind[1], byrow = T)
	ind <- cut(x, xb,labels=FALSE) # labels=FALSE R ver.
	ind1 <- ind[2:length(ind)]
	ind <- ind[1:(length(ind) - 1)]
	ii <- c(1:length(ind))
	i <- ifelse(is.na(ind) & !is.na(ind1), ii, NA)
	i <- i[!is.na(i)]
	i1 <- ifelse(!is.na(ind) & is.na(ind1), ii, NA)
	i1 <- i1[!is.na(i1)]
	i2 <- c(1:length(ind))
	i2 <- i2[is.na(ind)]
	x2 <- y2 <- x3 <- y3 <- matrix(NA, nrow(xx), 3)
	x2[, 1] <- x[i]
	x2[, 2] <- x[i + 1]
	y2[, 1] <- y[i]
	y2[, 2] <- y[i + 1]
	x3[, 1] <- x[i1]
	x3[, 2] <- x[i1 + 1]
	y3[, 1] <- y[i1]
	y3[, 2] <- y[i1 + 1]
	t1 <- (xb[1] - x2[, 1])/(x2[, 2] - x2[, 1])
	t2 <- (xb[2] - x3[, 1])/(x3[, 2] - x3[, 1])
	x1 <- y1 <- matrix(NA, nrow(xx), 3)
	x1[, 1] <- x2[, 1] + t1 * (x2[, 2] - x2[, 1])
	x1[, 2] <- x3[, 1] + t2 * (x3[, 2] - x3[, 1])
	y1[, 1] <- y2[, 1] + t1 * (y2[, 2] - y2[, 1])
	y1[, 2] <- y3[, 1] + t2 * (y3[, 2] - y3[, 1])
	y[i2] <- NA
	x[i2] <- NA
	x[i] <- x1[, 1]
	y[i] <- y1[, 1]
	x[i1 + 1] <- x1[, 2]
	y[i1 + 1] <- y1[, 2]
	ind <- cut(y, c(-999999, yb, 999999),labels=FALSE) # labels=FALSE R ver.
	ind1 <- c(1:length(ind))
	ind1 <- ind1[ind != 2]
	x[ind1] <- NA
	y[ind1] <- NA
	return(list(x = x, y = y, x1 = c(x1[, 1]), y1 = c(y1[, 1])))
}

