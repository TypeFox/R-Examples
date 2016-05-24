#' Manipulate grid and z-values for contouring.
#' 
#' Manipulate grid and z-values for contouring.
#' 
#' 
#' @param grd grid
#' @param z z-value
#' @param mask mask
#' @param set set
#' @param col.names column names, defaults to \code{lat} and \code{lon}
#' @return Returns a list of: \item{grd}{grid} \item{z}{value over grid}
#' @note Used in \code{geocontour}-functions.
#' @keywords manip
#' @export Set.grd.and.z
Set.grd.and.z <-
function(grd, z, mask, set = NA, col.names = c("lon", "lat"))
{
	# z is a name of a column in the dataframe grd.
	if(is.data.frame(grd) && is.character(z)) z <- grd[, z]
	if(is.data.frame(grd) && nrow(grd) == length(z)) {
		i1 <- match(col.names[1], names(grd))
		i2 <- match(col.names[2], names(grd))
		xgr <- sort(unique(grd[, i1]))
		ygr <- sort(unique(grd[, i2]))
		xgr.1 <- c(matrix(xgr, length(xgr), length(ygr)))
		ygr.1 <- c(t(matrix(ygr, length(ygr), length(xgr))))
		xgr.data <- data.frame(x = xgr.1, y = ygr.1)
		names(xgr.data) <- col.names
		xgr.data$z <- rep(set, nrow(xgr.data))
		index <- paste(xgr.data[, 1], xgr.data[, 2], sep = "-")
		index1 <- paste(grd[, i1], grd[, i2], sep = "-")
		j <- match(index1, index)
		xgr.data$z[j] <- z
		grd1 <- list(xgr, ygr)
		names(grd1) <- col.names
		return(list(grd = grd1, z = xgr.data$z))
	}
	# grd is a list like returned by pointkriging
	if(is.list(grd) && !is.data.frame(grd)) {
		i1 <- match(col.names[1], names(grd))
		i2 <- match(col.names[2], names(grd))
		xgr <- grd[[i1]]
		ygr <- grd[[i2]]
		if(length(xgr) * length(ygr) != length(z)) {
			cat("Incorrect length on z")
			return(invisible())
		}
		xgr.1 <- c(matrix(xgr, length(xgr), length(ygr)))
		ygr.1 <- c(t(matrix(ygr, length(ygr), length(xgr))))
		xgr.data <- data.frame(x = xgr.1, y = ygr.1)
		names(xgr.data) <- col.names
		xgr.data$z <- z
		grd1 <- list(xgr, ygr)
		names(grd1) <- col.names
		return(list(grd = grd1, z = xgr.data$z))
	}
}

