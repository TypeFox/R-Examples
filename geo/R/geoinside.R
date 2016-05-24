#' Finds a subset of a given set of data which is inside a given region.
#' 
#' Finds a subset of a given set of data which is inside or outside a given
#' region and returns a submatrix of those values, boolean vector of indexes or
#' index vector.
#' 
#' 
#' @param data a dataframe, should include vectors \$lat and \$lon, but will
#' except other names, see col.names.
#' @param reg The region we want to determine whether the data is inside of,
#' should include vectors with same names as data.
#' @param option Allows you to determine on what format you recieve the output:
#' <s-example> option = 1: returns the submatrix of data which is inside reg.
#' option = 2: returns the submatrix of data which is outside reg.  option = 3:
#' returns a boolean vector saying whether a given index is inside reg.  option
#' = 4: returns a boolean vector saying whether a give index is outside reg.
#' option = 5: returns a vector of indexes to data, data[return[i],] is the
#' i-th point in data inside reg.  option = 6: returns a vector of indexes to
#' data, data[return[i],] is the i-th point in data outside reg. </s-example>
#' @param col.names Default col.names = c("lat","lon"), determines the names of
#' the base vectors of the space we are viewing.  May be replaced by for
#' instance col.names = c("x","y")
#' @param na.rm If true values where \$lat or \$lon are NA are removed.
#' Default is false.
#' @param robust If true a robust search is done, if false the function runs
#' faster.  Default is true. Robust = T will not work if the regions edges
#' overlap each other, if region is sensibly defined this will not happen and
#' robust =F should be used.
#' @return Returns and output vector or matrix, see option.
#' @section Side Effects: None.
#' @seealso \code{\link{geoplot}}, \code{\link{geolocator}}.
#' @examples
#' 
#' \dontrun{   seafishing <- geoinside(fishing,island,option=2)  
#'    # Removes those datapoints from fishing where fishing took
#'    # place inside Iceland (misspells).
#' 
#'    grd <- list(lat=c(64,64,63,63),lon=c(-23,-22,-22,-23))
#'    ins.lat.64.63.lon.23.22 <- geoinside(fishing,grd,robust =T)  
#'    # Extracts those points from fishing where fishing
#'    # took place inside the given box.
#' 
#' 
#'       
#'    #######################################################
#'    # Example                                             # 
#'    #######################################################
#' 
#'    par(mfrow=c(2,1))
#'    stations<-data.frame(lat=stodvar$lat,lon=stodvar$lon)
#'    stations<-stations[!is.na(stations$lon),]            
#'    stations<-stations[!is.na(stations$lat),]
#' 
#'    geoplot(grid=F)
#'    geopoints(stations,pch=".",col=25)
#'    title(main="Before geoinside")
#' 
#'    sea.stations <- geoinside(stations,island,option=2)
#'    # Removes those datapoints from stations where
#'    # measurments took place inside Iceland (misspells).        
#'  
#'    geoplot(grid=F) 
#'    geopoints(sea.stations,pch=".",col=25)       
#'    title(main="After geoinside")
#' 
#'    #######################################################
#'   
#' }
#' @export geoinside
geoinside <-
function(data, reg, option = 1, col.names = c("lat", "lon"), na.rm = T, robust
	 = F)
{
	if(!is.data.frame(data)) {
		i <- match(col.names, names(data))
		data <- data.frame(data[[i[1]]], data[[i[2]]])
		names(data) <- col.names
	}
	i <- match(col.names, names(data))
	index <- rep(NA, nrow(data))
	j <- rep(T, nrow(data))
	tmp <- data
	if(na.rm) {
		j <- !is.na(data[, i[1]]) & !is.na(data[, i[2]])
		data <- data[j,  ]
	}
	i1 <- match(col.names, names(reg))
	regx <- reg[[i1[1]]]
	regy <- reg[[i1[2]]]
	n <- length(regx)
	k <- (regx[1] != regx[n] || regy[1] != regy[n]) && length(regx) != 2
	if(k && !is.na(k)) {
		regx <- c(regx, regx[1])
		regy <- c(regy, regy[1])
	}
	reg <- list(x = regx, y = regy)
	if(length(reg$x) == 2)
		reg <- list(x = c(reg$x[1], reg$x[2], reg$x[2], reg$x[1], reg$
			x[1]), y = c(reg$y[1], reg$y[1], reg$y[2], reg$y[2],
			reg$y[1]))
	data <- list(x = data[[i[1]]], y = data[[i[2]]])
	border <- adapt(reg$y, reg$x, projection = "none")
	inside <- rep(0, length(data$x))
	# Robust method using trigonometric functions.  
	if(robust) {
		a <- a1 <- rep(0, length(reg$x))
		inside <- .C("marghc", PACKAGE = "geo", 
			as.double(data$x),
			as.double(data$y),
			as.integer(length(data$y)),
			as.double(border$x),
			as.double(border$y),
			as.integer(length(border$y)),
			as.integer(border$lxv),
			as.integer(length(border$lxv)),
			as.integer(inside),
			as.double(a),
			as.double(a1))
		inside <- inside[[9]]
	}
	else {
		# Faster method.  
		tmpinside <- rep(0, length(border$lxv))
		inside <- .C("geomarghc", PACKAGE = "geo", 
			as.double(data$x),
			as.double(data$y),
			as.integer(length(data$y)),
			as.double(border$x),
			as.double(border$y),
			as.integer(border$lxv),
			as.integer(length(border$lxv)),
			as.integer(inside),
			as.integer(tmpinside))
		inside <- inside[[8]]
	}
	index[j] <- inside
	inside <- index
	ind <- c(1:length(inside))
	ind <- ind[inside > 0 & !is.na(inside)]
	if(option == 1) {
		tmp <- tmp[ind,  ]
		return(tmp)
	}
	else if(option == 2) {
		tmp <- tmp[ - ind,  ]
		return(tmp)
	}
	else if(option == 3)
		return(inside)
	else if(option == 4)
		return(1 - inside)
	else if(option == 5) {
		ind <- c(1:length(inside))
		ind <- ind[inside == 0]
		return(ind)
	}
	else if(option == 6) {
		ind <- c(1:length(inside))
		ind <- ind[inside != 0]
		return(ind)
	}
	else return(ind)
}

