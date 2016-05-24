#' Split a polygon into a list
#' 
#' Split a 'polygon with NAs' into a list of true polygons before calculating
#' its area.
#' 
#' 
#' @param data A dataframe of coordinates in latitude and longitude
#' @return Returns a list of true polygons as dataframes, each with coordinates
#' \item{lat}{Latitude} \item{lon}{Longitude}
#' @note Needs further elaboration.
#' @seealso Called by \code{\link{geoarea}}.
#' @keywords manip
#' @export geo.Split.poly
geo.Split.poly <-
function(data)
{
	while(is.na(data[nrow(data), "lat"])) data <- data[ - nrow(data),  ]
	while(is.na(data[1, "lat"])) data <- data[-1,  ]
	n <- nrow(data)
	if(any(is.na(data$lat))) {
		tmp <- list()
		i <- 1:nrow(data)
		i <- i[is.na(data$lat)]
		i1 <- c(1, i + 1)
		i2 <- c(i - 1, nrow(data))
		for(i in 1:length(i1)) {
			tmp[[i]] <- data[i1[i]:i2[i],  ]
			n <- nrow(tmp[[i]])
			if(tmp[[i]]$lat[1] != tmp[[i]]$lat[n] || tmp[[i]]$
				lon[1] != tmp[[i]]$lon[n])
				tmp[[i]] <- rbind(tmp[[i]], tmp[[i]][1,  ])
			if(nrow(tmp[[i]]) < 4) {
				cat("minimum of 3 points needed to define area"
					)
				print(tmp[[i]][1:(nrow(tmp[[i]]) - 1)])
			}
		}
	}
	else {
		if(data$lat[1] != data$lat[n] || data$lon[1] != data$lon[n])
			data <- rbind(data, data[1,  ])
		tmp <- list(c1 = data)
	}
	return(tmp)
}

