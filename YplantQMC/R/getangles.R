#'Read angles from a .p file.
#'
#'Reads leaf angles (angle, orientation or azimuth) from a plant file (a Yplant
#'input file with extension \code{.p}, known as a \code{pfile}).
#'
#'If the leaf angle is returned (An.3, the default), all angles are converted
#'so that they are between 0 and 90 degrees. A warning is printed when any
#'angle > 360, which may indicate problems in the data (this is uncommon).
#'
#'Other angles may be read, see \code{\link{modifypfile}} for a list of the
#'angles in a pfile.
#'
#'@param plant An object of class \code{plant3d}, or the name of a
#'\code{pfile}.
#'@param whichangle The name of the angle, in quotes (see
#'\code{\link{modifypfile}})
#'@return A vector of angles (in degrees).
#'@author Remko Duursma
#'@keywords misc
#'@examples
#'
#'
#'\dontrun{
#'# Two options:
#'# Get leaf angles from a pfile
#'ang <- getangles("someplant.p")
#'
#'# Or from a constructed plant:
#'myplant <- constructplant("someplant.p","someleaf.l")
#'ang <- getangles(myplant)
#'}
#'
#'
#'
#'
getangles <- function(plant, whichangle="An.3"){
	
	if(class(plant) == "plant3d"){
		if(plant$inputformat == "P"){
			pdata <- plant$pdata 
			inputformat <- "P"
		} else {
			pdata <- plant$qdata
			inputformat <- "Q"
		}
	} else if(file.exists(plant)){
		pdata <- readp(plant)
		inputformat <- "P"
	}
	
	if(inputformat == "P")
		angles <- pdata[pdata$Lt >= 1, whichangle]
	else
		angles <- pdata$An.3
	
	# Fix angles > 90
	if(whichangle == "An.3"){
	angles[angles > 90 & angles < 180] <- angles[angles > 90 & angles < 180] - 90
	angles[angles > 180 & angles < 270] <- 270 - angles[angles > 180 & angles < 270]
	angles[angles > 270 & angles < 360] <- angles[angles > 270 & angles < 360] - 270
	
	angles[angles < 0] <- -angles[angles < 0]
	}
	
	if(max(angles) > 360)warning("Angle >360deg in file", plant$pfile)
	
return(angles)
}
