#' An Infrastructure Proxy Function
#' 
#' This function takes latitude and longitude coordinates from a data frame
#' and downloads images from map servers for those coordinates.  
#' It will return the file sizes of those map images which can
#' be used as a useful proxy for the level of infrastructure 
#' development.  Map servers currently used are Google, Bing
#' and Openstreetmap.  
#' @param dataFrame The name of your data frame
#' @param server The map server you want to use.  Options are google, bing, and open.  Defaults to google
#' @param zoom The zoom level you want to use.  Defaults to 14.  Valid values for google are 0 to 22; values for bing are 1 to 19; values for open are 0 to 19
#' @keywords infrastructure maps google bing openstreetmap
#' @export
#' @examples
#' # In this example, we'll build a very low resolution (two degrees per grid cell) grid of Japan, 
#' # retrieve the map values from Bing at zoom level 14, and render a simple map. 
#' 
#' # First, define the latitude and longitude bounds of Japan and create a data frame. 
#' # (This uses a very restricted version of the territory of Japan for the sake of speed.)
#' 
#' lats <- rep(rev(seq(32,44, by=2)), each=7)
#' lngs <- rep(seq(130,142, by=2), 7)
#' 
#' coords <- data.frame(lats,lngs)
#' 
#' infraBing <- getInfra(dataFrame = coords, zoom=14, server="bing")
#' 
#' # Now that we have the values in the newly-created infraBing 
#' # variable we can plot them.  First, we need to define our colour 
#' # gamut, running from 0 to the maximum map value retrieved.  
#' 
#' grey <- gray(0:max(infraBing) / max(infraBing))
#' 
#' # Now we can plot.  
#' 
#' png(filename="japanInfraBing.png", width=480, height=480, units="px")
#' plot(lngs, lats, col=grey[infraBing], pch=15, cex=14)
#' dev.off()

getInfra <- function(dataFrame=NULL, server="google", zoom=14) { 

	dataFrameMissing = 0
		
	if (missing(dataFrame)) { 
		
		print ("No data frame declared.", quote=FALSE)
		print ("Syntax is: getInfra(dataFrame, server, zoom)", quote=FALSE)
		print ("Exiting.", quote=FALSE)
		dataFrameMissing = 1
		
	}
	
	if (dataFrameMissing == 0) { 
	
		if (missing(server)) { 

			print ("Server was not declared.  Options are 'google', 'bing' or 'open'.  Defaulting to google.", quote=FALSE)

		}

		if (missing(zoom)) { 

			print ("Zoom was not declared.  Defaulting to 14.  For google, valid values are 0 (whole world) to 22 (street level); 1 to 19 for bing; 0 to 19 for open.", quote=FALSE) 

		}

	
	}
	
	serverFound = 0
	
	if (server == "google") { 
	
		serverFound = 1
		
		lineCounter = 1
		
		googleFileSize <- NULL
		
		progressBar <- txtProgressBar(min = 0, max = length(dataFrame$lat), style = 3)
		
		while(lineCounter < (length(dataFrame$lat) + 1)) { 
		
			degLat = dataFrame$lat[lineCounter]
			degLng = dataFrame$lng[lineCounter]
			zValue = zoom

			mapSize = 2^zValue
			radLat = degLat * (pi / 180)

			pixelY = 0.5*(log((1+sin(radLat))/(1-sin(radLat))))

			pixelY = pixelY / (pi / (mapSize / 2))

			pixelY = (mapSize / 2) - pixelY
			pixelY = round(pixelY)

			pixelX = (mapSize / 2) + (degLng * ((mapSize / 2) / 180))
			pixelX = round(pixelX)

			urlToGet = paste("http://mt0.google.com/vt/lyrs=m@0&hl=en&x=", pixelX, "&y=", pixelY, "&z=", zValue, sep="")
			download.file(url=urlToGet, destfile="googleFileContents.txt", quiet = TRUE, mode = "w", cacheOK = TRUE, extra = getOption("download.file.extra"))
			fileSize = file.info("googleFileContents.txt")["size"]$size

			googleFileSize = c(googleFileSize , fileSize) 
			
			lineCounter = lineCounter + 1
			
			setTxtProgressBar(progressBar, lineCounter)
		
		}
		
		return(googleFileSize)
		file.remove("googleFileContents.txt")
		
		close(progressBar)
	
	}
	
	if (server == "bing") { 
	
		serverFound = 1
		
		lineCounter = 1
		
		bingFileSize <- NULL
		
		progressBar <- txtProgressBar(min = 0, max = length(dataFrame$lat), style = 3)
		
		while(lineCounter < (length(dataFrame$lat) + 1)) { 
		
			inputLat = dataFrame$lat[lineCounter]
			inputLng = dataFrame$lng[lineCounter]
			
			minLat = -85.05112878
			maxLat = 85.05112878
			
			imageWidthAndHeight = bitwShiftL(256, zoom)
			
			# first, convert lat lng to pixel x y
			# make sure they fit in the MS bounds
			
			if (inputLat < minLat) { 
			
				inputLat = minLat
			
			}
			
			if (inputLat > maxLat) { 
			
				inputLat = maxLat
			
			}
			
			firstX = (inputLng + 180) / 360
			sinLat = sin(inputLat * pi / 180) 
			firstY = 0.5 - log((1 + sinLat) / (1 - sinLat)) / (4 * pi)
			
			pixelX = round(firstX * imageWidthAndHeight)
			pixelY = round(firstY * imageWidthAndHeight)
			
			tileX = pixelX / 256
			tileY = pixelY / 256
			
			quadKey = NULL
			
			for(i in zoom:1){
			
				digit = 0
				mask = bitwShiftL(1,i-1)
				
				if (bitwAnd(tileX, mask) != 0) { 
				
					digit = digit + 1
				
				}
				
				if (bitwAnd(tileY, mask) != 0) { 
								
					digit = digit + 2
				
				}
								
				quadKey = paste(quadKey, digit, sep="")
			
			}
			
			urlToGet = paste("http://ak.dynamic.t0.tiles.virtualearth.net/comp/ch/", quadKey, "?mkt=en-gb&it=G,VE,BX,L,LA", sep="")
			
			download.file(url=urlToGet, destfile="bingFileContents.txt", quiet = TRUE, mode = "w", cacheOK = TRUE, extra = getOption("download.file.extra"))
			fileSize = file.info("bingFileContents.txt")["size"]$size
			
			bingFileSize = c(bingFileSize, fileSize)
			
			lineCounter = lineCounter + 1
			setTxtProgressBar(progressBar, lineCounter)
		
		}
		
		return(bingFileSize)
		file.remove("bingFileContents.txt")
				
		close(progressBar)
	
	}
	
	if (server == "open") { 
	
		serverFound = 1
		
		lineCounter = 1
		
		openMapFileSize <- NULL
		
		progressBar <- txtProgressBar(min = 0, max = length(dataFrame$lat), style = 3)
		
		while(lineCounter < (length(dataFrame$lat) + 1)) { 
		
			inputLat = dataFrame$lat[lineCounter]
			inputLng = dataFrame$lng[lineCounter]
			
			numberOfTiles = bitwShiftL(1, zoom)
			
			inputLngMod = inputLng + 180
			pixX = (numberOfTiles / 360) * inputLngMod
			
			inputLatMod = inputLat * pi / 180
			n <- 2 ^ zoom
			
			pixY = floor((1 - log(tan(inputLatMod) + (1 / cos(inputLatMod))) / pi) / 2 * n)
			
			pixX = round(pixX)
			pixY = round(pixY)
			
			urlToGet = paste("http://c.tile.openstreetmap.org/", zoom, "/", pixX, "/", pixY, ".png", sep="")
			openFileNotFound = 0
		
			tryCatch(download.file(url=urlToGet, destfile="openFileContents.txt", quiet = TRUE, mode = "w", cacheOK = TRUE, extra = getOption("download.file.extra")),
		
			warning = function(w) {
				
				Sys.sleep(10)
				
				print("Unable to download image from openstreetmap.  Retrying after 10 seconds...", quote=FALSE) 
		
				tryCatch(download.file(url=urlToGet, destfile="openFileContents.txt", quiet = TRUE, mode = "w", cacheOK = TRUE, extra = getOption("download.file.extra")),
		
				warning = function(w2) {
				
					openFileNotFound = 1
		
					print("Unsuccessful.  This location will be coded as NA.", quote=FALSE)  
		
				}
			
				)
		
				}              
		
			)

			if (openFileNotFound == 0) {
			
				fileSize = file.info("openFileContents.txt")["size"]$size
				
			}
			
			if (openFileNotFound == 1) {
			
				fileSize = "NA"; 
			
			}
			
			openMapFileSize = c(openMapFileSize, fileSize)
			
			lineCounter = lineCounter + 1
			setTxtProgressBar(progressBar, lineCounter)
		
		}
		
		return(openMapFileSize)
		file.remove("openFileContents.txt")
				
		close(progressBar)
	
	}

	
	if (serverFound == 0) { 
	
		print (paste("Server name (", server, ") not recognised.  Options are 'google', 'bing', 'open'.  Exiting."), quote=FALSE)
	
	}
	
	

}
