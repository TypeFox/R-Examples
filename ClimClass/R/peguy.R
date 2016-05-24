# TODO: Add comment
# 
# Author: ecor
###############################################################################

NULL
#' @title Peguy Climograph
#' 
#' 
#' @description Representation of Peguy Climograph from monthly weather data (Mean Temperature, Precipitation)
#' 
#' 
#' 
#' 
#' @param data input dataset with climatological monthly weather data
#' @param TemperatureTriangleCoords Temperature coordinates for triangle vertices in the Peguy Climograph. Default coordinates are expressed in Celsius Degrees.
#' @param PrecipitationTriangleCoords Precipitation coordinates for triangle vertices in the Peguy Climograph. Default coordinates are expressed in millimeters.
#' @param xlab,ylab xy axis labels
#' @param lambda.label numeric value used to locate climate attribute labels
#' @param climate.label string vector containing  climate attributes. Default is \code{c("Temperate", "Cold", "Arid", "Hot")}. Alternatively it can be translated into any other languange.   
#' @param xyField column names of \code{data} for the x and y variables used in the Peguy Climate Diagram.
#' @param pointsField column name of \code{data} containing the fields to be represented with different point colors. Default is \code{"month"}.
#' @param StationsField column name of \code{data} containing the fields  with station ID names. Default is \code{"station"}. 
#' @param color.scale character scale indicating a use of a specific color scale. Default is \code{"monthly"}. 
#' @param ...  further arguments
#'
#' @author Emanuele Cordano
#' @import ggplot2 reshape2
#' 
#' @references Peguy, C.P. (1970) Precis de climatologie, ed. Masson, Paris.
#' @importFrom grDevices rainbow
#' @export
#' 
#' @examples
#' 
#' library(stringr)
#' data(Trent_climate)
#### #TrentinoClimateDf <- melt(clima_81_10,id=names(clima_81_10[[1]]))
##### #names(TrentinoClimateDf)[names(TrentinoClimateDf)=="L1"] <- "station"
#' 
#' 
#' TrentinoClimateDf <- do.call(rbind,clima_81_10)
#' names <- rownames(TrentinoClimateDf)
#' TrentinoClimateDf$station <- unlist(lapply(X=str_split(names,pattern="[.]"),FUN=function(x) {x[1]}))
#'  
#' 
#' 
#' data <- TrentinoClimateDf[TrentinoClimateDf$station %in% unique(TrentinoClimateDf$station)[1:3],]
#' p <- peguy(data=data)
#' 
#' 

# http://sape.inf.usi.ch/quick-reference/ggplot2/shape
# 
peguy <- function(data=NULL,TemperatureTriangleCoords=c(0,23.4,15),PrecipitationTriangleCoords=c(0,40,200),ylab="Precipitation[mm]",xlab="Mean Temperature [degC]",
		lambda.label=1.75,climate.label=c("Temperate", "Cool", "Arid", "Hot"),xyField=c("Tn","P"),pointsField="month",StationsField="station",color.scale="monthly",
		...) {
	    mpoints <- NA
		textlabel <- NA 
		idstation <- NA 
		if (!is.null(data)) {
			
			data$TemperatureTriangleCoords <- data[,xyField[1]]
			data$PrecipitationTriangleCoords <- data[,xyField[2]]
			if (pointsField %in% names(data)) data$mpoints <- sprintf("%02d",data[,pointsField])
			if (StationsField %in% names(data)) data$idstation <- data[,StationsField]
			
		} 
#  http://sape.inf.usi.ch/quick-reference/ggplot2/geom_polygon	
#	d=data.frame(x=c(1,2,2, 3,4,4), y=c(1,1,2, 2,2,3), t=c('a', 'a', 'a',  'b', 'b', 'b'), r=c(1,2,3, 4,5,6))
#	ggplot() +
#			geom_polygon(data=d, mapping=aes(x=x, y=y, group=t)) +P
#			geom_point(data=d, aes(x=x, y=y, color=t)) +
#			geom_text(data=d, aes(x=x, y=y, label=r), hjust=0, vjust=1, size=4) +
#			opts(title="geom_polygon", plot.title=theme_text(size=40, vjust=1.5))	
	
	
	
#	
	
	
# http://www.ocf.berkeley.edu/~mikeck/?p=407
#
#
	
## LOOK at geom_path !!!! 

	dfcoords <- data.frame(TemperatureTriangleCoords=TemperatureTriangleCoords,PrecipitationTriangleCoords=PrecipitationTriangleCoords)
	


	
	dfcoords[nrow(dfcoords)+1,] <- dfcoords[1,]
	dfcoords$textlabel <- climate.label
	textcoords <- dfcoords
	for (i in 1:(ncol(textcoords)-1)) {
		textcoords[1,i] <- mean(dfcoords[-1,i])
		
	}
	
	if (length(lambda.label)==(nrow(textcoords)-1)) lambda.label <- c(0,lambda.label)
	if (length(lambda.label)<nrow(textcoords)) lambda.label <- array(lambda.label,nrow(textcoords))
	
	
	for (r in 2:nrow(textcoords)) {
		cond <- which(!(names(textcoords) %in% "textlabel"))
		
		textcoords[r,cond] <- (1-lambda.label[r])*dfcoords[r,cond]+lambda.label[r]*textcoords[1,cond]
		
	}
	
	#
	# INSERIRE LE SCRITTE NEL climogramma!!!!!
	#
	
	out <- ggplot()+geom_path(data=dfcoords,aes(x=TemperatureTriangleCoords,y=PrecipitationTriangleCoords))+xlab(xlab)+ylab(ylab)
	out <- out+geom_text(data=textcoords,aes(x=TemperatureTriangleCoords,y=PrecipitationTriangleCoords,label=textlabel), hjust=0, vjust=1, size=4) 
	
	if (!is.null(data)) {
		
		if (("mpoints" %in% names(data)) & ("idstation" %in% names(data))) out <- out+geom_point(data=data,mapping=aes(x=TemperatureTriangleCoords,y=PrecipitationTriangleCoords,colour=mpoints,shape=idstation))
		if (color.scale=="monthly") {
			lenscale=12
			first=4
			scale=rainbow(lenscale)[lenscale:1][c(first:lenscale,1:(first-1))]
			out <- out+scale_colour_manual(values = scale)
		} else if (is.numeric(color.scale[1])) {
			
			scale=color.scale
			out <- out+scale_colour_manual(values = scale)
			
		}
		
	}
	
	
	
	return(out)
	
	
	
	
	
	
	
}

