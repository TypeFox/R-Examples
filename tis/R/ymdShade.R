ymdShade <- function(ymds, col = grey(0.8), border = FALSE, xpd = FALSE,
                     xrange = NULL, openShade = TRUE, ...){
  dotArgs <- list(...)
  locations <- ymdXy(ymds, openShade = openShade, xrange = xrange)
  polyArgs <- list(x = locations$x,
                   y = locations$y,
                   col = col,
                   border = border,
                   xpd = xpd)
  if(length(dotArgs) > 0)
    polyArgs <- updateList(polyArgs, dotArgs)
  do.call("polygon", polyArgs)
  if(!is.null(locations$vLine))
    abline(v = locations$vLine, lwd = 1.5, col = col)
}

ymdXy <- function(ymds, xrange = NULL, openShade = TRUE){
  ## Returns a list of x and y coordinates that can be fed to polygon() to
  ## draw shadings on the current plot.  If openShade is FALSE and the
  ## last row of nberDates() has a Start but End is NA, and the Start is
  ## within the range of the current plot, the returned list will also include
  ## a vLine element that gives the x coordinate of the last Start. 
  
  if(as.vector(dev.cur()) < 2) ## there is no active graphics device
	return(list(x = NULL, y = NULL, vLine = NULL))
  
  usr <- par("usr")
  if(is.null(xrange))
    xrange <- usr[1:2]
  yrange <- usr[3:4]
  if(par("ylog")) 
    yrange <- 10^yrange

  ## ymds should be a matrix of dates in yyyymmdd form with columns
  ## named "Start" and "End".
  ymds[] <- ymd(ymds)
  if(!is.matrix(ymds)) stop("ymds is not a matrix")
  if(ncol(ymds) != 2) stop("ymds does not have two columns")
  if(is.null(colnames(ymds))){
    if(mean(ymds[,1] - ymds[,2], na.rm = T) <= 0)
      colnames(ymds) <- c("Start", "End")
    else
      colnames(ymds) <- c("End", "Start")
  }
  colnames(ymds) <- capitalize(tolower(colnames(ymds)))
  if(any(sort(colnames(ymds)) != c("End", "Start")))
    stop("ymds must have columns named 'Start' and 'End'")
  naSpots <- is.na(ymds)
  x <- ymds
  x[!naSpots]   <- time(jul(x[!naSpots]))

  ## handle case where latest recession end is NA
  finalRow <- x[nrow(x),]
  if(is.na(finalRow["End"]) && (finalRow["Start"] < xrange[2]) && !openShade){
    x <- x[1:(nrow(x)-1),]
    vLine <- finalRow["Start"]
  }
  else vLine <- NULL
  
  firstRow <- sum(x[,"End"] <= xrange[1], na.rm = TRUE) + 1
  lastRow  <- sum(x[,"Start"] <= xrange[2], na.rm = TRUE)
  if((firstRow > nrow(x)) || lastRow == 0){
	## entire plot is out of the range covered by x
	return(list(x = NULL, y = NULL, vLine = NULL))
  }
  x <- x[firstRow:lastRow,,drop = FALSE]
  n <- nrow(x)
  x[1, "Start"] <- max(xrange[1], x[1, "Start"])
  if(is.na(x[n, "End"]) || x[n, "End"] > xrange[2])
	x[n, "End"] <- xrange[2]

  poly.x <- poly.y <- numeric(0)
  for(i in 1:n){
	poly.x <- c(poly.x, rep(x[i, "Start"], 2), rep(x[i, "End"], 2))
	poly.y <- c(poly.y, yrange, rev(yrange))
	if(i < n){
	  poly.x <- c(poly.x, NA)
	  poly.y <- c(poly.y, NA)
	}
  }
  list(x = poly.x, y = poly.y, vLine = vLine)
}
