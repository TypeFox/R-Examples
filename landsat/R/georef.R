georef <-
function(target, tofix, maxdist = 1000, startx = 0, starty = 0)
{
    ## find best-match shift of tofix to target
    ## returns coefficients, but geoshift() must be used to 
    ## actually adjust the matrix/dataframe/SpatialGridDataFrame

    # minimum-path fitting of tofix matrix to target matrix

    target <- as.matrix(target)
    tofix <- as.matrix(tofix)

    # if startx or starty != 0, shift matrices initially
    # this is intended to get away from local minima
    if( startx != 0 | starty != 0) {
	padx <- pady <- max(abs(startx), abs(starty))
	target <- geoshift(target, padx, pady, 0, 0)
	tofix <- geoshift(tofix, padx, pady, startx, starty)
    }

    if(!all(dim(target) == dim(tofix))) stop("target and tofix must be the same size.\n")

    # initial configuration
    thisx <- thisy <- 0
    currrmse <- sqrt(sum((as.vector(target) - as.vector(tofix))^2, na.rm=TRUE) / (sum(!is.na(as.vector(target)) & !is.na(as.vector(tofix)))))
    prevrmse <- currrmse + 1
    maxx <- maxy <- 1
    newx <- newy <- 0
    initrmse <- currrmse

    while(currrmse < prevrmse) {
#        cat(newx, ", ", newy, ": ", currrmse, "\n")

        results <- matrix(NA, nrow=9, ncol=3)
        colnames(results) <- c("x", "y", "RMSE")

        target2 <- geoshift(target, maxx, maxy, 0,  0)
        target2 <- as.vector(target2)

        currrow <- 1

        for(x in seq(newx-1, newx+1, by=1)) {
            for(y in seq(newy-1, newy+1, by=1)) {
                tofix2 <- geoshift(tofix, maxx, maxy, x, y)
                tofix2 <- as.vector(tofix2)
                results[currrow, 1:2] <- c(x, y)
                results[currrow, 3] <- sqrt(sum((target2 - tofix2)^2, na.rm=TRUE) / (sum(!is.na(target2) & !is.na(tofix2))))
                currrow <- currrow + 1
            }
        }

        prevrmse <- currrmse
        currrmse <- min(results[, "RMSE"])

        newx <- results[results[,"RMSE"] == currrmse, "x"]
        newy <- results[results[,"RMSE"] == currrmse, "y"]
        maxx <- max(abs(newx-1), abs(newx+1))
        maxy <- max(abs(newy-1), abs(newy+1))
	
	# check to see if the loop should stop anyway
	if(abs(newx) > maxdist | abs(newy) > maxdist) currrmse <- 9999
    
    }
    list(shiftx=newx, shifty=newy, initrmse=initrmse, currrmse=currrmse)
}

