plot.influenceSSN <-
function(x, color.palette = NULL, nclasses = NULL,
	inflcol = "_resid_", breaktype = "quantile", brks = NULL, pch = 19, ...)
{
    par.orig <- par(no.readonly = TRUE)
    if(class(x) != "influenceSSN") return("Not a influenceSSN object")
    layout(matrix(1:2, nrow = 1), widths = c(4,1))
    par(mar = c(5,5,3,0))
    plot(x$ssn.object@bbox[1,],x$ssn.object@bbox[2,], type = "n",
         xlab = "x-coordinate", ylab = "y-coordinate",
         main = paste("Influence Diagnostic = ",inflcol), cex.main = .9, ...)
    for(i in 1:length(x$ssn.object@lines))
    	for(j in 1:length(x$ssn.object@lines[[i]]))
            lines((x$ssn.object@lines[[i]]@Lines[[j]]@coords), ...)
    data <- x$ssn.object@obspoints@SSNPoints[[1]]@point.data
    if(is.null(nclasses)) nclasses <- 10
    lower.breaks <- matrix(0, nrow = nclasses, ncol = 1)
    upper.breaks <- matrix(0, nrow = nclasses, ncol = 1)
    if(breaktype == "quantile") {
        brks <- quantile(data[,inflcol],
                         probs = (1:(nclasses-1))/nclasses, na.rm = T)
    	lower.breaks <- c(min(data[,inflcol], na.rm = T), brks)
    	upper.breaks <- c(brks, max(data[,inflcol], na.rm = T))
    }
    if(breaktype == "even") {
    	brks <- min(data[,inflcol], na.rm=T) +
            (max(data[,inflcol], na.rm=T) - min(data[,inflcol], na.rm=T)) *
      		(1:(nclasses-1))/nclasses
    	lower.breaks <- c(min(data[,inflcol], na.rm = T), brks)
    	upper.breaks <- c(brks, max(data[,inflcol], na.rm = T))
    }
    if(breaktype == "user") {
        if(is.null(brks)) return("Must specify brks if breaktype = user")
        minD <- min(data[,inflcol], na.rm=TRUE)
        maxD <- max(data[,inflcol], na.rm=TRUE)
        if(minD < min(brks)) brks <- c(brks, minD)
        if(maxD > max(brks)) brks <- c(brks, maxD)
        brks <- sort(unique(unlist(brks)))
        nclasses <- length(brks) - 1
        lower.breaks <- brks[1:nclasses]
        upper.breaks <- brks[2:(nclasses+1)]
    }
    if(length(color.palette) == 0)
	color.palette <- rainbow(nclasses, start = .66, end = .99)
    for (j in 1:nclasses){
    	jmax <- upper.breaks[j]
        jmin <- lower.breaks[j]
        indj <- data[,inflcol] >= jmin & data[,inflcol] <= jmax
      	points(x$ssn.object@obspoints@SSNPoints[[1]]@point.coords[indj,,drop = FALSE],
               col = color.palette[j], pch = pch, ...)
    }

    dec.dig <- 2
    left <- as.character(as.numeric(as.integer(lower.breaks*10^dec.dig))/10^dec.dig)
    rght <- as.character(as.numeric(as.integer(upper.breaks*10^dec.dig))/10^dec.dig)
    leglabs <- paste(left," to ",rght)
    par(mar = c(0,0,0,0))
    plot(c(0,0), c(1,1), type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab ="",
         bty = "n")
    legend(x = -1, y = 1.1, legend = leglabs, bty = "n",
           pch = rep(pch, times = length(leglabs)),
           col = color.palette, cex = .8)
    par(par.orig)
    return(invisible(data.frame(lower.breaks = lower.breaks, upper.breaks = upper.breaks)))
}

