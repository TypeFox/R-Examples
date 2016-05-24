plot.glmssn.predict <-
function(x, VariableName = NULL, VarPlot = "Both",
  color.palette = rainbow(nclasses, start = .66, end = .99),
  nclasses = 10, breaktype = "quantile", dec.dig = 2,
  SEcex.min = 0.5, SEcex.max = 2, brks = NULL, add = FALSE, ...)
{
    par.orig <- par(no.readonly = TRUE)
    if(class(x) != "glmssn.predict") return("Not a glmssn.predict object")
    if(!any(VarPlot == c("Both", "Predictions", "Standard Errors")))
	return("VarPlot must be one of Both, Predictions, or Standard Errors")
    if(!any(breaktype == c("quantile", "even", "user")))
	return("breaktype must be one of quantile, even, or user")

    if(is.null(as.list(match.call()[-1])$pch)) {
        plch = 19
    } else plch <- as.list(match.call()[-1])$pch
    if(is.null(VariableName)) zcol <- x$args$zcol
    if(!is.null(VariableName)) {
	VarPlot <- "Predictions"
	zcol <- VariableName
    }
    zSEcol <- paste(zcol, ".predSE", sep = "")
    zcolplot <- zcol
    if(VarPlot == "Standard Errors") zcolplot <- zSEcol
    if(add == FALSE) {
        layout(matrix(1:2, nrow = 1), widths = c(4,1))
        par(mar = c(5,5,3,0))
        plot(x$ssn.object@bbox[1,],x$ssn.object@bbox[2,], type = "n",
             xlab = "x-coordinate", ylab = "y-coordinate",
             main = paste("Prediction Variable = ", zcol, ":  Plotting", VarPlot),
             cex.main = .9, ...)
        for(i in 1:length(x$ssn.object@lines))
            for(j in 1:length(x$ssn.object@lines[[i]]))
                lines((x$ssn.object@lines[[i]]@Lines[[j]]@coords), ...)
    }
    predpointsID <- x$args$predpointsID
    for(i in 1:length((x$ssn.object@predpoints@SSNPoints)))
        if(x$ssn.object@predpoints@ID[i] == predpointsID){
            datap <- x$ssn.object@predpoints@SSNPoints[[i]]@point.data
            pcoord <- x$ssn.object@predpoints@SSNPoints[[i]]@point.coords
        }
    if(is.null(nclasses)) nclasses <- 10
    if(breaktype == "quantile") {
	brks <- quantile(datap[,zcolplot],
                         probs = (1:(nclasses-1))/nclasses, na.rm = T)
	lower.breaks <- c(min(datap[,zcolplot], na.rm = T), brks)
	upper.breaks <- c(brks, max(datap[,zcolplot], na.rm = T))
    }
    if(breaktype == "even") {
	brks <- min(datap[,zcolplot]) +
            (max(datap[,zcolplot]) - min(datap[,zcolplot])) *
		(1:(nclasses-1))/nclasses
	lower.breaks <- c(min(datap[,zcolplot], na.rm = T), brks)
	upper.breaks <- c(brks, max(datap[,zcolplot], na.rm = T))
    }
    if(breaktype == "user") {
            if(is.null(brks)) return("Must specify brks if breaktype = user")
            minD <- min(datap[,zcolplot], na.rm=TRUE)
            maxD <- max(datap[,zcolplot], na.rm=TRUE)
            brks <- as.vector(unlist(brks))
            if(minD < min(brks)) brks <- c(brks, minD)
            if(maxD > max(brks)) brks <- c(brks, maxD)
            brks <- sort(unique(unlist(brks)))
            nclasses <- length(brks) - 1
            lower.breaks <- brks[1:nclasses]
            upper.breaks <- brks[2:(nclasses+1)]
            if(nclasses > length(color.palette)) {
                cat("Warning: not enough colours in color.palette, adding random colours as required\n")
                newColours <- setdiff(colors(),color.palette)
                color.palette <- c(color.palette,sample(newColours,nclasses-length(color.palette)))
            }
            if(nclasses < length(color.palette)) {
                cat("Warning: too many colours specified in color.palette. Dropping extra ones\n")
                color.palette <- color.palette[1:nclasses]
            }
        }
    if(VarPlot == "Both") SErange <- max(datap[,zSEcol]) - min(datap[,zSEcol])
    if(add == TRUE) {
	par(new = TRUE)
  	layout(matrix(1:2, nrow = 1), widths = c(4,1))
  	par(mar = c(5,5,3,0))
	par(mfg = c(1,1))
	plot(x$ssn.object@bbox[1,],x$ssn.object@bbox[2,], type = "n", bty = "n",
             xlab = "", ylab = "",...)
    }
    for (j in 1:nclasses){
	jmax <- upper.breaks[j]
	jmin <- lower.breaks[j]
	indj <- datap[,zcolplot] >= jmin & datap[,zcolplot] <= jmax
	if(VarPlot == "Both") {
            points(pcoord[indj,], col = color.palette[j], pch = 19,
                   cex = SEcex.max - (SEcex.max - SEcex.min)*
                   (datap[indj,zSEcol] - min(datap[,zSEcol]))/SErange, ...)
	} else points(pcoord[indj,], col = color.palette[j], pch = 19, ...)
    }
    if(add == FALSE) {
	left <- as.character(as.numeric(
                                        as.integer(lower.breaks*10^dec.dig))/10^dec.dig)
	rght <- as.character(as.numeric(
                                        as.integer(upper.breaks*10^dec.dig))/10^dec.dig)
	leglabs <- paste(left,"to",rght)
	par(mar = c(5,0.1,3,0))
	plot(c(0,0), c(1,1), type = "n", xaxt = "n", yaxt = "n", xlab = "",
             ylab ="", bty = "n")
        title=NULL
        if(VarPlot=="Both") title="Predictions"
	legend("bottomleft", legend = leglabs, bty = "n",
               pch = rep(plch, times = length(leglabs)),
               col = color.palette, title=title, cex = .8)
        if(VarPlot=="Both") {
            title = "Standard Errors"
            cexVals <- SEcex.max - (SEcex.max - SEcex.min)*(datap[,zSEcol] - min(datap[,zSEcol]))/SErange
            cexVals <- quantile(cexVals,seq(0,1,by=0.2))
            cexLab <-  signif(quantile(datap[,zSEcol],seq(1,0,by=-0.2)),3)
            legend("topleft", legend = cexLab, bty = "n",
                   pch = 1,pt.cex=cexVals,col=1, title=title, cex = .8)
        }
	par(par.orig)
	return(invisible(data.frame(lower.breaks = lower.breaks, upper.breaks = upper.breaks)))
    }
    if(add == TRUE) {
	par(mar = c(0,0,0,0))
	par(mfg=c(1,2))
	plot(c(0,0), c(1,1), type = "n", xaxt = "n", yaxt = "n",
             xlab = "", ylab ="", bty = "n")
	lp <- as.character(as.numeric(as.integer(
                                                 min(datap[,zcolplot])*10^dec.dig))/10^dec.dig)
	up <- as.character(as.numeric(as.integer(
                                                 max(datap[,zcolplot])*10^dec.dig))/10^dec.dig)
	text(0, 1.3, paste("lowest pred =", lp), cex = .8)
	text(0, .7, paste("highest pred =", up), cex = .8)
	par(par.orig)
    }

}

