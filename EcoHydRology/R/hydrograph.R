hydrograph <-
function (input = matrix(ncol = 2, nrow = 2), streamflow = input[, 
    2], timeSeries = input[, 1], streamflow2 = NULL, precip = NULL, 
    begin = 1, endindex = length(streamflow), P.units = "", S.units = "", 
    S1.col = "black", S2.col = "red", stream.label = "Streamflow", 
	streamflow3=NULL,streamflow4=NULL, precip2=NULL) 
{
    if (is.null(streamflow2) & (ncol(input) > 3)) 
        streamflow2 <- input[, 4]
    if (is.null(precip) & (ncol(input) > 2)) {
        precip <- input[, 2]
        streamflow <- input[, 3]
    }
    if (!is.null(precip))  {
        par(mar = c(3, 5, 1, 4))
        barplot(precip[begin:endindex], yaxt = "n", space = NULL, 
            ylim = rev(c(0, 4 * max(na.omit(precip[begin:endindex])))), 
            xaxt = "n")
        axis(side = 3, pos = 0, tck = 0,xaxt = "n")
        axis(side = 4, at = seq(0, floor(max(na.omit(precip[begin:endindex])) + 
            1), length = (1 + ifelse(floor(max(na.omit(precip[begin:endindex])) + 
            1) < 10, floor(max(na.omit(precip[begin:endindex])) + 1), 
            4))), labels = as.integer(seq(0, floor(max(na.omit(precip[begin:endindex])) + 
            1), length = (1 + ifelse(floor(max(na.omit(precip[begin:endindex])) + 
            1) < 10, floor(max(na.omit(precip[begin:endindex])) + 1), 
            4)))))
        if (P.units=="") {
			mtext(paste("Precipitation", P.units), 4, line = 2, cex = 0.9, adj = 1)
		} else  mtext(paste("Precipitation (", P.units, ")", sep=""), 4, line = 2, cex = 0.9, adj = 1)
        par(new = T)
    }
if (!is.null(precip2)){
barplot(precip2[begin:endindex], yaxt = "n", space = NULL, col="pink",
            ylim = rev(c(0, 4 * max(na.omit(precip[begin:endindex])))), 
            xaxt = "n")
        par(new = T)
}
    
	plot(streamflow[begin:endindex], col = S1.col, type = "l", 
		lwd = 1, ylab = stream.label, xaxt = "n", xlab = "date", 
		ylim = c(0, 1.2 * max(na.omit(streamflow[begin:endindex]), na.omit(streamflow2[begin:endindex]))), 
		axes = FALSE)
	#mtext (expression(paste("                              ", " (" , m^3/s, ")", sep="")), 2,3)
	if (S.units=="m3/s" | S.units=="m3s"){
		mtext (expression(paste(" (" , m^3/s, ")", sep="")), 2,1.5)
	} else if (S.units=="ft3/s" | S.units=="ft3s") {
		mtext (expression(paste(" (" , ft^3/s, ")", sep="")), 2,1.5)
	} else if (S.units!="") mtext (paste(" (" , S.units, ")", sep=""), 2,1.5)
	lines(streamflow2[begin:endindex], col = S2.col, lwd = 1, 
        lty = 2, xaxt = "n")
if (!is.null(streamflow3)){
lines(streamflow3[begin:endindex], col = "blue", lwd = 1,   ##potential for more streamflows
lty = 3, xaxt = "n")
}
if (!is.null(streamflow4)){
lines(streamflow4[begin:endindex], col = "green", lwd = 1,   ##potential for more streamflows
lty = 4, xaxt = "n")
}
	axis(side = 1, at = seq(1, (endindex - begin + 1), length = 14), 
		pos = 0, labels = format(timeSeries[seq(begin, endindex, 
		length = 14)], "%d-%b-%y"))
	
    axis(side = 2, pos = 0)
}
