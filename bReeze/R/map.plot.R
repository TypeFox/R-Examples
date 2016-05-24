map.plot <-
function(mast, type=c("satellite", "terrain", "hybrid", "roadmap"), zoom, label, ...) {
### plotting map or satellite image of mast location

	stopifnot(requireNamespace("RgoogleMaps", quietly=TRUE))
	
	if(class(mast)!="mast") stop(substitute(mast), " is no mast object")
	if(missing(type)) type <- "satellite"
	type <- match.arg(type)
	if(is.null(mast$location)) stop("No location found")
	lat <- mast$location[1]
	lon <- mast$location[2]
	if(missing(zoom)) zoom <- 15
	if(missing(label)) label <- paste(lat, lon, sep=",")

	plot.param <- list(...)
	if(any(names(plot.param)=="pch")) pch <- plot.param$pch
	else pch <- 8
	if(any(names(plot.param)=="col")) col <- plot.param$col
	else col <- "#E41A1C"
	if(any(names(plot.param)=="cex")) cex <- plot.param$cex
	else cex <- 1.5
	if(any(names(plot.param)=="col.lab")) col.lab <- plot.param$col.lab
	else col.lab <- col
	if(any(names(plot.param)=="cex.lab")) cex.lab <- plot.param$cex.lab
	else cex.lab <- 1
	if(any(names(plot.param)=="pos.lab")) pos.lab <- plot.param$pos.lab
	else pos.lab <- 4
	
	tmp.file <- gsub("[^0-9]", "", substr(Sys.time(), 1, 19))
	tmpmap <- RgoogleMaps::GetMap(center=c(lat, lon), zoom=zoom, destfile=file.path(tempdir(), paste0("map", tmp.file, ".png")), maptype=type, format="png32", verbose=0)
	RgoogleMaps::PlotOnStaticMap(tmpmap, lat=lat, lon=lon, destfile=file.path(tempdir(), paste0("map", tmp.file, ".png")), cex=cex, pch=pch ,col=col, NEWMAP=FALSE)
	if(!is.na(label)) RgoogleMaps::TextOnStaticMap(tmpmap, lat=lat, lon=lon, labels=label, cex=cex.lab, col=col.lab, pos=pos.lab, add=TRUE)
	
	unlink(file.path(tempdir(), paste0("map", tmp.file, ".png.rda")), recursive=FALSE, force=FALSE)
	unlink(file.path(tempdir(), paste0("map", tmp.file, ".png")), recursive=FALSE, force=FALSE)
}
