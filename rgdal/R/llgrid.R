llgridlines = function(obj, easts, norths, ndiscr = 20, 
		lty = 2, offset=0.5, side="WS", 
                llcrs = "+proj=longlat +datum=WGS84",
		plotLines = TRUE, plotLabels = TRUE, ...) {
	obj_ll <- spTransform(obj, CRS(llcrs))
	if (missing(easts))
		easts = pretty(bbox(obj_ll)[1,])
	if (missing(norths))
		norths = pretty(bbox(obj_ll)[2,])
	grd <- gridlines(obj_ll, easts = easts, norths = norths, ndiscr = ndiscr)
	grd_x <- spTransform(grd, CRS(proj4string(obj)))
	if (plotLines)
		plot(grd_x, add = TRUE, lty = lty, ...)
        if (packageVersion("sp") >= "0.9.84") {
	    grdat_ll <- gridat(obj_ll, easts = easts, norths = norths,
                side=side, offset = offset)
        } else {
	    grdat_ll <- gridat(obj_ll, easts = easts, norths = norths, offset = offset)
        }
	grdat_x <- spTransform(grdat_ll, CRS(proj4string(obj)))
	if (plotLabels)
		text(coordinates(grdat_x), labels=parse(text=grdat_x$labels),
  			pos=grdat_x$pos, offset=grdat_x$offset, ...)
}
