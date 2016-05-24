tracksPanel = function(x, y, sp.layout, ...) {
    sppanel(sp.layout, panel.number())
	panel.xyplot(x, y, ...)
}

segPanel = function(x, y, subscripts, ..., x0, y0, x1, y1, 
		arrows, length, col, sp.layout) {
    sppanel(sp.layout, panel.number())
	if (arrows)
		panel.arrows(x0[subscripts], y0[subscripts], 
			x1[subscripts], y1[subscripts], length = length, 
			col = col[subscripts], ...)
	else
		panel.segments(x0[subscripts], y0[subscripts], 
			x1[subscripts], y1[subscripts], 
			col = col[subscripts], ...)
}

stplotTracksCollection = function(obj, ..., by, groups,
		scales = list(draw = FALSE), segments = TRUE, attr = NULL,
		ncuts = length(col.regions), col.regions = bpy.colors(), cuts,
		xlab = NULL, ylab = NULL, arrows = FALSE, length = 0.1,
		xlim = bbexpand(bbox(obj)[1,], 0.04), 
		ylim = bbexpand(bbox(obj)[2,], 0.04),
		sp.layout = NULL) {
	sp = obj@tracksCollection[[1]]@tracks[[1]]@sp
	scales = longlat.scales(sp, scales, xlim, ylim)
	args = list(..., asp = mapasp(sp, xlim, ylim), scales = scales, 
		xlab = xlab, ylab = ylab, arrows = arrows, length = length,
		xlim = xlim, ylim = ylim, sp.layout = sp.layout)
	if (!is.null(attr)) {
		df = as(obj, "segments")
		args$x0 = df$x0
		args$y0 = df$y0
		args$x1 = df$x1
		args$y1 = df$y1
		# compute color:
		z = df[[attr]]
		attr = na.omit(z)
		if (missing(cuts))
			cuts = seq(min(attr), max(attr), length.out = ncuts)
        if (ncuts != length(col.regions)) {
            cols = round(1 + (length(col.regions) - 1) * (0:(ncuts -
                1))/(ncuts - 1))
            fill = col.regions[cols]
        } else
            fill = col.regions
		grps = cut(as.matrix(z), cuts, dig.lab = 4, include.lowest = TRUE)
		args$col = fill[grps]
		# set colorkey:
		args$legend = list(right = list(fun = draw.colorkey,
                args = list(key = list(col = col.regions, at = cuts), 
				draw = FALSE)))
		if (is.null(args$panel))
			args$panel = "segPanel"
		cn = c("x0", "y0")
	} else {
		if (is.null(args$panel))
			args$panel = "tracksPanel"
		df = as(obj, "data.frame")
		cn = coordnames(obj)
		args$type = "l"
	}
	if (!missing(by))
		args$x = as.formula(paste(cn[2], "~", cn[1], "|", by))
	else
		args$x = as.formula(paste(cn[2], cn[1], sep = " ~ "))
	if (!missing(groups))
		args$groups = df[[groups]]
	args$data = df
	do.call(xyplot, args)
}
setMethod("stplot", "TracksCollection", stplotTracksCollection)
