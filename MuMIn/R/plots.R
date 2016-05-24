`plot.model.selection` <-
function(x,
	ylab = NULL, xlab = NULL,
	labels = attr(x, "terms"), labAsExpr = FALSE,
	col = c("SlateGray", "SlateGray2"), col2 = "white",
	border = par("col"),
	par.lab = NULL, par.vlab = NULL,
	axes = TRUE, ann = TRUE,
	...) {
	
	if (is.null(xlab)) xlab <- NA  
	if (is.null(ylab)) ylab <- expression("Cumulative Akaike weight" ~~(omega))

	op <- par(..., no.readonly = TRUE)
	on.exit(par(op))
	
	cumweight <- cumsum(weight <- Weights(x))
	stdweight <- weight / max(weight)
	
	n <- nrow(x)
	m <- length(attr(x, "terms"))
	plot.new()
	plot.window(xlim = c(0, m), ylim = c(1, 0), xaxs = "i", yaxs = "i")

	pal <- if(is.na(col2)) rbind(col) else 
		vapply(col, function(x) grDevices::rgb(grDevices::colorRamp(c(col2, x))(stdweight),
			maxColorValue = 255), character(n))
	npal <- ncol(pal)
		
	for(i in 1L:m)
		rect(i - 1, c(0, cumweight), i, c(cumweight, 1),
			col = ifelse(is.na(x[, i]), NA, pal[, 1L + ((i - 1L) %% npal)]),
			border = border)

	if(ann) {
		labCommonArg <- list(col = par("col.axis"), font = par("font.axis"), cex = par("cex.axis"))
		if(labAsExpr) {
			labels <- gsub(":", "%*%", labels, perl = TRUE)
				labels <- gsub("\\B_?(\\d+)(?![\\w\\._])", "[\\1]", labels, perl = TRUE)
			labels <- parse(text = labels)
		}	
		arg <- c(list(side = 3L, padj = 0.5, line = 1, las = 2), labCommonArg)
		for(i in names(par.lab)) arg[i] <- par.lab[i]
		
		if(is.expression(labels)) {
			if(length(labels) != m) stop("length of 'labels' is not equal to number of terms")
			for(i in 1L:m) do.call("mtext", c(list(text = as.expression(labels[[i]]), at = i - 0.5), arg))
		} else if (!is.null(labels) && !is.na(labels)) {
			if(length(labels) != m) stop("length of 'labels' is not equal to number of terms")
			do.call("mtext", c(list(text = labels, at = 1L:m - 0.5), arg))
		}
	   
		arg <- c(list(side = 4, las = 2, line = 1, adj = 1), labCommonArg)
		for(i in names(par.vlab)) arg[i] <- par.vlab[i]
		ss <- weight > -(1.2 * strheight("I", cex = arg$cex))
		arg[['at']] <- (c(0, cumweight[-n]) + cumweight)[ss] / 2
		arg[['text']] <- rownames(x)[ss]
		arg$line <- arg$line + max(strwidth(arg[['text']], cex = arg$cex,
			units = "in")) / par("mai")[4L] * par("mar")[4L]
		do.call(mtext, arg)

		title(ylab = ylab, xlab = xlab)
	}
	if(axes) {
		axis(2L, col = border, col.ticks = border)
		box(col = border)
	}
	invisible(x)
}
