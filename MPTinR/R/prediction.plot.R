
prediction.plot <- function(results, model.filename, dataset = 1, absolute = TRUE, spacing = 2, axis.labels = NULL, ylim, model.type = c("easy", "eqn", "eqn2"),
args.plot = list(), args.rect = list(), args.box = list(), args.points = list(), args.labels = list(),
numbers = c("individual", "continuous"),
pos.numbers = c("plot", "axis"),
args.numbers = list(), args.abline = list(), abline) {
	tree <- .get.mpt.model(model.filename, model.type)
	categories.per.type <- vapply(tree, length, 0)
	#browser()
	n.trees <- length(categories.per.type)
	# get y-values
	if (class(results[["data"]][["observed"]]) != "list") {
		data <- results[["data"]][["observed"]] 
		exp <- results[["data"]][["predicted"]]
	} else {
		if (dataset == "aggregated") {
			data <- results[["data"]][["observed"]][["aggregated"]]
			exp <- results[["data"]][["predicted"]][["aggregated"]]
		} else {
			data <- results[["data"]][["observed"]][["individual"]][dataset,]
			exp <- results[["data"]][["predicted"]][["individual"]][dataset,]
		}
	}
	if (absolute) y.values <- data - exp
	else {
		y.values <- data
		y.values[data != 0] <- 2 * data[data != 0] * (log(data[data != 0]) - log(exp[data != 0]))
	}
	# y.values <- data^2 * if (sign) sign(data) else 1 
	# get x-values
	x.bool <- unlist(mapply(function(x, y, spacing) rep(c(y, FALSE), c(x, spacing)), categories.per.type, seq(categories.per.type),spacing = spacing, SIMPLIFY = FALSE))
	x.values <- which(x.bool != 0)
	#get limits
	xlim <- c(1, max(x.values))
	if (missing(ylim)) ylim <- c(min(y.values), max(y.values))
	# plot
	def.args.plot <- list(xlab = "", ylab = "", main = "")
	args.to.plot <- c(x = list(x.values), y = list(y.values), ylim = list(ylim), xaxt = "n", type = "n", args.plot, def.args.plot[!names(def.args.plot) %in% names(args.plot)])
	do.call("plot", args = args.to.plot )
	# rect
	mean.trees <- vector("numeric", n.trees)
	def.args.rect <- list(col = "grey", border = "transparent", density = 30, angle = 45)
	args.to.rect <- c(args.rect, def.args.rect[!names(def.args.rect) %in% names(args.rect)])
	for (c in seq_along(categories.per.type)) {
		tmp.args.to.rect <- c(xleft  = min(which(x.bool == c)) - 1, ybottom = 1.5 * ylim[1], xright = max(which(x.bool == c)) + 1,  ytop = 1.5 * ylim[2], args.to.rect)
		do.call("rect", tmp.args.to.rect)
		mean.trees[c] <- mean(which(x.bool == c))
	}
	#box
	do.call("box", args.box)
	#abline
	if (missing(abline)) if (pos.numbers[1] == "axis") abline <- TRUE
	else abline <- FALSE
	def.args.abline <- list(col = "darkgrey")
	args.to.abline <- c(v = list(x.values), args.abline, def.args.abline[!names(def.args.abline) %in% names(args.abline)])
	if (abline) do.call("abline", args.to.abline)
	# points
	def.args.points <- list(pch = 1, cex = 2.25)
	args.to.points <- c(x = list(x.values), y = list(y.values), args.points, def.args.points[!names(def.args.points) %in% names(args.points)])
	do.call("points", args.to.points)
	def.args.axis <- list(line = -1)
	pos.numbers <- match.arg(pos.numbers, c("plot", "axis"), several.ok = TRUE)
	# numbers
	if (!is.null(numbers)) {
		numbers <- match.arg(numbers, c("individual", "continuous"), several.ok = TRUE)
		if (numbers[1] == "individual") numbers.pch <- unlist(lapply(rle(x.bool[x.bool != 0])[["lengths"]], seq_len))
		if (numbers[1] == "continuous") numbers.pch <- seq_along(y.values)
		if (pos.numbers[1] == "plot") {
			def.args.text <- list(labels = as.character(numbers.pch), cex = 0.7)
			args.to.text <- c(x = list(x.values), y = list(y.values), args.numbers, def.args.text[!names(def.args.text) %in% names(args.numbers)])
			do.call("text", args.to.text)
		}
		if (pos.numbers[1] == "axis") {
			def.args.axis2 <- list(labels = numbers.pch, cex.axis = 0.6, mgp = c(0,0.3,0))
			args.to.axis2 <- c(side = 1, at = list(x.values), args.numbers, def.args.axis2[!names(def.args.axis2) %in% names(args.numbers)])
			do.call("axis", args.to.axis2)
			def.args.axis <- list(line = 1)
		}
	}
	#axis
	if (is.null(axis.labels)) axis.labels <- paste("Tree", seq_len(n.trees))
	args.to.axis <- c(side = 1, at = list(mean.trees), labels = list(axis.labels), tick = FALSE, args.labels, def.args.axis[!names(def.args.axis) %in% names(args.labels)])
	do.call("axis", args.to.axis)
	invisible(list(x = x.values, y = y.values))
}
