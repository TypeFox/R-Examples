## Modified 1/25/2012 11:43AM by M. friendly
# -- added fill_col argument, specifying a function to be used to fill the tiles
# -- added xscale, yscale arguments to show the marginal frequencies at top & right
# -- added line_col to change the color of the diagonal line
## Modified 3/24/2012 11:38AM by M. friendly
# -- fixed buglet with yscale=TRUE and reverse_y=FALSE

"agreementplot" <- function (x, ...)
	UseMethod ("agreementplot")

"agreementplot.formula" <-
		function (formula, data = NULL, ..., subset)
{
	m <- match.call(expand.dots = FALSE)
	edata <- eval(m$data, parent.frame())
	if (inherits(edata, "ftable") || inherits(edata, "table")) {
		data <- as.table(data)
		varnames <- attr(terms(formula), "term.labels")
		if (all(varnames != "."))
			data <- margin.table(data, match(varnames, names(dimnames(data))))
		agreementplot(data, ...)
	}
	else {
		if (is.matrix(edata))
			m$data <- as.data.frame(data)
		m$... <- NULL
		m[[1]] <- as.name("model.frame")
		mf <- eval(m, parent.frame())
		if (length(formula) == 2) {
			by <- mf
			y <- NULL
		}
		else {
			i <- attr(attr(mf, "terms"), "response")
			by <- mf[-i]
			y <- mf[[i]]
		}
		by <- lapply(by, factor)
		x <- if (is.null(y))
					do.call("table", by)
				else if (NCOL(y) == 1)
					tapply(y, by, sum)
				else {
					z <- lapply(as.data.frame(y), tapply, by, sum)
					array(unlist(z), dim = c(dim(z[[1]]), length(z)),
							dimnames = c(dimnames(z[[1]]),
									list(names(z))))
				}
		x[is.na(x)] <- 0
		agreementplot(x, ...)
	}
}

"agreementplot.default" <- function(x,
				reverse_y = TRUE,
				main = NULL,
				weights = c(1, 1 - 1 / (ncol(x) - 1)^2),
				margins = par("mar"),
				newpage = TRUE,
				pop = TRUE,
				xlab = names(dimnames(x))[2],
				ylab = names(dimnames(x))[1],
				xlab_rot = 0, xlab_just = "center",
				ylab_rot = 90, ylab_just = "center",
				fill_col = function(j) gray((1 - (weights[j]) ^ 2) ^ 0.5),
				line_col = "red",
				xscale=TRUE, yscale = TRUE,
                                return_grob = FALSE,
                                prefix = "",
				...)
{
	if (length(dim(x)) > 2)
		stop("Function implemented for two-way tables only!")
	if (ncol(x) != nrow(x))
		stop("Dimensions must have equal length!")

	nc <- ncol(x)

	## compute relative frequencies
	n <- sum(x)
	colFreqs <- colSums(x) / n
	rowFreqs <- rowSums(x) / n

	## open viewport
	if (newpage) grid.newpage()
	pushViewport(plotViewport(margins,
                                  name = paste(prefix,"agreementplot")))
	pushViewport(viewport(width = unit(1, "snpc"),
                              height = unit(1, "snpc")))

	if(!is.null(main))
		grid.text(main, y = unit(1.1, "npc"),
				gp = gpar(fontsize = 25))

	## axis labels
	grid.text(xlab, y = -0.12, gp = gpar(fontsize = 20))
	grid.text(ylab, x = -0.1, gp = gpar(fontsize = 20), rot = 90)

	grid.rect(gp = gpar(fill = "transparent"))

	xc <- c(0, cumsum(colFreqs))
	yc <- c(0, cumsum(rowFreqs))

	my.text <- if(reverse_y)
				function(y, ...) grid.text(y = y, ...)
			else
				function(y, ...) grid.text(y = 1 - y, ...)

	my.rect <- if(reverse_y)
				function(xleft, ybottom, xright, ytop, ...)
					grid.rect(x = xleft, y = ybottom, width = xright - xleft,
							height = ytop - ybottom, just = c("left","bottom"), ...)
			else
				function(xleft, ybottom, xright, ytop, ...)
					grid.rect(x = xleft, y = 1 - ybottom, width = xright - xleft,
							height = ytop - ybottom, just = c("left","top"), ...)

	A <- matrix(0, length(weights), nc)
	for (i in 1:nc) {
		## x - axis
		grid.text(dimnames(x)[[2]][i],
				x = xc[i] + (xc[i + 1] - xc[i]) / 2,
				y = - 0.04, check.overlap = TRUE, rot = xlab_rot, just = xlab_just, ...)

		## y - axis
		my.text(dimnames(x)[[1]][i],
				y = yc[i] + (yc[i + 1] - yc[i]) / 2,
				x = - 0.03, check.overlap = TRUE, rot = ylab_rot, just = ylab_just, ...)

		## expected rectangle
		my.rect(xc[i], yc[i], xc[i + 1], yc[i + 1])

		## observed rectangle
		y0 <- c(0, cumsum(x[i,])) / sum(x[i,])
		x0 <- c(0, cumsum(x[,i])) / sum(x[,i])

		rec <- function (col, dens, lb, tr)
			my.rect(xc[i] + (xc[i + 1] - xc[i]) * x0[lb],
					yc[i] + (yc[i + 1] - yc[i]) * y0[lb],
					xc[i] + (xc[i + 1] - xc[i]) * x0[tr],
					yc[i] + (yc[i + 1] - yc[i]) * y0[tr],
					gp = gpar(fill = fill_col(j), col = col, rot = 135)
			)

		for (j in length(weights):1) {
			lb <- max(1, i - j + 1)
			tr <- 1 + min(nc, i + j - 1)
			A[j, i] <- sum(x[lb:(tr-1),i]) * sum(x[i, lb:(tr-1)])
			rec("white", NULL, lb, tr) ## erase background
			rec("black", if (weights[j] < 1) weights[j] * 20 else NULL, lb, tr)
		}

		## correct A[j,i] -> not done by Friendly==Bug?
		for (j in length(weights):1)
			if (j > 1) A[j, i] <- A[j, i] - A[j - 1, i]
	}
	if (reverse_y)
		grid.lines(c(0, 1), c(0, 1), gp = gpar(col = line_col, linetype = "longdash"))
	else
		grid.lines(c(0, 1), c(1, 0), gp = gpar(col = line_col, linetype = "longdash"))

	if (xscale) {
		cx <- xc[-(nc+1)] + diff(xc)/2
		grid.text(colSums(x),
				x = cx,
				y =  1.03, rot = xlab_rot, just = xlab_just, ...)
		grid.xaxis(at = xc, label = FALSE, main=FALSE,
				gp = gpar(fontsize=10), draw = TRUE, vp = NULL)
	}
	if (yscale) {
		cy <- yc[-(nc+1)] + diff(yc)/2
		my.text(rowSums(x),
				x = 1.04,
				y = cy, rot = 0, just = ylab_just, ...)
		grid.yaxis(at = if(reverse_y) yc else 1-yc,
				FALSE, main=FALSE,
				gp = gpar(fontsize=10), draw = TRUE, vp = NULL)
	}

	if (pop) popViewport(2) else upViewport(2)

	## Statistics - Returned invisibly
	ads <- crossprod(diag(x))
	ar  <- n * n * crossprod(colFreqs, rowFreqs)
        if (return_grob)
            invisible(structure(list(
                Bangdiwala = ads / ar,
                Bangdiwala_Weighted = (sum(weights * A)) /  ar,
                weights = weights),
                                grob = grid.grab()
                                )
                      )
        else
            invisible(list(
                Bangdiwala = ads / ar,
                Bangdiwala_Weighted = (sum(weights * A)) /  ar,
                weights = weights)
                      )
}


