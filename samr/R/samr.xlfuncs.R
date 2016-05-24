samr.options <- list(debug=TRUE, #whether to turn on debugging or not
	err.file=ifelse(.Platform$OS.type=='windows', 'C:/samrtrace.txt', 'C:/samrtrace.txt'),
	image.file=ifelse(.Platform$OS.type=='windows', 'C:/samrimage.Rdata','samrimage.Rdata')
)
#
# Our error handler
#
samr.xl.error.trace <- function() {
	err.message <- geterrmessage()
	if (!is.null(samr.options$image.file)) {
		save.image(samr.options$image.file)
	}
	if (!is.null(samr.options$err.file)) {
		sink(samr.options$err.file)
		print(err.message)
		traceback()
		sink()
	}
	winDialog(type = "ok", message = err.message)
}
##
## Upon loading, if we are in a windows environment, we use
#   the windows
## dialog mechanism to display errors. Useful for debugging
#   COM apps
##
.onLoad <- function(lib, pkg) {
	if (.Platform$OS.type == "windows") {
		# options(error=function() winDialog(type='ok',
		#   message=geterrmessage()))
		options(error = samr.xl.error.trace)
	}
}
##
## Upon unload, we set things back the way they were...
##
.onUnload <- function(libpath) {
	if (.Platform$OS.type == "windows") {
		options(error = NULL)
	}
}
samr.xl.build.data <- function(x, y, geneid, genenames, 
	logged2, resp.type) {
	censoring.status = NULL
	eigengene.number = NULL
	if (resp.type == samr.const.survival.response) {
		junk = samr.xl.parse.survival(y)
		y = junk$tim
		censoring.status = junk$censoring.status
	}
	if (resp.type == samr.const.patterndiscovery.response) {
		eigengene.number = samr.xl.parse.patterndiscovery(y)
		if (eigengene.number < 1 | eigengene.number > ncol(x)) {
			stop("pattern discovery specified; response line should be of the form\n`eigengenek', where k is an integer in 1...n (n is the number of samples)")
		}
	}
	return(list(x = x, y = y, censoring.status = censoring.status, 
		eigengene.number = eigengene.number, geneid = geneid, 
		genenames = genenames, logged2 = logged2))
}
samr.xl.compute.plot.xy <- function(samr.obj, del, 
	min.foldchange = 0, plot = FALSE) {
	## make observed-expected plot
	##
	## takes foldchange into account too
	## also returns info needed to make plot in Excel
	## xval, yval and col= color for each point (1=black, 2=red
	#   for top right, 3=green for bottom left)
	LARGE = 1e+10
	b <- detec.slab(samr.obj, del, min.foldchange)
	bb <- c(b$pup, b$plow)
	b1 = LARGE
	b0 = -LARGE
	if (!is.null(b$pup)) {
		b1 <- min(samr.obj$tt[b$pup])
	}
	if (!is.null(b$plow)) {
		b0 <- max(samr.obj$tt[b$plow])
	}
	c1 <- (1:samr.obj$n)[sort(samr.obj$tt) >= b1]
	c0 <- (1:samr.obj$n)[sort(samr.obj$tt) <= b0]
	c2 <- c(c0, c1)
	foldchange.cond.up = rep(T, length(samr.obj$evo))
	foldchange.cond.lo = rep(T, length(samr.obj$evo))
	if (!is.null(samr.obj$foldchange[1]) & (min.foldchange > 
		0)) {
		foldchange.cond.up = samr.obj$foldchange >= min.foldchange
		foldchange.cond.lo = samr.obj$foldchange <= 1/min.foldchange
	}
	col = rep(samr.const.black.color, length(samr.obj$evo))
	col[b$plow] = samr.const.green.color
	col[b$pup] = samr.const.red.color
	if (!is.null(samr.obj$foldchange[1]) & (min.foldchange > 
		0)) {
		col[!foldchange.cond.lo & !foldchange.cond.up] = samr.const.black.color
	}
	col.ordered = col[order(samr.obj$tt)]
	if (plot) {
		ylims <- range(samr.obj$tt)
		xlims <- range(samr.obj$evo)
		plot(samr.obj$evo, sort(samr.obj$tt), xlab = "expected score", 
			ylab = "observed score", ylim = ylims, xlim = xlims, 
			type = "n")
		points(samr.obj$evo, sort(samr.obj$tt), col = col.ordered)
		abline(0, 1)
		abline(del, 1, lty = 2)
		abline(-del, 1, lty = 2)
	}
	sorted.tt = sort(samr.obj$tt)
	neg = list(x = samr.obj$evo[col.ordered == samr.const.green.color], 
		y = sorted.tt[col.ordered == samr.const.green.color])
	pos = list(x = samr.obj$evo[col.ordered == samr.const.red.color], 
		y = sorted.tt[col.ordered == samr.const.red.color])
	rest = list(x = samr.obj$evo[col.ordered != samr.const.green.color & 
		col.ordered != samr.const.red.color], y = sorted.tt[col.ordered != 
		samr.const.green.color & col.ordered != samr.const.red.color])
	return(list(pos = pos, neg = neg, rest = rest, xRange = range(samr.obj$evo)))
}
samr.xl.parse.survival = function(y) {
	remove.leading.spaces = function(x) {
		xx = x
		n = nchar(x)
		if (substring(x, 1, 1) == " ") {
			j = 1
			while (substring(x, j, j) == " ") {
				j = j + 1
			}
			xx = substring(x, j, n)
		}
		return(xx)
	}
	remove.trailing.spaces = function(x) {
		xx = x
		n = nchar(x)
		if (substring(x, n, n) == " ") {
			j = n
			while (substring(x, j, j) == " ") {
				j = j - 1
			}
			xx = substring(x, 1, j)
		}
		return(xx)
	}
	n = length(y)
	tim = rep(NA, n)
	censoring.status = rep(NA, n)
	for (i in 1:n) {
		str = y[i]
		nn = nchar(str)
		j = 1
		while (substring(str, j, j) != "(" & j < nn) {
			j = j + 1
		}
		jj = j + 1
		while (substring(str, jj, jj) != "," & jj < nn) {
			jj = jj + 1
		}
		tim[i] = substring(str, j + 1, jj - 1)
		j = jj + 1
		while (substring(str, j, j) != ")" & j <= nn) {
			j = j + 1
		}
		censoring.status[i] = as.numeric(substring(str, jj + 
			1, j - 1))
	}
	if (sum(is.na(tim)) > 0 | sum(is.na(censoring.status)) > 
		0) {
		stop(" Format error in  survival times")
	}
	for (i in 1:n) {
		tim[i] = remove.leading.spaces(tim[i])
		tim[i] = remove.trailing.spaces(tim[i])
		censoring.status[i] = remove.leading.spaces(censoring.status[i])
		censoring.status[i] = remove.trailing.spaces(censoring.status[i])
	}
	tim = as.numeric(tim)
	censoring.status = as.numeric(censoring.status)
	return(list(tim = tim, censoring.status = censoring.status))
}
samr.xl.parse.patterndiscovery = function(y) {
	# extract eigengene number if pattern discovery specified
	y = y[1]
	if (substring(y, 1, 9) != "eigengene") {
		stop("pattern discovery specified; response line should be of the form\n`eigengenek', where k is an integer in 1...n (n is the number of samples)")
	}
	if (nchar(y) < 10) {
		stop("pattern discovery specified; response line should be of the form\n`eigengenek', where k is an integer in 1...n (n is the number of samples)")
	}
	eigengene.number = round(as.numeric(substring(y, 10, nchar(y))))
	return(eigengene.number)
}
samr.xl.impute.data <- function(data.obj) {
	require(impute)
	data.obj$x = impute.knn(data.obj$x, k = samr.xl.var.knn.neighbors)
	data.obj
} 
