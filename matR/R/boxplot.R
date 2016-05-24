
#-----------------------------------------------------------------------------------------
#  boxplots of biom objects.
#-----------------------------------------------------------------------------------------

boxplot.biom <- function(
	x, y=NULL, 
	..., 
	map=NULL,
	columns=TRUE) {

	x <- x [, columns]
	if (!missing(y)) {
		y <- y [, columns]
		if (!all.equal (colnames (x), colnames (y)))
			stop ("\'x\' and \'y\' must have identical column ids")
		}
	arg <- list(...)

####  for first plot (analogously for second):
####    drop "y.*" args
####    strip "x.*" from arg names
####    apply metadata substitution to "names"
####    drop "y.*" mappings
####    strip "x.* from mapping names
####    apply mapping

	x.arg <- arg [substr(names(arg),1,2) != "y."]
	names(x.arg) <- ifelse(
		substr(names(x.arg),1,2) == "x.",
		substring(names(x.arg),3),
		names(x.arg))
	y.arg <- arg [substr(names(arg),1,2) != "x."]
	names(y.arg) <- ifelse(
		substr(names(y.arg),1,2) == "y.",
		substring(names(y.arg),3),
		names(y.arg))

	x.arg$names <- subColumn(x.arg$names, x)
	y.arg$names <- subColumn(y.arg$names, y)

	map <- as.list (map)
	x.map <- map [substr(names(map),1,2) != "y."]
	names(x.map) <- ifelse(
		substr(names(x.map),1,2) == "x.",
		substring(names(x.map),3),
		names(x.map))
	x.map <- unlist (x.map)

	y.map <- map [substr(names(map),1,2) != "x."]
	names(y.map) <- ifelse (
		substr(names(y.map),1,2) == "y.",
		substring(names(y.map),3),
		names(y.map))
	y.map <- unlist (y.map)

	x.arg [names (x.map)] <- parMap (x, x.map, x.arg)
	y.arg [names (y.map)] <- parMap (y, y.map, y.arg)

	par0 <- list(
		x = as.matrix (x, TRUE),
		main = NULL,						# no title
		names = colnames(x),				# x-axis annotations
#		notch = TRUE,						# box notch -- often looks bad & produces warnings
		las	= 2, 			 				# label orientation (x: vertical; y: horizontal)
		cex.axis = 0.5, 					# axis annotation size (smallish)
		font.axis = 3,						# axis annotation font style (italic)
		outpch = 19,						# outlier character (filled circle)
		outcex = 0.5)						# outlier size (small)
	par <- resolve (x.arg, par0)

	zz <- list()
	if (missing(y)) {
		zz <- do.call (graphics::boxplot, par)
		zz$call <- match.call()
		return(invisible(zz))
		}
	screen.n <- split.screen(c(2,1))

	screen (screen.n [1])
	zz$x <- do.call (graphics::boxplot, par)

	screen (screen.n [2])
	par0$x <- y
	par <- resolve (y.arg, par0)
	zz$y <- do.call (graphics::boxplot, par)

	close.screen (screen.n)
	zz$call <- match.call()
	return(invisible(zz))
	}
