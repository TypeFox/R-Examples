#'@title sequential HSV colors
#'
#'@description
#'This functions allows to get a sequence of colors in an HSV model with
#'optional pre-especified numbers for saturation, value, and alpha. It is a
#'very flexible function to play with different combinations of saturation,
#'value, and alpha.
#'
#'@details
#'The idea bechind this function is to explore a sequence of colors given some
#'fixed numbers of saturation, valur or alpha for an HSV color model. The
#'argument \code{what} will be taken to generate the sequence in the given
#'\code{percentage} increment steps. In addition, we can specify a number for
#'\code{s, v, alpha}. For example, if \code{what="value"}, we can fix the
#'saturation in \code{s=0.8}, obtaining a sequence of colors with different
#'values but with the same level of saturation.
#'
#'The argument \code{fun} allows to apply a transformation to the generated
#'sequence. By default \code{fun="linear"}, no transformation is applied. If
#'\code{fun="sqrt"}, the square root of the generated sequence will be taken.
#'If \code{fun="log"}, the logarithmic of the generated sequence will be taken.
#'
#'@param color an R color name or a color in hexadeciaml notation
#'@param percentage numeric value indicating the increment steps of the
#'sequence in percentage
#'@param what character string indicating what parameter to taki into account
#'to generate the sequence. Possible values are \code{"saturation"},
#'\code{"value"}, and \code{alpha}
#'@param s optional decimal value (between 0 and 1) to fix the color saturation
#'@param v optional decimal value (between 0 and 1) to fix the color value
#'@param alpha optional decimal value (between 0 and 1) to fix the color alpha
#'transparency
#'@param fun character string indicating the applied transformation to the
#'generated sequence. Possible values are \code{"linear"}, \code{"sqrt"}, and
#'\code{"log"}
#'@param plot logical value indicating whether to plot the sequence
#'@param verbose logical value indicating whether to return the color names of
#'the sequence
#'@author Gaston Sanchez
#'@seealso \code{\link{pizza}}
#'@export
#'@examples
#'
#' # sequence for 'orange'
#' sequential("orange")
#' 
#' # sequence for 'orange' with fun='sqrt' transformation
#' sequential("orange", fun = "sqrt")
#' 
#' # sequence for 'orange' with fun='log' transformation
#' sequential("orange", fun = "log")
#' 
#' # sequential sequence for value with fix saturation s=0.7 and fun='log'
#' sequential("orange", what = "value", s = 0.7, fun = "log")
#' 
#' # sequential sequence for saturation, with fix value s=0.8, alpha=0.5, percentage 10, and fun='log'
#' sequential("orange", 10, what = "value", s = 0.7, alpha = 0.5, fun = "log")
#'
sequential <-
function(color, percentage=5, what="saturation", 
	s=NULL, v=NULL, alpha=NULL, fun="linear", plot=TRUE, verbose=TRUE)
{	
	# convert to HSV
	col_hsv = rgb2hsv(col2rgb(color))[,1]
	# transparency
	if (is.null(alpha))
		alpha = 1
	if (substr(color, 1, 1) == "#" && nchar(color) == 9)
		alpha = substr(color, 8, 9)
	# get hue, saturation, and value
	hue = col_hsv[1]
	if (is.null(s)) s = col_hsv[2]
	if (is.null(v)) v = col_hsv[3]
	# sequence function
	getseq = switch(fun, 
		linear = seq(0, 1, by=percentage/100),
		sqrt = sqrt(seq(0, 1, by=percentage/100)),
		log = log1p(seq(0, 1, by=percentage/100)),
		log10 = log10(seq(0, 1, by=percentage/100))
		)
	# what type of sequence?
	if (what == "saturation") {
		sat = getseq
		fixed = paste("v=", round(v,2), " and alpha=", alpha, sep="")
		if (is.numeric(alpha))
			seq_col = hsv(hue, s=sat, v=v, alpha=alpha)
		if (is.character(alpha)) {
			seq_col = hsv(hue, s=sat, v=v)
			seq_col = paste(seq_col, alpha, sep="")
		}
	}
	if (what == "value") {
		val = getseq
		fixed = paste("s=", round(s,2), " and alpha=", alpha, sep="")
		if (is.numeric(alpha))
			seq_col = hsv(hue, s=s, v=val, alpha=alpha)
		if (is.character(alpha)) {
			seq_col = hsv(hue, s=s, v=val)
			seq_col = paste(seq_col, alpha, sep="")
		}
	}
	if (what == "alpha") {
		alpha = getseq
		fixed = paste("s=", round(s,2), " and v=", round(v,2), sep="")
		seq_col = hsv(hue, s=s, v=v, alpha=alpha)
	}
	# if plot TRUE
	if (plot)
	{
		n = length(seq(0, 1, by=percentage/100))
		fx = unlist(fixed)
		#dev.new()
		plot(0, 0, type="n", xlim=c(0,1), ylim=c(0,1), axes=FALSE, xlab="", ylab="")
		rect(0:(n-1)/n, 0, 1:n/n, 1, col=seq_col, border="lightgray")
		mtext(seq_col, side=1, at=0.5:(n)/n, cex=0.8, las=2)
		title(paste("Sequential colors based on ", what, "\n with fixed ", fx, sep=""),
			cex.main=0.9)
	}
	# result
	if (verbose)
		seq_col
}

