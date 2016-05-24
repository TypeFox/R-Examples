#' Add subplot
#' 
#' Add subplot.
#' 
#' 
#' @param fun Graphcical function to call
#' @param x,y Coordinates
#' @param size Subplot size (inches?)
#' @param vadj Vertical adjustment
#' @param hadj Horizontal adjustment
#' @param pars Parameters to set?
#' @return Returns parameter list invisibly, side effect subplot added to
#' current plot.
#' @note Needs elaboration.
#' @seealso Called by \code{\link{geosubplot}}, calls
#' \code{\link{selectedpar}}.
#' @keywords aplot
#' @export subplot
subplot <-
function(fun, x, y, size = c(1, 1), vadj = 0.5, hadj = 0.5, pars)
{
	if(missing(fun))
		stop("missing argument \"fun\"")
	if(missing(pars)) {
		# set graphical parameters
		opars <- selectedpar()#(no.readonly = TRUE)
		# save old parameters
		par(err = -1)
		fin <- par()$fin
		# dimensions of figure, inches
		#
		#     mai doesn't deal with pty='s', etc
		#     mai <- par()$mai	# bottom, left, top, right margins, in inches
		pltin <- par("plt") * fin[c(1, 1, 2, 2)]
		mai <- c(pltin[3], pltin[1], fin[2] - pltin[4], fin[1] - pltin[
			2])
		#
		#
		usr <- par()$usr
		# limits in user units : xmin,xmax,ymin,ymax
                # uin paramter does not exist in R.  
		uin <- par()$pin/(c(usr[2]-usr[1],usr[4]-usr[3]))#par()$uin
		# inches per user units , x then y.
		if(missing(x)) if(missing(size)) {
				cat("Using function \"locator(2)\" to place opposite corners of subplot\n"
					)
				x <- locator(2)
			}
			else {
				cat("Using function \"locator(1)\" to place subplot\n"
					)
				x <- locator(1)
			}
		if(!is.null(x$x) && !is.null(x$y)) {
			y <- x$y
			x <- x$x
		}
		if(length(x) == 2 && length(y) == 2) {
			# then x,y represent corners of plot
			# reparameterize to lower left corner, size
			x <- sort(x)
			y <- sort(y)
			size[1] <- (x[2] - x[1]) * uin[1]
			size[2] <- (y[2] - y[1]) * uin[2]
			x <- x[1]
			y <- y[1]
			hadj <- 0
			vadj <- 0
		}
		if(length(x) != 1 || length(y) != 1)
			stop("length of x and y must both be same: 1 or 2")
		# convert x, y to inches from edges of plot, xi and yi
		xi <- mai[2] + (x - usr[1]) * uin[1]
		yi <- mai[1] + (y - usr[3]) * uin[2]
		hoff <- size[1] * hadj
		voff <- size[2] * vadj
		newmai <- c(yi - voff, xi - hoff, fin[2] - yi - size[2] + voff,
			fin[1] - xi - size[1] + hoff)
		newmex <- sqrt(max(size)/min(fin))
		if(any(newmai < 0))
			stop("subplot out of bounds")
		par(pty = "m", mex = newmex, mai = newmai)
	}
	else opars <- par(pars)
	# don't set new graphical parameters
	opars$new <- F
	#sometimes axes-less image plots will make it stick
	par(new = T)
	# don't erase current plot
	on.exit(par(opars))
	# don't erase current plot
	eval(fun, sys.parent(1))
	invisible(par()[names(opars)])
}

