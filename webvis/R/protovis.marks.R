###############################################################################
#
# webvis: An R package to create web visualizations.
# author: Shane Conway <shane.conway@gmail.com>
#
# This is released under a BSD license.
#
# Documentation was created using roxygen:
# roxygenize('webvis', roxygen.dir='webvis', copy.package=FALSE, unlink.target=FALSE)
#
###############################################################################

#' Add a line to the visualization.
#'
#' \code{pv.line} Adds a line plot to the visualization
#'
#' @param bottom.name The name of the field in the supplied data frame or vector. 
#' @param left.name The name of the field in the supplied data frame or vector. 
#' @param bottom.scale The name of the field in the supplied data frame or vector. 
#' @param left.scale The name of the field in the supplied data frame or vector. 
#' @param ... The parameters from pv.chart 
#' @return A wv object.
#' @keywords hplot
#' @author Shane Conway \email{shane.conway@@gmail.com}
#' @references
#' \url{http://vis.stanford.edu/protovis/}
#' @seealso \code{\link{pv.chart}} that creates the webvis object.
#' @examples
#' #pv.line()
#' #pv.line(anchor="top")	
#' #pv.line(data=c(1, 1.2, 1.7, 1.5, .7, .5, .2), bottom.name="y", left.name="x", bottom.scale="linear.y.y", left.scale="linear.x.x", line.width=5, render=TRUE)
#' 
#' # line example 1
#' pv.line(data=c(1, 1.2, 1.7, 1.5, .7, .5, .2), render=FALSE)
#' pv.line(data=data.frame(z=c(1, 1.2, 1.7, 1.5, .7, .5, .2)), bottom.name="z", render=TRUE)
#' pv.line(data=c(1, 1.2, 1.7, 1.5, .7, .5, .2), bottom.name="y", left.name="x", bottom.scale=NULL, left.scale="linear.x.x", render=TRUE)
#' 
#' # line example 1 (using layers)
#' wv <- pv.panel(width=150, height=150)
#' render.webvis(wv + pv.line(wv=wv, data=c(1, 1.2, 1.7, 1.5, .7, .5, .2), bottom.name="y", left.name="x", bottom.scale="linear.y.y", left.scale="linear.x.x"))
#' 
#' # line example 2 (need to make sure that it doesn't go over the edge
#' wv <- pv.panel(width=150, height=150)
#' wv <- wv + (pv.line(wv=wv, data=c(1, 1.2, 1.7, 1.5, .7, .5, .2), bottom.name="y", left.name="x", bottom.scale="linear.y.y", left.scale="linear.x.x")
#'      + pv.dot())
#' render.webvis(wv)
#' 
#' # line example 4
#' wv <- pv.panel(width=150, height=150)
#' render.webvis(wv + pv.line(wv=wv, data=c(1, 1.2, 1.7, 1.5, .7, .5, .2), top.name="y", left.name="x", top.scale="linear.y.y", left.scale="linear.x.x"))
#' 
#' # line example 5
#' wv <- pv.panel(width=150, height=150)
#' render.webvis(wv + pv.line(wv=wv, data=c(1, 1.2, 1.7, 1.5, .7, .5, .2), top.name="y", right.name="x", top.scale="linear.y.y", right.scale="linear.x.x"))
#' 
#' # line example 7
#' wv <- pv.panel(width=150, height=150)
#' render.webvis(wv + pv.line(wv=wv, data=c(1, 1.2, 1.7, 1.5, .7, .5, .2), bottom.name="y", left.name="x", bottom.scale="linear.y.y", left.scale="linear.x.x", interpolate="step-after"))
#' 
pv.line <- function(bottom.name="y", left.name="x", bottom.scale=paste("linear", bottom.name, "y", sep="."), left.scale="linear.x.x", ...) {
	vis <- pv.chart(type="Line", bottom.name=bottom.name, left.name=left.name, bottom.scale=bottom.scale, left.scale=left.scale, ...)
	vis
}

#' Add a bar to the visualization.
#'
#' \code{pv.bar} Adds a bar plot to the visualization
#'
#' @param height.name The name of the field in the supplied data frame or vector. 
#' @param left.name The name of the field in the supplied data frame or vector. 
#' @param height.scale The scale of the field in the supplied data frame or vector. 
#' @param left.scale The scale of the field in the supplied data frame or vector. 
#' @param bottom The bottom of the bar, with respect to the panel.
#' @param width The width of each bar.
#' @param width.name The name of the field in the data for the width.
#' @param xmax The max scale for the x-axis.
#' @param scale.min The minimum y-value for scaling.
#' @param ... The parameters from pv.chart 
#' @return A wv object.
#' @keywords hplot
#' @author Shane Conway \email{shane.conway@@gmail.com}
#' @references
#' \url{http://vis.stanford.edu/protovis/}
#' @seealso \code{\link{pv.chart}} that creates the webvis object.
#' @examples
#' plot.webvis(x=c(1, 2, 1.5, 3, 1.2), type="bar", scale.min=0)
#' pv.bar(data=c(1, 1.2, 1.7, 1.5, .7, .5, .2), height.name="y", left.name="x", height.scale="linear.y.y", left.scale="linear.x.x", bottom=0, width=25, render=TRUE)
#' 
#' pv.bar(data=c(1, 1.2, 1.7, 1.5, .7, .5, .2), width=20, render=TRUE)
#' pv.bar(data=data.frame(z=c(1, 1.2, 1.7, 1.5, .7, .5, .2)), height.name="z", height.scale="linear.z.y", width=20, render=TRUE)
#' pv.bar(data=c(1, 1.2, 1.7, 1.5, .7, .5, .2), height.name="y", left.name="x", height.scale=NULL, left.scale="linear.x.x", render=TRUE)
#' 
#' # bar example 1 
#' wv <- pv.panel(width=150, height=150)
#' render.webvis(wv + pv.bar(wv=wv, data=c(1, 1.2, 1.7, 1.5, .7), height.name="y", left.name="x", height.scale="linear.y.y", left.scale="linear.x.x"))
#' 
#' # bar example 2 (doesn't work properly)
#' d <- data.frame(y=c(1, 1.2, 1.7, 1.5, .7), z=c(0, 0.5, 0.9, 0.2, 0.7))
#' d <- cbind(d, k=d$y-d$z)
#' wv <- pv.panel(width=150, height=150)
#' render.webvis(wv + pv.bar(wv=wv, data=d, height=20, height.name=NULL, bottom=NULL, bottom.name="x", width=NULL, width.name="k", left.name="z", left.scale="linear.z.y", bottom.scale="linear.x.x", width.scale="linear.k.x"))
#' 
#' # bar example 3
#' wv <- pv.panel(width=150, height=150)
#' render.webvis(wv + pv.bar(wv=wv, data=c(1, 1.2, 1.7, 1.5, .7), height.name="y", left.name="x", bottom=NULL, top=0, height.scale="linear.y.y", left.scale="linear.x.x"))
#' 
pv.bar <- function(height.name="y", left.name="x", height.scale="linear.y.y", left.scale="linear.x.x", bottom=0, width=NULL, width.name=NULL, xmax=NULL, scale.min=0, ...) {
	args <- list(...)
	panel.width <- args[names(args)=="wv"]$wv$width
	n <- args[names(args)=="data"]$data
	if(is.null(width) && is.null(width.name)) if(is.data.frame(n)) width <- (panel.width/nrow(n))/1.2 else width <- (panel.width/length(n))/1.2
	if(is.null(xmax) && esse(panel.width)) xmax <- panel.width - width
	vis <- pv.chart(type="Bar", height.name=height.name, left.name=left.name, height.scale=height.scale, left.scale=left.scale, bottom=bottom, width=width, width.name=width.name, xmax=xmax, scale.min=scale.min, ...)
	vis
}


#' Add a bar to the visualization.
#'
#' \code{pv.area} Adds a bar plot to the visualization
#'
#' @param height.name The name of the field in the supplied data frame or vector. 
#' @param left.name The name of the field in the supplied data frame or vector. 
#' @param height.scale The scale of the field in the supplied data frame or vector. 
#' @param left.scale The scale of the field in the supplied data frame or vector. 
#' @param bottom The distance from the bottom of the panel.
#' @param scale.min The minimum value for the y-axis scale.
#' @param ... The parameters from pv.chart 
#' @return A wv object.
#' @keywords hplot
#' @author Shane Conway \email{shane.conway@@gmail.com}
#' @references
#' \url{http://vis.stanford.edu/protovis/}
#' @seealso \code{\link{pv.chart}} that creates the webvis object.
#' @examples
#' # line example 1 (using layers)
#' wv <- pv.panel(width=150, height=150)
#' render.webvis(wv + pv.area(wv=wv, data=c(1, 1.2, 1.7, 1.5, .7, .5, .2), bottom=0, height.name="y", left.name="x", height.scale="linear.y.y", left.scale="linear.x.x"))
#' 
#' wv <- pv.panel(width=150, height=150)
#' render.webvis(wv + pv.area(wv=wv, data=c(1, 1.2, 1.7, 1.5, .7, .5, .2)))
#' 
#' # line example 2 (need to make sure that it doesn't go over the edge
#' wv <- pv.panel(width=150, height=150)
#' wv <- wv + (pv.area(wv=wv, data=c(1, 1.2, 1.7, 1.5, .7, .5, .2), height.name="y", left.name="x", height.scale="linear.y.y", left.scale="linear.x.x"))
#' render.webvis(wv)
#' 
#' # line example 4
#' wv <- pv.panel(width=150, height=150)
#' render.webvis(wv + pv.area(wv=wv, data=c(1, 1.2, 1.7, 1.5, .7, .5, .2), top.name="y", left.name="x", top.scale="linear.y.y", left.scale="linear.x.x"))
#' 
#' # line example 5
#' wv <- pv.panel(width=150, height=150)
#' render.webvis(wv + pv.area(wv=wv, data=c(1, 1.2, 1.7, 1.5, .7, .5, .2), bottom=0, top.name="y", right.name="x", top.scale="linear.y.y", right.scale="linear.x.x"))
#' 
#' # line example 7
#' wv <- pv.panel(width=150, height=150)
#' render.webvis(wv + pv.area(wv=wv, data=c(1, 1.2, 1.7, 1.5, .7, .5, .2), bottom.name="y", left.name="x", bottom.scale="linear.y.y", left.scale="linear.x.x", interpolate="step-after"))
#' 
#' #pv.area(data=c(1, 1.2, 1.7, 1.5, .7, .5, .2), height.name="y", left.name="x", height.scale="linear.y.y", left.scale="linear.x.x", ymin=10, ymax=100, bottom=0, render=TRUE)
#' 
pv.area <- function(bottom=0, height.name="y", left.name="x", height.scale=paste("linear", height.name, "y", sep="."), left.scale=paste("linear", left.name, "x", sep="."), scale.min=0, ...) {
	vis <- pv.chart(type="Area", bottom=bottom, height.name=height.name, left.name=left.name, height.scale=height.scale, left.scale=left.scale, scale.min=scale.min, ...)
	vis
}


#' Add a wedge to the visualization (for pie charts, etc).
#'
#' \code{pv.wedge} Add a wedge to the visualization (for pie charts, etc).
#'
#' @param wv A webvis object (defaults to an empty panel). 
#' @param left Where the pie chart will be centered w.r.t. the left side of the panel. 
#' @param bottom Where the pie chart will be centered w.r.t. the bottom of the panel.
#' @param angle.name The name of the data field which contains the angle of each wedge.
#' @param inner.radius The inner radius of the chart. 
#' @param outer.radius The outer radius of the chart (defaults to the width/height). 
#' @param outer.radius.name The name of a data field for varying radius for each wedge. 
#' @param ... The parameters from pv.chart 
#' @return A wv object.
#' @keywords hplot
#' @author Shane Conway \email{shane.conway@@gmail.com}
#' @references
#' \url{http://code.google.com/p/protovis-js/wiki/PvWedge}
#' @seealso \code{\link{pv.chart}}, a more low-level charting function.
#' @examples
#' pv.wedge(data=c(1, 1.2, 1.7, 1.5, .7, .5, .2), render=TRUE)
#' pv.wedge(data=c(1, 1.2, 1.7, 1.5, .7, .5, .2), outer.radius=70, angle.name="y", render=TRUE)
#' pv.wedge(data=c(1, 1.2, 1.7, 1.5, .7, .5, .2), inner.radius=50, outer.radius=70, angle.name="y", render=TRUE)
#' pv.wedge(data=data.frame(y=c(1, 1.2, 1.7, 1.5, .7, .5, .2), rad=15*(1:7)), inner.radius=50, outer.radius.name="rad", angle.name="y", render=TRUE)
pv.wedge <- function(wv=pv.panel(), left=NULL, bottom=NULL, angle.name="y", inner.radius=NULL, outer.radius=NULL, outer.radius.name=NULL, ...) {
	panel.width <- wv$width
	panel.height <- wv$height
	if(is.null(left)) left <- panel.width/2
	if(is.null(bottom)) bottom <- panel.height/2
	if(is.null(outer.radius) && is.null(outer.radius.name)) outer.radius <- min(panel.width, panel.height)/2
	vis <- pv.chart(type="Wedge", left=left, bottom=bottom, angle.name=angle.name, inner.radius=inner.radius, outer.radius=outer.radius, outer.radius.name=outer.radius.name, ..., normalize=TRUE)
	vis
}



#' Add a dot to the visualization.
#'
#' \code{pv.dot} Add a dot to the visualization.
#'
#' @param bottom.name The name of the field in the supplied data frame or vector. 
#' @param left.name The name of the field in the supplied data frame or vector. 
#' @param bottom.scale The name of the field in the supplied data frame or vector. 
#' @param left.scale The name of the field in the supplied data frame or vector. 
#' @param ... The parameters from pv.chart 
#' @return A wv object.
#' @keywords hplot
#' @author Shane Conway \email{shane.conway@@gmail.com}
#' @references
#' \url{http://code.google.com/p/protovis-js/wiki/PvWedge}
#' @seealso \code{\link{pv.chart}}, a more low-level charting function.
#' @examples
#' wv <- pv.panel(width=150, height=150)
#' render.webvis(wv + pv.dot(wv=wv, data=c(1, 1.2, 1.7, 1.5, .7, .5, .2), ymax=140, xmax=140, xmin=10, ymin=10, scale.min=0))
pv.dot <- function(left.name="x", bottom.name="y", bottom.scale=paste("linear", bottom.name, "y", sep="."), left.scale=paste("linear", left.name, "x", sep="."), ...) {
	vis <- pv.chart(type="Dot", left.name=left.name, bottom.name=bottom.name, bottom.scale=bottom.scale, left.scale=left.scale, ...)
	vis
}

#' Add an image to the visualization.
#'
#' \code{pv.image} Add an image to the visualization.
#'
#' @param url The url for the image to be displayed. 
#' @param bottom.name The name of the field in the supplied data frame or vector. 
#' @param left.name The name of the field in the supplied data frame or vector. 
#' @param bottom.scale The name of the field in the supplied data frame or vector. 
#' @param left.scale The name of the field in the supplied data frame or vector. 
#' @param ... The parameters from pv.chart 
#' @return A wv object.
#' @keywords hplot
#' @author Shane Conway \email{shane.conway@@gmail.com}
#' @references
#' \url{http://code.google.com/p/protovis-js/wiki/PvWedge}
#' @seealso \code{\link{pv.chart}}, a more low-level charting function.
#' @examples
#' wv <- pv.panel(width=150, height=150)
#' render.webvis(wv + pv.image(url="http://vis.stanford.edu/protovis/ex/stanford.png", left.name=NULL, bottom.name=NULL))
#' pv.image(url="http://vis.stanford.edu/protovis/ex/stanford.png", width=100, height=100, left.name=NULL, bottom.name=NULL, render=TRUE)
#' pv.image(data=c(1, 1.2, 1.7, 1.5, .7, .5, .2), url="http://vis.stanford.edu/protovis/ex/stanford.png", width=50, height=50, render=TRUE)
pv.image <- function(url=NULL, left.name="x", bottom.name="y", bottom.scale=paste("linear", bottom.name, "y", sep="."), left.scale=paste("linear", left.name, "x", sep="."), ...) {
	vis <- pv.chart(type="Image", url=url, left.name=left.name, bottom.name=bottom.name, bottom.scale=bottom.scale, left.scale=left.scale, ...)
	vis
}


#' Add an rule to the visualization.
#'
#' \code{pv.rule} Add an rule to the visualization.
#'
#' @param axis can be either "x" or "y"
#' @param data The data to be used in the function.
#' @param left an exact position w.r.t. the left side of the panel
#' @param bottom an exact position w.r.t. the bottom side of the panel
#' @param bottom.name The name of the field in the supplied data frame or vector. 
#' @param left.name The name of the field in the supplied data frame or vector. 
#' @param bottom.scale The name of the field in the supplied data frame or vector. 
#' @param left.scale The name of the field in the supplied data frame or vector. 
#' @param ... The parameters from pv.chart 
#' @return A wv object.
#' @keywords hplot
#' @author Shane Conway \email{shane.conway@@gmail.com}
#' @references
#' \url{http://code.google.com/p/protovis-js/wiki/PvWedge}
#' @seealso \code{\link{pv.chart}}, a more low-level charting function.
#' @examples
#' data <- data.frame(y=c(1, 1.2, 1.7, 1.5, .7, .5, .2))
#' wv <- pv.panel(width=300, height=200, left=20, top=20, right=20, bottom=20)
#' wv <- wv + (pv.rule(wv=wv, data=data, axis="y", stroke.style="rgba(128,128,128,.2)", scale.min=0) + pv.label(anchor="right", text.name="y"))
#' wv <- wv + (pv.rule(wv=wv, data=1:10, axis="x", stroke.style="rgba(128,128,128,.2)") + pv.label(anchor="bottom", text.name="x"))
#' wv <- wv + pv.line(wv=wv, data=data, scale.min=0)
#' render.webvis(wv)
pv.rule <- function(axis=NULL, data, left, left.name, left.scale, bottom, bottom.name, bottom.scale, ...) {
	if(!is.null(axis)) {
		if(axis=="x" && !esse(left) && !esse(left.name)) {
			left.name <- "x"
			left.scale <- "linear.x.x"
		} else if(axis=="y" && !esse(bottom) && !esse(bottom.name)) {
			bottom.name <- "y"
			bottom.scale <- "linear.y.y"
		} 
	}
	vis <- pv.chart(type="Rule", data=data, left=left, left.name=left.name, left.scale=left.scale, bottom=bottom, bottom.name=bottom.name, bottom.scale=bottom.scale, ...)
	vis
}

#' Add a label to the visualization.
#'
#' \code{pv.label} Add a label to the visualization.
#'
#' @param text The text for the label.
#' @param text.name The name of the column from supplied data.
#' @param ... The parameters from pv.chart 
#' @return A wv object.
#' @keywords hplot
#' @author Shane Conway \email{shane.conway@@gmail.com}
#' @references
#' \url{http://code.google.com/p/protovis-js/wiki/PvWedge}
#' @seealso \code{\link{pv.chart}}, a more low-level charting function.
#' @examples
#' wv <- pv.panel(width=300, height=200, left=20, top=20, right=20, bottom=20)
#' wv <- wv + (pv.rule(wv=wv, data=1:10, axis="y", stroke.style="rgba(128,128,128,.2)") + pv.label(anchor="right", text.name="y"))
#' wv <- wv + (pv.rule(wv=wv, data=1:10, axis="x", stroke.style="rgba(128,128,128,.2)") + pv.label(anchor="bottom", text.name="x"))
#' wv <- wv + pv.dot(wv=wv, data=c(1, 1.2, 1.7, 1.5, .7, .5, .2), scale.min=0)
#' render.webvis(wv)
pv.label <- function(text, text.name, ...) {
	vis <- pv.chart(type="Label", text=text, text.name=text.name, ...)
	vis
}
