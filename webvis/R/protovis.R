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

#' Provides web graphics for R by wrapping visualization API's.  Currently supports Protovis.
#'
#' \tabular{ll}{
#' Package: \tab webvis\cr
#' Type: \tab Package\cr
#' Version: \tab 0.0.1\cr
#' Date: \tab 2010-03-22\cr
#' License: \tab BSD (>= 2)\cr
#' LazyLoad: \tab no\cr
#' }
#'
#' Uses Protovis to provide web graphics for R.  Package is still under active development and shouldn't be considered
#' stable until version 0.1.  Currently uses a web browser to process JavaScript, although future version will process 
#' JavaScript directly and return the SVG output.
#' 
#' @name webvis-package
#' @aliases webvis
#' @docType package
#' @title Web graphics for R.
#' @author Shane Conway \email{shane.conway@@gmail.com}
#' @references
#' \url{http://vis.stanford.edu/protovis}
#' @keywords package
NULL


#' Add a new layers to a visualization.
#'
#' \code{+.webvis} Add a new layers to a visualization.
#'
#' @param parent The root node.
#' @param child The leaf node.
#' @return A webvis object.
#' @keywords hplot
#' @author Shane Conway \email{shane.conway@@gmail.com}
#' @references
#' \url{http://vis.stanford.edu/protovis/}
`+.webvis` <- function (parent, child) {
	# check that the parent is a "webvis" object; if not, use normal + operation
	i <- length(parent$branch)
	# this section was testing to see if I could inherit the wv object to the child object,
	# but it doesn't work if there are nested statements (e.g. c + (a + b))
	#f <- substitute(child)
	#f[["wv2"]] <- parent
	#child <- eval(f) 
	parent$branch[[i+1]] <- child 
	parent
}
#`%+%` <- function (parent, child) {
#	`+.webvis`(parent=parent, child=child)
#}

#' Create a new webvis object to store each layer of the visualization.
#'
#' \code{new.webvis} Create a new webvis object to store each layer of the visualization.
#'
#' @param width  The width in pixels.
#' @param height The height in pixels.
#' @param name The name of the visualization.
#' @param description A description of the visualization.
#' @param dataset A dataset associated with the visualization.
#' @param root The root node of the visualization (the primary root should be a panel).
#' @param branch A node layer underneath the root visualization.
#' @param render The render command for the given visualization.
#' @return A webvis object.
#' @keywords hplot
#' @author Shane Conway \email{shane.conway@@gmail.com}
#' @references
#' \url{http://vis.stanford.edu/protovis/}
#' @examples
#' new.webvis()
new.webvis <- function(name="vis", description=NULL, width=300, height=200, dataset=NULL, root=NULL, branch=list(), render=NULL) {
	wv <- list(name=name,
		description=description, 
		width=width,
		height=height,
		data=dataset,
		root=root,
		branch=branch,
		render=render)
    class(wv) <- "webvis"
	return(wv)
}

#' Convert webvis to HTML.
#'
#' \code{webvisToHTML} Convert webvis to HTML.
#'
#' @param wv A webvis.flat object (from the unfold.webvis function()).
#' @param div.id The div tag id.
#' @param html.wrap Whether to wrap the visualization in other supplied HTML.
#' @param title The title of the HTML page.
#' @param head The HTML above the webvis.
#' @param tail The HTML below the webvis.
#' @param protovis.path The path to the protovis javascript.
#' @return The HTML output
#' @keywords hplot
#' @author Shane Conway \email{shane.conway@@gmail.com}
#' @references
#' \url{http://vis.stanford.edu/protovis/}
#' @examples
#' webvisToHTML(wv=unfold.webvis(new.webvis() + pv.line(data=c(1, 1.2, 1.7, 1.5, .7, .5, .2), bottom.name="y", left.name="x", bottom.scale="linear.y.y", left.scale="linear.x.x", line.width=5, render=FALSE)))
webvisToHTML <- function(wv, div.id="id", html.wrap=TRUE, title=NULL, head=getHead(title=title, protovis.path=protovis.path), tail=getTail(), protovis.path) {
	if(!esse(protovis.path))
		if(!exists("PROTOVIS.PATH"))
			protovis.path <- "http://protovis-js.googlecode.com/svn/trunk/protovis-d3.1.js"
		else protovis.path <- PROTOVIS.PATH
	if(!is.webvis.flat(wv)) stop("webvisToHTML requires a webvis.flat object; need to call unfold.webvis first.")
	wv.html <- c(collapse("<center><div id='", div.id, "'>"), 
		"<script type='text/javascript+protovis'>",
		as.character(wv),
		"</script></div></center>")
	if(html.wrap) wv.html <- c(head, wv.html, tail)
	wv.html
}

#' A protovis mark parameter.
#'
#' \code{pv.param} A protovis mark parameter.
#'
#' @param name The name of mark parameter.
#' @param data The data used in the parameter settings.
#' @param data.name The name of the field in the dataset.
#' @param value An explicit value for the parameter.
#' @param scale Whether the value or data should be scaled.  Can be "linear", "log", or ...
#' @param scale.min The minimum scaled value (or defaults to zero) in pixels.
#' @param scale.max The maximum scaled value (or defaults to the height/width of the visualization) in pixels.
#' @param xmin The minimum x value for the scaled output.
#' @param xmax The maximum x value for the scaled output.
#' @param ymin The minimum y value for the scaled output.
#' @param ymax The maximum y value for the scaled output.
#' @param default The default value for the parameter.
#' @param quote Whether character values should be quoted.
#' @return A webvis.param object.
#' @keywords hplot
#' @author Shane Conway \email{shane.conway@@gmail.com}
#' @references
#' \url{http://vis.stanford.edu/protovis/}
#' @examples
#' pv.param(name="data", value="d")
pv.param <- function(name, data=NULL, data.name=NULL, value=NULL, scale=NULL, scale.min=NULL, scale.max=NULL, xmin=NULL, xmax=NULL, ymin=NULL, ymax=NULL, default=NULL, quote=TRUE) {
	if(missing(name)) stop("'name' is a required field for a webvis.param")
	if(!esse(data.name) && !esse(value)) return()
	param <- list(name=name, data=data, data.name=data.name, value=value, scale=scale, scale.min=scale.min, scale.max=scale.max, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, quote=quote)
	class(param) <- "webvis.param"
	param
}

#' A scaling function for protovis.
#'
#' \code{pv.scale} A scaling function for protovis.
#'
#' @param type The type is a "." separated string which determines whether the scaling is linear, log, or ..., and what parameter should be scaled.
#' @param width The width of the panel.
#' @param height The height of the panel.
#' @param data The data for the panel.
#' @param data.name The name of the variable to be scaled.
#' @param scale.min The minimum scaled value (or defaults to zero) in pixels.
#' @param scale.max The maximum scaled value (or defaults to the height/width of the visualization) in pixels.
#' @param xmin The minimum x value for the scaled output.
#' @param xmax The maximum x value for the scaled output.
#' @param ymin The minimum y value for the scaled output.
#' @param ymax The maximum y value for the scaled output.
#' @return The HTML output
#' @keywords hplot
#' @author Shane Conway \email{shane.conway@@gmail.com}
#' @references
#' \url{http://vis.stanford.edu/protovis/}
#' @examples
#' pv.scale(type="linear.value.y", width=200, height=200, data=data.frame(value=c(1:5)))
pv.scale <- function(type, width, height, data=NULL, data.name=NULL, scale.min=NULL, scale.max=NULL, xmin=NULL, xmax=NULL, ymin=NULL, ymax=NULL) {
	type <- unlist(strsplit(type, ".", fixed=TRUE))
	if(length(type) != 3) stop("scale type must be of format type.datarange.scale.range (e.g. linear.y.y)")
	if(is.null(data.name)) data.name <- type[2]
	range.max <- if(type[3] == "y") { if(is.null(ymax)) height else ymax } else { if(is.null(xmax)) width else xmax }
	range.min <- if(type[3] == "y") { if(is.null(ymin)) 0 else ymin } else { if(is.null(xmin)) 0 else xmin }
	type <- type[1]
	collapse("pv.Scale.", type, "(", 
			if(is.null(scale.min)) min(data[,data.name]) else scale.min, ", ", 
			if(is.null(scale.max)) max(data[,data.name]) else scale.max, 
			").range(", range.min, ",", range.max, ")")
}

#' Takes a parameter and a webvis object and parses them.
#'
#' \code{pv.parse} Takes a parameter and a webvis object and parses them.
#'
#' @param param A webvis param object from pv.param().
#' @param wv A webvis object.
#' @param data A dataset.
#' @return The HTML output
#' @keywords hplot
#' @author Shane Conway \email{shane.conway@@gmail.com}
#' @references
#' \url{http://vis.stanford.edu/protovis/}
#' @examples
#' pv.parse(pv.param(name="text", data.name="y"))
#' pv.parse(pv.param(name="text", data.name="y"), data=data.frame(y=1:5))
pv.parse <- function(param, wv, data) {
	if(!class(param) == "webvis.param") stop(paste("Function pv.parse expects a webvis.param input but received", class(param), "instead"))
	if(esse(data))
		if(!field.exists("x", data)) 
			data$x <- 1:length(data$y)
	if(!is.null(param$value)) if(all(param$value == "d")) param$value <- data
	if(is.null(param$data) && esse(data)) param$data <- data
	if(!is.null(param$data) && field.exists(field=param$data.name, data=param$data)) { 
		return(collapse(".", param$name, "(function(d) ", 
				if(!is.null(param$scale)) pv.scale(type=param$scale, width=wv$width, height=wv$height, data=param$data, scale.min=param$scale.min, scale.max=param$scale.max, xmin=param$xmin, xmax=param$xmax, ymin=param$ymin, ymax=param$ymax) else "", 
					"(d.", param$data.name, ")", ")"))
	} else if(!is.null(param$data.name)) {
		return(collapse(".", param$name, "(function(d) ", 
						"(d.", param$data.name, ")", ")"))
	} else if(!is.null(param$value)) {
		return(collapse(".", param$name, "(", pv.data(param$value, quote=param$quote), ")")) 
	}  
	if(!is.null(param$default)) {
		collapse(".", param$name, "(", pv.data(param$default), ")") 
	} else {
		collapse(".", param$name)
	}
}

#' A protovis panal.
#'
#' \code{pv.panel} A protovis panal.
#'
#' @param wv A webvis param object from pv.param().
#' @param data A webvis object.
#' @param width The width of the panel (in pixels).
#' @param height The height of the panel (in pixels).
#' @param left Where the panel should start with respect to the left of the window.
#' @param right Where the panel should start with respect to the right of the window.
#' @param bottom Where the panel should start with respect to the bottom of the window.
#' @param top Where the panel should start with respect to the top of the window.
#' @param anchor Whether the panel should be anchored to the parent object.
#' @return The HTML output
#' @keywords hplot
#' @author Shane Conway \email{shane.conway@@gmail.com}
#' @references
#' \url{http://vis.stanford.edu/protovis/}
#' @examples 
#' pv.panel() # the default panel size
#' pv.panel(width=NULL, height=NULL, anchor="bottom") # a panel with no settings anchored to a parent object
pv.panel <- function(wv=NULL, data, width=300, height=200, left, right, bottom, top, anchor=NULL) {
	if(!esse(wv)) { wv <- new.webvis(width=width, height=height) }
	params <- list((if(esse(data)) pv.param(name="data", value="d") else NULL), 
		pv.param(name="width", value=width),
		pv.param(name="height", value=height),
		pv.param(name="left", value=left),
		pv.param(name="right", value=right),
		pv.param(name="bottom", value=bottom),
		pv.param(name="top", value=top))
	missing.param <- (unlist(lapply(params, is.null)))
	params <- params[which(!missing.param)]
	panel <- pv.mark(wv=wv, data=data, type="Panel", params, anchor=anchor)
	if(is.null(wv$root)) wv$root <- panel
	else wv <- wv + panel
	wv
}

#' Add a dataset as a variable to the visualization.
#'
#' \code{pv.dataset} Add a dataset as a variable to the visualization.
#'
#' @param data The dataset to be used in the graphic.
#' @param name The name of the dataset.
#' @return A string of the relevant javascript.
#' @keywords hplot
#' @author Shane Conway \email{shane.conway@@gmail.com}
#' @references
#' \url{http://vis.stanford.edu/protovis/}
#' @examples
#' pv.dataset(data=data.frame(wheat=1:10), name="wheat")
pv.dataset <- function(data, name) {
	paste("var", name, "=", pv.data(data))
}

#' Generic function for all Protovis mark types.
#'
#' \code{pv.mark} Generic function for all Protovis mark types.  This function can be used
#'   to create any kind of Protovis object regardless of whether it has been exposed separately.
#'
#' @param wv A webvis object
#' @param type Can be "Line", "Bar", etc. (see Protovis API)
#' @param data A dataset for plotting.
#' @param ... Any number of pv.param objects.
#' @param anchor If anchoring to another object.
#' @return A webvis object.
#' @keywords hplot
#' @author Shane Conway \email{shane.conway@@gmail.com}
#' @references
#' \url{http://vis.stanford.edu/protovis/}
#' @examples
#' data <- data.frame(y=1:5)
#' pv.mark(wv=new.webvis(), type="Line", data=data, 
#' 		pv.param(name="data", value=data), 
#' 		pv.param(name="bottom", data.name="y", scale="linear.y.y"))
#' pv.mark(type="Label", ...=pv.param(name="text", data.name="y"))
#' pv.parse(pv.param(name="text", data.name="y"), data=data.frame(y=1:5))
pv.mark <- function(wv=NULL, type, data=NULL, ..., anchor=NULL) {
	if(!esse(data)) data <- NULL
	args <- if(length(list(...)) > 0) { if(is.webvis.param(list(...)[[1]])) list(...) else list(...)[[1]] } else list()
	vis <- list(type=collapse("pv.", type),
			parameters=collapse(
					unlist(lapply(args, function(x, data, wv) {
										#if(!is.null(x$data) && is.null(x$data.name)) pv.parse(x, wv=wv) else 
										pv.parse(param=x, wv=wv, data=data)
									}, data=data, wv=wv)),";"),
			anchor=anchor)
	vis
}

#' Used in pv.chart to correctly combine parameters from the input.
#'
#' \code{append.param} Used in pv.chart to correctly combine parameters from the input.  A "value" has a higher priority than
#'   a name in the data.
#'
#' @param paramlist A list of parameters, to be appended.
#' @param name The name of the Protovis parameter (e.g. "width").
#' @param value A specific value for the parameter (e.g. 200).
#' @param param.name The name of the field in the dataset (if not supplying specific value).
#' @param param.scale The scaling for the parameter.
#' @param scale.min The minimum scaled value (or defaults to zero) in pixels.
#' @param scale.max The maximum scaled value (or defaults to the height/width of the visualization) in pixels.
#' @param xmin The minimum x value for the scaled output.
#' @param xmax The maximum x value for the scaled output.
#' @param ymin The minimum y value for the scaled output.
#' @param ymax The maximum y value for the scaled output.
#' @return A list of parameters.
#' @keywords hplot
#' @author Shane Conway \email{shane.conway@@gmail.com}
#' @references
#' \url{http://vis.stanford.edu/protovis/}
append.param <- function(paramlist, name, value, param.name, param.scale, scale.min=NULL, scale.max=NULL, xmin, xmax, ymin, ymax) {
	if(esse(value)) {
		paramlist[[length(paramlist) + 1]] <- pv.param(name=name, value=if(name=="data") "d" else value) 
	} else if(esse(param.name)) {
		paramlist[[length(paramlist) + 1]] <- pv.param(name=name, data.name=param.name, scale=(if(esse(param.scale)) param.scale else NULL), scale.min=scale.min, scale.max=scale.max, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)
	}
	paramlist
}

#' Add a chart to the visualization.
#'
#' \code{pv.chart} Adds a chart to the visualization.  This is the core charting function which all other charting functions reference.
#' Function is generally not used directly, but is called by the higher level functions (e.g. pv.line, pv.dot).  
#'
#' @param type The type of chart (can be "Line", "Bar", "Dot", "Image", "Area") 
#' @param wv The webvis object containing the visualization. 
#' @param data The webvis object containing the visualization. 
#' @param bottom.name   A protovis mark parameter (see API documentation http://protovis-js.googlecode.com/svn/trunk/jsdoc/index.html).
#' @param bottom.scale A protovis mark parameter (see API documentation http://protovis-js.googlecode.com/svn/trunk/jsdoc/index.html).
#' @param top.name A protovis mark parameter (see API documentation http://protovis-js.googlecode.com/svn/trunk/jsdoc/index.html).
#' @param top.scale   A protovis mark parameter (see API documentation http://protovis-js.googlecode.com/svn/trunk/jsdoc/index.html).
#' @param left.name   A protovis mark parameter (see API documentation http://protovis-js.googlecode.com/svn/trunk/jsdoc/index.html).
#' @param left.scale A protovis mark parameter (see API documentation http://protovis-js.googlecode.com/svn/trunk/jsdoc/index.html).
#' @param right.name A protovis mark parameter (see API documentation http://protovis-js.googlecode.com/svn/trunk/jsdoc/index.html).
#' @param right.scale   A protovis mark parameter (see API documentation http://protovis-js.googlecode.com/svn/trunk/jsdoc/index.html).
#' @param height A protovis mark parameter (see API documentation http://protovis-js.googlecode.com/svn/trunk/jsdoc/index.html).
#' @param height.name   A protovis mark parameter (see API documentation http://protovis-js.googlecode.com/svn/trunk/jsdoc/index.html).
#' @param height.scale A protovis mark parameter (see API documentation http://protovis-js.googlecode.com/svn/trunk/jsdoc/index.html).
#' @param width   A protovis mark parameter (see API documentation http://protovis-js.googlecode.com/svn/trunk/jsdoc/index.html).
#' @param width.name A protovis mark parameter (see API documentation http://protovis-js.googlecode.com/svn/trunk/jsdoc/index.html).
#' @param width.scale   A protovis mark parameter (see API documentation http://protovis-js.googlecode.com/svn/trunk/jsdoc/index.html).
#' @param line.width.name   A protovis mark parameter (see API documentation http://protovis-js.googlecode.com/svn/trunk/jsdoc/index.html).
#' @param line.width.scale A protovis mark parameter (see API documentation http://protovis-js.googlecode.com/svn/trunk/jsdoc/index.html).
#' @param size A protovis mark parameter (see API documentation http://protovis-js.googlecode.com/svn/trunk/jsdoc/index.html).
#' @param size.name   A protovis mark parameter (see API documentation http://protovis-js.googlecode.com/svn/trunk/jsdoc/index.html).
#' @param size.scale A protovis mark parameter (see API documentation http://protovis-js.googlecode.com/svn/trunk/jsdoc/index.html).
#' @param shape   A protovis mark parameter (see API documentation http://protovis-js.googlecode.com/svn/trunk/jsdoc/index.html).
#' @param shape.name A protovis mark parameter (see API documentation http://protovis-js.googlecode.com/svn/trunk/jsdoc/index.html).
#' @param shape.scale   A protovis mark parameter (see API documentation http://protovis-js.googlecode.com/svn/trunk/jsdoc/index.html).
#' @param inner.radius A protovis mark parameter (see API documentation http://protovis-js.googlecode.com/svn/trunk/jsdoc/index.html).
#' @param inner.radius.name   A protovis mark parameter (see API documentation http://protovis-js.googlecode.com/svn/trunk/jsdoc/index.html).
#' @param inner.radius.scale A protovis mark parameter (see API documentation http://protovis-js.googlecode.com/svn/trunk/jsdoc/index.html).
#' @param outer.radius A protovis mark parameter (see API documentation http://protovis-js.googlecode.com/svn/trunk/jsdoc/index.html).
#' @param outer.radius.name   A protovis mark parameter (see API documentation http://protovis-js.googlecode.com/svn/trunk/jsdoc/index.html).
#' @param outer.radius.scale A protovis mark parameter (see API documentation http://protovis-js.googlecode.com/svn/trunk/jsdoc/index.html).
#' @param angle   A protovis mark parameter (see API documentation http://protovis-js.googlecode.com/svn/trunk/jsdoc/index.html).
#' @param angle.name A protovis mark parameter (see API documentation http://protovis-js.googlecode.com/svn/trunk/jsdoc/index.html).
#' @param angle.scale   A protovis mark parameter (see API documentation http://protovis-js.googlecode.com/svn/trunk/jsdoc/index.html).
#' @param start.angle   A protovis mark parameter (see API documentation http://protovis-js.googlecode.com/svn/trunk/jsdoc/index.html).
#' @param start.angle.name A protovis mark parameter (see API documentation http://protovis-js.googlecode.com/svn/trunk/jsdoc/index.html).
#' @param start.angle.scale   A protovis mark parameter (see API documentation http://protovis-js.googlecode.com/svn/trunk/jsdoc/index.html).
#' @param end.angle   A protovis mark parameter (see API documentation http://protovis-js.googlecode.com/svn/trunk/jsdoc/index.html).
#' @param end.angle.name A protovis mark parameter (see API documentation http://protovis-js.googlecode.com/svn/trunk/jsdoc/index.html).
#' @param end.angle.scale   A protovis mark parameter (see API documentation http://protovis-js.googlecode.com/svn/trunk/jsdoc/index.html).
#' @param text A protovis mark parameter (see API documentation http://protovis-js.googlecode.com/svn/trunk/jsdoc/index.html).
#' @param text.name   A protovis mark parameter (see API documentation http://protovis-js.googlecode.com/svn/trunk/jsdoc/index.html).
#' @param text.scale A protovis mark parameter (see API documentation http://protovis-js.googlecode.com/svn/trunk/jsdoc/index.html).
#' @param font A protovis mark parameter (see API documentation http://protovis-js.googlecode.com/svn/trunk/jsdoc/index.html).
#' @param text.style A protovis mark parameter (see API documentation http://protovis-js.googlecode.com/svn/trunk/jsdoc/index.html).
#' @param text.align A protovis mark parameter (see API documentation http://protovis-js.googlecode.com/svn/trunk/jsdoc/index.html).
#' @param text.baseline   A protovis mark parameter (see API documentation http://protovis-js.googlecode.com/svn/trunk/jsdoc/index.html).
#' @param text.margin   A protovis mark parameter (see API documentation http://protovis-js.googlecode.com/svn/trunk/jsdoc/index.html).
#' @param text.angle A protovis mark parameter (see API documentation http://protovis-js.googlecode.com/svn/trunk/jsdoc/index.html).
#' @param stroke.style.name   A protovis mark parameter (see API documentation http://protovis-js.googlecode.com/svn/trunk/jsdoc/index.html).
#' @param stroke.style.scale A protovis mark parameter (see API documentation http://protovis-js.googlecode.com/svn/trunk/jsdoc/index.html).
#' @param fill.style.name   A protovis mark parameter (see API documentation http://protovis-js.googlecode.com/svn/trunk/jsdoc/index.html).
#' @param fill.style.scale A protovis mark parameter (see API documentation http://protovis-js.googlecode.com/svn/trunk/jsdoc/index.html).
#' @param segmented   A protovis mark parameter (see API documentation http://protovis-js.googlecode.com/svn/trunk/jsdoc/index.html).
#' @param x.padding   A protovis mark parameter (see API documentation http://protovis-js.googlecode.com/svn/trunk/jsdoc/index.html).
#' @param y.padding   A protovis mark parameter (see API documentation http://protovis-js.googlecode.com/svn/trunk/jsdoc/index.html).
#' @param xmin A protovis mark parameter (see API documentation http://protovis-js.googlecode.com/svn/trunk/jsdoc/index.html).
#' @param xmax A protovis mark parameter (see API documentation http://protovis-js.googlecode.com/svn/trunk/jsdoc/index.html).
#' @param ymin A protovis mark parameter (see API documentation http://protovis-js.googlecode.com/svn/trunk/jsdoc/index.html).
#' @param ymax A protovis mark parameter (see API documentation http://protovis-js.googlecode.com/svn/trunk/jsdoc/index.html).
#' @param scale.min   A protovis mark parameter (see API documentation http://protovis-js.googlecode.com/svn/trunk/jsdoc/index.html).
#' @param scale.max   A protovis mark parameter (see API documentation http://protovis-js.googlecode.com/svn/trunk/jsdoc/index.html).
#' @param anchor A protovis mark parameter (see API documentation http://protovis-js.googlecode.com/svn/trunk/jsdoc/index.html).
#' @param render A protovis mark parameter (see API documentation http://protovis-js.googlecode.com/svn/trunk/jsdoc/index.html).
#' @param url   A protovis mark parameter (see API documentation http://protovis-js.googlecode.com/svn/trunk/jsdoc/index.html).
#' @param normalize   A protovis mark parameter (see API documentation http://protovis-js.googlecode.com/svn/trunk/jsdoc/index.html).
#' @param ...   A protovis mark parameter (see API documentation http://protovis-js.googlecode.com/svn/trunk/jsdoc/index.html).
#' @param bottom The width of the panel in pixels.
#' @param top The width of the panel in pixels.
#' @param left The width of the panel in pixels.
#' @param right The width of the panel in pixels.
#' @param line.width The width of the panel in pixels.
#' @param stroke.style The width of the panel in pixels.
#' @param fill.style The width of the panel in pixels.
#' @param interpolate The width of the panel in pixels.
#' @return A wv object.
#' @keywords hplot
#' @author Shane Conway \email{shane.conway@@gmail.com}
#' @references
#' \url{http://vis.stanford.edu/protovis/}
#' @seealso \code{\link{new.webvis}} that creates the webvis object.
pv.chart <- function(type, wv=NULL, data=NULL, bottom, bottom.name, bottom.scale, top, top.name, top.scale, left, left.name, left.scale, right, right.name, right.scale, 
		height, height.name, height.scale, width, width.name, width.scale, line.width, line.width.name, line.width.scale, 
		size, size.name, size.scale, shape, shape.name, shape.scale, inner.radius, inner.radius.name, inner.radius.scale, outer.radius, outer.radius.name, outer.radius.scale, 
		angle, angle.name, angle.scale, start.angle, start.angle.name, start.angle.scale, end.angle, end.angle.name, end.angle.scale, 
		text, text.name, text.scale, font, text.style, text.align, text.baseline, text.margin, text.angle,  
		stroke.style, stroke.style.name, stroke.style.scale, fill.style, fill.style.name, fill.style.scale, segmented, interpolate, 
		x.padding, y.padding, xmin=NULL, xmax=NULL, ymin=NULL, ymax=NULL, scale.min=NULL, scale.max=NULL, anchor=NULL, render=FALSE, 
		url, normalize=FALSE, ...) {
	if(!esse(wv)) { wv <- pv.panel() }
	if(esse(data)) {
		if(is.vector(data))
			data <- data.frame(y=data)
		if(!field.exists("x", data))
			data$x <- 1:length(data[,1])
		if(normalize)
			data[,angle.name] <- (data[,angle.name]/sum(data[,angle.name])) * 2 * pi
	}	
	# build the final parameter list
	paramlist <- list()
	paramlist <- append.param(paramlist=paramlist, name="data", value=data)
	paramlist <- append.param(paramlist=paramlist, name="bottom", value=bottom, param.name=bottom.name, param.scale=bottom.scale, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, scale.min=scale.min, scale.max=scale.max)
	paramlist <- append.param(paramlist=paramlist, name="top", value=top, param.name=top.name, param.scale=top.scale, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)
	paramlist <- append.param(paramlist=paramlist, name="right", value=right, param.name=right.name, param.scale=right.scale, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)
	paramlist <- append.param(paramlist=paramlist, name="left", value=left, param.name=left.name, param.scale=left.scale, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)
	paramlist <- append.param(paramlist=paramlist, name="height", value=height, param.name=height.name, param.scale=height.scale, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, scale.min=scale.min, scale.max=scale.max)
	paramlist <- append.param(paramlist=paramlist, name="lineWidth", value=line.width, param.name=line.width.name, param.scale=line.width.scale, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, scale.max=scale.max)
	paramlist <- append.param(paramlist=paramlist, name="width", value=width, param.name=width.name, param.scale=width.scale, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)
	paramlist <- append.param(paramlist=paramlist, name="strokeStyle", value=stroke.style, param.name=stroke.style.name, param.scale=stroke.style.scale, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)
	paramlist <- append.param(paramlist=paramlist, name="fillStyle", value=fill.style, param.name=fill.style.name, param.scale=fill.style.scale, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)
	paramlist <- append.param(paramlist=paramlist, name="interpolate", value=interpolate)
	paramlist <- append.param(paramlist=paramlist, name="segmented", value=segmented)
	paramlist <- append.param(paramlist=paramlist, name="shape", value=shape, param.name=shape.name, param.scale=shape.scale, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)
	paramlist <- append.param(paramlist=paramlist, name="size", value=size, param.name=size.name, param.scale=size.scale, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)
	paramlist <- append.param(paramlist=paramlist, name="innerRadius", value=inner.radius, param.name=inner.radius.name, param.scale=inner.radius.scale, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)
	paramlist <- append.param(paramlist=paramlist, name="outerRadius", value=outer.radius, param.name=outer.radius.name, param.scale=outer.radius.scale, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)
	paramlist <- append.param(paramlist=paramlist, name="angle", value=angle, param.name=angle.name, param.scale=angle.scale, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)
	paramlist <- append.param(paramlist=paramlist, name="startAngle", value=start.angle, param.name=start.angle.name, param.scale=start.angle.scale, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)
	paramlist <- append.param(paramlist=paramlist, name="endAngle", value=end.angle, param.name=end.angle.name, param.scale=end.angle.scale, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)
	paramlist <- append.param(paramlist=paramlist, name="text", value=text, param.name=text.name, param.scale=text.scale, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)
	paramlist <- append.param(paramlist=paramlist, name="font", value=font)
	paramlist <- append.param(paramlist=paramlist, name="textStyle", value=text.style)
	paramlist <- append.param(paramlist=paramlist, name="textAlign", value=text.align)
	paramlist <- append.param(paramlist=paramlist, name="textBaseline", value=text.baseline)
	paramlist <- append.param(paramlist=paramlist, name="textMargin", value=text.margin)
	paramlist <- append.param(paramlist=paramlist, name="textAngle", value=text.angle)
	paramlist <- append.param(paramlist=paramlist, name="url", value=url)
	# assemble the mark
	vis <- new.webvis(root=pv.mark(wv=wv, data=data, type=type, paramlist, anchor=anchor))
	if(render) { render.webvis(wv=(wv + vis)); return(wv + vis) } else vis
}

#' Unfolds the visualization tree structure into a flat form.
#'
#' \code{unfold.webvis} Unfolds the visualization tree structure into a flat form.
#'
#' @param wv The webvis object containing the visualization. 
#' @param name The name of the visualization (will show up as variables in the javascript).
#' @param parent If the node has a parent (particularly used when function is called recursively).
#' @return A wv.flat object.
#' @keywords hplot
#' @author Shane Conway \email{shane.conway@@gmail.com}
#' @references
#' \url{http://vis.stanford.edu/protovis/}
#' @seealso \code{\link{new.webvis}} that creates the webvis object.
unfold.webvis <- function(wv, name="vis", parent=NULL) {
	root <- paste("var", name, "=", if(is.null(parent)) paste("new ", wv$root$type, "()", sep="") else paste(parent, if(!is.null(wv$root$anchor)) paste(".anchor('", wv$root$anchor, "')", sep="") else "", ".add(", wv$root$type, ")", sep=""), wv$root$parameters)
	wv2 <- as.list(wv$branch)
	if(length(wv2)) {
		wv2 <- unlist(lapply(wv2, function(x) {
			x <- if(class(x) == "webvis") {
					unfold.webvis(x, name=paste(name, ceiling(runif(1, 1, 100)), sep=""), parent=name)
				 } else { 
					paste(paste(name, if(!is.null(x$anchor)) paste(".anchor('", x$anchor, "')", sep="") else "", ".add(", x$type, ")", sep=""), x$parameters, sep="") 
				 }
			return(x)
			}))
	}
	wv2 <- unlist(c(root, wv2, wv$render))
	class(wv2) <- "webvis.flat"
	wv2
}

#' Creates a theme to be used by a webvis chart object.
#'
#' \code{webvis.theme} Creates a theme to be used by a webvis chart object.  Not currently used anywhere.
#'
#' @param background The background color.
#' @param font The found used for any text.
#' @param colors The color spectrum.
#' @param grid The grid colors.
#' @return A webvis.theme object.
#' @author Shane Conway \email{shane.conway@@gmail.com}
#' @references
#' \url{http://vis.stanford.edu/protovis/}
webvis.theme <- function(background, font, colors, grid) {
	wv.t <- list(background=background,
			font=font, 
			colors=colors,
			grid=grid)
	class(wv.t) <- "webvis.theme"
	wv.t
}

#' Simplified plot function for web vis plots.
#'
#' \code{plot.webvis} Simplified plot function for web vis plots
#'
#' @param x Either the "x" axis data or all the data for the visualization (can be vector or dataset). 
#' @param y Optional, can specify the "y" axis data. 
#' @param type The type of plot.  Can be "bar", "line", "area", "pie", "dot", or "shape"
#' @param width The width of the panel in pixels.
#' @param height The width of the panel in pixels. 
#' @param add.grid Logical value for whether to add a grid. 
#' @param add.axes Whether to add x-y axes. 
#' @param scale.min Whether the y-axis should be scaled to zero or the minimum value in the data 
#' @param ... Other parmaeters for pv.chart.
#' @return Opens a plot in a browser window.
#' @keywords hplot
#' @author Shane Conway \email{shane.conway@@gmail.com}
#' @references
#' \url{http://vis.stanford.edu/protovis/}
#' @seealso \code{\link{new.webvis}} that creates the webvis object.
#' @examples
#' plot.webvis(x=c(1, 2, 1.5, 3, 1.2), type="line")
#' plot.webvis(x=c(1, 2, 1.5, 3, 1.2), type="area")
#' plot.webvis(c(1, 2, 1.5, 3, 1.2, 1.7, 2.5, 6, 5), add.grid=FALSE)
#' plot.webvis(c(1, 2, 1.5, 3, 1.2, 1.7, 2.5, 6, 5), type="area")
#' plot.webvis(c(1, 2, 1.5, 3, 1.2, 1.7, 2.5, 6, 5), type="line", scale.min=0)
#' plot.webvis(c(1, 2, 1.5, 3, 1.2, 1.7, 2.5, 6, 5), type="line", scale.min=NULL)
#' plot.webvis(x=10*rnorm(20), width=500, height=500, type="line")
#' plot.webvis(x=100*rnorm(20), y=100*rnorm(20), width=500, height=500, type="dot")
#' plot.webvis(x=c(1, 2, 1.5, 3, 1.2), type="pie")
#' plot.webvis(x=c(1, 2, 1.5, 3, 1.2), type="pie", inner.radius=80)
#' plot.webvis(x=1:5, y=c(1, 2, 1.5, 3, 1.2), type="area")
plot.webvis <- function(x, y=NULL, type="dot", width=300, height=200, add.grid=TRUE, add.axes=TRUE, scale.min=NULL, ...) {
	if(is.null(y) && is.vector(x)) {
		data <- data.frame(y=x, x=1:length(x))
	} else { data <- data.frame(x=x, y=y) }
	wv <- pv.panel(width=width+50, height=height+50, left=50, bottom=50, right=50, top=50)
	wv <- if(type=="bar") {
		wv + pv.bar(wv =wv, data=data, ...)
	} else if(type=="line") {
		wv + pv.line(wv=wv, data=data, scale.min=scale.min, ...)
	} else if(type=="dot") {
		wv + pv.dot(wv=wv, data=data, scale.min=scale.min, ...)
	} else if(type=="pie") {
		wv + pv.wedge(wv=wv, data=data, ...)
	} else if(type=="area") {
		wv + pv.area(wv=wv, data=data, scale.min=scale.min, ...)
	}
	if(type!="pie") {
		if(add.grid) {
			wv <- wv + (pv.rule(wv=wv, data=as.numeric(sprintf("%1.02f", seq(if(!esse(scale.min)) min(data$y) else scale.min, max(data$y), length.out=(height/20)))), axis="y", stroke.style="rgba(128,128,128,.2)") + pv.label(anchor="left", text.name="y"))
			wv <- wv + (pv.rule(wv=wv, data=data.frame(x=as.numeric(sprintf("%.2f", seq(min(data$x), max(data$x), length.out=(width/20))))), axis="x", stroke.style="rgba(128,128,128,.2)") + pv.label(anchor="bottom", text.name="x", text.angle=-pi/4, text.align="right"))
		}
		if(add.axes) {
			wv <- wv + pv.rule(wv=wv, bottom=0)
			wv <- wv + pv.rule(wv=wv, left=0)
		}
	}
	render.webvis(wv)
}

#' Converts R data into Protovis data.
#'
#' \code{pv.data} Converts R data into Protovis data.
#'
#' @param data An R data object. 
#' @param quote Whether characters should be quoted. 
#' @return A protovis data object.
#' @keywords hplot
#' @author Shane Conway \email{shane.conway@@gmail.com}
#' @references
#' \url{http://vis.stanford.edu/protovis/}
pv.data <- function(data, quote=FALSE) {
	if(all(data == "null")) return(data)
	if(is.character(data) && length(data)==1) if(quote) return(collapse("'", data, "'")) else data
	if(is.numeric(data) && length(data)==1) return(data)
	if(is.logical(data) && length(data)==1) return(tolower(as.character(data)))
	if(is.vector(data)) {
		data <- paste("[", 
				paste(as.matrix(data), collapse=", "), 
				"]", sep="")		
	} else if(class(data) %in% c("data.frame")) {
		data <- paste("[", paste(lapply(1:nrow(data), function(i) { k <- data[i,]; nm <- colnames(data); paste("{", paste(paste(nm, ":", k), collapse=", "), "}") }), collapse=","), "]")
	}
	return(data)
}

#' Create the final visualization from the webvis object.
#'
#' \code{render.webvis} Renders the visualization from the webvis object.
#'
#' @param wv The webvis object containing the visualization. 
#' @param vis.name The file name of the output HTML.
#' @param path The file path to the HTML file. 
#' @param file.name The file path to the HTML file. 
#' @param title The file path to the HTML file. 
#' @param protovis.path The file path to the HTML file.
#' @return The path to the output visualization.
#' @keywords hplot
#' @author Shane Conway \email{shane.conway@@gmail.com}
#' @references
#' \url{http://vis.stanford.edu/protovis/}
#' @seealso \code{\link{new.webvis}} that creates the webvis object.
render.webvis <- function(wv, vis.name=NULL, file.name, title="", protovis.path) {
	file.name <- if(exists("OUTPUT.PATH")) {
			if(!is.null(OUTPUT.PATH)) collapse(OUTPUT.PATH, vis.name, ".html") 
		} else { collapse(tempfile(), ".html") }
	if(!esse(protovis.path))
		if(!exists("PROTOVIS.PATH"))
			protovis.path <- "http://protovis-js.googlecode.com/svn/trunk/protovis-d3.1.js"
		else protovis.path <- PROTOVIS.PATH
	con=file(file.name, "w")
	if(!isOpen(con)) stop("unable to connect to output file")
	#check.webvis(wv)
	wv$render <- "vis.root.render();"
	wv.out <- webvisToHTML(unfold.webvis(wv), title=title, protovis.path=protovis.path)
	writeLines(wv.out, con=con)
	close(con)
	browseURL(url=file.name)
	return(file.name)
}
