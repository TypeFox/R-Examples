#' @title Convert 'ijroi' and 'ijzip' objects to spatstat spatial patterns
#' @description Converts \code{\link[=read.ijroi]{ijroi}} and \code{\link[=read.ijzip]{ijzip}} objects to a list of \link[spatstat]{spatstat} spatial patterns.
#' @param X \code{\link[=read.ijroi]{ijroi}} or \code{\link[=read.ijzip]{ijzip}} object to be converted.
#' @param window the \link[=owin]{window} for returned spatial patterns. Can be an \code{\link{owin}} object defining a common window for all returned patterns, a character string \code{'range'} leading to a common window based \code{\link{range}} of all returned patterns, or \code{NULL} (default) leading to separate windows for each pattern.
#' @param pattern.type a character string specifying the desired pattern type to be returned (\code{\link[=ppp]{"ppp"}}, \code{\link[=psp]{"psp"}} or \code{\link[=owin]{"owin"}}). Works only if \code{X} is an 'ijroi' object. Ignored otherwise. Defaults to an appropriate pattern type depending on the ROI type (see 'Details').
#' @param unitname Name of the unit of length for the resulting window(s) (see \code{\link[=owin]{owin}}).
#' @param scale A numeric value defining the scale of photograph in pixels / \code{unitname}. Defaults to 1.
#' @param return.type should the type of ROI object(s) be returned in addition to spatstat spatial patterns? Defaults to \code{FALSE}.
#' @param convert.only a character vector specifying the \code{strType} of ROI objects to be converted (see \code{\link{plot.ijroi}} for possible pattern types). Pattern types not mentioned will not be converted. Works only if \code{X} is an 'ijzip' object. Ignored otherwise.
#' @return Returns a list of \link[spatstat]{spatstat} patterns of approperiate type (see 'Details'). If \code{return.type = TRUE} returns a list with two levels specifying the spatstat pattern and the ROI type.
#' @details The function converts \code{\link[=read.ijroi]{ijroi}} and \code{\link[=read.ijzip]{ijzip}} objects to \link[spatstat]{spatstat} spatial patterns for further calculations with the objects. By default, areal types ("rect", "oval", "ELLIPSE", "polygon") are converted to \code{\link{owin}} objects. Line types ("line" (including "ARROW"), "freeline", "polyline", "angle", "freehand" (excluding "ELLIPSE")) are converted to \code{\link{psp}} objects and "point" types to \code{\link{ppp}} objects.
#' @author Mikko Vihtakari
#' @seealso \code{\link{read.ijroi}} \code{\link{read.ijzip}} 
#' @examples file <- file.path(system.file(package = "RImageJROI"), "extdata", "ijroi", "ijzip.zip")
#' x <- read.ijzip(file)
#' ij2spatstat(x)
#' @importFrom spatstat ppp psp owin as.psp as.ppp marks marks<-
#' @export

ij2spatstat <- function(X, window = NULL, pattern.type = NULL, unitname = NULL, scale = 1, return.type = FALSE, convert.only = NULL){
#convert.only = NULL; window = NULL; pattern.type = NULL; unitname = NULL; scale = 1; return.type = FALSE

conv.fun <- function(x, ...) {
  
  type <- x$strType
  scale.elements <- c("top", "left", "bottom", "right", "width", "height", "x1", "y1", "x2", "y2", "coords", "xrange", "yrange")
  x[names(x) %in% scale.elements] <- lapply(x[names(x) %in% scale.elements], function(k) k/scale)

  if(type == "polygon") {
    out <- owin(poly = list(x = x$coords[,1], y = x$coords[,2]), unitname = unitname)
    if(!is.null(pattern.type)){
    if(pattern.type %in% c("psp", "ppp")) {
      out <- as.psp(out, window = window)
        marks(out) <- factor(x$name)}
    if(pattern.type == "ppp") {
      out <- as.ppp(out, window = window)
    }}
  }
  
  if(type == "rect") {
    out <- owin(xrange = c(x$left, x$right), yrange = c(x$top, x$bottom), unitname = unitname)
    if(!is.null(pattern.type)){
    if(pattern.type %in% c("psp", "ppp")) {
      out <- as.psp(out, window = window)
        marks(out) <- factor(x$name)}
    if(pattern.type == "ppp") {
      out <- as.ppp(out, window = window)
    }}
  }
  
  if(type == "oval") {
    theta <- seq(0, 2*pi, len=360)
    out <- owin(poly = list(x = rev(x$left + x$width/2*(1 + sin(theta))), y = rev(x$top + (x$height)/2*(1 + cos(theta)))), unitname = unitname)
    if(!is.null(pattern.type)){
    if(pattern.type %in% c("psp", "ppp")) {
      out <- as.psp(out, window = window)
        marks(out) <- factor(x$name)}
    if(pattern.type == "ppp") {
      out <- as.ppp(out, window = window)
    }}
  }
  
  if(type == "line") {
    out <- psp(x0 = x$x1, x1 = x$x2, y0 = x$y1, y1 = x$y2, marks = factor(x$name), 
            window = if(is.null(window)) {owin(xrange = x$xrange, yrange = x$yrange)} else window)
    if(!is.null(pattern.type)){
    if(pattern.type == "owin") stop(paste0("ROI types ", "'", type, "'", " cannot be assigned as windows (see ?owin)"))
    if(pattern.type == "ppp") {
      out <- as.ppp(out)
    }}
  }
  
  if(type == "freeline") {
    xx <- x$coords[,1]
    yy <- x$coords[,2]
    out <- psp(x0 = xx[1:(length(xx)-1)], x1 = xx[2:length(xx)], y0 = yy[1:(length(yy)-1)], y1 = yy[2:length(yy)],
      marks = factor(rep(x$name, length(xx)-1)), 
      window = if(is.null(window)) {owin(xrange = x$xrange, yrange = x$yrange)} else window)
    if(!is.null(pattern.type)){
    if(pattern.type == "owin") stop(paste0("ROI types ", "'", type, "'", " cannot be assigned as windows (see ?owin)"))
    if(pattern.type == "ppp") {
      out <- as.ppp(out)
    }}
  }
  
  if(type == "polyline") {
    xx <- x$coords[,1]
    yy <- x$coords[,2]
    out <- psp(x0 = xx[1:(length(xx)-1)], x1 = xx[2:length(xx)], y0 = yy[1:(length(yy)-1)], y1 = yy[2:length(yy)],
      marks = factor(rep(x$name, length(xx)-1)), 
      window = if(is.null(window)) {owin(xrange = x$xrange, yrange = x$yrange)} else window)
    if(!is.null(pattern.type)){
    if(pattern.type == "owin") stop(paste0("ROI types ", "'", type, "'", " cannot be assigned as windows (see ?owin)"))
    if(pattern.type == "ppp") {
      out <- as.ppp(out)
    }}
  }
  
  if(type == "noRoi") stop("Nothing to convert. Remove 'noRoi' objects.")
  
  if(type == "freehand") {
    if(exists("strSubtype", where = x)) {
       if(x$strSubtype == "ELLIPSE") {
        centerX <- (x$x1 + x$x2)/2
        centerY <- (x$y1 + x$y2)/2
        theta <- seq(0, 2*pi, len=360)
        dx <- x$x2 - x$x1
        dy <- x$y2 - x$y1
        major <- sqrt(dx^2 + dy^2)
        minor <- major*x$aspectRatio
        a <- major/2
        b <- minor/2
        phi <- atan2(dy, dx)
        ellipX <- centerX + a*cos(theta)*cos(phi) - b*sin(theta)*sin(phi)
        ellipY <- centerY + a*cos(theta)*sin(phi) + b*sin(theta)*cos(phi)
        out <- owin(poly = list(x = ellipX, y = ellipY), unitname = unitname)
         if(!is.null(pattern.type)){
          if(pattern.type %in% c("psp", "ppp")) {
            out <- as.psp(out, window = window)
            marks(out) <- factor(x$name)
            }
              if(pattern.type == "ppp") {
              out <- as.ppp(out, window = window)
              }}
  }} else {
    if(type == "freehand") {
      out <- owin(poly = list(x = rev(x$coords[,1]), y = rev(x$coords[,2])), unitname = unitname)
      if(!is.null(pattern.type)){
          if(pattern.type %in% c("psp", "ppp")) {
            out <- as.psp(out, window = window)
            marks(out) <- factor(x$name)
            }
              if(pattern.type == "ppp") {
              out <- as.ppp(out, window = window)
              }}
  }}
  }
  
  if(type == "angle") {
    xx <- x$coords[,1]
    yy <- x$coords[,2]
    out <- psp(x0 = xx[1:(length(xx)-1)], x1 = xx[2:length(xx)], y0 = yy[1:(length(yy)-1)], y1 = yy[2:length(yy)],
      marks = factor(rep(x$name, length(xx)-1)), 
      window = if(is.null(window)) {owin(xrange = x$xrange, yrange = x$yrange)} else window)
    if(!is.null(pattern.type)){
    if(pattern.type == "owin") stop(paste0("ROI types ", "'", type, "'", " cannot be assigned as windows (see ?owin)"))
    if(pattern.type == "ppp") {
      out <- as.ppp(out)
    }}
  }
  
  if(type == "point") {
    out <- ppp(x = x$coords[,1], y = x$coords[,2], marks = ordered(seq_along(x$coords[,1])),
      window = if(is.null(window)) {owin(xrange = x$xrange, yrange = x$yrange)} else window)
    if(!is.null(pattern.type)){
    if(pattern.type == "owin") stop(paste0("ROI types ", "'", type, "'", " cannot be assigned as windows (see ?owin)"))
    if(pattern.type == "psp") {
    xx <- x$coords[,1]
    yy <- x$coords[,2]
    out <- psp(x0 = xx[1:(length(xx)-1)], x1 = xx[2:length(xx)], y0 = yy[1:(length(yy)-1)], y1 = yy[2:length(yy)],
      marks = factor(rep(x$name, length(xx)-1)), 
      window = if(is.null(window)) {owin(xrange = x$xrange, yrange = x$yrange)} else window)
    }}
  }
  
  if(return.type) {
    return(list(pp = out, type = type))} else {
  return(out)}
}
  
if(class(X) == "ijroi"){
  if(!is.null(window)) if(window == "range") window <- owin(xrange = X$xrange, yrange = X$yrange, unitname = unitname)
  tmp <- conv.fun(X, window = window, pattern.type = pattern.type, unitname = unitname, scale = scale, return.type = return.type)
  return(tmp)
}
if(class(X) == "ijzip"){
  if(!is.null(convert.only)){
    X <- X[unlist(lapply(X, function(k) k$strType %in% convert.only))]
    }
  if(!is.null(window)) if(window == "range") {
    Xrange <- range(unlist(lapply(X, function(k) k$xrange)))
    Yrange <- range(unlist(lapply(X, function(k) k$yrange)))
    window <- owin(xrange = Xrange, yrange = Yrange, unitname = unitname)
  }
  lapply(X, function(k) 
        conv.fun(k, window = window, unitname = unitname, scale = scale, return.type = return.type)
    )  
}
  
}
