# GGobi brushing (size, location, glyph, etc)

brushSize <- function(x, ...) UseMethod("brushSize")

brushSize.GGobi <- function(.data = 1, .gobi = ggobi_get(), units=0) {
  ## Passing negative dimensions means that we just get the current value back
  .GGobiCall("setBrushSize", as.integer(c(-1,-1)), .gobi = x) 
}

brushLocation <- function(x, ...) UseMethod("brushLocation")

brushLocation.GGobi <- function(.data=1, .gobi = ggobi_get(), units=0) {
 setBrushLocation.ggobi(as.integer(-1), as.integer(-1), .data, .gobi=.gobi)
}

"brushLocation<-" <- function(x, value) {
  x <- value[1]
  y <- value[2]
  
  tmp <- .GGobiCall("setBrushLocation", as.integer(value), .gobi = x)
  if(!is.null(tmp))
    names(tmp) <- c("x", "y")

  return(tmp) 
}

setBrushLocation.ggobi <- function(x, y, .data = 1, .gobi = ggobi_get(), update = TRUE, units = 0) {
 
}


"brushSize<-" <- function(x, value) UseMethod("brushSize")
"brushSize<-.GGobi" <- function(x, value) {
  w <- as.integer(value[1])
  h <- as.integer(value[2])
  
  tmp <- .GGobiCall("setBrushSize", as.integer(c(w,h)), .gobi = x) 
  if(!is.null(tmp))
    names(tmp) <- c("width", "height")

  return(tmp)
}

"brushColor<-" <- function(x, value) UseMethod("brushColor")
"brushColor<-.GGobi" <- function(x, value) {
  if(value < 1 || value > 65535 || is.na(value))
    stop("Color index must be positive and less than 65536.")
  
  .GGobiCall("setBrushColor", as.integer(id-1), .gobi=x) 
}

brushColor <- function(x, ...) UseMethod("brushColor")

brushColor.GGobi <- function(x) {
  .GGobiCall("getBrushColor", .gobi=x)
}

setBrushGlyph.ggobi <- function(type = -1, size = -1, .gobi = ggobi_get()) {

 if(missing(type) & missing(size))
  stop("Must specify a glyph size or type")

 if(is.character(type)) {
   type <- mapGlyphType(type)
 }

 .GGobiCall("setBrushGlyph", as.integer(c(type, size)), .gobi=.gobi)
 return(TRUE)
}

brushGlyph <- function(x, ...) UseMethod("brushGlyph")

brushGlyph.GGobi <- function(x) {
  x <- .GGobiCall("getBrushGlyph", .gobi = x)
  if(is.null(x))
    return(x)

  n <- getGlyphTypes.ggobi()

  names(x) <- c( names(n)[x[1] == n], "size")
  x
}

glyphTypes <- function() .GGobiCall("getGlyphTypes")
glyphSizes <- function() .GGobiCall("getGlyphSizes")


