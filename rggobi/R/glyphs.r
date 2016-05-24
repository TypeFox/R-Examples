# "Manual" brushing
# Set colours and glyphs
# ============================================================================ 

# Map glyph type.
# Map glyph character code to number
# 
# @arguments vector of glyph character codes
# @keyword dynamic 
# @keyword internal 
mapGlyphType <- function(types) {
	if (!is.character(types)) 
    return(as.integer(types))
  stop("Character glyph types not yet supported")
	#sys.types <- getGlyphTypes.ggobi()
	#ids <- match(types, names(sys.types))
	#if (any(is.na(ids))) {
	#	stop(paste("Invalid glyph name(s):",paste(types[is.na(ids)], collapse=", ")))
	#}
	#ids
}

glyph_type <- function(x) UseMethod("glyph_type", x)
"glyph_type<-" <- function(x, value) UseMethod("glyph_type<-", x)
glyph_size <- function(x) UseMethod("glyph_size", x)
"glyph_size<-" <- function(x, value) UseMethod("glyph_size<-", x)
glyph_color <- glyph_colour <- function(x) UseMethod("glyph_colour", x)
"glyph_color<-" <- "glyph_colour<-" <- function(x, value) UseMethod("glyph_colour<-", x)

# Get glyph type.
# Get glyph type.
# 
# Glyph type refers to the shape of the glyph, one of:
# 
# \itemize{
# 	\item a filled circle
# 	\item an empty circle
# 	\item a filled square
# 	\item an empty square
# 	\item a single pixel
# 	\item a plus sign
# 	\item a cross
# }
# 
# @alias glyph_type
# @arguments GGobiData
# @keyword dynamic 
# @seealso \code{\link{glyph_type<-}}
glyph_type.GGobiData <- function(x) {
	which <- 1:nrow(x)
	.GGobiCall("getCaseGlyphs", as.integer(which - 1), x, .gobi = attr(x, "ggobi"))$type + 1
}

# Set glyph type
# Set glyph type
# 
# @alias glyph_type<-
# @arguments GGobiData
# @arguments vector of new types
# @keyword dynamic 
# @seealso \code{\link{glyph_type}}
"glyph_type<-.GGobiData" <- function(x, value) {
	which <- 1:nrow(x)
	types <- mapGlyphType(value)
	types <- rep(types - 1, length=length(which))
	sizes <- rep(as.integer(-1), length(which))
	
	.GGobiCall("setCaseGlyphs", as.integer(types), as.integer(sizes), as.integer(which - 1), x, .gobi = attr(x, "ggobi"))
	x
}


# Get glyph size
# Get glyph size
# 
# Glyph size is an integer between 1 and 6.
# 
# @alias glyph_size
# @arguments GGobiData
# @keyword dynamic 
# @seealso \code{\link{glyph_size<-}}
glyph_size.GGobiData <- function(x) {
	which <- 1:nrow(x)
	.GGobiCall("getCaseGlyphs", as.integer(which - 1), x, .gobi = attr(x, "ggobi"))$size + 1
}

# Set glyph size
# Set glyph size
# 
# Glyph size is an integer between 1 and 6.
# 
# @alias glyph_size<-
# @arguments GGobiData
# @arguments vector of new sizes
# @keyword dynamic 
# @seealso \code{\link{glyph_size}}
"glyph_size<-.GGobiData" <- function(x, value) {
	which <- 1:nrow(x)
	sizes <- rep(value - 1, length=length(which))
	types <- rep(as.integer(-1), length(which))
	
	.GGobiCall("setCaseGlyphs", as.integer(types), as.integer(sizes), as.integer(which - 1), x, .gobi = attr(x, "ggobi"))
	x
}


# Get glyph colour
# Get glyph colour
# 
# @alias glyph_colour
# @arguments GGobiData
# @alias glyph_color
# @keyword dynamic 
# @seealso \code{\link{glyph_colour<-}}
glyph_colour.GGobiData <- function(x)  {
	which <- 1:nrow(x)
	.GGobiCall("getCaseColors", as.integer(which - 1), x, .gobi = attr(x, "ggobi")) + 1
}

# Set glyph colour
# Set glyph colour
# 
# @alias glyph_colour<-
# @alias glyph_color<-
# @arguments GGobiData
# @arguments vector of new colours
# @seealso \code{\link{glyph_colour}}
# @keyword dynamic 
"glyph_colour<-.GGobiData" <- function(x, value) {
	#colors <- resolveColors.ggobi(value)
	colors <- rep(as.integer(value), length=nrow(x))
	ind <- 0:(length(colors)-1)

	if(any(colors < 1) || any(colors > 65535) || is.na(colors))
		stop("All color indices must be positive and less than 65536.")

	.GGobiCall("setCaseColors", colors, as.integer(ind), x)
	x
}



# Sampling, selection, hiding, exclusion
# ============================================================================ 

# Get and set exclusion
# Get and set the exclusion status of points
#
# @keyword internal
# @alias excluded<-
excluded <- function(x) UseMethod("excluded", x)
"excluded<-" <- function(x, value) UseMethod("excluded<-", x)


# Get excluded status
# Get the exclusion status of points.
# 
# If a point is excluded it is not drawn.
# 
# @arguments ggobiDataget
# @keyword dynamic 
# @seealso \code{\link{excluded<-}}
excluded.GGobiData <- function(x) {
 .GGobiCall("getExcludedIndices", .data = x)
}

# Set excluded status
# Set the exclusion status of points.
# 
# If a point is excluded it is not drawn.
# 
# @arguments GGobiData
# @arguments logical vector 
# @keyword dynamic 
# @seealso \code{\link{excluded}}
"excluded<-.GGobiData" <- function(x, value) {
	.GGobiCall("setExcludedIndices", rep(as.logical(value), length=nrow(x)), .data = x)
	x
}

# Get selection status
# Returns logical vector indicating if each point is under the brush
# 
# @alias selected
# @arguments GGobiData
# @arguments logical vector 
# @keyword dynamic 
selected.GGobiData <- function(x) {
	sel <- .GGobiCall("getSelectedIndices", .data=x)

	if(is.null(sel))
		logical(nrow(x))
	else {
		!is.na(match(seq(1, length = nrow(x)), sel))
	}
}
selected <- function(x) UseMethod("selected", x)

# Get shadowed status
# Get the exclusion status of points.
# 
# If a point is shadowed it is drawn in a dark gray colour, 
# behind all non-shadowed points.  It cannot be selected.
# 
# @alias shadowed
# @arguments ggobiDataget
# @keyword dynamic 
# @seealso \code{\link{shadowed<-}}
shadowed.GGobiData <- function(x) {
	.GGobiCall("getCasesHidden", .data=x)
}
shadowed <- function(x) UseMethod("shadowed", x)


# Set shadowed status
# Set the exclusion status of points.
# 
# If a point is shadowed it is drawn in a dark gray colour, 
# behind all non-shadowed points.  It cannot be selected.
# 
# @alias shadowed<-
# @arguments GGobiData
# @arguments logical vector 
# @keyword dynamic 
# @seealso \code{\link{shadowed}}
"shadowed<-.GGobiData" <- function(x, value) {
	.GGobiCall("setCasesHidden", rep(as.logical(value), length=nrow(x)), as.integer((1:nrow(x))-1), .data = x)
	x
}
"shadowed<-" <- function(x, value) UseMethod("shadowed<-", x)
