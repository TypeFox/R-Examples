# Projection mode
# Functions to get and set projection mode
#
# In GGobi, the projection mode determines the type
# of plot that is displayed.  For example, a 1D ASH plot,
# or a 2D tour.  
# 
# You can see the list of available projection modes
# using the \code{\link{pmodes}} function.  This accepts
# either a GGobiDisplay, or the name of the display type.
# 
# @arguments GGobiDisplay object
# @keyword dynamic
# @alias pmode<-
#X g <- ggobi(mtcars)
#X d <- displays(g)[[1]]
#X pmode(d)
#X pmodes(d)
#X pmode(d) <- "1D Plot"
pmode <- function(x) .GGobiCall("getPModeName", x)
"pmode<-" <- function(x, value) {
  value <- match.arg(value, pmodes(x))
	.GGobiCall("setPMode", as.character(value), x)
	x
}

# Interaction mode
# Functions to get and set interaction mode
#
# In GGobi, the interaction mode determines the how
# you interact with a plot: brushing, identify etc.
# Each projection mode also has a default interaction
# mode that allows you to select variables and control
# other parameters of the display
# 
# You can see the list of available interaction modes
# using the \code{\link{imodes}} function.  This accepts
# either a GGobiDisplay, or the name of the display type.
# 
# @arguments GGobiDisplay object
# @keyword dynamic
# @alias imode<-
#X g <- ggobi(mtcars)
#X d <- displays(g)[[1]]
#X imode(d)
#X imodes(d)
#X imode(d) <- "Brush"
imode <- function(x) .GGobiCall("getIModeName", x)
"imode<-" <- function(x, value) {
  value <- match.arg(value, imodes(x))
	.GGobiCall("setIMode",  as.character(value), x)
	x
}


# List available projection/interaction modes
# These generic methods retrieves possible projection/interaction modes.
#
# The default methods retrieve a list modes for all display types,
# the character methods retrieve modes for a given display type, and
# the GGobiDisplay methods return possible modes for a given display
#
# @alias imodes
# @alias pmodes.default
# @alias pmodes.GGobiDisplay
# @alias imodes.GGobiDisplay
# @alias pmodes.character
# @alias imodes.character
# @keyword internal
pmodes <- function(x=NULL) UseMethod("pmodes", x)

pmodes.default <- function(x) {
	types <- ggobi_display_types()
	modes <- lapply(types, pmodes)
	unlist(unname(mapply(
		function(modes, type) {
			all <- append(modes, type)
			lookup <- rep(type, length=length(all))
			names(lookup) <- all
			lookup
		}, 
		modes, names(modes), SIMPLIFY=FALSE
	)))
}
pmodes.GGobiDisplay <- function(x) pmodes(class(x)[1])
pmodes.character <- function(x) .GGobiCall("getPModeNames", x)

imodes <- function(x) UseMethod("imodes", x)
imodes.GGobiDisplay <- function(x) imodes(class(x)[1])
imodes.character <- function(x) .GGobiCall("getIModeNames", x)
