setClass("Color", representation())
setClass("RGBColor", representation("Color", "numeric"), prototype = c(red=0,green=0,blue=0))
setClass("ColorScheme",
	representation(
		colors="list",
		background="Color",
		annotations="Color",
		criticalvalue = "integer",                         
		description = "character",
		type = "integer",
		system = "integer",
		name="character"
	)
)


# Retrieve R descriptions of GGobi color schemes
# Retrieve R descriptions of GGobi color schemes
# 
# This retrieves R objects representing a collection of GGobi color schemes 
# directly from GGobi itself. It can access either a particular GGobi 
# instance's  color scheme list or the session-wide schemes.   Color schemes 
# can also be read from their XML  representation, but this allows one to query
# the actual ones that are available to a GGobi  session or instance.
#
# @arguments the GGobi instance from which to query the color schemes, defaults to active
# @value A list of the color schemes in the GGobi list. Each element is of class \code{ColorScheme}.
# @ref \references{\url{http://www.ggobi.org}}
# @seealso \link{addColorScheme}, \link{setActiveColorScheme}
# @keyword color
getColorSchemes <- function(.gobi = NULL)  {
 .GGobiCall("getColorSchemes", .gobi = .gobi)
}  


# Adds a color scheme to a GGobi instance or the global color  scheme list.
# Adds a color scheme to a GGobi instance or the global color  scheme list.
#
# This takes an R color scheme description and adds an equivalent representation
# to either the color schemes available to a particular GGobi instance or to the global 
# list of color schemes within a GGobi session.
#
# @arguments the \code{ColorScheme object. }
# @arguments the name to identify/describe the scheme.
# @arguments a logical value indicating whether we should replace any existing entry with this name  with the new scheme (\code{TRUE}), or simply append the new scheme to the end of the list.
#@arguments a logical value that indicates whether to also make this the active color scheme. See \code{setActiveColorScheme}
# @arguments the GGobi instance to which to append the new color scheme, or if \code{NULL}, the color scheme is added to the session options and will be available to new GGobi instances created after this.
# @value The index/position in the appropriate list of color schemes to which the new scheme was added.
# @keyword color
addColorScheme <- function(scheme, name = scheme@name, overwrite = FALSE, set = FALSE, .gobi = NULL)   {
  ans <- .GGobiCall("addColorScheme", scheme, as.character(name), as.logical(overwrite), .gobi = .gobi)

  if (set) setActiveColorScheme(name, .gobi = .gobi)

  ans
}  

# Map names to color values in a ggobi instance.
# Map names to color values in a ggobi instance.
#
# This currently does not work but will when color schemes are enabled in GGobi.
#
# Colors in ggobi are represented by entries in a color table or matrix as 
# Red, Green and Blue values. To set a color, one can identify the entry in that 
# table either by number/index or by name. This function maps color names to the 
# corresponding index so that these values can be passed to the low-level C routines.  
# This function is typically called internally.
#
# @arguments names of colors to be resolved in the specified color scheme of the ggobi instance.
# @arguments the name of color scheme of interest in which to resolve the color names. 
# @arguments the ggobi instance whose color schemes are to be used to resolve the \code{scheme},  if \code{NULL}, the color schemes from the GGobi session options are used.
# @value An integer vector specifying the indices of the entries in the color scheme corresponding to the specified names.
# @keyword color
resolveColors.ggobi <- function(colNames, scheme) {
	if (is.numeric(colNames)) return(as.integer(colNames))

	if(missing(scheme) || as.character(scheme) == "") {
		scheme = getActiveColorScheme()		
	} else {
		scheme = getColorSchemes()[[as.character(scheme)]]
	}

	match(as.character(colNames), names(scheme@colors))
}
