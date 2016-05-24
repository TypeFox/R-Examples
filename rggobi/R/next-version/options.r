

# Get ggobi options
# @keyword internal 
goptions <- function(x) UseMethod("goptions", x)

# Get ggobi options
# @keyword internal 
goptions.GGobiDisplay <- function(x) {
  .GGobiCall("getDisplayOptions", x)
}
goptions.GGobi <- function(x) {
  goptions.GGobiDisplay(NULL)
}

# @keyword internal 
"goptions<-" <- function(x, value) UseMethod("goptions<-", x)

# Set defaults for display options
#
"goptions<-.GGobi" <- function(x, value) {
	goptions.GGobiDisplay(NULL) <- value
}

# Changes the settings of the display options for a display.
# These options correspond to options available in the display menus in GGobi.
# 
# Possible options:
# 
#  * "Show Points"
#  * "Show axes"
#  * "Show tour axes"
#  * "Show axes labels"
#  * "Directed edges"
#  * "Undirected edges"
#  * "Arrowheads"
#  * "Show whiskers" 
# 
# @arguments GGobiDisplay to modify
# @arguments list of new settings
# @keyword manip 
"goptions<-.GGobiDisplay" <- function(x, value) {
	old <- goptions(x)
	cur <- update(old, value)
	.GGobiCall("setDisplayOptions", cur, x)

	x
}
