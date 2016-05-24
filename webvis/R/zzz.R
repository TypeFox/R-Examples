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

.onLoad <- function(lib, pkg) {
	print("Welcome to the webvis pacakge.")
	PROTOVIS.PATH <<- as.character(Sys.getenv("PROTOVIS_PATH"))
	if(PROTOVIS.PATH == "") PROTOVIS.PATH <<- "http://protovis-js.googlecode.com/svn/trunk/protovis-d3.1.js"
	OUTPUT.PATH <<- as.character(Sys.getenv("WEBVIS_PATH"))
	if(OUTPUT.PATH == "") OUTPUT.PATH <<- tempdir()
}
