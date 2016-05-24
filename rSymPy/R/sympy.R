
jythonStart <- function(jython.jar) {
	.jinit(jython.jar)
	assign(".Jython", .jnew("org.python.util.PythonInterpreter"), .GlobalEnv)
	invisible(.Jython)
}

sympyStart <- function() {

	# like system.file but on Windows uses \ in path rather than /
	system.file. <- function(...) {
		s <- system.file(...)
		if (.Platform$OS == "windows") gsub("/", "\\", s, fixed = TRUE) else s
	}

    assign(".Jython", rJython( modules = system.file( "Lib", package = "rSymPy" ) ), .GlobalEnv)

	.Jython$exec("import sys")
	.Jython$exec("from sympy import *")

}

sympy <- function(..., retclass = c("character", "Sym", "NULL"), debug = FALSE) {
	if (!exists(".Jython", .GlobalEnv)) sympyStart()
    retclass <- match.arg(retclass)
	if (retclass != "NULL") {
		.Jython$exec(paste("__Rsympy=", ...))
		if (debug) .Jython$exec("print __Rsympy") 
		Rsympy <- .Jython$get("__Rsympy")
		out <- if (!is.null(Rsympy)) .jstrVal(Rsympy)
        if (!is.null(out) && retclass == "Sym") structure(out, class = "Sym")
		else out
	} else .Jython$exec(paste(...))
}



