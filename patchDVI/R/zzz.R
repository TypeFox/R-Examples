JSweave <- function(file, ...) {
  # Check for commands
  cmds <- Sys.which(c("uplatex", "dvipdfmx"))
  if (any(nchar(cmds) == 0)) {
    warning(gettextf("Vignette '%s' requires 'uplatex' and 'dvipdfmx'", basename(file)))
    return()
  }
  
  # Need two runs to make TOC.  Skip dvipdfm on the first run,
  # skip Sweave on the second run
  SweaveDVIPDFM(file, latex = "uplatex", dvipdfm = "echo", 
    	          encoding = "UTF-8")
  tex <- sub("[.][RrSs](nw|tex)$", ".tex", file)
  SweaveDVIPDFM(tex, latex = "uplatex", dvipdfm = "dvipdfmx", 
		  encoding = "UTF-8")
}

.onLoad <- function(libname, pkgname) {
  vignetteEngine("JSweave", weave = JSweave, tangle = Stangle)
}
