isJGR <- function () {
	## Search for .jgr.works on the whoile search path, starting from GlobalEnv
	if (exists(".jgr.works", envir = .GlobalEnv, inherits = TRUE)) {
		get(".jgr.works", envir = .GlobalEnv, inherits = TRUE)
	} else FALSE  # JGR is probably not (correctly) installed
}
