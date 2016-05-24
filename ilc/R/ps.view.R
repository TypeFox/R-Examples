ps.view <-
function(file, path = NULL) {
	# Change this if necessary:
	gs <- names(which(Sys.which(c("gsview32", "gsview64")) != ""))[1]
	if (is.na(gs)) warning("Gsview application is not found.
        Install the application and/or update system variable PATH")
  else {
    str <- paste(gs, paste(c(path, file), collapse='/'), sep=' ')
    system(str, wait=F, invisible=F)
  }
}
