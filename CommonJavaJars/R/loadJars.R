loadJars <- function(jars, java="J5") {
	jarsFullname <- c()	
	classes <- system.file("java", package = "CommonJavaJars", lib.loc = NULL)
	files <- list.files(classes, full.names = FALSE)
	if (java=="J5") {
		files <- grep("J6", files, TRUE, value = TRUE, invert = TRUE)
	}
	for (j in jars) {
		# Always take the newest jar per default:
		x <- sort(grep(j, files, TRUE, value = TRUE), decreasing = TRUE)
		if (length(x)==0) {
			stop(paste("No jar that matches \"",j,"\" could be found.",sep=""))
		}
		jarsFullname <- c(jarsFullname, x[1])
	}
	
	rJava::.jpackage("CommonJavaJars", jars=jarsFullname)
	return(invisible(jarsFullname))
}