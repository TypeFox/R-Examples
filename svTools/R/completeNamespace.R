completeNamespace <- function (line)
{	
	## export
	if (regexpr("^[[:space:]]*export[[:space:]]*\\(", line) > 0) {
		ex <- gsub(".*[(,][[:space:]]*", "", line)
### TODO: parse the source files for functions
	}
	
	## import
	if (regexpr("^[[:space:]]*import[[:space:]]*\\(", line) > 0) {
		im <- gsub(".*[(,][[:space:]]*", "", line)
		allpacks <- pkgInstalled(pattern = im)[, c("Package", "Title")] 
		return(list(data = allpacks, type = "package"))
	}
	
	## importFrom
	if (regexpr("^[[:space:]]*importFrom[[:space:]]*\\([^,]*$", line) > 0) {
		im <- gsub(".*[(][[:space:]]*", "", line)
		allpacks <- pkgInstalled(pattern = im)[, c("Package", "Title")] 
		return(list(data = allpacks, type = "package"))
	}
}
