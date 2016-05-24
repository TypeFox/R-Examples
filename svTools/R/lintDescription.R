lintDescription <- function (descfile, txt = readLines(descfile))
{  
	txt <- unlist(strsplit(txt, "\\\n"))
	resetErrors(file = descfile)   
	
	addErr <- function (line, message, file = descfile)
		addError(line = line, message = message, file = file, type = "error")
		
	addWarn <- function (line, message, file = descfile)
		addError(line = line, message = message, file = file, type = "warning")
  
	## Check mandatory fields
	for (mandatory in c("Package", "Version", "License", "Description", "Title",
		"Author", "Maintainer"))
		if (!any(regexpr(sprintf("^%s", mandatory), txt) > 0))
			addErr(line = 1,
				message = sprintf("field `%s` is mandatory", mandatory))
	
	## Check the fields
	fields <- txt[regexpr("^[^:]+:", txt) > 0]
	fields <- gsub("[[:space:]]*:.*$", "", fields)
	okFields <- fields %in% .descriptionFields[, 1]
	if (!all(okFields)) {
		wrongFields <- fields[!okFields]
		lapply(wrongFields, function (x) {
			rx.out <- regexpr(sprintf("^%s *:", x), txt)
			line <- which(rx.out != -1)
			addWarn(line = line, message = sprintf("unknown field : `%s`", x))
		})
	}
	
	## Check the package name
	package <- grep("^Package[[:space:]]*:", txt)
	if (length(package)) {
		if (length(package) > 1)
			addErr(line = package[2], message = "multiple `Package` declarations")
		packageName <- sub("(^[^:]*: *| )", "", txt[package[1]])
		if (regexpr("^[a-zA-Z][\\.0-9a-zA-Z]*$", packageName) < 1)
			addErr(line = package, message = "wrong package name")
	}
  
	## Check the version
	version <- grep("^Version[[:space:]]*:", txt)
	if (length(version)) {
		if (length(version) > 1)
			addErr(line = version[2], message = "multiple `Version` declarations")
		versionNumber <- sub("(^[^:]*:| )", "", txt[version[1]])
### TODO: handle translation packages
		if (regexpr("[^0-9\\.-]", "", versionNumber) > 0)
			addWarn(line = version,
				message = "wrong format for the version number")
		nfields <- length(unlist(strsplit(versionNumber, "[-\\.]"))) 
		if (nfields < 2)
			addWarn(line = version,
				message = "wrong version number, need at least two fields")
		if (nfields > 3)
			addWarn(line = version,
				message = "wrong version number, too many fields")
	}
  
	## Check maintainer
	maintainer <- grep("^Maintainer[[:space:]]*:", txt)
	if (length(maintainer)) {
		if (length(maintainer) > 1)
			addErr(line = maintainer[2], message = "multiple `Maintainer` declarations")
		maintainerLine <- sub("^[^:]*: *", "", txt[maintainer[1]])
		if (regexpr("[\\.,]$", maintainerLine) > 0)
			addErr(line = maintainer[1],
				message = "the `Maintainer field should not end with a period or commas")
    }
    if (length(unlist(strsplit(maintainerLine, "@"))) != 2) {
		addErr(line = maintainer[1],
			message = "only one email adress in `Maintainer` field")
    }
    email <- gsub("(^[^<]*<|>[^>]*$)", "", maintainerLine)
    if (regexpr("[^@]+@[^@]+", email) < 1 || regexpr("[[:space:]]", email) > 0)
		addErr(line = maintainer[1], 
			message = paste("wrong `Maintainer` email adress: '", email, "'", sep = ""))
  
	## Check date                              
	date <- grep("^Date[[:space:]]*:", txt)
	if (length(date)) {
		if (length(date) > 1)
			addErr(line = date[2], message = "multiple `Date` declarations")
		dateLine <- gsub("(^[^:]*:| )", "", txt[date[1]]) 
		if (regexpr("^[0-9]{4}-[0-9]{1,2}-[0-9]{1,2}$", dateLine) < 1)
			addErr(line = date[1],
				message = "the `Date`field should be in format yyyy-mm-dd")
	}
	
### TODO: check OS_type that can only be unix or windows
  
	## Check the dependencies
### FIXME: all the stuff below comes from tools, I need to figure out what to do with it
	toolsNS <- getNamespace("tools")
	db <- toolsNS$.read_description(descfile)
	depends <- toolsNS$.get_requires_from_package_db(db, "Depends")
	imports <- toolsNS$.get_requires_from_package_db(db, "Imports")
	suggests <- toolsNS$.get_requires_from_package_db(db, "Suggests")
	reqs <- unique(c(depends, imports,
		if (!identical(as.logical(Sys.getenv("_R_CHECK_FORCE_SUGGESTS_")), 
		FALSE)) suggests))
	installed <- character(0)
	for (lib in .libPaths()) {
		pkgs <- list.files(lib)
		pkgs <- pkgs[file.access(file.path(lib, pkgs, "DESCRIPTION"), 4) == 0]
		installed <- c(pkgs, installed)
	}
	installed <- sub("_.*", "", installed)
	reqs <- reqs[!reqs %in% installed]
	stdPkgNames <- toolsNS$.get_standard_package_names()
	m <- reqs %in% stdPkgNames$stubs
	if (length(reqs[!m])) 
		addErr(line = grep("^(Depends|Suggests|Enhances)", txt), 
			message = paste("package `", reqs[!m],
			"` required but not installed", sep = ""))
	if (length(reqs[m])) 
		addErr(line = grep("^(Depends|Suggests|Enhances)", txt), 
			message = paste("package `", reqs[m],
			"` required but stub", sep = ""))
	
	return(getErrors(file = descfile))	
}
