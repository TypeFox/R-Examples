completeDescription <- function (file, row, column, text = readLines(file), 
author = getOption("svTools.description.author"))
{	
	if (missing(text)) {
		n <- if (missing(row)) -1 else row
		rl <- readLines(file, n = n)
		row <- length(rl)
		if (missing(column)) column <- nchar(rl[row])
	} else {
		rl <- unlist(strsplit(text, "\\\n"))
		row <- length(rl)
		column <- nchar(rl[row])
	}
	rl[row] <- substring(rl[row], 1, column)
	lastLine <- rl[row]
	
	if (regexpr("^( +|[^:]+:)", lastLine) > 0) {
		## Extract the last field 
		lastField <- tail(which(regexpr("^[^:]+:", rl) > 0), 1)
		field <- gsub("(:.*$|[[:space:]]+)", "", rl[lastField])
		
		## Complete package names 
		if (field %in% c("Depends", "Suggests", "Enhances", "Imports")) {
			start <- gsub(".*[,[:space:]]", "", lastLine) 
			packages <- pkgInstalled(pattern = start)[,
				c("Package", "Title"), drop = FALSE]
			return(list(data = packages, token = start, ok = 1,
				type = "package"))
		} 
		
		## Use the "svTools.description.author" option to complete
		if (field %in% c("Author", "Maintainer")) {
			if (!is.null(author)) {
				return(list(ok = 1, data = cbind(author, ""),
					token = gsub(".*: *", "", lastLine), type = "other"))
			} else return(list(ok = 0))
		}
		
		## Possible licenses
### TODO: add 'see LICENSE' if the file exists (or make sure it exists?!)
		if (field == "License") {
			possibleLicenses <- rbind(      
				c("GPL-2",        'The "GNU General Public License" version 2'),
				c("GPL-3",        'The "GNU General Public License" version 3'),
				c("LGPL-2",       'The "GNU Library General Public License" version 2'),
				c("LGPL-2.1",     'The "GNU Lesser General Public License" version 2.1'),
				c("LGPL-3",       'The "GNU Lesser General Public License" version 3'),
				c("AGPL-3",       'The "GNU Affero General Public License" version 3'),
				c("Artistic-1.0", 'The "Artistic License" version 1.0'),
				c("Artistic-2.0", 'The "Artistic License" version 2.0'))
			return(list(ok = TRUE, data = possibleLicenses,
				token = gsub(".*: *", "", lastLine), type = "other"))
		}        
		
		## Propose today's date
		if (field == "Date") {
			data <- cbind(format(Sys.time(), "%Y-%m-%d"), "Today")
			return(list(ok = TRUE, data = data,
				token = gsub(".*: *", "", lastLine), type = "other"))
		}
		
		## Fields that are supposed to accept only yes/no values
		if (field %in% c("LazyLoad", "LazyData", "ZipData")) {
			data <- rbind(c("yes", ""), c("no", ""))
			return(list(ok = TRUE, data = data,
				token = gsub(".*: *", "", lastLine), type = "other"))
		}
		
		## Encoding... only propose most current ones, or a more exhaustive list?
		if (field == "Encoding") {
			data <- rbind(c("latin1" , ""), c("latin2" , ""), c("UTF-8"  , ""))
			return(list(ok = TRUE, data = data,
				token = gsub(".*: *", "", lastLine), type = "other"))
		}
		
		## Package type
		if (field == "Type") {
			data <- rbind(c("Package", "Usual package"),
				c("Translation", "Translation package"),
				c("Frontend", "Frontend package"))
			return(list(ok = TRUE, data = data,
				token = gsub(".*: *", "", lastLine), type = "other"))
		}
		
		## Give up
		return(list(ok = FALSE))
		
 	} else if (regexpr("[^[:alpha:]]", lastLine) > 0) {
		return(list(ok = FALSE))
	} else {
		keep <- (regexpr(lastLine, .descriptionFields[, 1]) > 0  |
			regexpr(lastLine, .descriptionFields[, 3]) > 0)
		data <- as.matrix(.descriptionFields[keep, c(1, 3), drop = FALSE])
		data[, 1] <- paste(data[, 1], ": ", sep = "")
		return(list(data = data, ok = TRUE, token = lastLine, type = "fields"))
	}
}

