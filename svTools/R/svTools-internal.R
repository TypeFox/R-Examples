svTools.env <- new.env()

.onAttach <- function (libname, pkgname)
	addError(emptyError())

## TODO: eliminate the descriptionFields data!
.descriptionFields <- data.frame(
	field = c("Package", "Version", "License", "Description", "Title", "Author",
		"Maintainer", "Date", "Depends", "URL", "SystemRequirements", "Imports",
		"Collate", "LazyLoad", "LazyData", "ZipData", "Encoding", "Type",
		"LinkingTo", "Suggests", "Enhances", "Language", "OS_type",
		"BugReports", "Classification/ACM", "Classification/JEL",
		"Classification/MSC"),
	optional = c(rep("Mandatory", 7), rep("Optional", 20)),
	description = c("Name of the package", "Version number (x.y-z)", "License",
		"Comprehensive description of the package", "Short (1 line) title",
		"Author(s) of the package", "Maintainer (John Doe <john@doe.com>)",
		"Date (yyyy-mm-dd)", "Dependencies of this package", "Homepage",
		"System Requirements", "Packages that are imported", "Order of files",
		"Use the lazy loading mechanism", "Use the lazy loading of data",
		"Zip the data", "Encoding of this file",
		"Type of the package (Package|Frontend|Translation)",
		"The C code of this package links to ...",
		"Suggested (optional) dependencies", "Package(s) enhanced by this one",
		"Language used in this package", "Type of OS (unix|windows)",
		"Where to report bugs?",
		"Classification according to the Association of Computing Machinery",
		"Classification according to the Journal of Economic Litterature",
		"Classification according to the American Mathematical Society")
)

.looksLikeAFunction <- function (p)
{
	# Sometimes, p is not subsettable => use try here
	p1 <- try(p[[1]], silent = TRUE)
	if (inherits(p1, "try-error")) return(FALSE)
	if (length(p1) != 1) return(FALSE)
	if (!as.character(p1) %in% c("<-", "<<-", "=")) return(FALSE)
	if (length(p) <= 2) return(FALSE) 
	if (is.null(p[[3]])) return(FALSE)
	if (length(p[[3]]) == 1) return(FALSE)
	asc <- as.character(p[[3]][[1]])
	if (length(asc) > 1 || asc != "function") return(FALSE)
	return(TRUE)
}

.looksLikeAnIf <- function (p)
{
	# Sometimes, p is not subsettable => use try here
	p1 <- try(p[[1]], silent = TRUE)
	if(length(p1) != 1) return(FALSE)
	return(as.character(p1) == "if")
}

.getIfSrcRef <- function (p)
{
	x <- attr(p, "srcref")
	if (is.null(x)) {
		x <- attr(p[[3]], "srcref")
		if (length(p) == 4)
			x <- append(x, attr(p[[4]], "srcref"))
	}
	y <- lapply(x, as.integer)
	srcref <- c(head(y, 1)[[1]][1:2], tail(y, 1)[[1]][3:4])
	data.frame(srcref1 = srcref[1], srcref2 = srcref[2], srcref3 = srcref[3], 
	  srcref4 = srcref[4], stringsAsFactors = FALSE)
}

.addIfNode <- function (value = TRUE, env = env, parent, nextnode)
{	
	data <- env[["data"]]
	srcref <- attr(nextnode, "srcref")
	if (is.null(srcref)) {
		if (!.looksLikeAnIf(nextnode)) {
			return(parent)
		} else {          
			srcref <- attr(nextnode[[3]], "srcref")
			if (length(nextnode) == 4)
				srcref <- append(srcref, attr(nextnode[[3]], "srcref"))
		}
	}
	id <- max(data$id) + 1
	lap.out <- lapply(srcref, as.integer)
	srcref <- t(c(head(lap.out, 1)[[1]][1:2], tail(lap.out, 1)[[1]][3:4])) 
	colnames(srcref) <- paste("srcref", 1:4, sep = "") 
	mode <- paste("if", value, sep = ":")
	description <- mode
	env[["data"]] <- rbind(env[["data"]], data.frame(id = id, srcref,
		description = description, mode = mode, parent = parent))
	return(id)
}

.as.characterSrcRef <- function (x, useSource = TRUE, encoding = "unknown")
{
    srcfile <- attr(x, "srcfile")
    
	getSrcFileLines <- function (srcfile, first, last, encoding = "unknown") {
		if (first > last) return(character(0))
		lines <- tail(readLines(srcfile, n = last, warn = FALSE, encoding = encoding), -(first - 1))
		return(lines)
	}
	
	if (isTRUE(useSource))
        lines <- try(getSrcFileLines(srcfile, x[1], x[3], encoding = encoding), TRUE)
    if (!isTRUE(useSource) || inherits(lines, "try-error")) {
        lines <- paste("<srcref: file \"", srcfile$filename,
            "\" chars ", x[1], ":", x[2], " to ", x[3], ":",
            x[4], ">", sep = "")
	} else {
		if (length(lines) < x[3] - x[1] + 1) x[4] <- .Machine$integer.max
        lines[length(lines)] <- substring(lines[length(lines)], 1, x[4])
        lines[1] <- substring(lines[1], x[2])
    }
    return(lines)
}

### Note: this is apparently not used elsewhere?!
.dump <- function (data, id = 0, level = 0)
{
	offset <- paste(rep("\t", level), collapse = "") 
	ids <- data$id[data$parent == id]
	if (length(ids)) {
		for (i in 1:length(ids)) {
			description <- data$description[ids[i]]
			if (description != "{") cat(offset, description, "\n")
			.dump(data, id = ids[i], level = level + 1)
		}
	}
}

### Get arguments of the function defined after a given point
.argsFunAfter <- function (file, row, all.args = FALSE, p.out = parse(file),
target = NULL)
{
	positions <- sapply(attr(p.out, "srcref"), function (x) as.numeric(x)[1])
	target <- if (!is.null(target)) target else which(positions >= row)[1]
	funstart <- positions[target]
	  
	## Find the arguments that have already been documented
	definedPars <- try(if (!isTRUE(all.args) && funstart > 1) { 
		rl <- readLines(file, n = funstart - 1)
		test <- regexpr("#'", rl) > 0 
		definedPars <- if (tail(test, 1)) {
			pos <- which(!rev(test))[1] - 1
			last <- if(is.na(pos)) funstart else funstart - pos
			oxBlock <- rl[last:(funstart - 1)]
			oxBlock <- oxBlock[regexpr("^.*@param[[:space:]]+", oxBlock) > 0]
			oxBlock <- gsub("^.*@param[[:space:]]+", "", oxBlock)
			oxBlock <- gsub("[[:space:]].*$", "", oxBlock)
			oxBlock <- oxBlock[regexpr(".", oxBlock) > 0]
			oxBlock
		} 
	}, silent = TRUE)
	if (inherits(definedPars, "try-error")) definedPars <- NULL
	
	## Find the arguments of the function
	chunk <- p.out[[target]]
	env <- new.env()
	if (length(chunk) == 3) {
		if (chunk[[1]] == "<-" || chunk[[1]] == "=") {    
			eval(chunk, envir = env)
			contents <- ls(env)
			if (length(contents)) {
				object <- env[[contents]]
				if (inherits(object, "function")) {
					allpars <- names(formals(object))
					return(setdiff(allpars, definedPars))
				}
			}
		}
	}
	return(invisible(NULL))
}

### Note: not used elsewhere?
.isAfter <- function (p.out, row, col)
{
	srcref <- attr(p.out, "srcref")
	rows <- sapply(srcref, "[", 3)
	cols <- sapply(srcref, "[", 4)
	return(row > rows | (rows == row & col > cols))
}

.isBefore <- function (p.out, row, col)
{
	srcref <- attr(p.out, "srcref")
	rows <- sapply(srcref, "[", 1)
	cols <- sapply(srcref, "[", 2)
	return(row < rows | (rows == row & col < cols))
}

.isInside <- function (p.out, row, col)
{
	srcref <- attr(p.out, "srcref")
	startRows <- sapply(srcref, "[", 1)
	startCols <- sapply(srcref, "[", 2)
	endRows <- sapply(srcref, "[", 3)
	endCols <- sapply(srcref, "[", 4)
	return(((row == startRows & col >= startCols) | (row >= startRows)) &
		((row == endRows & col <= endCols) | (row <= endRows)))
}

.isFunction <- function (chunk)
{
	env <- new.env()
	out <- try({
		eval(chunk, envir = env)
		name <- ls(env) 
	}, silent = TRUE)
	testFun <- !inherits(out, "try-error") && length(name) == 1 &&
		inherits(env[[name]], "function") 
	if (testFun) attr(testFun, "fun") <- name
	return(testFun)
}
