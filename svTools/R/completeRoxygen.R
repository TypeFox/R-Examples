completeRoxygen <- function (line = "#'")
{	
	roxygenTags <- rbind( 
		c("author",              "\\author"                 , "Author of the function"),
		c("aliases",             "\\alias, ..."             , ""),
		c("concept",             "\\concept"                , ""),
		c("examples",            "\\examples"               , ""),
		c("keywords",            "\\keyword, ..."           , ""),
		c("method",              "\\method"                 , ""),
		c("name",                "\\name"                   , ""),
		c("note",                "\\note"                   , ""),
		c("param",               "\\arguments{\\item, ...}" , ""),
		c("references",          "\\references"             , ""),
		c("return",              "\\value"                  , ""),
		c("seealso",             "\\seealso"                , ""),
		c("title",               "\\title"                  , ""),
		c("TODO",                ""                         , ""),
		c("usage",               "\\usage"                  , ""),
		c("callGraph",           ""                         , "Create a call graph of the default depth, excluding primitive functions"),
		c("callGraphPrimitives", ""                         , "Create a call graph of the default depth, including primitive functions"),
		c("callGraphDepth",      ""                         , "Change the depth of the callgraph from the default of 2"),
		c("include",             ""                         , "See ?make.collate.roclet"),
		c("export",              "export"                   , ""),
		c("exportClass",         "exportClass"              , ""),
		c("exportMethod",        "exportMethod"             , ""),
		c("exportPattern",       "exportPattern"            , ""),
		c("S3method",            "S3method"                 , ""),
		c("import",              "import"                   , ""),
		c("importFrom",          "importFrom"               , ""),
		c("importClassesFrom",   "importClassesFrom"        , ""),
		c("importMethodsFrom",   "importMethodsFrom"        , "")) 
	
	if (line == "#'") {
		template <- " @%s "
		completions <- roxygenTags
		token <- ""
	} else if (line == "#' ") {
		template <- "@%s "
		completions <- roxygenTags
		token <- ""
	} else {
		template <- "%s "
		tag <- gsub("^#' *@", "", line)    
		matchingKeywords <- unique(c(grep(tag, roxygenTags[, 1], ignore.case = TRUE), 
			grep(tag, roxygenTags[, 3], ignore.case = TRUE)))
		completions <- if (!length(matchingKeywords)) roxygenTags else
			roxygenTags[matchingKeywords, , drop = FALSE]
		token <- tag
	}
	return(list(token = token, completions = sprintf(template, completions[, 1]),
		tooltip = completions[, 3])) 
}

completeRoxygenParam <- function (file, row, line = "#' @param ")
{
	potential <- paste(.argsFunAfter(file, row, all.args = FALSE),
		" ", sep = "")
	line <- gsub("^#' *@param", "", line)
	if (regexpr("^ +$", line) > 0)
		return(list(token = "", completions = potential))
	
	start <- gsub("^[[:space:]]+", "", line)
	if (regexpr("[[:space:]]+", start) > 0)
		return(list(token = "", completions = character(0)))
	
	completions <- grep(start, potential, value = TRUE)
	if (length(completions)) {
		return(list(token = start, completions = completions))
	} else {
		return(list(token = "", completions = character(0)))
	}
}

generateRoxygenTemplate <- function (file, row, column, 
author = getOption("svTools.roxygen.author"),
type = c("verbatim", "supperabbrev"))
{	
	p.out <- parse(file)
	where <- if (any(inside <- .isInside(p.out, row, column))) {
		which(inside)
	} else if (any(before <- .isBefore(p.out, row, column))) { 
		which(before)[1]
	} else length(p.out)
	
	isfun <- .isFunction(p.out[[where]])
	if(!isfun) return(list(ok = 0))
	funname <- attr(isfun, "fun")
	
	startPos <- as.numeric(attr(p.out, "srcref")[[where]][1:2])
	arguments <- .argsFunAfter(file = file, all.args =  TRUE, 
		p.out = p.out, row = startPos[1], target = where)
	
	template <- "#' ${1:Title (short) }\n#' \n#' ${2:Description (2-3 lines)}\n#' @export"
	if (length(arguments)) {
		template <- paste(template, paste("#' @param ", arguments, " ${",
			2 + 1:length(arguments), ": define `", arguments, "` }", sep = "",
			collapse = "\n"), sep = "\n") 
	}
	index <- length(arguments) + 3
	template <- paste(template, paste("#' @return ${", index,
		": What does the function return}", sep = ""), sep = "\n")
	index <- index + 1
	if (!is.null(author)) {
		author <- gsub("([^@])@([^@])", "\\1@@\\2", author)
		template <- paste(template, paste("#' @author ${", index, ":", author,
			"}", sep = ""), sep = "\n") 
	}
	index <- index + 1
	template <- paste(template, paste("#' @callGraph\n#' @examples\n#' ${",
		index, ":# executable code for `", funname, "`}\n", sep = ""), sep = "\n") 
	
	type <- match.arg(type)
	
	## Remove the super abbrev. stuff
	if (type == "verbatim") {
		template <- gsub("(?s)\\$\\{[[:digit:]]+: *([^}]+)\\}", "\\1", template,
			perl = TRUE)
		template <- paste(template, "\n", sep = "")
	}
	
	return(list(template = template, row = attr(p.out, "srcref")[[where]][1], ok = 1))
}
