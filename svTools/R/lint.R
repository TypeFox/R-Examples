### A function that lints provided code and returns a data.frame, a flat text
### output or a rjson object (for Komodo)
lint <- function (file, text = NULL, filename = NULL,
encoding = getOption("encoding"), type = c("data.frame", "flat", "rjson"),
sep = "+++")
{
	if (missing(file)) {
		if (is.null(text) || !is.character(text))
			stop("If you do not provide 'file', you must provide R code in 'text'")
		## Place the code in a temporary file
		f <- tempfile()
		on.exit(unlink(f))
		cat(text, sep = "\n", file = f)
	} else f <- file
	type <- match.arg(type)
	## Run .lint() on this file
	res <- .lint(f, encoding = encoding)
	## Is it something to return?
	if (nrow(res) == 0) return("") else {
		## For type == data.frame, change nothing
		if (type == "data.frame") return(res)
		## I prefer to get warning|error+++(filename)+++line+++col+++message
		res <- res[, c(5, 1:4)]
		if (is.null(filename)) {
			res <- res[, -2]  # Eliminate filename
		} else {
			## Replace file by filename
			res$file <- rep(filename, nrow(res))
		}
		if (type == "rjson") {
			## Print a rjson object version of the data
			cat(toRjson(res), sep = "")
		} else {
			## Print a flat version of the results
			cat(apply(res, 1, paste, collapse = sep), sep = "\n")
		}
		return(invisible(res))
	}
}

### Wrapper for the checkUsage function in codetools. 
### Romain Francois <francoisromain@free.fr>
.lint <- function (file, encoding = getOption("encoding"))
{	
	if (is.character(file) && regexpr('^rwd:', file) > 0)
		file <- sub('^rwd:', getwd(), file)
	file <- tools::file_path_as_absolute(file)
	
	old.op <- options(encoding = encoding)
	on.exit(options(old.op))
	
	resetErrors(file = file)
	
	## First parse for errors
	p.out <- tryParse(file, action = addError, encoding = encoding)
	if (inherits(p.out, "data.frame")) return(getErrors(file = file))
	if (length(p.out) == 0) return(emptyError())
	
	## Hack to retrieve information from codetools
	chkres <- new.env()
	chkres$findings <- NULL
	report <- function (x)
		assign("findings", c(chkres$findings, x), envir = chkres)
	
	addErrorFile <- function (line, msg)
		addError(line = line, message = gsub("(\\\n|^: )", "", msg),
			file = file, type = "warning")
	
	finding <- function (txt, p, i, rx, rx2) {
		param <- sub(rx, "\\1", txt) 
		param <- gsub("\\.", "\\\\.", param) 
		exprs <- p[[i]][[3]][[3]]
		srcref <- do.call(rbind, lapply(attr(exprs, "srcref"), as.integer))  
		regex <- sprintf(rx2, param)
		for (j in 1:length(exprs)) {
			## I don't understand this, since retrieve of source code is easier
			#src <- .as.characterSrcRef(attr(exprs, "srcref")[[j]],
			#	useSource = TRUE, encoding = encoding)
			#matchingLines <- grep(regex, src)
			matchingLines <- grep(regex, attr(exprs, "srcref")[[j]])
			if (length(matchingLines))
			  	return(matchingLines + as.integer(srcref[j, 1]) - 1)
		}
	}
	
	## OK, tested with codetools 0.2-8
	findParamChangedByAssign <- function (txt, p, i)
		return(finding(txt, p, i, 
			rx = "^.*: parameter .(.*). changed by assignment.*\\\n", 
			rx2 = "[^.a-zA-Z0-9_]*%s[[:space:]]*(=|<-|<<-)"))
	
	## OK, tested with codetools 0.2-8
	findUnusedLocalAssign <- function (txt, p, i)
		return(finding(txt, p, i, 
			rx = "^.*: local variable .(.*). assigned but may not be used.*\\\n", 
			rx2 = "^[^.a-zA-Z0-9_(,]*%s[[:space:]]*(=|<-|<<-)"))
	
	## TODO: check this!
	findNoGlobalDef <- function (txt, p, i)
		return(finding(txt, p, i, 
			rx = "^.*: no visible global function definition for .(.*)..*\\\n", 
			rx2 = "[^.a-zA-Z0-9_]*%s[[:space:]]*\\("))
	
	## TODO: check this!
	findNoLocalDefAsFun <- function (txt, p, i)
		return(finding(txt, p, i, 
			rx = "^.*: local variable .(.*). used as function with no apparent local function definition.*\\\n", 
			rx2 = "[^.a-zA-Z0-9_]*%s[[:space:]]*\\("))
	 
	## TODO: check this!
	findNoBindingGlobalVar <- function (txt, p, i)
		return(finding(txt, p, i, 
			rx = "^.*: no visible binding for global variable .(.*)..*\\\n", 
			rx2 = "[^.a-zA-Z0-9_]*%s[^.a-zA-Z0-9_]*"))
	
	## OK, tested with codetools 0.2-8
	findMultipleLocalDef <- function (txt, p, i)
		return(finding(txt, p, i, 
			rx = "^.*: multiple local function definitions for .(.*). with different formal arguments.*\\\n", 
			rx2 = "[^.a-zA-Z0-9_]*%s[[:space:]]*(=|<-|<<-)[[:space:]]*function"))
	
	searchAndReport <- function (regex, fun) {
		if (length(test.match <- grep(regex, chkres$findings)))
			for (j in test.match) {
				out <- fun(chkres$findings[j], p.out, i)
				if (length(out))
					addErrorFile(out, sub(" \\([^)]+\\)\\\n$", "",
						chkres$findings[j]))
			}
	}
	
	for (i in 1:length(p.out)) {
		if (.looksLikeAFunction(p.out[[i]])) {
			env <- new.env()
			eval(p.out[[i]], envir = env)
			fname <- ls(env) 
			if (length(fname) == 1) {
				chkres$findings <- NULL
				checkUsage(env[[fname]], report = report, all = TRUE, name = "")
				#cat(chkres$findings, "\n")
				if (length(chkres$findings)) {
					searchAndReport("changed by assignment", findParamChangedByAssign) 
					searchAndReport("assigned but may not be used", findUnusedLocalAssign) 
					searchAndReport("no visible global function definition", findNoGlobalDef)
					searchAndReport("no apparent local function definition", findNoLocalDefAsFun)
					searchAndReport("no visible binding for global variable", findNoBindingGlobalVar)
					searchAndReport("multiple local function definitions", findMultipleLocalDef)  					
### TODO: this needs to be improved to deal with nested functions
#					if (length(test.assign <- grep(" may not be used", chkres$findings)))
#						for (j in test.assign)
#							addErrorFile(attr(p.out, "srcref")[[i]][1], chkres$findings[j])
				}
			}
		}
	}
	return(getErrors(file = file))
}
