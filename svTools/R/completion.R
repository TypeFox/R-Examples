completeCode <- function (...)
{
	utilsNS <- getNamespace("utils")
	out <- completion(..., min.length = 1)
	if (is.null(out)) {
		out <- matrix("", ncol = 4, nrow = 0) 
		token <- utilsNS$.guessTokenFromLine()
	} else {
		token <- attr(out, "token")
		if (is.null(token)) token <- utilsNS$.guessTokenFromLine()
		types <- rep("function", nrow(out))
		completions <- out[, 1]
		types[regexpr("= *$", completions) > 0] <- "argument"
		types[regexpr(":: *$", completions) > 0] <- "package"
		## Arguments first, then functions, then packages
		out <- cbind(out, types)[order(types), , drop = FALSE]
	}
	
	fun <- utilsNS$inFunction()
  	if (length(fun) && !is.na(fun)) {
		tooltip <- callTip(fun)
  	} else {
		tooltip <- NULL
		fun <- ""
  	}
	
	return(list(token = token, completions = out, fun = fun, tooltip = tooltip))
}

completePch <- function (line)
{	
	allPchs <- 1:25
	if (regexpr("pch *= *[^,)]*$", line) > 0) {
		start <- gsub("^.*= *", "", line)
		if (regexpr("^ *[0-9]+ *$", start) <= 0) { 
			pchs <- allPchs
		} else {
			int <- try(as.integer(start), silent = TRUE)
			if (inherits(int, "try-error")) pchs <- allPchs
			out <- as.integer(grep(start, allPchs, value = TRUE))
			pchs <- if (length(out)) out else allPchs
		}
		return(list(token = line, completions = pchs))
	}
}

completeCol <- function (line)
{	
	token <- sub("^.*=[[:space:]]*", "", line)
	start <- gsub("[[:space:]]+", "", token)
	
	if (regexpr("['\"]", start) > 0) {
		## Look at named colors
		start <- gsub("['\"]", "", start) 
		
		allColors <- colors()
		cols <- allColors[regexpr(start, allColors) > 0]
		if (!length(cols)) cols <- allColors
		rgb <- t(col2rgb(cols))
		cols <- paste('"', cols, '"', sep = "")
	} else {
		## Look at colors in the palette
		pal <- palette()
		if (nchar(start)) {
			cols <- 1:length(pal)
			cols <- cols[regexpr(start, cols) > 0]
			if (!length(cols)) cols <- 1:length(pal)
		} else cols <- 1:length(pal)
		rgb <- t(col2rgb(pal[cols]))
	}
	return(list(token = line, completions = cols, col.rgb = rgb))
}

completeLty <- function (line)
{	
	ltys <- c("blank", "dashed", "solid", "dotted", "longdash", "twodash" )
	
	token <- gsub("^.*=[[:space:]]*", "", line)
	start <- gsub("([[:space:]]+|'|\")" , "", token)
	
	matches <- ltys[regexpr(start, ltys) > 0]
	if (!length(matches)) matches <- ltys
	
	return(list(token = line, completions = ltys))
}
