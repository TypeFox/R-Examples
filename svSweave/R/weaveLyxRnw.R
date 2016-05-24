## TODO: see the tables package for better looking tables, including in LaTeX!

weaveLyxRnw <- function (file, driver = RweaveLatex(),
syntax = getOption("SweaveSyntax"), encoding = "UTF-8", width = 80,
useFancyQuotes = TRUE, logFile = file.path(tempdir(), ".lyxSweave.log"), ...)
{
	## Run in LyX as the Sweave -> LaTex (from Sweave|pdflatex|plain) converter:
	##> R -e svSweave::weaveLyxRnw(\"$$i\",log=\"/tmp/.lyxSweave.log\"[,driver=highlight::HighlightWeaveLatex()]) -q --vanilla

	## Switch encoding (we do work with UTF-8 by default)
	oenc <- options(encoding = encoding)
	on.exit(options(oenc))
	Sys.setlocale("LC_CTYPE", "UTF-8")

	## By default, use fancy quotes
	ofc <- options(useFancyQuotes = useFancyQuotes)
	on.exit(options(ofc), add = TRUE)
	
	## Set default width for text to a reasonable value
	owidth <- options(width = width)
	on.exit(options(owidth), add = TRUE)

	## Issue warnings immediately
	owarn <- options(warn = 1)
	on.exit(options(owarn), add = TRUE)

	## Process 'file'
	if (!file.exists(file)) {
		stop("You must provide the name of an existing .Rnw file to process!")
	} else {
		## Redirect output
		unlink(logFile)
		on.exit({
			sink(type = "message")
			sink()
			## Echo results
			try(cat(readLines(logFile), sep = "\n"), silent = TRUE)
		}, add = TRUE)
		con <- file(logFile, open = "wt")
		sink(con)
		sink(con, type = "message")
		cat("Weaving ", basename(file), " ...\n", sep = "")
		
		## Clean the R noweb file
		opts <- cleanLyxRnw(file, encoding = encoding)
		
		## Weave or Knit the file
		if (!is.list(opts) || !length(opts$kind))
			stop("no Sweave kind defined (must be Sweave or Knitr)")
		cat("Processing the document using ", opts$kind, "\n", sep = "")
		switch(opts$kind,
			Sweave = Sweave(file, driver = driver, syntax = syntax, ...),
			Knitr = knit(file, ...),
			stop("Wrong kind of Sweave document: '", opts$kind, "'")
		)	
	}
}

knitLyxRnw <- function (file, encoding = "UTF-8", width = 80,
useFancyQuotes = TRUE, logFile = file.path(tempdir(), ".lyxSweave.log"), ...)
{
	## Run in LyX as the Sweave -> LaTex (from Sweave|pdflatex|plain) converter:
	##> R -e svSweave::knitLyxRnw(\"$$i\",log=\"/tmp/.lyxSweave.log\") -q --vanilla

	## Switch encoding (we do work with UTF-8 by default)
	oenc <- options(encoding = encoding)
	on.exit(options(oenc))
	Sys.setlocale("LC_CTYPE", "UTF-8")

	## By default, use fancy quotes
	ofc <- options(useFancyQuotes = useFancyQuotes)
	on.exit(options(ofc), add = TRUE)
	
	## Set default width for text to a reasonable value
	owidth <- options(width = width)
	on.exit(options(owidth), add = TRUE)

	## Issue warnings immediately
	owarn <- options(warn = 1)
	on.exit(options(owarn), add = TRUE)

	## Process 'file'
	if (!file.exists(file)) {
		stop("You must provide the name of an existing .Rnw file to process!")
	} else {
		## Redirect output
		unlink(logFile)
		on.exit({
			sink(type = "message")
			sink()
			## Echo results
			try(cat(readLines(logFile), sep = "\n"), silent = TRUE)
		}, add = TRUE)
		con <- file(logFile, open = "wt")
		sink(con)
		sink(con, type = "message")
		cat("Weaving ", basename(file), " with knitr...\n", sep = "")
		
		## Clean the R noweb file
		cleanLyxRnw(file, encoding = encoding)
		
		## Weave the file using knitr
		knit(file, ...)
	}
}
