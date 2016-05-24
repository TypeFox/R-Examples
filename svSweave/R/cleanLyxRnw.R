cleanLyxRnw <- function (RnwCon, RnwCon2 = RnwCon, encoding = "UTF-8")
{
	## By default, it is supposed to be a Sweave document
	opts <- list(kind = "Sweave")
	
	## Run in LyX as the Sweave copier using something like:
	##> R -e svSweave::cleanLyxRnw(\"$$i\",\"$$o\") -q --vanilla --slave
	
	## Make sure default encoding is the same
	options(encoding = encoding)
	Sys.setlocale("LC_CTYPE", "UTF-8")
	
	## Read the data in the Rnw file
	Rnw <- readLines(RnwCon, encoding = encoding)

	## Is it a knitr document?
	if (any(grepl("\\%Sweave-kind=knitr", Rnw))) opts$kind <- "Knitr"
	
	## If the Rnw file is produced with LyX and SciViews Sweave module, chunks are
	## separated by \rchunk{<<[pars]>>= ... @}

	## Beginning of R-Chunks (rewrite into <<[pars]>>=)
	#starts <- grepl("^\\\\rchunk\\{<<.*>>=$", Rnw)
	#Rnw[starts] <- sub("^\\\\rchunk\\{", "", Rnw[starts])
	starts <- grepl("^<<.*>>=$", Rnw)

	## End of R-Chunks (rewrite as @)
	#ends <- grepl("^@[ \t]*\\}$", Rnw)
	#Rnw[ends] <- "@"
	ends <- Rnw == "@"
	
	parts <- cumsum(starts | ends)
	chunk <- parts %% 2 > 0 	# R chunks are odd parts

	## Do we need to change something?
	if (!any(chunk)) {
		writeLines(Rnw, RnwCon2)
		return(invisible(opts))
	}

	## Eliminate empty strings not followed by empty strings inside chunks
	Rnw[chunk & Rnw == "" & Rnw != c(Rnw[-1], " ")] <- NA
	isna <- is.na(Rnw)
	Rnw <- Rnw[!isna]
	chunk <- chunk[!isna]

	## Eliminate successive empty strings (keep only one)
	Rnw[chunk & Rnw == "" & Rnw == c(Rnw[-1], " ")] <- NA
	Rnw <- Rnw[!is.na(Rnw)]

	## Convert tabulations into four spaces inside chunks (8 spaces by default)
	Rnw[chunk] <- gsub("\t", "    ", Rnw[chunk])

	## Write the result to the new .Rnw file
	## TODO: test useBytes = TRUE!
	writeLines(Rnw, RnwCon2, useBytes = TRUE)
	
	## Invisibly return opts
	invisible(opts)
}
