### Tries to parse a file
### Romain Francois <francoisromain@free.fr>
tryParse <- function (file, action, encoding = getOption("encoding"))
{	
	if (is.character(file)) {
		filename <- file
		file <- file(filename, encoding = encoding)
		on.exit(close(file))
	} else {
		filename <- summary(file)$description
	}
	
	## On Windows, filename could be c:/dir/file.R and the ':' does interfere
	## with the mechanism used to retrieve the error by parseError()
	## Hence, we prefer to change working dir to the file dir and pass only
	## the file basename to parse!
	dir <- dirname(filename)
	basefile <- basename(filename)
	odir <- getwd()
	setwd(dir)
	on.exit(setwd(odir), add = TRUE)
	
	out <- try(parse(file, srcfile = srcfile(basefile)), silent = TRUE)
	if (inherits(out, "try-error")) {
		err <- parseError(out)
		err$file <- rep(filename, nrow(err))
		if (!missing(action)) action(err)
		return(invisible(err))
	} else return(invisible(out))  
}
