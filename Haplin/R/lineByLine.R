lineByLine <- function(infile, outfile, linefunc = identity, choose.lines = NULL, choose.columns = NULL, col.sep = " ", ask = TRUE, blank.lines.skip = TRUE, verbose = TRUE, ...){
##
## Reads file line by line, modifies each line using the argument linefunc, and then writes to outfile.
##
## By default, linefunc returns its argument.
## Default sets ask = TRUE. No file is overwritten unless specified by user.
## choose.lines and choose.columns enable the selection of certain lines and/or columns of the infile. 
## Both are set to NULL as default, which means that all lines and/or columns are selected. If not set to NULL, must be numeric vectors with all values > 0 or all values < 0.
## Values > 0 mean that these lines and/or columns are chosen, whereas values < 0 mean that the line modifications are carried out for all lines and/or columns except for those listed in these arguments.
## col.sep specifies the separator used to split the columns in the file. By default, col.sep = " ". To split at all types of spaces or blank characters, set col.sep = "[[:space:]]" or col.sep = "[[:blank:]]".
## ask is a logical variable. If "TRUE", lineByLine will ask before overwriting an already existing 'outfile'.
## blank.lines.skip ignores blank lines in the input if TRUE.
## If TRUE, verbose displays the line number for each iteration, in addition to output from linefunc.
##
#
## Error in lineByLine if outfile == infile
if(file.exists(outfile) & outfile == infile) stop("\"outfile\" is equal to \"infile\"", call. = F)
#
## Misc errors regarding choose.lines and choose.columns
if(any(duplicated(choose.lines))) stop("\"choose.lines\" contains duplicated values", call. = F)
if(!is.null(choose.lines) && !is.numeric(choose.lines)) stop("\"choose.lines\" must have the value \"NULL\" or be numeric", call. = F)
if(is.numeric(choose.lines) && !(all(choose.lines <= 0) | all(choose.lines >= 0))) stop("Invalid line number(s). \"choose.lines\" cannot consist of both positive and negative values", call. = F)
if(!is.null(choose.columns) && !is.numeric(choose.columns)) stop("\"choose.columns\" must have the value \"NULL\" or be numeric", call. = F)
if(is.numeric(choose.columns) && !(all(choose.columns <= 0) | all(choose.columns >= 0))) stop("Invalid column number(s). \"choose.columns\" cannot consist of both positive and negative values", call. = F)
#
## Open infile for reading and writing
.infile <- file(description = infile, open = "r")
#
## Make sure the connection is closed when function exits
on.exit(close(.infile))
#
## If outfile exists & ask == TRUE, query user
if(file.exists(outfile) & ask){
	.answer <- readline(paste('Overwrite ', outfile, '? (y/n)', sep = ""))
	if(.answer != "y"){
		cat("Stopped without overwriting file\n")
		return(invisible())
	}
}
unlink(outfile)
#
## (Re)open outfile for writing
.outfile <- file(description = outfile, open = "w")
#
## Make sure the connections are closed when function exits
on.exit(close(.outfile), add = TRUE)
#
## Loop over lines
.i <- 0
#
.k <- 0 # Equals line number if line does not have enough columns 
#
.s <- 0 # Number of lines modified
#
repeat{
	.i <- .i + 1
	#
	## Break off at a given number of lines
	if(!is.null(choose.lines) && all(choose.lines >= 0) && (.i > max(choose.lines) + 0.1)) break
	#
	## Read a single line. NOTE: Reads numeric as character
	.line <- readLines(.infile, n = 1) 
	#
	## Break off at end of file
	if(length(.line) == 0) break 
	#
	## Skip blank lines if blank.lines.skip == TRUE. Otherwise a warning is given
	if(!nzchar(.line) & blank.lines.skip){
		.i <- .i-1
		next
	}
	if(!nzchar(.line) & !blank.lines.skip) warning(paste("Line ", .i, " is empty", sep = ""), call. = F)
	#
	## Choose lines 
	if(!is.null(choose.lines)){
		if(all(choose.lines >= 0) && !is.element(.i, choose.lines)) next
		else if(all(choose.lines <= 0) && is.element(.i, abs(choose.lines))) next		
	}	
	#
	## Split line into elements
	if(!(identical(linefunc, identity) & is.null(choose.columns))){
		.fixed <- TRUE
		if(col.sep != " " | col.sep != "\t") .fixed <- FALSE
		.line <- strsplit(.line, split = col.sep, fixed = .fixed)[[1]]
	}
	#
	## Choose columns, allowing reordering
	if(!is.null(choose.columns)){
		if(sum(abs(choose.columns) <= length(.line)) != length(choose.columns))	.k <- .i
		if(all(choose.columns >= 0)) .line <- .line[choose.columns[which(choose.columns <= length(.line) & choose.columns > 0)]]
		else if(all(choose.columns <= 0)) .line <- .line[choose.columns[which(abs(choose.columns) <= length(.line) & abs(choose.columns) > 0)]]
	}	
	#
	## Display line number (and output from linefunc) 
 	if(verbose) cat(.i, " --- ", sep = "")
	#
	.s <- .s + 1
	#
	## Convert line
	if("verbose" %in% names(formals(linefunc))) .line <- linefunc(x = .line, verbose = verbose, ...)
	else .line <- linefunc(x = .line, ...)
	#
	## Newline
	if(verbose) cat("\n")
	#
	## Display invalid column numbers
	if((.k == .i) & verbose) cat("Invalid column number(s).\n")
	#
	## Paste elements, separated by space
	if(!(identical(linefunc, identity) & is.null(choose.columns))) .line <- paste(.line, collapse = " ")
	#
	## Skip blank lines in order to delete the requested rows/lines
	if(!nzchar(.line) & blank.lines.skip) next
	#
	## Write to new file
	writeLines(.line, .outfile) 
}# end repeat
#
## Display number of lines read and converted
if(!is.null(choose.lines)){
	if(!sum(abs(choose.lines) <= .i-1) == length(choose.lines)) warning(paste("Invalid line number(s). Read and converted information for", .s, "line(s).\n", sep = " "), call. = F)
	else cat("Read and converted information for", .s, "line(s).\n", sep = " ")
} else cat("Read and converted information for", .s, "line(s).\n", sep = " ")
#
## Warning if choose.column contains invalid column number(s)
if(.k != 0) warning("\"choose.columns\" contains invalid column number(s).", call. = F)
#
## Return number of lines read
return(invisible(.i - 1))
}

