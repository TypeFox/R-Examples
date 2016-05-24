rbindFiles <- function(infiles, outfile, col.sep, header = FALSE, ask = TRUE, verbose = FALSE, add.file.number = FALSE, blank.lines.skip = FALSE){
##
## Binds two or more files, by rows. (Reads each file line by line, then writes to outfile.)
##
## infiles is a character vector of the files to be bound.
## outfile is a character string giving the name and path of the new file.
## col.sep specifies the separator used to split the columns in the files. To split at all types of spaces or blank characters, set col.sep = "[[:space:]]" or col.sep = "[[:blank:]]".
## header is a logical variable implying if the first line in each file contains the names of the variables. Equals FALSE by default, i.e., no headers assumed.
## ask is a logical variable. If "TRUE", rbindFiles will ask before overwriting an already existing 'outfile'.
## By default, the logical variable verbose equals "FALSE". If "TRUE", output is displayed for each iteration.
## If add.file.number equals TRUE, an extra first column will be added to the outfile, consisting of the file numbers for each line.
## blank.lines.skip ignores blank lines in the input if TRUE.
##
##
#
.file.numbers <- length(infiles)
#
## Misc errors
if(!all(file.exists(infiles))) stop("\"infiles\" contains invalid file name(s)", call. = F)
if(file.exists(outfile) & any(outfile == infiles)) stop("\"outfile\" is equal to one of the files given in \"infiles\"", call. = F)
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
## If .file.numbers equals 1, add.file.number == FALSE and blank.lines.skip == FALSE, copy file directly, display output and exit program
if(.file.numbers == 1 & !add.file.number & !blank.lines.skip){
	file.copy(infiles, outfile, overwrite = TRUE, copy.mode = TRUE)
	cat(paste(infiles, "copied and saved to file", outfile, "\n", sep = " "))
	return(invisible())
}
#
## (Re)open outfile for writing
.outfile <- file(description = outfile, open = "w")
#
## Make sure the connection is closed when function exits
on.exit(close(.outfile))
#
## Create temporary file
.temp.outfile <- tempfile(pattern = "temporaryOutfile", fileext = ".txt")
on.exit(unlink(.temp.outfile), add = TRUE)	
#
## Misc initializations
.s <- 0
.column.length <- rep(-1, 2) # To check if all files have the same number of columns in first two rows
.header <- NULL
.fixed <- TRUE
if(col.sep != " " | col.sep != "\t") .fixed <- FALSE
#
## New function
rbindModify <- function(x, add.file.number = FALSE, file.number = NULL, verbose = FALSE){
	#
	.xnew <- x
	#
	## Adding new column
	if(add.file.number) .xnew = c(file.number,.xnew)
	#
	## Display 10 first elements, without newline
	if(verbose & !length(.xnew) == 0) cat(paste(.xnew[1:min(10, length(.xnew))], collapse = " "))
	#
	## Return converted vector
	return(.xnew)
}
#
for(i in 1:.file.numbers){
	#
	## Close connection
	if(i != 1) close(.infile)
	#
	## Open infiles for reading and writing
	.infile <- file(description = infiles[i], open = "r")
	#
	## Make sure connection is closed when function exits
	if(i == 1) on.exit(close(.infile), add = TRUE)
	#
	.line <- rep(NA, 2)
	.k <- 0
	#
	repeat{
		## Read a single line. NOTE: Reads numeric as character
		.line.temp <- readLines(.infile, n = 1) 
		#
		## Break off at end of file
		if(length(.line.temp) == 0) break 
		#
		.line[.k+1] <- .line.temp
		#
		## Skip blank lines if blank.lines.skip == TRUE
		if(!nzchar(.line.temp) & blank.lines.skip) next
		#
		.k <- .k + 1
		#
		## Check the two first lines only
		if(.k == 2) break
	}
	#
	## Error if first file is empty. Warnings if the following files are empty
	if(.k == 0 && all(is.na(.line))){
		if(i == 1) stop("The first file is empty", call. = F)
		else{
			warning(paste("File number", i, "is empty", sep = " "), call. = F)
			next
		}	
	}
	#
	## Error if first file contains only blank lines and blank.lines.skip == TRUE. Warnings if the following files contain only blank lines
	if(.k == 0 && !all(is.na(.line) && blank.lines.skip)){
		if(i == 1) stop("The first file contains only blank lines", call. = F)
		else{
			warning(paste("File number", i, "contains only blank lines", sep = " "), call. = F)
			next
		}	
	}
	#	
	## Split lines into elements
	.line <- strsplit(.line, split = col.sep, fixed = .fixed)
	#
	## Do the files have the same number of columns? (Check first two lines of each file.)
	if(i == 1){ 
		.column.length[1] <- length(.line[[1]])
		if(.k == 2) .column.length[2] <- length(.line[[2]])
	} else if((length(.line[[1]]) != .column.length[1]) | ((.k == 2 & .column.length[2] != -1) && (length(.line[[2]]) != .column.length[2]))) warning(paste("File number", i, "does not have the same number of columns as the first file. This may be due to an incorrect col.sep argument. As a consequence, \"outfile\" may be incorrect", sep = " "), call. = F)
	# 
	## Do the files have equal header?
	if(header){
		if(i == 1) .header <- .line[[1]]
		else if(!identical(.line[[1]],.header)) warning(paste("File number", i, "does not have the same header as the first file. As a result, \"outfile\" may be incorrect. ", sep = " "), call. = F)
	}
	#
	## Appending files
	cat(paste("Appending", infiles[i]), "\n", sep = " ")	
	#
	if(!header){
		if(!verbose & !add.file.number & !blank.lines.skip) file.append(outfile, infiles[i])
		else{
			suppressWarnings(lineByLine(infile = infiles[i], outfile = .temp.outfile, linefunc = rbindModify, col.sep = col.sep, ask = FALSE, verbose = verbose, blank.lines.skip = blank.lines.skip, add.file.number = add.file.number, file.number = i))
			file.append(outfile, .temp.outfile)
		}
	}	
	if(header){		
		if(i == 1) suppressWarnings(lineByLine(infile = infiles[i], outfile = .temp.outfile, linefunc = rbindModify, col.sep = col.sep, ask = FALSE, verbose = verbose, blank.lines.skip = blank.lines.skip, add.file.number = add.file.number, file.number = i))
		else suppressWarnings(lineByLine(infile = infiles[i], outfile = .temp.outfile, linefunc = rbindModify, choose.lines = c(-1), col.sep = col.sep, ask = FALSE, verbose = verbose, blank.lines.skip = blank.lines.skip, add.file.number = add.file.number, file.number = i))
		file.append(outfile, .temp.outfile)			
	}
	#
	.s <- .s + 1
}	
#
## Display output
cat(paste("Combined rows from", .s, "files, saved to file", outfile, "\n", sep = " "))
#
## Return empty
return(invisible())
#
}
