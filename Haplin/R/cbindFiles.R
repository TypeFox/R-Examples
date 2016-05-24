cbindFiles <- function(infiles, outfile, col.sep, ask = TRUE, verbose = TRUE){
##
## Combines two or more files, column-wise (side-by-side). 
## (Reads each file line by line, pastes corresponding lines, then writes to outfile.)
##
## infiles is a character vector of names (and paths) of the files to combine.
## outfile is a character string giving the name and path of the resulting file.
## col.sep specifies the column separator which will be inserted between files.
## ask is a logical variable. If "TRUE", convertPed will ask before overwriting an already existing 'outfile'.
## By default, the logical variable verbose equals "TRUE", which means that output is displayed for each iteration.
##
#
.file.numbers <- length(infiles)
#
## Misc errors
if(!all(file.exists(infiles))) stop("\"infiles\" contains invalid file names", call. = F)
if(any(outfile == infiles)) stop("\"outfile\" is equal to one of the files given in \"infiles\"", call. = F)
#
## Open infiles for reading and writing
.infiles <- vector("list", .file.numbers) 
for(i in 1:.file.numbers) .infiles[[i]] <- file(description = infiles[i], open = "r")
#
## Make sure all connections are closed when function exits
on.exit(lapply(.infiles, close))
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
## If .file.numbers equals 1, copy file directly, display output and exit program
if(.file.numbers == 1){
	file.copy(infiles, outfile, overwrite = TRUE, copy.mode = TRUE)
	cat(paste(infiles, "copied and saved to file", outfile, "\n", sep = " "))
	return(invisible())
}	
#
## (Re)open outfile for writing
.outfile <- file(description = outfile, open = "w")
#
## Make sure the connection is closed when function exits
on.exit(close(.outfile), add = TRUE)
#
## Loop over lines
.k <- 0
repeat{
	.k <- .k + 1
	#
	.breakoff <- rep(FALSE, .file.numbers) # Equals TRUE if end of file is reached
	#
	## Read a single line from each file. NOTE: Reads numeric as character
	.line <- sapply(.infiles, readLines, n = 1)
	.breakoff <- (sapply(.line, length) == 0)
	.line[which(.breakoff)] <- NULL
	.line <- paste(.line, collapse = col.sep)
	#
	## Break off if all files are empty
	if(all(.breakoff)) break
	#
	## Display output 
	if(verbose) cat(.k, " --- ", sep = "")
	if(verbose & any(.breakoff)) cat(paste("The line is empty for the following file(s):", paste(infiles[which(.breakoff)], collapse = ", "), sep = ""))
	#
	## Write to new file
	writeLines(.line, .outfile)
	#
	## Newline
	if(verbose) cat("\n")
}# end repeat
#
## Display output
cat(paste("Combined lines from", .file.numbers, "files, saved to file", outfile, "\n", sep = " "))
#
## Return empty
return(invisible())
#
}

