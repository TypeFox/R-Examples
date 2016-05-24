convertPed <- function(ped.infile, map.infile, ped.outfile, map.outfile, create.unique.id = FALSE, convert, snp.select = NULL, choose.lines = NULL, col.sep = " ", ask = TRUE, blank.lines.skip = TRUE, verbose = TRUE){
##
## Recodes a ped file with several options: 
## 	- Possibly creating unique individual IDs
## 	- Possibly converting the SNP alleles between letters (A,C,G,T) and numbers (1,2,3,4) 
## 	- Possibly making a selection of SNPs.
##
## ped.infile is a character string giving the name and path of the standard ped file to be modified.
## map.infile is a character string giving the name and path of the to-be-modified standard map file. Optional if snp.select = NULL.
## ped.outfile is a character string of the name and path of the converted ped file.
## map.outfile is a character string giving the name and path of the modified map file.
## snp.select is a numeric vector of SNP numbers or a character vector of the SNP identifiers (RS codes). Default is NULL, which means that all SNPs are selected. 
## choose.lines enables the selection of certain lines of the infile. 
## Default is NULL, which means that all lines are selected. If not set to NULL, must be numeric vectors with values > 0.
## col.sep specifies the separator used to split the columns in the file. By default, col.sep = " ". To split at all types of space or blank characters, set col.sep = "[[:space:]]" or col.sep = "[[:blank:]]".
## ask is a logical variable. If "TRUE", convertPed will ask before overwriting an already existing 'outfile'.
## blank.lines.skip is a logical varible. If "TRUE", convertPed ignores blank lines in ped.infile and map.infile.
## create.unique.id. If "TRUE", the function creates a unique individual ID.
## convert. The option "ACGT_to_1234" recodes the SNP alleles from A,C,G,T to 1,2,3,4, whereas "1234_to_ACGT" converts from 1,2,3,4 to A,C,G,T. If "no_recode", no conversion occurs.
## verbose. Default is "TRUE", which means that the line number is displayed for each iteration, i.e. each line read and modified, in addition to the first ten columns of the converted line.
##
#
#
.snp.select <- snp.select
#
if(!is.null(.snp.select) & missing(map.infile)) stop("\"map.infile\" is missing. Only optional if \"snp.select\" is NULL", call.=F)
#
if(!missing(map.infile)){
	## Error if map.outfile is equal to map.infile
	if(file.exists(map.outfile) && map.outfile == map.infile) stop("\"map.outfile\" is equal to \"map.infile\"", call. = F)
	#
	## Read map file
	.map.infile = read.table(map.infile, header = TRUE, stringsAsFactors = FALSE, blank.lines.skip = blank.lines.skip)
} else if(!missing(map.outfile)) warning("\"map.outfile\" is not needed when \"snp.select\" is equal to NULL and \"map.infile\" is missing", call. = F)

#
## Find number of columns in ped.infile
.ped.infile <- file(description = ped.infile, open = "r")
.line <- readLines(.ped.infile, n = 1)
.fixed <- TRUE
if(col.sep != " " | col.sep != "\t") .fixed <- FALSE
.columns.ped.infile <- length(strsplit(.line, split = col.sep, fixed = .fixed)[[1]])
close(.ped.infile)
#
## Checking if the map and ped file have the same number of SNPs
if(!missing(map.infile) && .columns.ped.infile != 2*nrow(.map.infile)+6) stop(paste("The map and ped file do not have the same number of SNPs. This may be due to an incorrect col.sep argument or a missing header in the map file. Number of columns in ped file: ", .columns.ped.infile, ". Number of rows in map file: ", nrow(.map.infile), sep = ""), call. = F)
#
## Find the column positions in the ped file corresponding to the vector .snp.select using the function snpPos
if(!is.null(.snp.select)) .snp.select.num <- snpPos(snp.select = .snp.select, map.file = map.infile, blank.lines.skip = blank.lines.skip)
else .snp.select.num <- 7:.columns.ped.infile ## If snp.select is set to NULL (default), all SNPs are selected
#
## Error if choose.lines has duplicated values
if(any(duplicated(choose.lines))) stop("\"choose.lines\" contains duplicated values", call. = F)
#
## Create new map file if map.outfile does not exist. Overwrite existing file if requested.
if(!missing(map.infile) && file.exists(map.outfile) & ask){
	.answer <- readline(paste('Overwrite ', map.outfile, '? (y/n)', sep = ""))
	if(.answer != "y"){
		cat("Stopped without overwriting file\n")
		return(invisible())
	}
}
#
## Convert data frame containing map data
if(is.numeric(.snp.select)) .snp.select <- .map.infile[.snp.select,2]
if(!is.null(.snp.select)) .converted.map <- .map.infile[match(.snp.select, .map.infile[,2]),]
else if(!missing(map.infile)) .converted.map <- .map.infile
if(!missing(map.infile)){
	write.table(.converted.map, map.outfile, row.names = F, col.names = T, quote = F)
	#
	## Display number of SNPs in map file
	cat(paste("Created map file for", length(.snp.select.num)/2, "SNP(s) \n", sep = " "))
}
# 
## Line-by-line conversion of ped-file. Message to user
cat("Starting line-by-line conversion of ped file \n")
#
## Convert ped file using the function lineByLine
lineByLine(infile = ped.infile, outfile = ped.outfile, linefunc = lineConvert, choose.lines = choose.lines, col.sep = col.sep, ask = ask, blank.lines.skip = blank.lines.skip, verbose = verbose, create.unique.id = create.unique.id, convert = convert, snp.select.num = .snp.select.num)
#
## Return empty
return(invisible())
#
}
