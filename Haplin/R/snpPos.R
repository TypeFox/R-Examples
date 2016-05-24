snpPos <- function(snp.select, map.file, blank.lines.skip = TRUE){
##
## Returns the corresponding column numbers of SNP identifiers or SNP numbers in a standard ped file, sorted in the same order as snp.select. 
##
## snp.select is a character vector of the SNP identifiers (RS codes) or a numeric vector of SNP numbers. 
## map.file is a character string giving the name and path of the standard map file to be used. 
## If blank.lines.skip is "TRUE", blank lines are ignored reading the map file.
##
#
## Read map file
.map.file = read.table(map.file, header = TRUE, stringsAsFactors = FALSE, blank.lines.skip = blank.lines.skip)
#
## Error if snp.select has duplicated values
if(any(duplicated(snp.select))) stop('\"snp.select\" contains duplicated values', call. = F)
#
## Find the corresponding columns in the ped file
if(is.character(snp.select)){
	.pos.temp <- match(snp.select, .map.file[,2])
	if(any(is.na(.pos.temp))) stop('"snp.select" contains SNP names(s) that cannot be found in the map file', call. = F)
} else if(is.numeric(snp.select)){ 
	if(max(snp.select) > nrow(.map.file) | min(snp.select)*(-1) >= 0) stop('"snp.select" contains invalid SNP number(s)', call. = F)
	.pos.temp <- snp.select
} else stop('"snp.select" must be a numeric vector of SNP numbers or a character vector of SNP names', call. = F)
#
.pos <- rep(0,2*length(.pos.temp))
for(i in 1:length(.pos.temp)){
	.pos[2*i-1] <- 2*.pos.temp[i] -1
	.pos[2*i] <- 2*.pos.temp[i]
}
.pos <- .pos + 6
#
## Return column numbers
return(.pos)
}
