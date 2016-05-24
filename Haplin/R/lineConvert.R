lineConvert <- function(x, create.unique.id = FALSE, convert, snp.select.num = 7:length(x), verbose = TRUE){
##
## Possibly creates a unique individual ID, converts genetic allele coding from A,C,G,T to 1,2,3,4 or from 1,2,3,4 to A,C,G,T 
## and extracts a selection of SNPs (reordering among the SNPs possible), performing operations on a single line (one individual) in a ped file.
## The function is used by lineByLine or convertPed to update the entire ped file, i.e. recoding information
## for every individual in the ped file.
##
## x is a character vector, representing a single line of the file, split at spaces.
## The function returns the modified line, given as a character vector.
##
## create.unique.id  == TRUE will create unique individual IDs. 
## convert must take one and only one of the values "ACGT_to_1234", "1234_to_ACGT" or "no_recode". None of the options are default.
## snp.select.num is a numeric vector of the columns in the ped file corresponding to the selected SNPs. By default, all SNPs are selected without reordering among the SNPs.  
## The first ten columns of the converted line are displayed if verbose equals "TRUE".
##
## Note: assumes a standard ped file as input
##
#
.xnew <- x
#
## Error if convert is not specified correctly
.options.convert <- c("ACGT_to_1234", "1234_to_ACGT", "no_recode")
if(!convert %in% .options.convert | length(convert) != 1) stop("argument \"convert\" is not specified correctly", call. = F)
#
## Error if snp.select.num contains number(s) corresponding to undefined columns in the ped file
if(!all(is.element(snp.select.num, 7:length(x)))) stop("\"snp.select.num\" contains number(s) corresponding to undefined columns in the ped file", call. = F)
#
## Create new id variables by pasting family and individual IDs
if(create.unique.id){
	if(.xnew[2] != "0") .xnew[2] <- paste(.xnew[1], .xnew[2], sep = "_")
	if(.xnew[3] != "0") .xnew[3] <- paste(.xnew[1], .xnew[3], sep = "_")
	if(.xnew[4] != "0") .xnew[4] <- paste(.xnew[1], .xnew[4], sep = "_")
}
#
## Select SNPs
.xnew <- .xnew[c(1:6, snp.select.num)]
#
## Convert genetic coding from A,C,G,T to 1,2,3,4 for selected SNPs
if(convert == "ACGT_to_1234"){
	.xnew[-(1:6)][.xnew[-(1:6)] == "A"] <- "1"
	.xnew[-(1:6)][.xnew[-(1:6)] == "C"] <- "2"
	.xnew[-(1:6)][.xnew[-(1:6)] == "G"] <- "3"
	.xnew[-(1:6)][.xnew[-(1:6)] == "T"] <- "4"
}
#
## Convert genetic coding from 1,2,3,4 to A,C,G,T for selected SNPs
if(convert == "1234_to_ACGT"){
	.xnew[-(1:6)][.xnew[-(1:6)] == "1"] <- "A"
	.xnew[-(1:6)][.xnew[-(1:6)] == "2"] <- "C"
	.xnew[-(1:6)][.xnew[-(1:6)] == "3"] <- "G"
	.xnew[-(1:6)][.xnew[-(1:6)] == "4"] <- "T"
}
#
## Display 10 first elements, without newline
if(verbose) cat(paste(.xnew[1:min(10, length(.xnew))], collapse = " "))
#
## Return converted vector
return(.xnew)
}

