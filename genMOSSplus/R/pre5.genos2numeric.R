pre5.genos2numeric <- function(file.ped, dir.ped, file.dat, dir.dat=dir.ped, dir.out, num.nonsnp.col=2, num.nonsnp.last.col=1, letter.encoding=TRUE, ped.has.ext=TRUE, dat.has.ext=TRUE, remove.bad.genos=FALSE, save.ids.name="") {
# Categorizes genotype data into 3 levels, 1, 2, 3.
# Genos with two different Alleles are encoded as "2". 
# Other genotypes are encoded as "1" or "3", where most frequent geno is "1". 
# No missing values allowed, must be done after imputation.
# Geno values should use letters A, T, C, G if letter.encoding=TRUE. 
 
# Example run:
# 
# For MaCH1 output format, with last column as disease status:
# pre5.genos2numeric(dir.ped="/home/briollaislab/olga/curr/data/pipe/f06_combined", dir.out="/home/briollaislab/olga/curr/data/pipe/f07_numeric", file.ped="CGEM_Breast_23.txt", file.dat="CGEM_Breast_chr23.removed.dat", num.nonsnp.col=2, num.nonsnp.last.col=1)
#
# For MaCH1 output format:
# pre5.genos2numeric(dir.ped="/home/briollaislab/olga/curr/data/mach1in/result", dir.out="/home/briollaislab/olga/curr/data/mach1out", file.ped="last25CASE.200.mlgeno", num.nonsnp.col=2, num.nonsnp.last.col=0)
#
# For PLINK format:
# pre5.genos2numeric(dir.ped="/home/briollaislab/olga/curr/data/output11", file.ped="shortie.ped", dir.out="/home/briollaislab/olga/curr/data/output11", num.nonsnp.col=6, num.nonsnp.last.col=0)
#
#
# file.ped: the name of file with genotypes, after imputation. 
#   Entries should be either tab or space separated.
# dir.ped: directory where file.ped can be found.
# file.dat: the .dat file, should be tab separated, and no header.
# dir.dat: directory where file.dat can be found. Defaults to dir.ped.
# dir.out: output directory to which resulting file should be saved.
#   the file will be named "Num.<file.ped>".
# num.nonsnp.col: the number of leading columns that do not correspond to geno values.
#   Ex. for MaCH1 input file format there are 5 non-snp columns; 
#       for MaCH1 output format .mlgeno it's 2;
#       for Plink it's 6.
# num.nonsnp.last.col: the number of last columns that do not correspond to geno values.
#   Ex. If last column is the disease status (0s and 1s), then set this variable to 1.
#       If 2 last columns correspond to confounding variables, set the variable to 2.
# letter.encoding: whether or not the ecoding used for Alleles is letters (A, C, T, G).
#   if True, then does additional check for Alleles corresponding to the letters, and
#   prints out warning messages if other symbols appear instead. 
# ped.has.ext: whether or not file.ped name has a filename extension (ex. ".ped", ".txt")
#   This is necessary for naming the output file.
# dat.has.ext: whether or not file.dat name has a filename extension (ex. ".dat", ".txt")
# remove.bad.genos: do you want to remove a geno if at least one of its values is not
#   valid (ex. "2" when only letters are expected, or "NA", etc). 
#   Warning: set this to TRUE only if the CASE and CONTROLs have been merged into the file.ped,
#     (otherwise we do not want to remove some SNPs from CASE but not from CONTROL
#      and generate two different .dat files)
# save.ids.name: If not empty, then will save IDs of patients into another file with this name.
#   Since dataset is generally split across many files, one chromosome each, the patient IDs
#   should be the same across these files, thus it is enough to extract the patient ID ONCE,
#   when running this code on the smallest chromosome. For runs on all other chromosomes,
#   leave save.ids.name="" to save time and avoid redundant work. 
#   Could name output file as "patients.fam".
#
# Note: in case of any bad values in the file.ped (ex. "NA", "0/0", "0", "1 1", etc),
# the output file Num_<file.ped> will still be produced, with '2' encoded by default
# in the place of bad input values, if remove.bad.genos=FALSE. Warning messages will be printed.
# If remove.bad.genos=TRUE, then these SNPs will be entirely removed, along with their 
# names in the .dat file.
#
# Result:
# <file.ped>_num<ending.ped> - in dir.out directory, the resultant binary file: 
#      the SNP columns + last columns (but no user IDs will be recorded).
#      where <ending.ped> is the filename extension of file.ped.
# <file.dat>_num.dat - in dir.out directory, the corresponding .dat file, will be different
#      from original <file.dat> if remove.bad.genos=TRUE. 
# <save.ids.name> - the patient IDs, if save.ids.name is not empty "".
#
# Returns:
# <file.ped>_num<ending.ped> filename - the name of the output file.

#TODO: remove this line:
#source("machout2id.R")
#source("get.ext.R")

if(missing(file.ped)) stop("Name of the .ped file must be provided.")
if(missing(dir.ped)) stop("Name of input directory for .ped file must be provided.")
if(missing(file.dat)) stop("Name of the .dat file must be provided")
if(missing(dir.out)) stop("Name of output directory must be provided.")


# Check that .dat file is provided, if we are asked to remove bad SNPs:
if(file.dat == "") {
	print("Since you wish to remove bad SNPs (if they exist), you must provide the .dat file")
	return(NULL)
}

# Check that file actually exists.
full.dat <- paste(dir.dat, file.dat, sep="/")
if(!file.exists(full.dat)) {
	print(paste("Error: the .dat file: ", full.dat, " does not exist.", sep=""))
	return(NULL)
}

dat <- read.table(full.dat, sep="\t", stringsAsFactors=FALSE)
if(ncol(dat) <= 1) {
	print(paste("Error: the .dat file: ", full.dat, " is either empty or not tab separated.", sep=""))
	return(NULL)
}

# If the tab separator is not working, then try to load file using space separator.
D <- read.table(paste(dir.ped, file.ped, sep="/"), sep="\t", stringsAsFactors=FALSE)
if(ncol(D) <=1) {
	D <- read.table(paste(dir.ped, file.ped, sep="/"), sep=" ", stringsAsFactors=FALSE)
	if(ncol(D) <=1) {
	        print(paste("File ", dir.ped, "/", file.ped, " does not seem to have proper format", sep=""))
	        return(NULL)
	}
}

ncols <- ncol(D)
# Remove the first num.nonsnp.col columns from the data, and the last num.nonsnp.last.col.
#frontD <- D[, 1:num.nonsnp.col]
endD <- NULL
if(num.nonsnp.last.col > 0)
	endD <- D[, (ncols-num.nonsnp.last.col+1):ncols]
D <- D[, (num.nonsnp.col+1):(ncols-num.nonsnp.last.col)]

nrows <- nrow(D)			
ncols <- ncol(D)
newD <- matrix(2, nrows, ncols)
genos.to.keep <- rep(TRUE, times=ncols)
j <- 1

# Iterate over all the columns. 
# For each column, determine the unique set of geno values that repeat.
while(j <= ncols) {

	a <- unique(D[,j])
	if(length(a) > 4 && !letter.encoding) { # if too many geno types are present in one column
		print(paste("Warning: file ", file.ped, " on SNP ", j, " has more than 4 different geno values:", sep=""))
		print(a)
	}

	# Sort the unique array, s.t. most frequent genos appear first.
	counts <- rep(0, length(a))
	for(k in 1:length(a)) 
		counts[k] <- length(which(D[,j] == a[k]))
	ordering <- sort(counts, decreasing=T, index.return=T)$ix
	a <- a[ordering]

	# Process the genos
	next.val <- 1
	for (k in 1:length(a)) { # 1 <= k <= 4
		# Determine whether the value is a palindrome or not
		status <- pali(a[k], letter.encoding)
		if(status == 0) { # in case if unlikely geno values are encountered
			print(paste("Warning: file ", file.ped, " on SNP ", j, " has a geno with illegal value: ", a[k], sep=""))
			if(remove.bad.genos == TRUE) 
				genos.to.keep[j] <- FALSE
		}
		if(status == -2) { # in case if bad Allele value is detected if letters are expected
			print(paste("Warning: file ", file.ped, " on SNP ", j, " has a geno with unexpected value: ", a[k], sep=""))
			if(remove.bad.genos == TRUE) 
				genos.to.keep[j] <- FALSE
		}
		if(status == 1) { # if the value is a palindrome, assign it 1 and then 3.
			num <- which(D[,j] == a[k])
			newD[num, j] <- rep(next.val, length(num))
			next.val <- 3
		}
		# if value is not palindrome, do nothing: the encoding of "2" is already default of newD array.
	}
	j <- j + 1
}

ending.dat <- get.ext(file.dat, dat.has.ext)
out.dat.file <- paste(dir.out, "/", ending.dat$part1, "_num", ending.dat$ext, sep="")

# Only keep genos that are good. Also crop and save the .dat file
if(remove.bad.genos == TRUE) {
	newD <- newD[, genos.to.keep]
	dat <- dat[genos.to.keep, ]
	write.table(dat, file=out.dat.file, col.names=F, row.names=F, quote=F, sep="\t")
	num.genos.removed <- sum(!genos.to.keep)
	print(paste("Removing: ", num.genos.removed, "/", ncols, " = ", (num.genos.removed/ncols * 100), "% SNPs.", sep=""))
} else {
	# Save the dat file as is.
	write.table(dat, file=out.dat.file, col.names=F, row.names=F, quote=F, sep="\t")
}

# Append the last columns to the file:
if(num.nonsnp.last.col > 0)
	newD <- cbind(newD, endD)

ending.ped <- get.ext(file.ped, ped.has.ext)
out.file <- paste(ending.ped$part1, "_num", ending.ped$ext, sep="")
out.file.full <- paste(dir.out, out.file, sep="/")
write.table(newD, file=out.file.full, col.names=F, row.names=F, quote=F, sep="\t")

# Save the patient IDs if needed:
if(save.ids.name != "") {
	machout2id(file.case=file.ped, dir.file=dir.ped, name.out=save.ids.name, dir.out=dir.out)
}

return(out.file)

}


pali <- function(word, letter.encoding) {
# Checks the geno 'word':
# Returns 1 if two Alleles are the same, ignoring the symbols in the middle, (ex "AA", "A/A", "A A", "C C", "GG"...)
# Returns -1 if the two Alleles are different, (ex "G/A", "A/G", "TC", "G T"...)
# Returns 0 if the Allele is invalid, (ex one letter long, NA, "NA", "0 0", "1 1", "0/0")
# Returns -2 if letter.encoding=T and either of the Alleles is not a valid letter (not one of "A", "C", "T", or "G").
	len <- nchar(word)
	if(len <= 1)
		return(0)
	if(is.na(word) || word == "NA")
		return(0)
		
	first <- substr(word, 1, 1)
	last <- substr(word, len, len)

	if(first == "0" || first == "1" || first == "N")
		return(0)

	# Check for palindrome:
	# - if we're not expecting letters, or
	# - if we are expecting letters and these letters are valid
	valids <- c("A", "T", "C", "G")
	if(!letter.encoding || (!is.na(match(first, valids)) && !is.na(match(last, valids)))) {

		if(first == last) 
			return(1)
		else
			return(-1)		
	} else 
		return(-2)
	
}
