pre6.merge.genos <- function(dir.file, dir.dat=dir.file, dir.out=dir.file, file.out="CGEM_Breast_complete.txt", dat.out="CGEM_Breast_complete.dat", prefix.file, prefix.dat, key.file="", key.dat="", ending.file=".txt", ending.dat=".dat", num.nonsnp.col=0, num.nonsnp.last.col=1, weak.check=FALSE, plan=FALSE) {
# Puts together all the genos files and their corresponding .dat files for all chromosomes.
# The files should have last column as the disease status, and the number of
# individuals (rows) must match across all files.
# Also the files are expected to have no leading non-snp columns. If they exist, they
# will be removed.  
# The dat files are expected to have the SNP names in their second column.
# If the first column of .dat file is 'M', then it will be replaced by the chomosome
# number of the file name (the number that follows prefix.dat).
#
# This function tries to make sure that the geno files and dat files correspond. 
#
# Example:
#
# merge.genos(dir.file="/home/briollaislab/olga/curr/data/step1out", dir.dat="/home/briollaislab/olga/curr/data/mach1out", dir.out="/home/briollaislab/olga/curr/data/step1out/merged", prefix.file="CGEM_Breast.", prefix.dat="CGEM_Breast_dat.", key.file="", key.dat=".dat", plan=FALSE, weak.check=FALSE)
#
# dir.file - directory containing files with geno information.
# 	- file names in this directory must have the last column as the disease status.
# dir.dat - directory containing .dat files. Should be a list of geno IDs, one ID per line, no header.
#       - note, the extension of this file name can be anything (not necessarily ".dat").
#       - Defaults to same directory as dir.genos.
# dir.out - directory where the two output files will go. 
#       - Defaults to same directory as dir.genos.
# file.out - the name of the output file which will contain the combined geno information
#       and the last column will be the disease status.
# dat.out - the name of the output file which will contain all the corresponding SNP values. 
# prefix.file - the string that appears at the beginning of all the geno input file names.
#       - The file names are expected to begin with <prefix.file>, and then be immediately followed
#         by chromosome number, for example, in <dir.file> directory files named like :
#              "cgem_breast.21.pure.txt"
#              "cgem_breast.5.pure.txt"
#              "cgem_breast.24_and_25.txt" 
#          must have prefix="cgem_breast."
# prefix.dat - the string that appears at the beginning of all the .dat file names.
#       - Similarly to <prefix.file>, it must be immediately followed by the chromosome number.
# key.file - some key string that uniquely identifies the desired file names for the genos within
#       the directory. This is to give more flexibility, as prefix might not be enough.
#               "cgem_breast.12.pure.txt" for geno files, and 
#               "cgem_breast.12.dat" for .dat files, then prefix is the same in both cases,
#               so the program will not be able to tell which one is geno and which is dat.
#               Thus to disambiguate, you could additionally state that all the geno files
#               also contain "pure" or "txt", so set key.file=".pure" and key.dat=".dat".
# key.dat - some key string that uniquely identifies the desired .dat files. Similar to key.file.
# ending.file - the string with which all the geno filenames end.
# ending.dat - the string with which all the .dat filenames end.
# num.nonsnp.col - number of first columns that represent non-SNP values.
#       If this code is run after pre6.discretize(), then there will be no leading
#       non-SNP columns. However if you wish to skip discretization and first merge all files,
#       then you would likely have leading columns that do not correspond to SNP values.
#       These columns will be deleted.
# num.nonsnp.last.col - number of last columns that are non-SNP values.
#       This variable is expected to be 1 - for the last column that has disease status.
#       If more confounding variables are included, state how many last columns are corresponding
#       to all these non-geno variables as the last columns. 
#       These last columns will appear at the end of the ouput file.
# weak.check - since this function will try to check correspondence of the number of genos in
#       the genos file to the .dat file, the function would expect there to be the same number
#       of genos and .dat files. If you wish to by-pass these checks, set weak.check=TRUE,
#       in which case only the total final number of the resultant geno and .dat files will be checked
#       for consistency, and only a warning message will be printed if there's a problem.
# plan - if this option is TRUE, then this function will "do" nothing, but
#       will simply print which files it plans to combine in which order, 
#       since combination step itself might take time for large files.
#
# Returns: the full file name of the combined genos file.
#
# Note: this function makes use of LINUX commands: 'paste', 'cat', and 'wc'.

# TODO: remove later:
#source("get.data.dims.R")
#source("get.file.name.R")
#source("get.chrom.num.R")

if(missing(dir.file)) stop("Name of input directory must be provided.")
if(missing(prefix.file)) stop("Prefix of the genos input file name must be provided.")
if(missing(prefix.dat)) stop("Prefix of the dat file name must be provided.")


# Obtain file names, sorted in order of chrom number, and
# the corresponding chrom numbers
info.genos <- extract.chrom.dir(dir.file=dir.file, prefix=prefix.file, key=key.file, ending=ending.file)
files.genos <- info.genos$files
chroms.genos <- info.genos$chroms

info.dat <- extract.chrom.dir(dir.file=dir.dat, prefix=prefix.dat, key=key.dat, ending=ending.dat)
files.dat <- info.dat$files
chroms.dat <- info.dat$chroms

# save the number of snp values in each file, since its re-calculation will take
# time for large dataset files.
NUM.GENOS <- rep(0, times=length(files.dat))
# save the number of patients/rows in each of the files
NUM.ROWS <- rep(0, times=length(files.dat))

if(weak.check == FALSE) {

	# Check that each genos file has a corresponding dat file.
	if(length(files.genos) != length(files.dat) || !all(chroms.genos == chroms.dat)){
		print(paste("Aborting. The number of geno files = ", length(files.genos), " does not correspond to number of .dat files = ", length(files.dat), sep=""))
		print("genos:")
		print(files.genos)
		print("dat:")
		print(files.dat)
		return(NULL)
	}
	
	# Check that each geno file has exactly the same genos 
	# as listed in its corresponding .dat file.
	# Also save the number of users/rows for each of the files

	i <- 1
	while(i <= length(files.genos)) {
		print(paste("Counting dim: ", dir.file, "/",  files.genos[i], sep=""))
		num.dims <- get.data.dims(paste(dir.file, files.genos[i], sep="/"))
		num.snps.genos <- num.dims$ncols - num.nonsnp.last.col
		NUM.GENOS[i] <- num.snps.genos	# save for future use
		NUM.ROWS[i] <- num.dims$nrows

		num.snps.dat <- get.data.dims(paste(dir.dat, files.dat[i], sep="/"))$nrows
		print(paste("num cols: ", (num.snps.genos-num.nonsnp.col), " and in .dat: ", num.snps.dat, sep=""))

		if((num.snps.genos-num.nonsnp.col) != num.snps.dat) {
			print(paste("Error: the geno file ", files.genos[i], " has ", (num.snps.genos-num.nonsnp.col), " SNP columns; whereas its matching dat file ", files.dat[i], " has ", num.snps.dat, " entries.", sep=""))
			return(NULL)
		}

		i <- i + 1
	}
	
	# Do the check for #rows:
	num.rows.uniq <- unique(NUM.ROWS)
	if(length(num.rows.uniq) != 1) {
		print("Error: not all the geno files have the same number of rows:")
		j <- 1
		while(j <= length(num.rows.uniq)) {
			tmp.row <- num.rows.uniq[j]
			print(" ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ ")
			print(paste("NUM = ", tmp.row, sep=""))
			print(files.genos[which(NUM.ROWS==tmp.row)])
			j <- j + 1
		}
		return(NULL)
	}
}

out.name.genos <- paste(dir.out, file.out, sep="/")
out.name.dat <- paste(dir.out, dat.out, sep="/")

if(plan == FALSE) {
	# Cut the first and last non-SNP columns out of every GENO file, and save into tmp.
	# For one of the files, save the last columns only.
	# Then paste the tmps together, along with the last columns at the end.
	# Remove all the tmp files.
	tmp.names <- paste(dir.out, "/tmp.", files.genos, sep="")
	disease.name <- paste(tmp.names[1], ".tmp", sep="") # name of file containing the disease status
	count.total.cols <- 0 # For printing purposes only, to know how many cols to expect in final file
	i <- 1
	while(i <= length(files.genos)) {
		curr.file <- paste(dir.file, files.genos[i], sep="/")
		if(weak.check == FALSE) 
			col.genos <- NUM.GENOS[i]
		else {
			print(paste("Counting dim: ", curr.file, sep=""))
			num.dims <- get.data.dims(curr.file)
			col.genos <- num.dims$ncols - num.nonsnp.last.col
			NUM.ROWS[i] <- num.dims$nrows
			print(paste("num cols: ", (col.genos-num.nonsnp.col), sep=""))
			count.total.cols <- count.total.cols + col.genos
		}
		try(system(paste("cut -f ", (num.nonsnp.col+1), "-", col.genos, " ", curr.file, " > ", tmp.names[i], sep="")))
		if(i==1) { # Once copy the last set of columns into disease name
			try(system(paste("cut -f 1-", col.genos, " --complement ", curr.file, " > ", disease.name, sep="")))
		}
		i <- i + 1
	}

	if(weak.check == FALSE) {
		count.total.cols <- sum(NUM.GENOS)
	} else {
		# Here check that number of rows is consistent:
		num.rows.uniq <- unique(NUM.ROWS)
		if(length(num.rows.uniq) != 1) {
			print("Error: not all the geno files have the same number of rows:")
			j <- 1
			while(j <= length(num.rows.uniq)) {
				tmp.row <- num.rows.uniq[j]
				print(" ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ ")
				print(paste("NUM = ", tmp.row, sep=""))
				print(files.genos[which(NUM.ROWS==tmp.row)])
				j <- j + 1
			}
			return(NULL)
		}
	}
	count.total.cols <- count.total.cols + num.nonsnp.last.col - (length(files.genos)*num.nonsnp.col) # Add columns for the disease status, and subtract out the initial columns of all files
	print(paste("Total number of columns would be: ", count.total.cols, sep=""))

	# Paste all files along with disease status.
	to.paste <- paste(tmp.names, collapse=" ")
	to.paste <- paste("paste ", to.paste, " ", disease.name, " > ", out.name.genos, sep="")
	print(to.paste)
	try(system(to.paste))

	file.remove(tmp.names)
	file.remove(disease.name)

	# ************************************************************** #
	# Put together all the DAT files.
	i <- 2
	DATS <- read.table(paste(dir.dat, files.dat[1], sep="/"), header=FALSE, sep="\t", stringsAsFactors=FALSE)
	DATS[,1] <- rep(chroms.dat[1], nrow(DATS))

	while(i <= length(files.dat)) {
		dats <- read.table(paste(dir.dat, files.dat[i], sep="/"), header=FALSE, sep="\t", stringsAsFactors=FALSE)
        	dats[,1] <- rep(chroms.dat[i], nrow(dats))
		DATS <- rbind(DATS, dats)
		i <- i+1
	}
	write.table(DATS, file=out.name.dat, col.names=FALSE, row.names=FALSE, quote=FALSE)
	
	# Copy all .fam files from the input directory over to output directory:
	get.file.copy(dir.in=dir.file, dir.out=dir.out, ending=".fam", verbal=FALSE)

	return(out.name.genos)

} else { # in case of plan, simply print out the file match
	print("*****************************************")
	print(paste("Matching GENOS from dir: ", dir.file, sep=""))
	print(files.genos)
        print("*****************************************")
	print(paste("With DAT from directory: ", dir.dat, sep=""))
	print(files.dat)
        print("*****************************************")
	print("Would save GENO result into: ")
	print(out.name.genos)
	print("and DAT result into: ")
	print(out.name.dat)
	return(NULL)
}

}


# **************************************************************************
# Helper function to return the list of files in directory dir.file, whose names
# begin with prefix, and the name contains 'key' substring, and this list will be
# sorted in the order of the chromosome number.
# The file names in dir.file should be formatted as:
# <prefix><number><whatever> - where prefix is given; number is following right
#    after prefix; and whatever can be any string containing 'key'. 
#    Note: key can be an empty string - "". 
# Returns:
#    out$files - the list of file names, no path
#    out$chroms - the list of chromosome numbers.
# **************************************************************************
extract.chrom.dir <- function(dir.file, prefix, key, ending) {

all.files <- get.file.name(dir=dir.file, prefix=prefix, key=key, ending=ending)

if(length(all.files) == 0)
        return()

chroms.all <- get.chrom.num(all.files, prefix=prefix)

sorted <- sort(chroms.all, index.return=TRUE)
ordering<- sorted$ix

# Sort the file names in the numeric order
all.files <- all.files[ordering]

return(list(files=all.files, chroms=sorted$x))

}






