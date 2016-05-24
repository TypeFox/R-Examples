pre4.combine.case.control.batch <- function(dir.file, dir.out=dir.file, prefix.case, prefix.control, prefix.out, key.case="", key.control="", ending.case=".mlgeno", ending.control=".mlgeno", separ=" ") {
#        
# Combines CASE and CONTROL files into one file, and appends disease status as the last column.
# The disease status is encoded as 1 for CASE and 0 for CONTROL.
# If you run genos2numeric() AFTER this function, don't forget to specify that
#   there is 1 ending non-SNP column.
# Example run:
#
#
# dir.file: the directory where CASE and CONTROL input files can be found.
# dir.out: the directory to which output files should go
# prefix.case: the beginning of the file name of CASE file (up until the chromosome number)
# prefix.control: the beginning of the file name of CONTROL file
# prefix.out: the beginning of the string that should be used for the output file name.
# key.case: any keyword in the name of the CASE file that distinguishes it from others
# key.control: any keyword in the name of the CONTROL file
# ending.case: the extension of the CASE file
# ending.control: the extension of the CONTROL file
# spear: the separator used in CASE and CONTROL files
#
# Outputs:
#
# - <dir.out>/<prefix.out><chromnum>.txt - the file containing both CASE and CONTROL values,
#     with the disease status as the last column.
#


if(missing(dir.file)) stop("Name of input directory must be provided.")
if(missing(prefix.out)) stop("Prefix of the output file name must be provided.")
if(missing(prefix.case)) stop("Prefix of the CASE file name must be provided.")
if(missing(prefix.control)) stop("Prefix of the CONTROL file name must be provided.")

# TODO: remove this line:
#source("pre4.combine.case.control.R")
#source("get.file.name.R")
#source("get.chrom.num.R")


# *******************************************
# 1. Obtain all .dat, CASE, and CONTROL files
all.case <- get.file.name(dir=dir.file, prefix=prefix.case, key=key.case, ending=ending.case)
all.control <- get.file.name(dir=dir.file, prefix=prefix.control, key=key.control, ending=ending.control)

if(length(all.case) == 0 || length(all.control) == 0)
	return()


# *******************************************
# 3. Match CASE and CONTROL and run the pre4.combine.case.control()
 
chroms.case <- get.chrom.num(all.case, prefix=prefix.case)
chroms.control <- get.chrom.num(all.control, prefix=prefix.control)

chroms.common <- intersect(chroms.case, chroms.control)

i <- 1
while (i <= length(chroms.common)) {
	curr.chrom <- chroms.common[i]

	curr.case <- all.case[match(curr.chrom, chroms.case)]
	curr.control <- all.control[match(curr.chrom, chroms.control)]

	curr.name <- paste(prefix.out, curr.chrom, ".txt", sep="")

	print(paste("Combining: ", curr.case, " + ", curr.control, sep=""))

	pre4.combine.case.control(case.file=curr.case, control.file=curr.control, dir.file=dir.file, name.out=curr.name, dir.out=dir.out, separ=separ)

        i <- i + 1
}

}






