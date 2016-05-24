genos.clean.batch <- function(dir.ped, dir.dat=dir.ped, dir.out, prefix.ped, prefix.dat, key.ped="", key.dat="", ending.ped=".txt", ending.dat=".dat", num.nonsnp.col=2, num.nonsnp.last.col=1, letter.encoding=TRUE, save.ids.name="patients.fam") {
#
# For all the .ped files in the directory dir.ped, 
# remove all SNPs that contain at least one illegal value.
# Leaves the format of genotype data as is (i.e. no conversion to 3 levels).
# Similar to pre5.genos2numeric().  
# Geno values should use letters A, T, C, G if letter.encoding=TRUE.
#
# Example:
# genos.clean.batch(dir.ped="/home/briollaislab/olga/curr/data/pipe/f06_combined", dir.out="/home/briollaislab/olga/curr/data/pipe/f07_numeric", prefix.ped="CGEM_Breast_", prefix.dat="CGEM_Breast_chr")
#
# dir.ped: directory where file.ped can be found.
# dir.dat: directory where file.dat can be found. Defaults to dir.ped.
# dir.out: output directory to which resulting file should be saved.
# prefix.ped: the beginning of the file name for the pedegree file (up until chrom number)
# prefix.dat: beginning of the file name for .dat file
# key.ped: any keyword in the name of the pedegree file that distinguishes it from others
# key.dat: any keyword in the name of the .dat file that distinguishes it from others
# ending.case: MUST be the extension of the pedegree file, including the dot ".".
# ending.dat: MUST be the extension of the .dat file, including the dot ".".
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
# save.ids.name: the name of the file to which all patient IDs should be saved. 
#
# Outputs:
#
# <file.ped>_clean<ending.ped> - in dir.out directory, the resultant binary file:
#      the SNP columns + last columns (but no user IDs will be recorded).
# <file.dat>_clean.dat - in dir.out directory, the corresponding .dat file, will be different
#      from original <file.dat> if remove.bad.genos=TRUE.
# <save.ids.name> - the file containing patient IDs.

if(missing(dir.ped)) stop("Name of input directory with pedegree files must be provided.")
if(missing(dir.out)) stop("Name of output directory must be provided")
if(missing(prefix.ped)) stop("Prefix of the pedegree file name must be provided.")
if(missing(prefix.dat)) stop("Prefix of the .dat file name must be provided.")

# TODO: remove this line:
#source("genos.clean.R")
#source("get.file.name.R")
#source("get.chrom.num.R")


# *******************************************
# 1. Obtain all ped and .dat files
all.ped <- get.file.name(dir=dir.ped, prefix=prefix.ped, key=key.ped, ending=ending.ped)
all.dat <- get.file.name(dir=dir.dat, prefix=prefix.dat, key=key.dat, ending=ending.dat)

if(length(all.ped) == 0 || length(all.dat) == 0)
	return()


# *******************************************
# 3. Match ped and .dat and run the pre5.genos2numeric()
 
chroms.ped <- get.chrom.num(all.ped, prefix=prefix.ped)
chroms.dat <- get.chrom.num(all.dat, prefix=prefix.dat)

chroms.common <- sort(intersect(chroms.ped, chroms.dat))

i <- 1
while (i <= length(chroms.common)) {
	curr.chrom <- chroms.common[i]

	curr.ped <- all.ped[match(curr.chrom, chroms.ped)]
	curr.dat <- all.dat[match(curr.chrom, chroms.dat)]

	print(paste("Checking: ", curr.ped, " + ", curr.dat, sep=""))

	# Save the patients only for last chromosome (hopefully it's a small one).
	if(i != length(chroms.common))
		genos.clean(file.ped=curr.ped, ending.ped=ending.ped, dir.ped=dir.ped, file.dat=curr.dat, ending.dat=ending.dat, dir.dat=dir.dat, dir.out=dir.out, num.nonsnp.col=num.nonsnp.col, num.nonsnp.last.col=num.nonsnp.last.col, letter.encoding=letter.encoding, save.ids.name="")
	else
		genos.clean(file.ped=curr.ped, ending.ped=ending.ped, dir.ped=dir.ped, file.dat=curr.dat, ending.dat=ending.dat, dir.dat=dir.dat, dir.out=dir.out, num.nonsnp.col=num.nonsnp.col, num.nonsnp.last.col=num.nonsnp.last.col, letter.encoding=letter.encoding, save.ids.name=save.ids.name)

        i <- i + 1
}

}






