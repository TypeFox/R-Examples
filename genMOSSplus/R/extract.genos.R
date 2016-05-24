extract.genos <- function(mach.prefix, dir.prefix, file.dat, dir.dat, file.legend, dir.legend, dir.out="", verbose=TRUE) { 
# Extracts the SNPs that appear in file.dat from the resultant <mach.prefix>.mlgeno file,
# creating new file named Filled_<file.dat>.ped in dir.dat directory.
# This is meant to be processing right AFTER calling MaCH1 algorithm with hapmap
# using .legend and .phase files that contain more SNPs than file.dat requires.
#
# Sample call:
# extract.genos("some200", dir.prefix="/home/briollaislab/olga/curr/data/mach1out", file.dat="Clean_genotypes_CASE_chr22.dat", dir.dat="/home/briollaislab/olga/curr/data/mach1in", file.legend="genotypes_chr22_CEU_r22_nr.b36_fwd_legend.txt", dir.legend="/home/briollaislab/olga/curr/data/mach1in/phases", dir.out="/home/briollaislab/olga/curr/data/mach1post")
#
# mach.prefix: the prefix given to MaCH1 followed by the num.subjects parameter,
#         s.t. MaCH1 creates <prefix><num.subjects>.mlgeno and other output files.
#         So when calling call.mach(prefix="some", num.subjects=200),
#         then in this function you should use: mach.prefix="some200".
# dir.prefix: directory where the <mach.prefix>.mlgeno can be found
#
# file.dat: the data file as required for MaCH1, of the format:
#         M SNP1
#         M SNP2
#         The name is expected to begin with "Clean_".
# dir.dat: directory where file.dat can be found
#
# file.legend: file that lists all the markers that appear in the phased haplotypes.
#         This file should NOT be zipped.
# dir.legend: directory where file.legend is located.
# dir.out: directory to which output files should go. Default of "" means dir.dat.
# verbose: Set it to FALSE if you wish to suppress most output. 
#
# Directory names should not end with "/"
#
# Example run:
#
# extract.genos(
#
#
# Outputs: 
#
# - Filled_<file.dat>.ped - the modified <mach.prefix>.mlgeno file containing only the
#    SNPs that appear both in file.dat and file.legend, saved in dir.out directory;
#    the <file.dat> will not have "Clean_" prefix anymore.
# - Filled_<file.dat> - the modified .dat file containing the list of SNPs that
#    appear in both file.dat and file.legend; these SNPs correspond to the
#    columns in Filled_<file.dat>.ped. Note: <file.dat> will not have "Clean_" prefix.
#

if(dir.out == "")
	dir.out <- dir.dat

Dat <- read.table(paste(dir.dat, file.dat, sep="/"), header=F, sep="\t", stringsAsFactors=FALSE)
numD <- nrow(Dat)
D <- Dat[, 2]

L <- read.table(paste(dir.legend, file.legend, sep="/"), header=T, sep="\t", stringsAsFactors=FALSE)

M <- read.table(paste(dir.prefix, "/", mach.prefix, ".mlgeno", sep=""), header=F, sep=" ", stringsAsFactors=FALSE)

print(paste("Number of SNPs in .dat file: ", numD, sep=""))
IndexesAns <- array(0, c(numD))
numIA <- 1 # up to how many values IndexesAns has been filled.
IndexesD <- array(0, c(numD))
numID <- 1 # up to how many values IndexesD have been filled.

i <- 1

while(i <= numD) {
	inL <- match(D[i], L[,1])
	if (is.na(inL)) {
		if(verbose)
			print(paste("SNP ", D[i], " from .dat file is not present in .legend file.", sep=""))
	}
	else {
		IndexesAns[numIA] <- inL
		numIA <- numIA + 1
		IndexesD[numID] <- i
		numID <- numID + 1
	}

	i <- i + 1
}

#print(paste("Total number of missing SNPs: ", numAbsent, "/", numD, " = ", 100*(numAbsent/numD), "%", sep=""))

# **********************************************************************************

if(numIA < 2) {
	print("Warning: NO SNPs in .dat correspond to any SNPs in .legend.")
	print(" Make sure that the files correspond:")
	print(paste(".dat: ", file.dat, sep=""))
	print(paste(".legend: ", file.legend, sep=""))

} else {
	# In case of any missing indexes, shorten the IA and ID arrays:
	IndexesAns <- IndexesAns[1:(numIA-1)]
	IndexesD <- IndexesD[1:(numID-1)]	

	out.Dname <- out.fname <- substr(file.dat, 7, nchar(file.dat))	# remove 'Clean_"
	out.Dname <- paste(dir.out, "/Filled_", out.Dname, sep="")       # prepend dir and 'Reshaped_" start
	out.Mname <- paste(substr(out.Dname, 1, nchar(out.Dname)-3), "ped", sep="")

	# Extract the IA indicies from .mlgeno file and from .dat files:
	AD <- D[IndexesD]
	AM <- M[, c(1, IndexesAns+2)]
	odds <- seq(1, nrow(AM)*2, by=2)		# odd numbers: 1, 3, 5, ...
	AM[,1] <- unlist(strsplit(AM[,1], "->"))[odds]	# Split the "id->id" on "->" and
							# take out the first value
	write.table(AD, file=out.Dname, col.names=F, row.names=F, quote=F, sep="\t")
	write.table(AM, file=out.Mname, col.names=F, row.names=F, quote=F, sep="\t")

}

print("done!")

}

