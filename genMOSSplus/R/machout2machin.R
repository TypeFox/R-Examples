machout2machin <- function(prefix, dir.prefix, dat, name.ped, dir.dat, dir.out, num.nonsnp.cols=5) {
#
# After MaCH1 is run with Hapmap, restore the computed values back into
# the MaCH input format. Will create a .ped file .
# Note: MaCH output will have genos in format "A/A" (with slash separator)
#    Before calling, make sure that ped file has the same style of separator. 
# Assumes that same patients appear in both .ped and prefix.mlgeno.
# 
# Sample run:
#
#
# Simple example run:
#
# machout2machin(prefix="a", dir.prefix="/home/briollaislab/olga/curr/data/mach1in/afterhap/try", dat="c.dat", name.ped="c.ped", dir.dat="/home/briollaislab/olga/curr/data/mach1in/afterhap/try", dir.out="/home/briollaislab/olga/curr/data/mach1in/afterhap/try/out")
#
# dat: the original .dat file in format as required by MaCH1 input.
# dir.dat: directory where .dat and .ped files can be found. 
# prefix: the prefix given to MaCH1 (second step of its computation).
# name.ped: the original pedegree file that was used to produce results beginning
#    with "prefix". Must contain genos of same format as prefix.mlgeno (i.e. slash separated);
#    spacing between entries is expected to be a tab.
# dir.prefix: directory where the MaCH1 output with the prefix can be found.
# num.nonsnp.cols: The number of non-SNP columns in the ped file.
#
# Returns the name of the output file, no path.
#
# Outputs result in dir.out directory:
#   - <name.ped>.conv.ped - the converted version; corresponds to the given "dat" file. 
#

# Read prefix.mlinfo and mlgeno 
mlinfo <- read.table(paste(dir.prefix, "/", prefix, ".mlinfo", sep=""), header=T, sep="\t", stringsAsFactors=FALSE)
mlgeno <- read.table(paste(dir.prefix, "/", prefix, ".mlgeno", sep=""), header=F, sep=" ", stringsAsFactors=FALSE)
ped <- read.table(paste(dir.dat, "/", name.ped, sep=""), header=F, sep="\t", stringsAsFactors=FALSE)

if(dat != "") {

	Dat <- read.table(paste(dir.dat, dat, sep="/"), header=F, sep="\t", stringsAsFactors=FALSE)
	numD <- nrow(Dat)
	D <- Dat[, 2]

	# **************************************************
	# 1. Deal with columns (SNP ids)
	# **************************************************
	IndexesAns <- array(0, c(numD))
	numIA <- 1 # up to how many values IndexesAns has been filled.
	IndexesPed <- array(0, c(numD)) # at which id in D SNPs appear in mlinfo.

	i <- 1
	while(i <= numD) {
        	inL <- match(D[i], mlinfo[,1])
	        if (!is.na(inL)) {
                	IndexesAns[numIA] <- inL
			IndexesPed[numIA] <- i
                	numIA <- numIA + 1
		}
        	i <- i + 1
	}
	# In case of any missing indexes, shorten the IA and IP arrays:
        IndexesAns <- IndexesAns[1:(numIA-1)]
	IndexesPed <- IndexesPed[1:(numIA-1)]

	# **************************************************
	# 2. Deal with rows (align patients)
	# **************************************************

	# From mlgeno extract patient name (first half of the string from 1st column):
	patient.name <- unlist(strsplit(mlgeno[,1], split="->"))[seq(1,nrow(mlgeno)*2,by=2)]
	# Sort both mlgeno and .ped files s.t. the ordering of patients is the same:
	sort.mlgeno.id <- sort(patient.name, index.return=T)$ix
	sort.ped.id <- sort(ped[,1], index.return=T)$ix

	mlgeno <- mlgeno[sort.mlgeno.id,]
	ped <- ped[sort.ped.id,]

	# TODO: Check that patients are the same
	# 

	# **************************************************
	# 3. Fill up the new .ped file information:
	# **************************************************
	
	# copy imputed information for SNPs that appeared in phase file.
	# Note: ped file has 5 leadning non-SNP columns, whereas mlgeno has only 2.
	ped[IndexesPed+num.nonsnp.cols] <- mlgeno[,(IndexesAns+2)]
	

	name.short <- paste(substr(name.ped, 1, nchar(name.ped)-4), ".conv.ped", sep="")
	name.out <- paste(dir.out, "/", name.short, sep="")
	write.table(ped, file=name.out, quote=FALSE, sep="\t", row.names=F, col.names=F)

	return(name.short)
}

return("")

}
