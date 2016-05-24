pre3.call.mach.batch <- function(dir.file, dir.ref="", dir.out, prefix.dat, prefix.ped, prefix.phase="", prefix.legend=prefix.phase, prefix.out="result", key.dat="", key.ped="", key.phase="", key.legend="", ending.dat=".dat", ending.ped=".ped", ending.phase=".phase", ending.legend="legend.txt", chrom.num, num.iters=2, num.subjects=200, step2.subjects=50, empty="0/0", resample=FALSE, mach.loc="/software/mach1") {

# Calls MaCH1 program on ONE chomosome ONLY, specified by 'chrom.num'.
# From approapriate directories, extracts the names of the files corresponding to given
# chromosome number, and runs MaCH1.
# MaCH1 can be run in 2 different ways:
#  1. with hapmap, and
#  2. without hapmap.
# This program first runs MaCH1 on file.ped with hapmap to fill in missing values for
# those SNPs that exist in the reference file; and then MaCH1 is run on the result without
# hapmap to fill in all the remaining missing values.
# If no reference file prefixes prefix.phase and prefix.legend are provided, then the program runs
# MaCH1 without hapmap only.
#
# Also, the MaCH1 algorithm requires 2 steps to be performed.
# (a) The first step of MaCH1 will be run on num.subjects randomly chosen from the set.
#     The file with randomly chosen individuals will be saved as file.ped.<num.subjects>.ped
#     in dir.file directory. If the file already exists for this num.subjects, the old file
#     will be used if resample=F.
#     If resample=T then old files will be ignored, and new sampling will take place.
#     This is useful if step1 runs well, but step2 crashes, then
#     re-calling this function will not waste time on re-running step1 over again.
# (b) 1. The second step will be run on all the subjects, when run with hapmap.
#     2. The second step without Hapmap takes exponentially long wrt number of subjects processed.
#        Thus the second step without hapmap will be run on bunches of subjects, step2.subjects at a time.
#
# Example run:
#
# pre3.call.mach.batch(dir.file="/home/briollaislab/olga/curr/data/pipe/f03_removed", dir.ref="/home/briollaislab/olga/curr/data/pipe/f04_ref", dir.out="/home/briollaislab/olga/curr/data/pipe/f05_machout", prefix.dat="CGEM_Breast_chr", prefix.ped="Reshape_CGEM_Breast_chr", prefix.phase="genotypes_chr", prefix.out="CGEM_Breast_CONTROL", key.ped="CONTROL", chrom.num=19)
#
# dir.file: directory where file.dat and file.ped can be found.
# dir.ref: directory where ref.phase and ref.legend can be found.
# dir.out: directory to which output files should be saved.
#
# prefix.dat: the prefix of the name of the data file as required for MaCH1, of the format:
#         M SNP1
#         M SNP2
#         Note: prefix is the string up until the chromosome number, so for file named
#         'CGEM_chr_12_reshape.dat', prefix='CGEM_chr_'
# prefix.ped: the prefix of the name of pedegree data file.
# prefix.phase: the prefix of the name of the the name of the reference file. 
#         The file must have no missing values, can be obtained from websites like:
#         http://hapmap.ncbi.nlm.nih.gov/downloads/phasing/2007-08_rel22/phased/
#         or similar/updated versions.
#         No zip. File must be normal and readable by R.
# prefix.legend: the prefix of the name of legend file for the phase file, 
#         obtained from same website. No zip. The prefix defaults to prefix.phase.
# prefix.out: the beginning string of how the output file should be named.
#
# key.dat: any keyword in the name of the data file that helps distinguish it from 
#         other files in the dir.file directory.
# key.ped: any keyword in the name of the pedegree file.
# key.phase: any keyword in the name of the phrase file.
# key.legend: any keyword in the name of the legend file.
#
# ending.dat: a string with which the name of the dat file ends.
# ending.ped: a string with which the name of the pedegree file ends.
# ending.phase: a string with which the name of the phase file ends.
# ending.legend: a string with which the name of the legend file ends.
# 
# chrom.num: Chromosome number for which processing should be done. 
# num.iters: how many iterations MaCH1 should make in its first step to
#         estimate its model parameters.
# num.subjects: how many individuals in the sample should be used for model building
#         by the first step of MaCH1. The random subset of inidividuals will be chosen
#         by this program. Recommended number of subjects is 200-500.
#         Value <= 0 corresponds to using ALL the subjects in the dataset.
# step2.subjects: how many individuals should be processed at a time during the
#         second step of MaCH computation.
#         Value <= 0 will use ALL the subjects in the dataset.
#
# Outputs:
#
#        

if(missing(dir.file)) stop("Name of input directory for .dat and .ped files must be provided.")
if(missing(dir.out)) stop("Name of output directory must be provided.")
if(missing(prefix.dat)) stop("Prefix of the .dat file name must be provided.")
if(missing(prefix.ped)) stop("Prefix of the .ped file name must be provided.")
if(missing(chrom.num)) stop("Chromosome number must be provided.")

# TODO: remove this line:
#source("pre3.call.mach.R")
#source("get.file.name.R")

# *******************************************
# 1. Obtain the .dat and .ped files
all.dat <- get.file.name(dir=dir.file, prefix=paste(prefix.dat, chrom.num, sep=""), key=key.dat, ending=ending.dat)
all.ped <- get.file.name(dir=dir.file, prefix=paste(prefix.ped, chrom.num, sep=""), key=key.ped, ending=ending.ped)

if(length(all.dat) == 0 || length(all.ped) == 0)
	return()

if(length(all.dat) != 1) {
	print(paste("Warning: more than one .dat file for chromosome ", chrom.num, " was detected, using the last:", sep=""))
	print(all.dat)
}

use.dat <- all.dat[length(all.dat)]	# use last .dat file (ignore all 2-digit chrom nums)
use.ped <- all.ped[length(all.ped)]	# use last .ped file (ignore all small splitted ones)


# *******************************************
# 2. Obtain the reference files: phase and legend, if any:

if(prefix.phase != "") {
	all.phase <- get.file.name(dir=dir.ref, prefix=paste(prefix.phase, chrom.num, sep=""), key=key.phase, ending=ending.phase)
	all.legend <- get.file.name(dir=dir.ref, prefix=paste(prefix.legend, chrom.num, sep=""), key=key.legend, ending=ending.legend)

	if(length(all.phase) == 0 || length(all.legend)==0)
		return()

	# Call MaCH with hapmap:	
	pre3.call.mach(file.dat=use.dat, file.ped=use.ped, dir.file=dir.file, ref.phase=all.phase[1], ref.legend=all.legend[1], dir.ref=dir.ref, dir.out=dir.out, out.prefix=prefix.out, chrom.num=chrom.num, num.iters=num.iters, num.subjects=num.subjects, step2.subjects=step2.subjects, empty=empty, resample=resample, mach.loc=mach.loc)

} else  { # Call MaCH without hapmap:

	pre3.call.mach(file.dat=use.dat, file.ped=use.ped, dir.file=dir.file, ref.phase="", ref.legend="", dir.ref="", dir.out=dir.out, out.prefix=prefix.out, chrom.num=chrom.num, num.iters=num.iters, num.subjects=num.subjects, step2.subjects=step2.subjects, empty=empty, resample=resample, mach.loc=mach.loc)
}


}






