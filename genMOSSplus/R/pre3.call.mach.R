pre3.call.mach <- function(file.dat, file.ped, dir.file, ref.phase="", ref.legend="", dir.ref="", dir.out, out.prefix="result", chrom.num="", num.iters=2, num.subjects=200, step2.subjects=50, empty="0/0", resample=FALSE, mach.loc="/software/mach1") {

# Calls MACH1 program on file.ped and file.dat.
# MaCH1 can be run in 2 different ways: 
#  1. with hapmap, and 
#  2. without hapmap.
# If ref.phase and ref.legend are provided, then this program runs MaCH1 on file.ped with hapmap,
# resulting in file with as many SNPs as is present in the reference files.
# If no reference files ref.phase and ref.legend are provided, then the program runs
# MaCH1 without hapmap. 
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
# pre3.call.mach(file.dat="CGEM_Breast_chr19.removed.dat", file.ped="Reshape_CGEM_Breast_chr19CASE.removed.ped", dir.file="/home/briollaislab/olga/curr/data/pipe/f03_removed", ref.phase="genotypes_chr19_CEU_r22_nr.b36_fwd.phase", ref.legend="genotypes_chr19_CEU_r22_nr.b36_fwd_legend.txt", dir.ref="/home/briollaislab/olga/curr/data/pipe/f04_ref", dir.out="/home/briollaislab/olga/curr/data/pipe/f05_machout", out.prefix="CGEM_Breast_", chrom.num="19", num.iters=1, num.subjects=50, step2.subjects=25) 
#
# file.dat: the data file as required for MaCH1, of the format:
#         M SNP1
#         M SNP2
# file.ped: the pedegree data file. Can be for either CASE or CONTROL or both.
# dir.file: directory where file.dat and file.ped can be found.
# ref.phase: the name of the reference file, must have no missing values,
#         can be obtained from websites like:
#         http://hapmap.ncbi.nlm.nih.gov/downloads/phasing/2007-08_rel22/phased/
#         or similar/updated versions. 
#         No zip. Must be a normal and readable by R file. 
# ref.legend: the legend file for file.phase, obtained from same website. No zip.
# dir.ref: directory where ref.phase and ref.legend can be found.
# dir.out: directory to which output files should be saved. 
# out.prefix: the beginning string of how the output file should be named. 
# chrom.num: optionally a string denoting the chromosome number, for better naming of
#         intermediate files.
# num.iters: how many iterations MaCH1 should make in its first step to
#         estimate its model parameters.
# num.subjects: how many individuals in the sample should be used for model building
#         by the first step of MaCH1. The random subset of inidividuals will be chosen
#         by this program. Recommended number of subjects is 200-500.
#         Value <= 0 corresponds to using ALL the subjects in the dataset.
# step2.subjects: how many individuals should be processed at a time during the
#         second step of MaCH computation if run without hapmap.
#         Value <= 0 will use ALL the subjects in the dataset.
# 
# Outputs: 
#


# TODO: remove this line:
#source("call.mach.R")
#source("call.mach.nohap.R")
#source("machout2machin.R")
#source("get.file.copy.R")

if(missing(file.dat)) stop("Name of the .dat file is required.")
if(missing(file.ped)) stop("Name of the .ped file is required.")
if(missing(dir.file)) stop("Name of the directory that contains .dat file is required.")
if(missing(dir.out)) stop("Name of the output directory is required.")

# 0. Define full path names of all subdirectories
dir.work <- paste(dir.out, "/working", sep="")
dir.main <- paste(dir.work, "/dir_", substr(file.ped, 1, nchar(file.ped)-4), "_", num.subjects, "_", num.iters, "_", step2.subjects, sep="")
#dir.s1 <- paste(dir.main, "/s1_withhap", sep="")
#dir.s2 <- paste(dir.main, "/s2_converted", sep="")
#dir.s3 <- paste(dir.main, "/s3_nohap", sep="")

#from.dir <- dir.file # will be replaced if hapmap is run.
#from.ped <- file.ped

# 0. Create subdirectory structure:
if(!file.exists(dir.work))
	dir.create(dir.work)
if(!file.exists(dir.main))
	dir.create(dir.main)
#if(!file.exists(dir.s3))
#        dir.create(dir.s3)


# Erase everything from the directory if fresh start is needed:
if(resample==TRUE) {
	file.remove(dir(dir.main, full.names=TRUE))
}


# 1. IF HAPMAP: Call MaCH1 WITH hapmap, and then convert the output to readable input format.
if(ref.phase != "" && ref.legend != "" && file.exists(paste(dir.ref, ref.phase, sep="/")) && file.exists(paste(dir.ref, ref.legend, sep="/"))) {

	#if(!file.exists(dir.s1))
	#	dir.create(dir.s1)
	#if(!file.exists(dir.s2))
	#	dir.create(dir.s2)

	# Erase everything from 2 directories if fresh start is needed:
#	if(resample==TRUE) {
#		file.remove(dir(dir.s1, full.names=TRUE))
#		file.remove(dir(dir.s2, full.names=TRUE))
#	}

	prefix1 <- paste(out.prefix, chrom.num, sep="")	

	call.mach(file.dat=file.dat, file.ped=file.ped, dir.file=dir.file, ref.phase=ref.phase, ref.legend=ref.legend, dir.ref=dir.ref, prefix=prefix1, dir.out=dir.main, num.iters=num.iters, num.subjects=num.subjects, resample=resample, mach.loc=mach.loc)

	get.file.copy(dir.in=dir.main, dir.out=dir.out, prefix=prefix1, ending="mlgeno")

	#machout2machin(prefix=prefix1, dir.prefix=dir.main, dat=file.dat, name.ped=file.ped, dir.dat=dir.file, dir.out=dir.out, num.nonsnp.cols=5)

#	from.dir <- dir.s2
} else {

	# 2. ELSE: Call Mach1 without hapmap


	# Erase everything from last directory:
#	if(resample==TRUE)
#		file.remove(dir(dir.main, full.names=TRUE)) #dir.s3

	call.mach.nohap(file.dat=file.dat, dir.dat=dir.file, file.ped=file.ped, dir.file=dir.file, dir.out=dir.out, dir.debug=dir.main, prefix=paste(out.prefix, chrom.num, sep=""), num.iters=num.iters, num.subjects=num.subjects, step2.subjects=step2.subjects, resample=resample, empty=empty, mach.loc=mach.loc) #dir.debug=dir.s3

}

get.file.copy(dir.in=dir.file, dir.out=dir.out, fname=file.dat)
	#if(!file.exists(paste(dir.out, file.dat, sep="/")))
	#	ans <- file.copy(paste(dir.file, file.dat, sep="/"), paste(dir.out, file.dat, sep="/"))

}

