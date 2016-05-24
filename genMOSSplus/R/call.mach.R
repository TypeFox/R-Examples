call.mach <- function(file.dat, file.ped="", dir.file, ref.phase, ref.legend, dir.ref, prefix="ans", dir.out=".", num.iters=2, num.subjects=200, resample=FALSE, hapmapformat=FALSE, mach.loc="/software/mach1") {
#
# Calls MaCH1 program with Hapmap on the file.ped and file.dat, making use of
# two reference files: ref.phase and ref.legend. 
#
# hapmapformat - refers to the use of old file format with .legend.txt.
#     If you use the 1000G haplotype dataset from MaCH's website, this 
#     parameter defaults to FALSE.
#
# The MaCH1 algorithm requires 2 steps to be performed.
# The first step of MaCH1 will be run on num.subjects randomly chosen from the set.
# The file with randomly chosen individuals will be saved as file.ped.<num.subjects>.ped
# in dir.file directory. If the file already exists for this num.subjects, the old file
# will be used if resample=F. 
# If resample=T then old files will be ignored, and new sampling will take place.  
# The second step will be run on all the subjects.
#
# Sample run:
#
# call.mach(file.dat="Clean_genotypes_CASE_chr22.dat", file.ped="Clean_genotypes_CASE_chr22.ped", dir.file="/home/briollaislab/olga/curr/data/mach1in", ref.phase="genotypes_chr22_CEU_r22_nr.b36_fwd.phase", ref.legend="genotypes_chr22_CEU_r22_nr.b36_fwd_legend.txt", dir.ref="/home/briollaislab/olga/curr/data/mach1in/phases", prefix="ans", dir.out="/home/briollaislab/olga/preprocess/mach/iter2_full_prog", num.subjects=25)
#
# file.dat: the data file as required for MaCH1, of the format:
#         M SNP1
#         M SNP2
# file.ped: the pedegree data file. Default value "" assumes that file.ped has the
#         same name as file.dat, only with extension ".ped".
# dir.file: the directory where file.dat and file.ped can be found
# ref.phase: file that includes a set of reference haplotypes,
#         as required for MaCH1 input if hapmap is used.
#         Make sure the file is NOT zipped (to unzip: "gzip -d file.gz")
# ref.legend: file that lists all the markers that appear in the phased haplotypes.
#         This file should NOT be zipped.
# dir.ref: directory where file.phase and file.legend are located.
# prefix: what prefix MaCH1 should use to output its files; if num.subjects > 0
#       then the num.subjects will be appended to the prefix. 
#       You may need to include the chromosome number and CASE/CONTROL status into the prefix.
# dir.out: the directory into which MaCH1 output should go.
# num.iters: how many iterations MaCH1 should make in its first step to 
#       estimate its model parameters. 
# num.subjects: how many individuals in the sample should be used for model building
#       by the first step of MaCH1. The random subset of inidividuals will be chosen
#       by this program. Recommended number of subjects is 200-500.
#       Value <= 0 corresponds to using ALL the subjects in the dataset.
# resample: whether or not to overwrite the existing file containing the
#       num.subjects entries produced by previous runs of this algorithm
#       with same file.dat, file.ped and num.subjects parameters.
#       By default, if the subjects have been sampled before, they are used.
# mach.loc: the location directory where "mach" executable can be found.
# 
#

# TODO: remove this line:
# source("rand.ints.R")

# *************************************************************************
# 1. Randomly choose a subset of num.subjects and create short versions of
#   file.ped in dir.file directory (will remove later).
#   Note: file.dat remains the same, as it lists SNPs and not individuals.
# *************************************************************************

out.ped <- paste(dir.file, file.ped, sep="/")
P <- read.table(out.ped, header=F, sep="\t", stringsAsFactors=FALSE)
numP <- nrow(P)
if (num.subjects > 0) {
	#prefix <- paste(prefix, ".", num.subjects, sep="")
	out.ped <- paste(dir.file, "/", substr(file.ped, 1, (nchar(file.ped)-4)), ".", num.subjects, ".ped", sep="")

	if (!file.exists(out.ped) || resample == TRUE) {
		chosen <- rand.ints(num.subjects, min=1, max=numP) #unique(as.integer(runif(num.subjects, min=1, max=numP)))
		write.table(P[chosen,], file=out.ped, col.names=F, row.names=F, quote=F, sep="\t")
		print(paste("Created the file: ", out.ped, sep=""))
	}
}

# *************************************************************************
# 2. Call MaCH1
# *************************************************************************

prefix1 <- paste(prefix, ".1", sep="")

step1 <- paste("time ", mach.loc, "/mach1 -d ", dir.file, "/", file.dat, " -p ", out.ped, " -s ", dir.ref, "/", ref.legend, " -h ", dir.ref, "/", ref.phase, " -r ", num.iters, " --autoflip --prefix ", dir.out, "/", prefix1, sep="")

step2 <- paste("time ", mach.loc, "/mach1 -d ", dir.file, "/", file.dat, " -p ", dir.file, "/", file.ped, " -s ", dir.ref, "/", ref.legend, " -h ", dir.ref, "/", ref.phase, " --crossover ", dir.out, "/", prefix1, ".rec --errormap ", dir.out, "/", prefix1, ".erate --greedy --autoflip --mle --mldetails --prefix ", dir.out, "/", prefix, sep="")

if(hapmapformat) {
	step1 <- paste(step1, " --hapmapFormat", sep="")
	step2 <- paste(step2, " --hapmapFormat", sep="")
}
#step1 <- paste("time ", mach.loc, "/mach1 -d ", dir.file, "/", file.dat, " -p ", out.ped, " -s ", dir.ref, "/", ref.legend, " -h ", dir.ref, "/", ref.phase, " --hapmapFormat --greedy -r ", num.iters, " --autoflip --prefix ", dir.out, "/", prefix1, sep="")

#step2 <- paste("time ", mach.loc, "/mach1 -d ", dir.file, "/", file.dat, " -p ", dir.file, "/", file.ped, " -s ", dir.ref, "/", ref.legend, " -h ", dir.ref, "/", ref.phase, " --hapmapFormat --crossover ", dir.out, "/", prefix1, ".rec --errormap ", dir.out, "/", prefix1, ".erate --greedy --autoflip --mle --mldetails --prefix ", dir.out, "/", prefix, sep="")

#step.call <- step1

step.call <- paste(step1, " ; ", step2)

print(step.call)
p <- proc.time()
system(step.call)
print(proc.time() - p)

}


