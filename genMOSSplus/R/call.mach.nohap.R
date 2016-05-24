call.mach.nohap <- function(file.dat, dir.dat, file.ped, dir.file, dir.out, dir.debug, prefix="ans", num.iters=2, num.subjects=200, step2.subjects=50, resample=FALSE, mach.loc="/software/mach1", empty="0/0") {
#
# Calls MaCH1 program with NO Hapmap on file.dat and on file.ped in groups of 
#    step2.subjects at a time for step 2 of MaCH1. Preferably all empty genos are
#    already removed. 
#
# The MaCH1 algorithm requires 2 steps to be performed.
# The first step of MaCH1 will be run on num.subjects randomly chosen from the set.
# The file with randomly chosen individuals will be saved as file.ped.<num.subjects>.ped
# in dir.file directory. If the file already exists for this num.subjects, the old file
# will be used if resample=F. 
# If resample=T then old files will be ignored, and new sampling will take place.  
# The step1 of MaCH will only be run if resample=T, or if the files that MaCH1 produces do not exist yet.
# Thus if step1 runs well, but step2 crashes, re-calling this function will not waste time
# on re-running step1 over again. 
#
# The second step without Hapmap takes exponentially long wrt number of subjects processed.
# Thus the second step will be run on bunches of subjects, step2.subjects at a time.
#
# Sample run:
#
# call.mach.nohap(file.dat="Clean_genotypes_CASE_chr22.noh.removed.dat", file.ped="Reshaped_genotypes_CONTROL_chr22.noh.removed.ped", dir.file="/home/briollaislab/olga/preprocess/mach/nohapmap/varyPple", prefix="test", dir.out="/home/briollaislab/olga/preprocess/mach/nohapmap/varyPple")
# 
# file.dat: the data file as required for MaCH1, of the format:
#         M SNP1
#         M SNP2
# file.ped: the pedegree data file. Default value "" assumes that file.ped has the
#         same name as file.dat, only with extension ".ped".
# dir.file: the directory where file.dat and file.ped can be found
# prefix: what prefix MaCH1 should use to output its files; if num.subjects > 0
#       then the num.subjects will be appended to the prefix.
# dir.out: the directory into which MaCH1 output should go.
# num.iters: how many iterations MaCH1 should make in its first step to 
#       estimate its model parameters. 
# num.subjects: how many individuals in the sample should be used for model building
#       by the first step of MaCH1. The random subset of inidividuals will be chosen
#       by this program. Recommended number of subjects is 200-500.
#       Value <= 0 corresponds to using ALL the subjects in the dataset.
# step2.subjects: how many individuals should be processed at a time during the
#       second step of MaCH computation. 
#       Value <= 0 will use ALL the subjects in the dataset.
# resample: whether or not to overwrite the existing file containing the
#       num.subjects entries produced by previous runs of this algorithm
#       with same file.dat, file.ped and num.subjects parameters. 
#       By default, if the subjects have been sampled before, they are used.
# mach.loc: the location directory where "mach" executable can be found.
# empty: how an empty entry for a SNP is represented in the file.ped
#
# Returns: the name of the resultant .mlgeno file, in format as outputed by MaCH step 2.
#
# ************************************************************** 
# OUTPUTS:
# 
# The output file that you NEED after it is all done will be in dir.out directory:
# <prefix>.<num.subjects>.<step2.subjects>.mlgeno   if step2.subjects > 0
# <prefix>.<num.subjects>.mlgeno 		if step2.subjects <= 0
#   - contains the full geno information for all the patients and all the SNPs.
#
# The rest of the output files are useful for debugging, and
# analyzing how the imputational process went, but is unnecessary otherwise 
# and can be deleted soon afterwards:
#
# The following helper files in "dir.debug" directory:
# <file.ped>.<num.subjects>.ped  - the randomly extracted "num.subjects" individuals
#       necessary for calling step1 of MaCH. 
# <file.ped>.group.<step2.subjects>.i.ped 
#       - contains the short .ped files with "step2.subjects" subjects each,
#       where "file.ped" does not have .ped extension, and "i" is the index
#       over the groups of subjects.
#
# And after performing Step1, in "dir.debug" directory:
# <prefix>.<num.subjects>.1.erate
# <prefix>.<num.subjects>.1.rec
#
# And after performing Step2, in "dir.debug" directory:
# <prefix>.<num.subjects>.<step2.subjects>.i.mldose
# 	- .mlgeno
#	- .mlinfo
#	- .mlprob
# 	- .mlqc
#	- .erate
#	- .rec
# The reason for naming output files of step2 this way:
# The name begins with user's defined prefix; 
# followed by the number of subjects that has been used for step1 
# (to encode the type of quality incorporated during step1);
# followed by the number of subjects considered at a time during Step2 execution;
# then 'i' is the consecutive number of the group/bundle of individuals in the sequence;
# ending with MaCH1 step2 default names.
#


# TODO: remove this:
#source("call.mach.nohap.step1.R")
#source("call.mach.nohap.step2.R")
#source("rand.ints.R")


# *************************************************************************
# 1. Randomly choose a subset of num.subjects and create short versions of
#   file.ped in dir.debug directory (to be removed later).
#   Note: file.dat remains the same, as it lists SNPs and not individuals.
# *************************************************************************

out.ped <- paste(dir.file, file.ped, sep="/")
P <- read.table(out.ped, header=F, sep="\t", stringsAsFactors=FALSE)
numP <- nrow(P)
if (num.subjects > 0) {
	prefix <- paste(prefix, ".", num.subjects, sep="")
	out.ped <- paste(dir.debug, "/", substr(file.ped, 1, (nchar(file.ped)-4)), ".", num.subjects, ".ped", sep="")

	if (!file.exists(out.ped) || resample == TRUE) {
		chosen <- rand.ints(num.subjects, min=1, max=numP) #unique(as.integer(runif(num.subjects, min=1, max=numP)))
		write.table(P[chosen,], file=out.ped, col.names=F, row.names=F, quote=F, sep="\t")
		print(paste("Created the file: ", out.ped, sep=""))
	}
}

# *************************************************************************
# 2. Call step1 of MaCH, if either we're forced by 'resample', or
# one of the .erate or .rec files do not exist in the dir.debug directory.
# *************************************************************************

prefix1 <- paste(prefix, ".1", sep="")

if(resample == TRUE || !file.exists(paste(dir.debug, "/", prefix1, ".erate", sep="")) || !file.exists(paste(dir.debug, "/", prefix1, ".rec", sep="")))
	call.mach.nohap.step1(file.dat=file.dat, dir.dat=dir.dat, out.ped=out.ped, num.iters=num.iters, prefix1=prefix1, dir.out=dir.debug, mach.loc=mach.loc)


# *************************************************************************
# 3. Split dataset into many files with 'step2.subjects' individuals each. 
# *************************************************************************

subjects.done <- 0

touse <- paste(dir.debug, "/", prefix, ".predicted.mlgeno", sep="")
toreturn <- paste(dir.out, "/", prefix, ".mlgeno", sep="")
if (step2.subjects > 0) {

	times <- floor(nrow(P) / step2.subjects)
	start.name <- paste(substr(file.ped, 1, (nchar(file.ped)-4)))
	#tocombine <- "cat "	# command to combine all .mlgeno files together	

	i <- 1
	while((subjects.done+step2.subjects) < nrow(P)) {
		smallP <- P[(subjects.done+1):(subjects.done + step2.subjects),]
		p.name <- paste(start.name, ".group.", step2.subjects, ".", i, ".ped", sep="")
		write.table(smallP, file=paste(dir.debug, p.name, sep="/"), col.names=F, row.names=F, quote=F, sep="\t")
		prefix2 <- paste(prefix, ".", step2.subjects, ".", i, sep="")
		curr.mlgeno <- paste(dir.debug, "/", prefix2, ".mlgeno", sep="")
		# Check if mlgeno already exists: if it exists AND resample=TRUE (i.e. we want to ignore old files and make new ones)
		# then we delete the old mlgeno, to run MaCH again to recompute it. 
		# Next, if mlgeno file doesn't exist, we call MaCH1 on it.
		# But if mlgeno does exist AND resample=F, then we don't want to touch it or re-run MaCH1, so we do nothing.
		if(file.exists(curr.mlgeno) && resample==TRUE) {
			file.remove(curr.mlgeno)
		}
		if(!file.exists(curr.mlgeno)) {
			call.mach.nohap.step2(file.dat, dir.dat, p.name, dir.debug, prefix1=prefix1, prefix=prefix2, dir.out=dir.debug, mach.loc=mach.loc)
		}
		#tocombine <- paste(tocombine, dir.debug, "/", prefix2, ".mlgeno ", sep="")
		subjects.done <- subjects.done + step2.subjects
		i <- i + 1
	}

	# Process remaining patients of the file 
	#if(nrow(P) %% step2.subjects != 0) {
	if(nrow(P) > subjects.done) {		
                smallP <- P[(subjects.done+1):nrow(P),]
                p.name <- paste(start.name, ".group.", step2.subjects, ".", i, ".ped", sep="")
                write.table(smallP, file=paste(dir.debug, p.name, sep="/"), col.names=F, row.names=F, quote=F, sep="\t")
                prefix2 <- paste(prefix, ".", step2.subjects, ".", i, sep="")

		curr.mlgeno <- paste(dir.debug, "/", prefix2, ".mlgeno", sep="")
		if(file.exists(curr.mlgeno) && resample==TRUE) {
			file.remove(curr.mlgeno)
		}
		if(!file.exists(curr.mlgeno)) {
	                call.mach.nohap.step2(file.dat, dir.dat, p.name, dir.debug, prefix1=prefix1, prefix=prefix2, dir.out=dir.debug, mach.loc=mach.loc)
		}
		#tocombine <- paste(tocombine, dir.debug, "/", prefix2, ".mlgeno ", sep="")
	}

	# Concatenate all the relevant .mlgeno files together
	#tocombine <- paste(tocombine, "> ", touse, sep="")
	#system(tocombine)

	# Concatenate all mlgeno files, in correct numerical order, into 'touse' file.

	system(paste("cd ", dir.debug, " ; ls -v *mlgeno > all.txt", sep=""))	
	system(paste("cd ", dir.debug, " ; cat $(ls -v *mlgeno) > ", touse, sep=""))

        # ****************************************************** #
        # Extract only the missing values from the imputed values;
        # since MaCH does imputation of ALL the SNP values, including those that are known.
        # ****************************************************** #
        fill.empties(ped.orig=paste(dir.file, file.ped, sep="/"), ped.new=touse, out.name=toreturn, empty=empty)
	
} else {
	# Process entire file, all output goes into dir.debug, and the necessary file gets copied to dir.out.
	fromfile <- paste(dir.debug, "/", prefix, ".mlgeno", sep="")
	if(file.exists(fromfile) && resample==TRUE) {
		file.remove(fromfile)
	}
	if(!file.exists(fromfile)) {
		call.mach.nohap.step2(file.dat, dir.dat, file.ped, dir.file, prefix1=prefix1, prefix=prefix, dir.out=dir.debug, mach.loc=mach.loc)
	}

	system(paste("cp", fromfile, touse, sep=" "))

	fill.empties(ped.orig=paste(dir.file, file.ped, sep="/"), ped.new=touse, out.name=toreturn, empty=empty)

}

return(toreturn)


}


fill.empties <- function(ped.orig, ped.new, out.name, empty="0/0", num.nonsnp.col.orig=5, num.nonsnp.col.new=2) {
# Replaces all values that are empty in ped.orig with corresponding values from ped.new.
# Note: ped.orig, ped.new, and out.name must be file names with FULL path.
# the ped.new is expected to be in MaCH output format.

orig <- read.table(ped.orig, header=F, sep=" ", stringsAsFactors=FALSE)
if(ncol(orig) <= 1)
	orig <- read.table(ped.orig, header=F, sep="\t", stringsAsFactors=FALSE)

neww <- read.table(ped.new, header=F, sep=" ", stringsAsFactors=FALSE)
if(ncol(neww) <= 1)
	neww <- read.table(ped.new, header=F, sep=" ", stringsAsFactors=FALSE)

# Make the rows consistent. First extract patient names of the new file
# Extract the values of first column and split them by '->'
patient.name <- unlist(strsplit(neww[,1], split="->"))
# If each entry has exactly one '->', then use the first half of each patient name as the id:
if(length(patient.name) == (nrow(neww)*2))
        patient.name <- patient.name[seq(1,nrow(neww)*2,by=2)]
else  # leave the name as is
        patient.name <- neww[,1]

# Build a mask for the original order of patients, at which indexes in new file these patients occur
i <- 1
M <- array(0, nrow(orig))
while(i <= nrow(orig)) {
	id <- which(patient.name == orig[i,1])
	if(!is.na(id))
		M[i] <- id	
	i <- i + 1
}

# Re-arrange rows of new file, s.t. it correponds to the original file
neww <- neww[M, ]

# Remove non-snp-columns from both files:
orig <- orig[,(num.nonsnp.col.orig+1):ncol(orig)]
new.front <- neww[,1:num.nonsnp.col.new]
neww <- neww[,(num.nonsnp.col.new+1):ncol(neww)]

# Replace all empties by new values:
Mempties <- orig==empty
orig[Mempties] <- neww[Mempties]
orig <- cbind(new.front, orig)

write.table(orig, file=out.name, col.names=FALSE, row.names=FALSE, quote=FALSE)

}





