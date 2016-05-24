tune1.subsets <- function(dir.dat, dir.ped, dir.ann, dir.pos.snp, dir.pos.ann, dir.pos.hap, dir.out, prefix.dat, prefix.ped, prefix.ann, prefix.pos.snp, prefix.pos.ann, prefix.pos.hap, key.dat="", key.ann="", key.pos.ann="", key.pos.hap="", ending.dat=".dat", ending.ped=".ped", ending.ann=".map", ending.pos.snp=".snps", ending.pos.ann="annotation.txt", ending.pos.hap=".hap.gz", pos.list.triple, ped.nonsnp=5, ann.header=FALSE, ann.snpcol=2, ann.poscol=4, ann.chrcol=0, pos.ann.header=TRUE, pos.ann.snpcol=5, pos.ann.poscol=2, pos.hap.nonsnp=2, out.name.subdir="seg1", out.prefix="subdata", rsq.thresh=0.5, num.iters=2, hapmapformat=FALSE, mach.loc="/software/mach1") {

# For chromosomes and small regions specified in pos.list.triple, runs MaCH1
# to get more detailed sampling of SNPs in the region, and prepares this
# subset of data to be processed by MOSS algorithm.

#pos.dir, pos.file, pos.header=FALSE, pos.snpcol=2, pos.poscol=4, pos.chrcol=1, dat.dir, dat.file, dat.header=FALSE, dat.snpcol=2, dat.ending=".dat", ped.dir=dat.dir, ped.file, ped.ending=".ped", ped.sep=" ", ped.read=FALSE, num.nonsnp.col=5, chrom.num, pos.list.double, out.dir, out.prefix) {

# Helper function that extracts SNPs (out of .dat and .ped files) 
# whose position is within limits defined by pos.list.double.
# Deals with only ONE chromosome. Either CASE or CONTROL at a time.
# Requires information about location, name and format of file that contains
# position info: pos.<info>; as well as information about .dat and .ped files.
# Files are assumed to be either space or tab delimited.
# Position file could be either .map, .legend.txt, .annotation.txt, or
# any other format that contains list of all SNPs of .dat file, and their
# positions within the chromosome. Please make sure it's unzipped and readable.
# Note: ped.read is a flag to indicate if you wish R to read in the .ped file
# (too time consuming for big .ped files), alternatively column extraction
# will be performed using Linux 'cut' and 'paste' functions.
# Example call:
#
# tune1.subsets(dir.dat="/home/briollaislab/olga/funtry/test_tune/ins/d03_removed", dir.ped="/home/briollaislab/olga/funtry/test_tune/ins/d03_removed", dir.ann="/home/briollaislab/olga/funtry/test_tune/ins/d01_plink", dir.pos.snp="/home/briollaislab/olga/funtry/test_tune/ins/d04_ref/snps", dir.pos.ann="/home/briollaislab/olga/funtry/test_tune/ins/d04_ref/anno", dir.pos.hap="/home/briollaislab/olga/funtry/test_tune/ins/d04_ref/hap", dir.out="/home/briollaislab/olga/funtry/test_tune", prefix.dat="genos_chr", prefix.ped="genotypes_", prefix.ann="genos_chr", prefix.pos.snp="snp", prefix.pos.ann="chr", prefix.pos.hap="20101123.", pos.list.triple=c(20, 13034962, 13521577, 20, 58174374, 60765836), out.name.subdir="seg1", out.prefix="subdata")
#
# Parameters:
#
# pos.dir: directory where position info can be found.
# pos.file: name of file containing position info.
# pos.header: whether or not the pos.file has a header.
#        .map = FALSE; .legend.txt = TRUE .annotation = TRUE
# pos.snpcol: which column of the file corresponds to SNP ids.
#        .map = 2;    .legend.txt = 1     .annotation = 5
# pos.poscol: which column of the file corresponds to position info.
#        .map = 4;    .legend.txt = 2     .annotation = 2
# pos.chrcol: which column of the file corresponds to chromosome number.
#        .map = 1;    .lengend.txt = 0 (doesn't exist default).
# dat.dir: directory where .dat file of the relevant dataset is located.
# dat.file: name of the .dat file.
# dat.header: whether or not the .dat file has a header.
#        .dat = FALSE; .mlinfo = TRUE
# dat.snpcol: which column in the .dat file contains SNP values.
#        .dat = 2;    .mlinfo = 1
# dat.ending: the end of .dat filename for output purposes (include the dot).
#        Ex: .dat, .mlinfo, .legend, .snps, etc.
#
# ped.dir: directory where relevant .ped file is located.
# ped.file: name of the .ped file.
# ped.ending: the end of .ped filename for output (include the dot).
# ped.read: TRUE if you wish R to read in the .ped file to process it.
#        FALSE if the file is too big to be read, and Linux commands will
#        be used to process the file's columns.
# num.nonsnp.col: how many leading non-SNP colums are present in .ped file.
#        .ped = 5;    .mlgeno = 2
# 
# chrom.num: which chromosome number  the .dat and .ped files correspond to.
# pos.list.double: a list of starting and ending positions within the chrom
#    c(pos1.start, pos1.end, pos2.start, pos2.end, ... posN.start, posN.end)
# 
# out.dir: directory to which output files should go.
# out.prefix: the beginning of the output filenames.
#

# TODO: remove these lines:
# source("subdir.create.R")
# source("get.hap.unzip.R")
# source("get.hap.zip.R")
#source("get.position.R")
#source("get.file.name.R")
#source("call.mach.R")
#source("clean.rsq.combine.R")

#source("pre5.genos2numeric.batch.R")
#source("pre5.genos2numeric.R")
#source("get.chrom.num.R")
#source("machout2id.R")
#source("get.ext.R")

# Check params:
if(missing(dir.dat)) stop("Name of input directory for .dat file must be provided.")
if(missing(dir.ped)) stop("Name of input directory for .ped must be provided.")
if(missing(dir.ann)) stop("Name of input directory with annotation file for the dataset must be provided.")
if(missing(dir.pos.snp)) stop("Name of input directory for hapmap SNP information must be provided.")
if(missing(dir.pos.ann)) stop("Name of input directory for hapmap position information must be provided.")
if(missing(dir.pos.hap)) stop("Name of input directory for hapmap data must be provided.")
if(missing(dir.out)) stop("Name of output directory must be provided.")

if(missing(prefix.dat)) stop("Prefix of the .dat file name must be provided.")
if(missing(prefix.ped)) stop("Prefix of the .ped file name must be provided.")
if(missing(prefix.ann)) stop("Prefix of the position information file name must be provided.")
if(missing(prefix.pos.snp)) stop("Prefix of the hapmap SNP file name must be provided.")
if(missing(prefix.pos.ann)) stop("Prefix of the hapmap annotation file name must be provided.")
if(missing(prefix.pos.hap)) stop("Prefix of the hapmap data file name must be provided.")

# Check validity of triplets list:
if(is.null(pos.list.triple)) stop("Paramter pos.list.triple is missing")
if((length(pos.list.triple)%%3) != 0) stop("Length of list with chromosome numbers and position intervals must be divisible by 3.")
i <- 1
while(i<= length(pos.list.triple)) {
	if(pos.list.triple[i] < 0) 
		stop("All values in the triplet list of chromosomes and position boundaries must be positive")
	if(i%%3 == 2 & pos.list.triple[i] > pos.list.triple[i+1])
		stop(paste("The triplet list seems to have start and end positions reversed at index ", i, " and ", i+1, ".", sep=""))
	i <- i + 1
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 0. Before creating subdirectory, check if 'position_info.txt'
#    file already exists. (If dne - go ahead and create subdirs.)
#    If it does, check that position triplets match. 
#      - match: ok, re-create subdirs anway.
#      - no match: print error message and stop.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
name.triplet <- paste(dir.out, out.name.subdir, "position_info.txt", sep="/")
if(file.exists(name.triplet)) {
	triplet.old <- read.table(name.triplet, header=F, sep=" ", stringsAsFactors=FALSE)
	if(!all(triplet.old == pos.list.triple))
		stop(paste("Subdirectory '", out.name.subdir, "' already exists for a different set of chromosomes/positions. Please change the name to new/different subdirectory.", sep=""))

} 

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1. Create subdirectory structure
#    dout$s1 = s01_trimmed
#    dout$s2 = s02_machout
#    dout$s3 = s03_combined
#    dout$s4 = s04_numeric
#    Write 'position.info.txt' file into s0.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dout <- subdir.create(dir.out=dir.out, out.name=out.name.subdir)
if(!file.exists(name.triplet))
	write.table(pos.list.triple, file=name.triplet, col.names=F, row.names=F, quote=F, sep=" ")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2. Shred the list of triplets into lists of doubles for each chromosome. 
#    By first converting list into matrix (first ROW lists chrom numbers)
#    Then obtain all chromosome numbers involved, and iterate over them.
#    Unzip the hapmap into its original folder.
#    Zip up the portion of hapmap file into s1_trimmed folder.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pos.matrix <- matrix(pos.list.triple, nrow=3)
chroms <- unique(pos.matrix[1,])

chrom.id <- 1
while(chrom.id <= length(chroms)) {
	chrom <- chroms[chrom.id]
	pos.list.double <- as.vector(pos.matrix[2:3, which(pos.matrix[1,]==chrom)])
	# Get all the hapmap files names:
	name.pos.hap <- get.file.name(dir=dir.pos.hap, prefix=paste(prefix.pos.hap, chrom, sep=""), key=key.pos.hap, ending=ending.pos.hap, onename=TRUE)
	name.pos.ann <- get.file.name(dir=dir.pos.ann, prefix=paste(prefix.pos.ann, chrom, sep=""), key=key.pos.ann, ending=ending.pos.ann, onename=TRUE)
	name.pos.snp <- get.file.name(dir=dir.pos.snp, prefix=paste(prefix.pos.snp, chrom, sep=""), ending=ending.pos.snp, onename=TRUE)

	print(paste("Unzipping: ", name.pos.hap, sep=""))
	unzipped.hap <- get.hap.unzip(ped.dir=dir.pos.hap, ped.file=name.pos.hap, num.nonsnp.col=pos.hap.nonsnp, chrom.num=chrom, out.prefix=prefix.pos.hap)

	print("Extracting and zipping... ")
	
	zipped.hap <- get.hap.zip(pos.dir=dir.pos.ann, pos.file=name.pos.ann, pos.header=pos.ann.header, pos.snpcol=pos.ann.snpcol, pos.poscol=pos.ann.poscol, dat.dir=dir.pos.snp, dat.file=name.pos.snp, dat.header=FALSE, dat.ending=ending.pos.snp, ped.dir=dir.pos.hap, ped.file=unzipped.hap, num.nonsnp.col=pos.hap.nonsnp, chrom.num=chrom, pos.list.double=pos.list.double, out.dir=dout$s1, out.prefix="hap.extract")

print("zipped:")
print(zipped.hap$ped)
print(zipped.hap$dat)

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# 3. Deal with CASE and CONTROL
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	case.type <- "CASE"
	case.id <- 1
	name.mlinfoCASE <- ""
	name.mlgenoCASE <- ""
	name.mlinfoCONTROL <- ""
	name.mlgenoCONTROL <- ""

	while(case.id <= 2) {
		if(case.id == 2) 
			case.type <- "CONTROL"

        	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        	# 4. Obtain all names of ped, dir and annotation files. 
		#    Extract desired positions from the ped file
        	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		print(case.type)
		name.ann <- get.file.name(dir=dir.ann, prefix=paste(prefix.ann, chrom, sep=""), key=key.ann, ending=ending.ann, onename=TRUE)
		if(ann.chrcol > 0) # There is one file for many chroms, so omit chrom number from file's name.
			name.ann <- get.file.name(dir=dir.ann, prefix=prefix.ann, key=key.ann, ending=ending.ann, onename=TRUE)

		name.dat <- get.file.name(dir=dir.dat, prefix=paste(prefix.dat, chrom, sep=""), ending=ending.dat, onename=TRUE)
		name.ped <- get.file.name(dir=dir.ped, prefix=paste(prefix.ped, chrom, sep=""), key=case.type, ending=ending.ped, onename=TRUE)

		name.extract <- get.position(pos.dir=dir.ann, pos.file=name.ann, pos.header=ann.header, pos.snpcol=ann.snpcol, pos.poscol=ann.poscol, pos.chrcol=ann.chrcol, dat.dir=dir.dat, dat.file=name.dat, dat.ending=ending.dat, ped.dir=dir.ped, ped.file=name.ped, ped.ending=ending.ped, ped.sep="\t", ped.read=FALSE, num.nonsnp.col=ped.nonsnp, chrom.num=chrom, pos.list.double=pos.list.double, out.dir=dout$s1, out.prefix=paste("data.extract", case.type, sep="")) 
print("Extracted .dat and .ped from data. Calling MaCH:")

                # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                # 5. Call MACH with hapmap into s2 folder.
                # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		mach.prefix <- paste("mach", chrom, case.type, sep="")
		mach.mlinfo <- paste(mach.prefix, ".mlinfo", sep="")
		mach.mlgeno <- paste(mach.prefix, ".mlgeno", sep="")

		if(case.id == 1) {
			name.mlinfoCASE <- mach.mlinfo
        		name.mlgenoCASE <- mach.mlgeno

		} else {
        		name.mlinfoCONTROL <- mach.mlinfo
        		name.mlgenoCONTROL <- mach.mlgeno
		}

		if(file.exists(paste(dout$s2, mach.mlinfo, sep="/")) & file.exists(paste(dout$s2, mach.mlgeno, sep="/"))) {
			print(paste("MaCH computed files for chromosome ", chrom, " ", case.type, "  aleady exist. Skipping", sep=""))
		} else {
			call.mach(file.dat=name.extract$dat, file.ped=name.extract$ped, dir.file=dout$s1, ref.phase=zipped.hap$ped, ref.legend=zipped.hap$dat, dir.ref=dout$s1, prefix=mach.prefix, dir.out=dout$s2, num.iters=num.iters, num.subjects=0, resample=FALSE, mach.loc=mach.loc) 
		}

                case.id <- case.id + 1
        } # end CASE/CONTROL iteration

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 6. Prune out bad estimates with RSQ<0.5 from 
	#    both CASE and CONTROL files, and combine them.
	#    Save results to s3 folder.
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	print("Pruning bad MaCH estimates for CASE and CONTROL.")

	name.pref <- "data"

	name.clean <- clean.rsq.combine(case.file=name.mlgenoCASE, control.file=name.mlgenoCONTROL, dir.file=dout$s2, case.info=name.mlinfoCASE, control.info=name.mlinfoCONTROL, out.name.start=paste(name.pref, chrom, sep=""), rsq.thresh=rsq.thresh, dir.out=dout$s3, separ=" ") 

	print(name.clean)

	chrom.id <- chrom.id + 1
} # end cromosome # iteration

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 7. Run pre5.genos2numeric.batch() to convert genos to 1,2,3.
#    into s4 folder.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

print("NUMERIZING...")	
pre5.genos2numeric.batch(dir.ped=dout$s3, dir.out=dout$s4, prefix.ped=name.pref, prefix.dat=name.pref, remove.bad.genos=TRUE, save.ids.name="patients.fam")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 8. Run pre6.discretize.batch() to convert genos to 0,1.
#    into s5 folder.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#pre6.discretize.batch(dir.file=dout$s4, dir.out=dout$s5, prefix.train=name.pref, train.output.append="binary", test.output.append="same")
 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 9. Run pre7.merge.genos() to merge all files,
#    into s0
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
name.final <- pre6.merge.genos(dir.file=dout$s4, dir.out=dout$s0, file.out=paste(out.prefix, ".txt", sep=""), dat.out=paste(out.prefix, ".dat", sep=""), prefix.file=name.pref, prefix.dat=name.pref, key.file="", key.dat="", weak.check=FALSE, plan=FALSE)



print("done")
return(name.final)

} # function end

