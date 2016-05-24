get.hap.zip <- function(pos.dir, pos.file, pos.header=TRUE, pos.snpcol=5, pos.poscol=2, dat.dir, dat.file, dat.header=FALSE, dat.snpcol=1, dat.ending=".snps", ped.dir, ped.file, num.nonsnp.col=2, chrom.num, pos.list.double, out.dir, out.prefix) {

# Zips up a section of the .hap/ped file (which has been unzipped and
# processed by get.hap.unzip() function). Output name will end with ".hap.gz".
# Helper function that extracts SNPs (out of .snps and .hap files) 
# whose position is within limits defined by pos.list.double.
# Deals with only ONE chromosome.
# Position/nnotation file is assumed to be either space or tab delimited.
#
# Example call:
#
# get.hap.zip(pos.dir="/home/briollaislab/olga/funtry/test_gethap/a1", pos.file="chr20.annotation.txt", dat.dir="/home/briollaislab/olga/funtry/test_gethap/s3", dat.file="chr20.snps", ped.dir="/home/briollaislab/olga/funtry/test_gethap/h2", ped.file="20101123.chr20.txt", chrom.num=20, pos.list.double=c(64139, 65900, 76700, 79299, 62894919, 62896149), out.dir="/home/briollaislab/olga/funtry/test_gethap/o4", out.prefix="reference")
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
# num.nonsnp.col: how many leading non-SNP colums are present in .ped file.
#        .hap.gz = 2;    .ped = 5;    .mlgeno = 2
# 
# chrom.num: which chromosome number  the .dat and .ped files correspond to.
# pos.list.double: a list of starting and ending positions within the chrom
#    c(pos1.start, pos1.end, pos2.start, pos2.end, ... posN.start, posN.end)
# 
# out.dir: directory to which output files should go.
# out.prefix: the beginning of the output filenames.
#
# Writes 2 files: .dat and .ped into out.dir, that only contain the SNPs 
#    whose position is within the pos.list.double boundaries.
#

# TODO: remove this line:
# source("get.position.R")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 0. Check if zip and corresponding dat files exist, do nothing.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ped.name.out <- paste(out.prefix, chrom.num, ".hap.gz", sep="")
ped.full.out <- paste(out.dir, ped.name.out, sep="/")

# NOTE: this is copied from get.position()
dat.name.out <- paste(out.prefix, "_", chrom.num, dat.ending, sep="")
dat.full.out <- paste(out.dir, dat.name.out, sep="/")

if(file.exists(ped.full.out) & file.exists(dat.full.out)) {
	print("Zipped file and its corresponding .dat exist. Skipping.")
	return(list(ped=ped.name.out, dat=dat.name.out))
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1. Call the get.position on the .ped file to get relevant data extracted.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
names <- get.position(pos.dir=pos.dir, pos.file=pos.file, pos.header=pos.header, pos.snpcol=pos.snpcol, pos.poscol=pos.poscol, pos.chrcol=0, dat.dir=dat.dir, dat.file=dat.file, dat.header=dat.header, dat.snpcol=dat.snpcol, dat.ending=dat.ending, ped.dir=ped.dir, ped.file=ped.file, ped.ending=".txt", ped.read=FALSE, num.nonsnp.col=num.nonsnp.col, chrom.num=chrom.num, pos.list.double=pos.list.double, out.dir=out.dir, out.prefix=out.prefix)

print(names)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2. Separate the first columns out.
#    Separate the remaining columns (all the data) out.
#    Remove all spaces from the data file.
#    Combine first columns with spaceless data.
#    Zip up the combined file.
#    Remove all intermediate steps.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

print("getting first columns...")

ped.name.cols <- paste(out.dir, "/cols.chr", chrom.num, ".txt", sep="")
cut.call <- paste("cut -d \" \" -f 1-", num.nonsnp.col, " ", out.dir, "/", names$ped, " > ", ped.name.cols,  sep="")
system(cut.call)
print("getting data cols...")

ped.name.data <- paste(out.dir, "/data.chr", chrom.num, ".txt", sep="")
cut.data <- paste("cut -d \" \" -f 1-", num.nonsnp.col, " --complement ", out.dir, "/", names$ped, " > ", ped.name.data,  sep="")
system(cut.data)
print("space removing from data...")

ped.name.delim <- paste(out.dir, "/undelim.chr", chrom.num, ".txt", sep="")
cut.delim <- paste("sed -e 's/ //g' ", ped.name.data, " > ", ped.name.delim,  sep="")
system(cut.delim)
print("combining first cols with delimeted data...")

ped.name.combo <- paste(out.dir, "/combo.chr", chrom.num, ".txt", sep="")
cut.combo <- paste("paste -d ' ' ", ped.name.cols, " ", ped.name.delim, " > ", ped.name.combo,  sep="")
system(cut.combo)
print("combined. Cleaning up intermediate files...")

#ped.name.out <- paste(out.prefix, chrom.num, ".hap.gz", sep="")
#ped.full.out <- paste(out.dir, ped.name.out, sep="/")
uncomp.call <- paste("gzip -c ", ped.name.combo, " > ", ped.full.out,  sep="")
system(uncomp.call)


system(paste("rm", ped.name.cols, names$ped, ped.name.data, ped.name.delim, ped.name.combo,  sep=" "))

print("done")

return(list(ped=ped.name.out, dat=names$dat))

}

