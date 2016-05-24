get.hap.unzip <- function(ped.dir, ped.file, num.nonsnp.col=2, chrom.num, out.dir=ped.dir, out.prefix) {

# Helper function that unzips and space delimits the data of the
# hapmap file.
# Deals with only ONE chromosome.
#
# Example call:
#
# get.hap.unzip(ped.dir="/home/briollaislab/olga/funtry/test_gethap/h2", ped.file="20101123.chr20.hap.gz", chrom.num=20, out.prefix="20101123.chr")
#
# Parameters:
#
# ped.dir: directory where relevant .ped file is located.
# ped.file: name of the hapmap file, ending with .hap.gz.
# num.nonsnp.col: how many leading non-SNP colums are present in .ped file.
#        .hap.gz = 2;    .ped = 5;    .mlgeno = 2
# 
# chrom.num: which chromosome number  the .dat and .ped files correspond to.
# 
# out.dir: directory to which output files should go.
# out.prefix: the beginning of the output filename.
#
# Writes the unzipped delimeted (not useable by MaCH) file
# into out.dir folder and returns its name (without path).

# Name of final output file.
ped.name.combo <- paste(out.prefix, chrom.num, ".txt", sep="")

# Check if it already exists:
if(file.exists(paste(out.dir, "/", ped.name.combo, sep=""))) {
	print("Unzipped prepared file already exists.")
	return(ped.name.combo)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1. Unzip the ped file.
#    Separate the first columns out, and if they were tab separated, make it space-separated.
#    Separate the last column (all the data) out.
#    Split all letters of the data file with spaces.
#    Combine first columns with delimeted data.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ped.name <- paste(out.dir, "/uncompressed.chr", chrom.num, ".hap", sep="")
uncomp.call <- paste("gzip -cd ", ped.dir, "/", ped.file, " > ", ped.name,  sep="")
system(uncomp.call)
print("Uncompressed. getting first cols...")

ped.name.cols <- paste(out.dir, "/cols.chr", chrom.num, ".txt", sep="")
cut.call <- paste("cut -d \" \" -f 1-", num.nonsnp.col, " ", ped.name, " > ", ped.name.cols,  sep="")
system(cut.call)
ped.name.cols.space <- paste(out.dir, "/colsspace.chr", chrom.num, ".txt", sep="")
cut.space <- paste("sed -e 's/\t/ /g' ", ped.name.cols, " > ", ped.name.cols.space,  sep="")
system(cut.space)
print("getting data...")

ped.name.data <- paste(out.dir, "/data.chr", chrom.num, ".txt", sep="")
cut.data <- paste("cut -d \" \" -f 1-", num.nonsnp.col, " --complement ", ped.name, " > ", ped.name.data,  sep="")
system(cut.data)
print("space delimiting the data...")

ped.name.delim <- paste(out.dir, "/delim.chr", chrom.num, ".txt", sep="")
cut.delim <- paste("sed -e 's/[A-Z]/& /g' ", ped.name.data, " > ", ped.name.delim,  sep="")
system(cut.delim)
print("combining first cols with delimeted data...")

cut.combo <- paste("paste -d ' ' ", ped.name.cols.space, " ", ped.name.delim, " > ", out.dir, "/", ped.name.combo,  sep="")
system(cut.combo)
print("combined. Cleaning up intermediate files...")

system(paste("rm", ped.name, ped.name.cols, ped.name.cols.space, ped.name.data, ped.name.delim, sep=" "))

print("done")

return(ped.name.combo)
}

