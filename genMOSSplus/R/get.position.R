get.position <- function(pos.dir, pos.file, pos.header=FALSE, pos.snpcol=2, pos.poscol=4, pos.chrcol=1, dat.dir, dat.file, dat.header=FALSE, dat.snpcol=2, dat.ending=".dat", ped.dir=dat.dir, ped.file, ped.ending=".ped", ped.sep=" ", ped.read=FALSE, num.nonsnp.col=5, chrom.num, pos.list.double, out.dir, out.prefix) {
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
# b <- get.position(pos.dir="/home/briollaislab/olga/funtry/test_getpos/d1", pos.file="genos_chr21.map", dat.dir="/home/briollaislab/olga/funtry/test_getpos/d3", dat.file="genos_chr21.removed.dat", ped.file="genotypes_21CASE.removed.ped", chrom.num=21, pos.list.double=c(14526772, 15465149, 20118479, 20979034), out.dir="/home/briollaislab/olga/funtry/test_getpos/s1", out.prefix="genosCASE")
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
# Writes 2 files: .dat and .ped into out.dir, that only contain the SNPs 
#    whose position is within the pos.list.double boundaries.
#


# Check params:
if((length(pos.list.double)%%2) != 0) stop("Length of list with position intervals must be even")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 0. check if files already exist, then do nothing.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
out.ped.name <- paste(out.prefix, "_", chrom.num, ped.ending, sep="")
out.ped <- paste(out.dir, out.ped.name, sep="/")

out.dat.name <- paste(out.prefix, "_", chrom.num, dat.ending, sep="")
out.dat <- paste(out.dir, out.dat.name, sep="/")

if(file.exists(out.ped) & file.exists(out.dat)) {
	print(paste("Extracted data files ", out.ped.name, " and ", out.dat.name, " already exist. Skipping", sep="")) 
	return(list(ped=out.ped.name, dat=out.dat.name))
}


POS <- read.table(paste(pos.dir, pos.file, sep="/"), header=pos.header, sep="\t", stringsAsFactors=FALSE)
if(ncol(POS) <= 1)
        POS <- read.table(paste(pos.dir, pos.file, sep="/"), header=pos.header, sep=" ", stringsAsFactors=FALSE)

dat.sep <- "\t"
D <- read.table(paste(dat.dir, dat.file, sep="/"), header=dat.header, sep="\t", stringsAsFactors=FALSE)
if(ncol(D) <= 1) {
	D <- read.table(paste(dat.dir, dat.file, sep="/"), header=dat.header, sep=" ", stringsAsFactors=FALSE)
	dat.sep <- " "
}

PED <- NULL
if(ped.read) {
	PED <- read.table(paste(ped.dir, ped.file, sep="/"), header=F, sep=ped.sep, stringsAsFactors=FALSE)
	if(ncol(PED) <= 1) {
		stop(paste("The ped file, '", ped.dir, "/", ped.file, "' seems not to be separated by '", ped.sep, "'.", sep=""))
	}
}

print("Read all files")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1. From pos file, obtain all SNP names of those SNPs whose position
#   satisfies the criteria of given list of position intervals.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

i = 1
M <- 0
while (i < length(pos.list.double)) {
	if (pos.chrcol > 0) # more than 1 chrom is in POS file
		m <- POS[which(POS[,pos.chrcol]==chrom.num & POS[,pos.poscol]>=pos.list.double[i] & POS[,pos.poscol]<=pos.list.double[i+1]), pos.snpcol]
	else
		m <- POS[which(POS[,pos.poscol]>=pos.list.double[i] & POS[,pos.poscol]<=pos.list.double[i+1]), pos.snpcol]

	if(i == 1) 
		M <- m
	else
		M <- append(M,m)
	i <- i + 2
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2. From dat file, obtain all SNP names that are  present in the above
#   list (intersect). 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if(length(M) < 1) 
	stop(paste("No SNPs in the specified range were found in position file: ", pos.dir, "/", pos.file, sep=""))

common.snps <- intersect(M, D[,dat.snpcol]) 

if(length(common.snps) < 1)
	stop("No SNPs in the specified range were found in common with .dat file")

print(paste("Obtained ", length(common.snps), " SNPs in .ann and .dat. Determining their indexes in .dat. This will take time...", sep=""))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3. Find row locations in .dat file where desired SNPs are located.
#   Then update .dat and .ped files accordingly.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

i <- 1
I <- vector(length=length(common.snps))
while(i <= length(common.snps)) {
       I[i] <- which(D[,dat.snpcol]==common.snps[i])
       i <- i + 1
}

#i <- 2
#I <- which(D[,dat.snpcol]==common.snps[1]) 
#while(i<= length(common.snps)) {
#	I <- append(I, which(D[,dat.snpcol]==common.snps[i])) 
#	i <- i + 1
#}

I <- sort(I)
print("Obtained. Saving .dat.")

newD <- D[I,]
if(dat.header) {
	dimnames(newD) <- list(NULL, unlist(dimnames(D)[2]))
}

# Save the results
#out.dat.name <- paste(out.prefix, "_", chrom.num, dat.ending, sep="")
#out.dat <- paste(out.dir, out.dat.name, sep="/")
write.table(newD, file=out.dat, col.names=dat.header, row.names=F, quote=F, sep=dat.sep)

# Deal with .ped file
#out.ped.name <- paste(out.prefix, "_", chrom.num, ped.ending, sep="")
#out.ped <- paste(out.dir, out.ped.name, sep="/")

if(ped.read) {
print("Using .ped's read file")
        newP <- cbind(PED[,1:num.nonsnp.col], PED[, I+num.nonsnp.col])
	write.table(newP, file=out.ped, col.names=F, row.names=F, quote=F, sep=ped.sep)
	# Remove variables:
	remove(PED)

} else {
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 4. To use LINUX commands on .ped file, prepare the I list:
#    Group the list into regions with start and end, ex:
#    I = [4, 5, 6, 8, 10, 11, 12]
#    To group like: 4-6, 8-8, 10-12:
#    G.start = [4, 8, 10]
#    G.end   = [6, 8, 12]   
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	G.start <- I[1]
	G.end <- NULL
	last.val <- I[1]
	i <- 2		# counter along I
	while(i <= length(I)) {
		if(I[i] != (last.val+1)) {
			G.end <- append(G.end, last.val)
			G.start <- append(G.start, I[i])
		}
		last.val <- I[i]
		i <- i + 1
	}
	G.end <- append(G.end, last.val) # last value
	print(paste("Using Linux to deal with .ped. There are ", length(G.start), " intervals of column ids that are grouped as follows: Group Starts:", sep=""))
	print(G.start)
	print("Group Ends:")
	print(G.end)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 5. For each of the groups of columns, 
#    write and execute Linux command to extract these cols from .ped file
#    Also extract first non-snp columns.
#    Then combine all these intermediate files together.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	print("Extracting groups from .ped: see files")

	g <- 1
	tmp.file <- paste(out.dir, "/tmp", chrom.num, "_0.txt", sep="")
	system(paste("cut -d \"", ped.sep ,"\" -f 1-", num.nonsnp.col, " ", ped.dir, "/", ped.file, " > ", tmp.file, sep=""))

	combo <- paste("paste -d \"", ped.sep, "\" ", tmp.file, sep="")
	todelete <- paste("rm ", tmp.file, sep="")

	while(g <= length(G.start)) {
		tmp.file <- paste(out.dir, "/tmp", chrom.num, "_", g, ".txt", sep="")
		# If only one column need extraction:
		if((G.end[g] - G.start[g]) == 0) 
			system(paste("cut -d \"", ped.sep ,"\" -f ", (G.start[g]+num.nonsnp.col), " ", ped.dir, "/", ped.file, " > ", tmp.file, sep=""))
		else
			 system(paste("cut -d \"", ped.sep ,"\" -f ", (G.start[g]+num.nonsnp.col), "-", (G.end[g]+num.nonsnp.col), " ", ped.dir, "/", ped.file, " > ", tmp.file, sep=""))

		combo <- paste(combo, tmp.file, sep=" ")
		todelete <- paste(todelete, tmp.file, sep=" ")
		g <- g + 1
	}
print("pasting peds together")
	combo <- paste(combo, " > ", out.ped, sep="") 
	system(combo)
	system(todelete)
}



return(list(ped=out.ped.name, dat=out.dat.name))
}

