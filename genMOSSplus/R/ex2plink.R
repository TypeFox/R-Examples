ex2plink <- function(dir.file, dir.out, file.name="genotypes_10_90.txt", annotation.name="Identifiers_comma.csv", out.prefix.ped="genotypes_", out.prefix.dat="genos_chr") {

# Converts the example dataset to PLINK format.
# This file is for demo purposes only. 
# You will need to modify it to go from your file format to PLINK.
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# The dataset stored in file.name is of the following format:
#
# Status	1	0	1 ...
# 1719214 AG      GG      AG    ...
# 2320341 TT      TT      TT   ...
# ...
# 
#   - Tab delimeted
#   - No header
#   - First row is the disease status
#   - First column is the list of Markers
#   - rows: geno information, no separator between alleles. 
#   - columns: individuals/patients/samples
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# The annotation.name file is of the following format:
# 
# Marker,RefSNP_ID,CHROMOSOME,CHROMOSOME_LOCATION     ...
# 1546,,1,2103664 ...
# 1996,rs1338382,1,2708522 ....
# 2841,"rs2887274,rs4369170",1,3504300 ...
# ...
# 
#   - Comma delimited (due to missing values)
#   - Has a header
#   - Col 1: Markers, most appear in Col 1 of file.name
#   - Col 2: RefSNP_ID:
#      * empty if missing
#      * one SNP ID
#      * two or 3 corresponding SNP IDs, comma separated, no spaces. 
#   - Col 3: chromosome number
#   - Col 4: physical locations
#   - First 4 columns are important, other columns will be ignored so do not matter.
#   - rows: correspond to all available SNP IDs
#
#
# ************************************************************
# ************************************************************
# ************************************************************ 
#
# This program will generate output files: 
#   for each chromosome 2 files: .ped and .dat 
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# The .ped output file will be of the format:
# 
# p1    p1    0       0       1       2       CC     NN     TC ...
# p2	p2    0       0       1       2       TT     AC     GG ...
# ...
#
#   - Tab separated
#   - No header
#   - 6 non-SNP leading columns
#   - Col 1 and Col 2: patient ID: some unique ID 
#   - Col 3 and Col 4: parents: mother/father: set to 0
#   - Col 5: gender, default to 1 (male)
#   - Col 6: disease status: 1 CONTROL and 2 CASE
#   - Col 7+: geno information, no separator between alleles.
#
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# The .map output file will be of the format:
#
# 19 rs32453434 0 5465475
# 19 rs6547434 0 23534543
# ...
#
#   - Space separated
#   - No header
#   - 4 columns:
#   - Col 1: Chromosome number (Col 3 from annotation file)
#   - Col 2: SNP ID or Marker if SNP is not known (Col 2 from annotation file, or Col1 if Col2="")
#   - Col 3: always 0
#   - Col 4: physical locations (Col 4 from annotation file)
#   - Number of rows is the number of SNPs used in the given chromosome. (= number of SNP columns of .ped)
#


# Read in the data file and annotation file
data.file <- read.table(paste(dir.file, file.name, sep="/"), sep="\t", header=FALSE, stringsAsFactors=FALSE)
ann.file<-read.table(paste(dir.file, annotation.name, sep="/"), sep=",", header=TRUE, stringsAsFactors=FALSE)

# Transpose the data.file, such that columns are SNPs, 
# and 1st column becomes disease status.
# and 1st row lists all the SNP Markers.
data.file <- t(data.file)

# Save the disease status and SNP Marker names separately
disease.status <- data.file[2:nrow(data.file), 1]
marker.names <- data.file[1,2:ncol(data.file)]

# Now set data.file to be pure data
data.file <- data.file[2:nrow(data.file), 2:ncol(data.file)]
ncols <- ncol(data.file)

# ************************************************************************ #
# Iterate over all the Markers of data file.
# For each marker, find its corresponding row in annotation file
# If a marker does not exist in annotation file, print error 
# (since we don't know chromosome number for it)

i <- 1 

# Array that keeps at which index in annotation file Marker was found.
ids.ann <- matrix(0, ncols, 1) 

# Since finding the indexes takes a long time, we can save them and
# use them instead of generating them every time.
index.name <- paste(dir.file, "indices.ann.txt", sep="/")

if(file.exists(index.name)){
	ids.ann <- read.table(index.name, header=FALSE, sep=" ", stringsAsFactors=FALSE)
	ids.ann <- unlist(ids.ann)
} else {
	# The following code shows how to generate that file with indices.
	print(paste("Processing ", ncols, " SNPs. This is slow...", sep=""))

	while(i <= ncols) {
		if(i%%1000==0)
			print(paste("i = ", i, sep=""))

		# Find index of current Marker in annotation file's 1st column
		id <- match(marker.names[i], ann.file[,1])
	
		# If the search failed, then we do not know anything about this marker
		if(is.na(id)) {
			print(paste("Warning: Marker ", data.file[1,i], " was not found in annotation file", sep=""))
		} else {
			ids.ann[i] <- id	# store the ID
		}

		i <- i + 1
	}
	# save the indexes
	write.table(ids.ann, file=index.name, sep=" ", col.names=FALSE, row.names=FALSE, quote=FALSE)
}

# ************************************************************************ #
# Now ids.ann contain annotation file IDs for each marker in data.file. 
# Get all the SNPs that are used and throw out the rest. 
# Set ann.file to contain all info from annotation file only for used SNPs,
#   ordered in the same way as SNPs are ordered in the data file.
# Get all chromosome numbers that are used (all.chroms) and sort them.
ann.file <- ann.file[ids.ann,1:4]

all.chroms <- unique(ann.file[, 3])

# Convert all chromosomes to numeric values (luckily for this dataset, all chroms are numeric,
# but if they were not, we would need to encode non-numeric values as numeric: 
# for example "X" as 23, "Y" as 24, etc).
all.chroms.sort <- sort(as.numeric(all.chroms))

# ************************************************************************ #
# For each chromosome, create 2 files: .ped and .map of the format described above.
# 

i <- 1
while(i <= length(all.chroms.sort)) {

	curr.chrom <- all.chroms.sort[i]
	
	# boolean has TRUE for all rows that correspond to current chromosome
	bool.chrom <- (ann.file[,3]==curr.chrom)

	# Data for this chromosome, its annotation, and its markers
	chrom.data <- data.file[,bool.chrom]
	chrom.ann <- ann.file[bool.chrom,]
	chrom.markers <- marker.names[bool.chrom]

	# Data should consist of Alleles separated by a slash,
	# whereas this dataset currently has no separator between Alleles
	chrom.data <- matrix(paste(substr(chrom.data,1,1), substr(chrom.data,2,2), sep="/"), nrow=nrow(chrom.data), byrow=F)

	# Prepare the .ped file format:
	# Col 1 and 2: invent some unique names for data rows
	# Col 3 and 4: remain 0s
	# Col 5: set to 1, as if all are males.
	# Col 6: disease status, originally we have 0-CONTROL and 1-CASE,
	#       now we re-encode it as 1-CONTROL and 2-CASE
	ped.file <- matrix(0, nrow(chrom.data), 6)
	ped.file[,1] <- paste("p", (1:nrow(chrom.data)), sep="")
	ped.file[,2] <- ped.file[,1]
	ped.file[,5] <- rep(1, nrow(chrom.data))
	ped.file[,6] <- as.numeric(disease.status) + 1

	ped.file <- cbind(ped.file, chrom.data) # combine first 6 cols with full data

	# Save .ped file:
	ped.name <- paste(dir.out, "/", out.prefix.ped, curr.chrom, ".ped", sep="")
	write.table(ped.file, file=ped.name, col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")

	# Save 2 .ped files: CASE and CONTROL separately
	#ped.name.case <- paste(dir.out, "/", out.prefix.ped, curr.chrom, "CASE.ped", sep="")
	#ped.name.control <- paste(dir.out, "/", out.prefix.ped, curr.chrom, "CONTROL.ped", sep="")
	#write.table(ped.file[(ped.file[,6]==2),], file=ped.name.case, col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")
	#write.table(ped.file[(ped.file[,6]==1),], file=ped.name.control, col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~	
	# Prepare the .map file format:
	# Col1: chrom number
	# Col2: SNP ID, or Marker if no SNP ID
	# Col3: 0
	# Col4: physical location, Col4 from annotation
	dat.file <- matrix(0, ncol(chrom.data), 4)
	dat.file[,1] <- rep(curr.chrom, ncol(chrom.data))

	# Iterate over all SNP IDs in annotation, extract the first SNP ID from
	# each row (since for any one entry there may be multiple SNP IDs, comma separated)
	# If there is no SNP ID for given entry, then use the Marker name
	id.splits <- strsplit(chrom.ann[,2], ",")

	j <- 1
	while(j <= ncol(chrom.data)) {
		dat.file[j,2] <- unlist(id.splits[j])[1]
		if(is.na(dat.file[j,2]))
			dat.file[j,2] <- chrom.markers[j]

		j <- j+1
	}
	dat.file[,4] <- chrom.ann[,4]

	# Save the .map file
	dat.name <- paste(dir.out, "/", out.prefix.dat, curr.chrom, ".map", sep="")
	write.table(dat.file, file=dat.name, col.names=FALSE, row.names=FALSE, quote=FALSE, sep=" ")

	print(paste("Chromosome ", curr.chrom, " written.", sep=""))
	i <- i + 1
}

}




