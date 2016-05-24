pre7.add.conf.var <- function(file.name, dir.file, file.fam, dir.fam=dir.file, file.conf, dir.conf=dir.file, file.out, fam.out=file.fam, dir.out) {
# Appends confounding variables listed in file.conf to the end of the file.name, 
# right before the disease status (last) column. The output will contain only the
# patients for which confounding variables exist (other patients will be omitted),
# so new family file will be written. 
# Note: file.name must be in binary format, with last column as disease status, 
#    tab separated
# Note: file.conf must have the following format:
#    patientID1	1	2 ...
#    patientID2	3	1 ...
#    patientID3	2	2 ...
#    ...
#    - Column 1: patient ID, exactly the same names should appear in file.fam; 
#               * order does not matter; 
#               * some patients may be missing;
#               * no new patients should appear in file.conf (if they don't exist in file.fam)
#    - Column 2: the confounding variable must have no more than 3 different values.
#    - Other columns are optional, may be included if there are more confounding variables (3 categories each)
#    - No header
#    - Tab separated
#    - No missings or NAs
#
# Example run:
# 
#
# file.name: name of the binary data file, whose last column is the disease status. 
# dir.file: directory where file.name can be found.
# file.fam: name of the family file: one column: one patient ID per line.
# dir.fam: directory where file.fam can be found.
# file.conf: The name of the file that contains confounding variable information in
#     format described above.
# dir.conf: directory where file.conf can be found.
# file.out: name of output file, which will contain all information of file.name, plus
#     confounding variables, only for the patients mentioned in file.conf.
# fam.out: name of the family output file.
# dir.out: directory to which file.out and fam.out should be saved.
#
# Result:
# <file.out> - in dir.out directory, the resultant binary file: 
# <fam.out> - in dir.out directory, the corresponding .fam file, will be different
#      from original <file.fam> if some patients were missing in file.conf. 
#

if(missing(file.name)) stop("Name of the binary data file must be provided.")
if(missing(dir.file)) stop("Name of input directory for binary data file must be provided.")
if(missing(file.fam)) stop("Name of the family file must be provided")
if(missing(file.conf)) stop("Name of the file with confounding variables must be provided")
if(missing(file.out)) stop("Name of the output file must be provided")
if(missing(dir.out)) stop("Name of output directory must be provided.")

# Read in the three files: binary data, family, and conf vars (New variables):
D <- read.table(paste(dir.file, file.name, sep="/"), sep="\t", stringsAsFactors=FALSE)
F <- read.table(paste(dir.fam, file.fam, sep="/"), sep=" ", stringsAsFactors=FALSE)
N <- read.table(paste(dir.conf, file.conf, sep="/"), sep="\t", stringsAsFactors=FALSE)

indx.file <- paste(dir.out, "/", file.conf, ".indx.txt", sep="")

MF <- array(0, nrow(F))

if(file.exists(indx.file)) {
	MF <- read.table(indx.file, sep="\t", stringsAsFactors=FALSE)
} else {
	# Iterate over the variables in N and create a mask MF that keeps for each position in F
	# the index of that patient ID in N (or 0 if it doesn't appear in N). 
	i <- 1
	while(i <= nrow(N)) {
		curr.indx <- which(F[,1] == N[i,1])[1]
		if(!is.na(curr.indx))
			MF[curr.indx] <- i
		i <- i + 1
	}

	write.table(MF, file=indx.file, col.names=FALSE,  row.names=FALSE, quote=FALSE, sep="\t")
}

# Now MF can be used to remove patients from D and F who do not appear in N:
boolMF <- (MF != 0)
D <- D[boolMF,]
F <- F[boolMF,]

# Create a short version of MF: without 0s - these indices will correspond to the new D and new F
MF <- MF[boolMF]

# Use that new MF to re-arrange the order of N s.t. it corresponds to D and F;
# Note: if N were to have any additional patients that are not in F, these will be ignored automatically
N <- N[MF,]

# Append the disease status (last column of D) to N, remove patient IDs (col 1 of N) 
print("dim D: ")
print(dim(D))
print("ncol D:")
print(ncol(D))
print("len D:")
print(length(D[,ncol(D)]))
#print("dim N: ")
#print(dim(N[,2:ncol(N)]))
print("len N:")
print(length(N[,2:ncol(N)]))


newN <- cbind(N[,2:ncol(N)], D[,ncol(D)])

# Read in the output file and append its info to the end of D
newD <- cbind(D[,1:(ncol(D)-1)], newN)
write.table(newD, file=paste(dir.out, file.out, sep="/"), sep=" ", col.names=FALSE, row.names=FALSE, quote=FALSE)
write.table(F, file=paste(dir.out, fam.out, sep="/"), sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)

}


