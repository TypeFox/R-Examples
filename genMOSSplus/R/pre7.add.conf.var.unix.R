pre7.add.conf.var.unix <- function(file.name, dir.file, file.fam, dir.fam=dir.file, file.conf, dir.conf=dir.file, file.out, fam.out=file.fam, dir.out) {
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


#TODO: remove this line:
#source("get.data.dims.R")
#source("get.file.copy.R")

if(missing(file.name)) stop("Name of the binary data file must be provided.")
if(missing(dir.file)) stop("Name of input directory for binary data file must be provided.")
if(missing(file.fam)) stop("Name of the family file must be provided")
if(missing(file.conf)) stop("Name of the file with confounding variables must be provided")
if(missing(file.out)) stop("Name of the output file must be provided")
if(missing(dir.out)) stop("Name of output directory must be provided.")

# Read in the three files: binary data, family, and conf vars (New variables):
#D <- read.table(paste(dir.file, file.name, sep="/"), sep=" ", stringsAsFactors=FALSE)

#if(ncol(D) <= 1) {
#print("try tab")
#	D <- read.table(paste(dir.file, file.name, sep="/"), sep="\t", stringsAsFactors=FALSE)
#}
#print("dim D init:")
#print(dim(D))

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

# Now MF can be used to remove patients from F who do not appear in N:
boolMF <- (MF != 0)
#D <- D[boolMF,]
F <- F[boolMF,]

# Use linux command to remove lines from ped file:
# First, identify the line numbers
lines.remove <- which(MF == 0)
file.ped.tmp <- paste(dir.out, "/", file.name, ".temp", sep="")
if(length(lines.remove) > 0) {
	print("removing: ")
	print(lines.remove)

	command.rem <- paste(lines.remove, "d", sep="", collapse=";")
	command.rem <- paste("sed -e '", command.rem, "' ", dir.file, "/", file.name, " > ", file.ped.tmp, sep="")
	print(command.rem)
	system(command.rem)
} else {
	file.copy(paste(dir.file, file.name, sep="/"), file.ped.tmp)
}

# Create a short version of MF: without 0s - these indices will correspond to the new D and new F
MF <- MF[boolMF]

# Use that new MF to re-arrange the order of N s.t. it corresponds to D and F;
# Note: if N were to have any additional patients that are not in F, these will be ignored automatically
N <- N[MF,]

# Append the disease status (last column of D) to N, remove patient IDs (col 1 of N), 
# and run Adrian's discretization step to make N binary.

# Find the number of the last column of the new .ped file,
# Then extract its last column (disease status)
file.dis <- paste(dir.out, "/", file.name, ".disease", sep="")
ped.dims <- get.data.dims(file.ped.tmp)
print(ped.dims)
command.dis.status <- paste("cut -f ", ped.dims$ncols, " ", file.ped.tmp, " > ", file.dis, sep="")
print(command.dis.status)
system(command.dis.status)
dis <- read.table(file.dis, header=FALSE, sep=" ", stringsAsFactors=FALSE)


newN <- cbind(N[,2:ncol(N)], dis)
new.name <- paste(dir.out, "tmpN.txt", sep="/")

write.table(newN, file=new.name, col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")
#new.name.bin <- pre6.discretize(file.train=new.name, dir.file=dir.out, test.output.append="useless")

# Append the output file's info to the end of new .ped D (without last column)

command.paste <- paste("cut -f 1-", (ped.dims$ncols-1), " ", file.ped.tmp, " | paste - ", new.name, " > ", dir.out, "/", file.out, sep="")
print(command.paste)
system(command.paste)

#newD <- cbind(D[,1:(ncol(D)-1)], binN)
#write.table(newD, file=paste(dir.out, file.out, sep="/"), sep=" ", col.names=FALSE, row.names=FALSE, quote=FALSE)
write.table(F, file=paste(dir.out, fam.out, sep="/"), sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)


# clean up temporary files
file.remove(new.name)
file.remove(file.ped.tmp)
file.remove(file.dis)


}


