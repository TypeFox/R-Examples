pre4.combine.case.control <- function(case.file, control.file, dir.file, name.out, dir.out=dir.file, separ=" ") {
# Combines CASE and CONTROL files into one file, and appends disease status as the last column.
# The disease status is encoded as 1 for CASE and 0 for CONTROL. 
# If you run genos2numeric() AFTER this function, don't forget to specify that 
#   there is 1 ending non-SNP column. 
# Example run:
#
# pre4.combine.case.control(dir.file="/home/briollaislab/olga/curr/data/mach1out/notabs", dir.out="/home/briollaislab/olga/curr/data/mach1out/", separ=" ", case.file="Num_last22CASE.200.mlgeno", control.file="Num_last22CONTROL.200.mlgeno", name.out="CGEM_Breast.22.txt")
#
# case.file: the name of CASE file
# control.file: the name of CONTROL file
# dir.file: the directory where CASE and CONTROL input files live.
# name.out: the desired name for the output file.
# dir.out: the directory to which output file should be written; 
# separ: the separator used in the CASE and CONTROL input files.
#
# Outputs: 
#
# - <dir.out>/<name.out> - the file containing both CASE and CONTROL values, 
#     with the disease status as the last column.
# - <dir.out>/<all.dat> - also will copy over all the files ending with ".dat" 
#     that exist in dir.file

if(missing(dir.file)) stop("Name of input directory must be provided.")
if(missing(name.out)) stop("The output file name must be provided.")
if(missing(case.file)) stop("The CASE file name must be provided.")
if(missing(control.file)) stop("The CONTROL file name must be provided.")

# TODO: remove
#source("get.file.copy.R")

Pcase <- read.table(paste(dir.file, case.file, sep="/"), header=F, sep=separ, stringsAsFactors=FALSE)
Pcontrol <- read.table(paste(dir.file, control.file, sep="/"), header=F, sep=separ, stringsAsFactors=FALSE)

if(ncol(Pcase) <= 1) {
	print(paste("Error: the CASE file's separator is not '", separ, "'.", sep=""))
	return(0)
}
if(ncol(Pcontrol) <= 1) {
        print(paste("Error: the CONTROL file's separator is not '", separ, "'.", sep=""))
        return(0)
}

# Combine Pcase and Pcontrol
if(ncol(Pcase) != ncol(Pcontrol)) {
	print(paste("Error: Number of columns does not match: ", case.file, " has ", ncol(Pcase), " and ", control.file, " has ", ncol(Pcontrol), sep=""))
	return(0)
}

P <- rbind(Pcase, Pcontrol)

# Append a column of 1s for CASE and 0s for CONTROL
Dcase <- rep(1, nrow(Pcase))
Dcontrol <- rep(0, nrow(Pcontrol))

D <- rbind(as.matrix(Dcase), as.matrix(Dcontrol))

PD <- cbind(P, D)

# Save the results
out <- paste(dir.out, "/", name.out, sep="")
write.table(PD, file=out, col.names=F, row.names=F, quote=F, sep="\t")

# Copy all .dat files from the input directory over to output directory:
get.file.copy(dir.in=dir.file, dir.out=dir.out, ending=".dat", verbal=FALSE)

}

