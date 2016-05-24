pre2.remove.genos <- function(file.dat, case.ped, control.ped, dir.dat, dir.out, dir.warning=dir.out, perc.snp=10, perc.patient=20, empty="0/0", num.nonsnp.col=5) {
# Remove columns (genos) that have too many missing values.
# All genos that have more than 'perc.snp' values missing in both 
#    case.ped AND control.ped files will be removed.
# All patients that have more than 'perc.patient' values missing will have their IDs 
#    written into "warning.<case.ped>.txt" files. 
# Output will be two clean versions of case.ped and control.ped files in dir.out directory,
#    and optionally the warning files in dir.warning directory.
#
# Example run:
#
# pre2.remove.genos(file.dat="CGEM_Breast_chr11.dat", case.ped="Reshape_CGEM_Breast_chr11CASE.ped", control.ped="Reshape_CGEM_Breast_chr11CONTROL.ped", dir.dat="/home/briollaislab/olga/curr/data/mach1in/hisformat", dir.out="/home/briollaislab/olga/curr/data/mach1in/", dir.warning="/home/briollaislab/olga/curr/data/mach1in/warnings", empty="0/0")
# pre2.remove.genos(file.dat="Clean_genotypes_CASE_chr21.noh.dat", case.ped="Reshaped_genotypes_CASE_chr21.noh.ped", control.ped="Reshaped_genotypes_CONTROL_chr21.noh.ped", dir.dat="/home/briollaislab/olga/curr/data/mach1in/old/extras", dir.out="/home/briollaislab/olga/curr/data/mach1in/old/extras/testremove", perc.snp=10, perc.patient=20)
#
# file.dat: the data file as required for MaCH1, of the format:
#         M SNP1
#         M SNP2
# case.ped: the pedegree data file with CASEs.
# control.ped: the pedegree data file with CONTROLs.
# dir.dat: directory where file.dat and file.ped can be found.
# dir.out: directory to which output files should be saved. 
#        Defaults to the same place as dir.dat.
# dir.warning: directory to which warnings about patients with too many missing SNPs
#        should go. Defaults to the same place as dir.out.
# perc.snp: the percentage (0-100%) of maximum empty values allowed for each geno (column).
# perc.patient: the percentage (0-100%) of empty values allowed for each patient (row).
# empty: is the representation of a missing SNP in the file ("0 0", "0/0", "1/1", "N N", etc)
# num.nonsnp.col: The number of columns that do not contain SNP values. The first columns of the
#        file represent non-SNP values (like patient ID, gender, etc). For MaCH1 input format,
#        the num.nonsnp.col=5, for PLINK it is =6 (due to extra disease status column).
#
# Outputs: 
#
# - <file.dat>.removed.dat - the .dat file containing only the SNPs that were not removed, 
#       will be placed in dir.out directory
# - <case.ped>.removed.ped - the CASE .ped file without columns that contain too many 
#       missing values based on the thresholds perc.snp; in dir.out directory
# - <control.ped>.removed.ped - the CONTROL .ped file without columns that contain too many
#       missing values based on the thresholds perc.snp; in dir.out directory
#
# - warning.<case.ped>.txt - file containing warning messages about patients that have
#       too many SNPs missing (based on perc.patients) in CASE.ped file, after the removal of
#       bad SNPs. 
# - warning.<control.ped>.txt - similar to warning.<case.ped>.txt, only for CONTROL file.
#

# TODO: remove this line:
#source("len.empty.vec.R")

if(missing(file.dat)) stop("Name of the .dat file is required.")
if(missing(case.ped)) stop("Name of the .ped file for CASE is required.")
if(missing(control.ped)) stop("Name of the .ped file for CONTROL is required.")
if(missing(dir.dat)) stop("Name of the directory that contains .dat file is required.")
if(missing(dir.out)) stop("Name of the output directory is required.")


D <- read.table(paste(dir.dat, file.dat, sep="/"), header=F, sep="\t", stringsAsFactors=FALSE)
if(ncol(D) <= 1)
	D <- read.table(paste(dir.dat, file.dat, sep="/"), header=F, sep=" ", stringsAsFactors=FALSE)

Pcase <- read.table(paste(dir.dat, case.ped, sep="/"), header=F, sep="\t", stringsAsFactors=FALSE)
if(ncol(Pcase) <= 1)
	Pcase <- read.table(paste(dir.dat, case.ped, sep="/"), header=F, sep=" ", stringsAsFactors=FALSE)

Pcontrol <- read.table(paste(dir.dat, control.ped, sep="/"), header=F, sep="\t", stringsAsFactors=FALSE)
if(ncol(Pcontrol) <= 1)
	Pcontrol <- read.table(paste(dir.dat, control.ped, sep="/"), header=F, sep=" ", stringsAsFactors=FALSE)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1. Deal with columns: remove SNPs that have few entries
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Combine Pcase and Pcontrol to compute statistics over both
if(ncol(Pcase) != ncol(Pcontrol)) {
	print(paste("Number of columns does not match: ", case.ped, " has ", ncol(Pcase), " and ", control.ped, " has ", ncol(Pcontrol), sep=""))
	return(0)
}

# Make sure .dat file is the correct one: check that #SNPs in .dat and .ped is the same:
if((ncol(Pcase)-num.nonsnp.col) != nrow(D)) {
	print(paste("The dat file ", file.dat, " is inconsistent with .ped ", case.ped, " file: #SNPs in .dat = ", nrow(D), " and in .ped = ", (ncol(Pcase)-num.nonsnp.col), sep=""))
	return(0)
}

P <- rbind(Pcase, Pcontrol)

# Compute the sum of empty values in each column, and express it as a percentage
sums <- apply(P[, (num.nonsnp.col+1):ncol(P)], 2, 'len.empty.vec', empty) / nrow(P) * 100
mask <- sums <= perc.snp

# For both .ped files, keep the first columns as well as masked ones.
Pcase <- Pcase[, c(matrix(TRUE, 1, num.nonsnp.col), mask)]
Pcontrol <- Pcontrol[, c(matrix(TRUE, 1, num.nonsnp.col), mask)]
D <- D[mask,]
nremoved <- length(which(mask==FALSE))
norig <- ncol(Pcase) - num.nonsnp.col + nremoved
print(paste(case.ped, " #SNPs removed: ", nremoved, "/", norig, " = ", (nremoved/norig*100), "%",sep=""))

# Save the results
out.dat <- paste(dir.out, "/", substr(file.dat, 1, (nchar(file.dat)-4)), ".removed.dat", sep="")
write.table(D, file=out.dat, col.names=F, row.names=F, quote=F, sep="\t")

out.ped <- paste(dir.out, "/", substr(case.ped, 1, (nchar(case.ped)-4)), ".removed.ped", sep="")
write.table(Pcase, file=out.ped, col.names=F, row.names=F, quote=F, sep="\t")

out.ped <- paste(dir.out, "/", substr(control.ped, 1, (nchar(control.ped)-4)), ".removed.ped", sep="")
write.table(Pcontrol, file=out.ped, col.names=F, row.names=F, quote=F, sep="\t")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2. Deal with rows: print warning messages for Patients with few entries
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if(norig - nremoved > 0) {

# Compute the sum of empty values in each row, and express it as a percentage
sums <- apply(Pcase[, (num.nonsnp.col+1):ncol(Pcase)], 1, 'len.empty.vec', empty) / ncol(Pcase) * 100
mask <- sums > perc.patient

out.warn <- paste(dir.warning, "/Warning_", substr(case.ped, 1, (nchar(case.ped)-4)), ".txt", sep="")
# Put together the warning message and write to file. 
if(!is.na(match(TRUE, mask))) {
	msg <- paste("Invalid for thresh ", perc.patient, "% after removal of empty SNPs > ", perc.snp, "%: patient ", Pcase[mask,1], " has ", sums[mask], "% SNPs missing", sep="")
	write(msg, file=out.warn, append=FALSE)
} else { # Remove any warning files if they exist for this case.ped
	if(file.exists(out.warn)) 
		file.remove(out.warn)
}

# Same for CONTROL file:
sums <- apply(Pcontrol[, (num.nonsnp.col+1):ncol(Pcontrol)], 1, 'len.empty.vec', empty) / ncol(Pcontrol) * 100
mask <- sums > perc.patient
out.warn <- paste(dir.warning, "/Warning_", substr(control.ped, 1, (nchar(control.ped)-4)), ".txt", sep="")
if(!is.na(match(TRUE, mask))) {
	msg <- paste("Invalid for thresh ", perc.patient, "% after removal of empty SNPs > ", perc.snp, "%: patient ", Pcontrol[mask,1], " has ", sums[mask], "% SNPs missing", sep="")
	write(msg, file=out.warn, append=FALSE)
}else { 
        if(file.exists(out.warn))
                file.remove(out.warn)
}
}
}

