clean.rsq.combine <- function(case.file, control.file, dir.file, case.info, control.info, out.name.start, rsq.thresh=0.5, dir.out=dir.file, separ=" ") {
#
# Removes all SNPs whose RSQ is lower than the threshold.
# Combines CASE and CONTROL files into one file, ppends disease status as the last column.
# The disease status is encoded as 1 for CASE and 0 for CONTROL. 
# The info file is expected to be MaCH's output .mlinfo: 7 columns,
# first col is SNP name, and last col is RSQ.

# If you run genos2numeric() AFTER this function, don't forget to specify that 
#   there is 1 ending non-SNP column. 
# Example run:
#
# clean.rsq.combine(case.file="machCASE.mlgeno",  control.file="machCONTROL.mlgeno", dir.file="/home/briollaislab/olga/funtry/test_tune/seg1/working/s02_machout", case.info="machCASE.mlinfo", control.info="machCONTROL.mlinfo", out.name.start="data20", rsq.thresh=0.5, dir.out="/home/briollaislab/olga/funtry/test_tune/seg1/working/s03_combined") 
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
if(missing(out.name.start)) stop("The output file name must be provided.")
if(missing(case.file)) stop("The CASE file name must be provided.")
if(missing(control.file)) stop("The CONTROL file name must be provided.")
if(missing(case.info)) stop("The CASE .mlinfo file name must be provided.")
if(missing(control.info)) stop("The CONTROL .mlinfo file name must be provided.")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 0. Read in all 4 files
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

Dcase <- read.table(paste(dir.file, case.info, sep="/"), header=T, sep="\t", stringsAsFactors=FALSE)
if(ncol(Dcase) <= 1)
        Dcase <- read.table(paste(dir.file, case.info, sep="/"), header=F, sep=" ", stringsAsFactors=FALSE)

Dcontrol <- read.table(paste(dir.file, control.info, sep="/"), header=T, sep="\t", stringsAsFactors=FALSE)
if(ncol(Dcontrol) <= 1)
        Dcontrol <- read.table(paste(dir.file, control.info, sep="/"), header=F, sep=" ", stringsAsFactors=FALSE)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1. Find all SNPs whose RSQ is bigger than threshold.
#    Create new .dat file with those SNPs in MaCH input format.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
I <- intersect(which(Dcase[,7]>=rsq.thresh), which(Dcontrol[,7]>=rsq.thresh)) 

newD <- cbind("M", Dcase[I,1])
name.dat <- paste(out.name.start, ".dat", sep="")
write.table(newD, file=paste(dir.out, name.dat, sep="/"), col.names=F, row.names=F, quote=F, sep="\t")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2. Extract the required SNPs from both .mlgeno files: 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
newPcase <- cbind(Pcase[,1:2], Pcase[,(I+2)], 1)
newPcontrol <- cbind(Pcontrol[,1:2], Pcontrol[,(I+2)], 0)

name.ped <- paste(out.name.start, ".txt", sep="")
write.table(newPcase, file=paste(dir.out, name.ped, sep="/"), col.names=F, row.names=F, quote=F, sep=" ")
write.table(newPcontrol, file=paste(dir.out, name.ped, sep="/"), col.names=F, row.names=F, quote=F, append=TRUE, sep=" ")

return(list(ped=name.ped, dat=name.dat))

}

