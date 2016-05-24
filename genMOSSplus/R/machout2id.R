machout2id <- function(file.case, file.control="", dir.file, name.out, dir.out) {
#
# Extracts the patient IDs from the first column of the file, for both
# .ped (MaCH1 input) and .mlgeno (MaCH1 output) formats, 
# from both CASE and (optionally) CONTROL files and save them in name.out file. 
# Just run this once for the smallest chromosome. 
# 
# Sample run:
#
# machout2id(file.case="last23CASE.200.mlgeno", file.control="last23CONTROL.200.mlgeno", dir.file="/home/briollaislab/olga/curr/data/f04_mach1out_result", name.out="CGEM.fam", dir.out="/home/briollaislab/olga/curr/data/f08_step2in")
#
# file.case: the name of the file for CASE.
#    Each line is expected to begin either with:
#    - the patient ID (which should NOT contain '->' symbol), or
#    - the entry of format of MaCH1 output: 
#         id->id2, where id is the patient ID to be extracted.
#    The file should have no header.
# file.control: the name of CONTROL file. If there is no control file,
#    then only individuals from CASE file will be stored.
# dir.file: directory where the file.name can be found.
# dir.out: the directory to which output should go.
#
# Outputs result in dir.out directory:
#   - <name.out> - the list of patient IDs 

case <- read.table(paste(dir.file, file.case, sep="/"), header=FALSE, sep=" ", stringsAsFactors=FALSE)
if(ncol(case) <= 1)
	case <- read.table(paste(dir.file, file.case, sep="/"), header=FALSE, sep="\t", stringsAsFactors=FALSE)

# Extract the values of first column and split them by '->'
patient.name <- unlist(strsplit(case[,1], split="->"))

# If each entry has exactly one '->', then use the first half of each patient name as the id:
if(length(patient.name) == (nrow(case)*2))
	patient.name <- patient.name[seq(1,nrow(case)*2,by=2)]
else  # leave the name as is
	patient.name <- case[,1]

# Read in and add control if needed.
if(file.control!="") {
	control <- read.table(paste(dir.file, file.control, sep="/"), header=FALSE, sep=" ", stringsAsFactors=FALSE)
	if(ncol(control) <= 1)
		control <- read.table(paste(dir.file, file.control, sep="/"), header=FALSE, sep="\t", stringsAsFactors=FALSE)

	patient.name.control <- unlist(strsplit(control[,1], split="->"))
	if(length(patient.name.control) == (nrow(control)*2))
		patient.name.control <- patient.name.control[seq(1,nrow(control)*2,by=2)]
	else
		patient.name.control <- control[,1]

	patient.name <- t(cbind(t(patient.name), t(patient.name.control)))
}

full.name.out <- paste(dir.out, name.out, sep="/")
write.table(patient.name, file=full.name.out, quote=FALSE, sep="\t", row.names=F, col.names=F)

}
