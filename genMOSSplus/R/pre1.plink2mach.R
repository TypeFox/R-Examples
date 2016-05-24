pre1.plink2mach <- function(file.ped="", file.map="", dir.in, dir.out) {
# Provided with Plink-format files file.ped and file.map in dir.in, 
# re-format it into MACH pedigree (file.ped) and data (file.dat) file fromats.
# and saves the reformatted files in dir.out.
#
# Note: if file.ped or file.map are empty strings (""), then this
#    file will not be processed (so that you can use this function 
#    to do ONLY PED files but not map, and vice versa) 
#
# Note: does NOT change unknown Allele values from "0" to "N", as MACH program can use either.
# Does NOT recode gender to "M" and "F", since MaCH1 doesn't care, but
# further file processing interprets "F" as "FALSE".

if(missing(dir.in)) stop("Input directory must be provided")
if(missing(dir.out)) stop("Output directory must be provided")

if(file.map != "") {
	# Deal with .map file first: "M  SNPID"
	map.file <- read.table(paste(dir.in, file.map, sep="/"), header=FALSE, sep="\t", stringsAsFactors=FALSE)
	if(ncol(map.file) <= 1) 
		map.file <- read.table(paste(dir.in, file.map, sep="/"), header=FALSE, sep=" ", stringsAsFactors=FALSE)

	m.code <- rep('M', nrow(map.file))
	map.result <- data.frame(m.code, map.file[,2])
	write.table(map.result, file=paste(dir.out, paste(substr(file.map, 1, nchar(file.map)-3), "dat", sep=""), sep="/"), col.names=F, row.names=F, quote=F, sep="\t")
}

if(file.ped != "") {
	# Now the .ped file: same as old ped file, without 6th column (disease status):
	ped.file <- read.table(paste(dir.in, file.ped, sep="/"), header=FALSE, sep="\t", stringsAsFactors=FALSE)
	if(ncol(ped.file) <= 1)
		ped.file <- read.table(paste(dir.in, file.ped, sep="/"), header=FALSE, sep=" ", stringsAsFactors=FALSE)

	status <- ped.file[,6]
	status.uniq <- unique(status)

	# For each allele, take first and last character and separate them by a slash "/", then format them back into matrix shape.
	#orig.alleles <- ped.file[, 7:ncol(ped.file)]
	#slash.alleles <- matrix(paste(substr(orig.alleles, 1,1), substr(orig.alleles, nchar(orig.alleles), nchar(orig.alleles)), sep="/"), nrow(orig.alleles), ncol(orig.alleles))
#print(orig.alleles)
#print("***********")
#return(orig.alleles)
#print(substr(orig.alleles, 1,1))	
#print("KKKKKKKKK")
#print(substr(orig.alleles, nchar(orig.alleles), nchar(orig.alleles)))

	ped.file <- ped.file[, c(1:5, 7:ncol(ped.file))]

#	ped.file <- cbind(ped.file[,1:5], slash.alleles)

	if(length(status.uniq) == 2) {
		# The larger value is CASE, smaller is the CONTROL
		# Swap them if default assignment was bad.
		id.case <- status.uniq[1]
		id.control <- status.uniq[2]
		if(id.case < id.control) {
			id.case <- status.uniq[2]
			id.control <- status.uniq[1]
		}
				
		ped.case <- ped.file[which(status==id.case),]
		ped.control <- ped.file[which(status==id.control),]

		write.table(ped.case, file=paste(dir.out, "/", substr(file.ped, 1, nchar(file.ped)-4), "CASE.ped", sep=""), col.names=F, row.names=F, quote=F, sep="\t")
		write.table(ped.control, file=paste(dir.out, "/", substr(file.ped, 1, nchar(file.ped)-4), "CONTROL.ped", sep=""), col.names=F, row.names=F, quote=F, sep="\t")
	
	} else {

		if(length(status.uniq) > 2) {
			print(paste("Warning: The disease status column has more than 2 values in file: ", file.ped, sep=""))
			print(status.uniq, collapse=" ")
		} else 
			print(paste("Warning: The disease status column has less than 2 values in file: ", file.ped, sep=""))

		# Save the result anyway into one file
		write.table(ped.file, file=paste(dir.out, file.ped, sep="/"), col.names=F, row.names=F, quote=F, sep="\t")
	}

}
}
