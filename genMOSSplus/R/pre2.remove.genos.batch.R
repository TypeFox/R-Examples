pre2.remove.genos.batch <- function(dir.dat, dir.ped=dir.dat, dir.out, dir.warning=dir.out, perc.snp=10, perc.patient=20, empty="0/0", num.nonsnp.col=5, prefix.dat, prefix.case, prefix.control, key.dat="", key.case="CASE", key.control="CONTROL", ending.dat=".dat", ending.case=".ped", ending.control=".ped") {
# Removes SNPs that contain more than perc.snp empty geno values,
# from all the corresponding .ped and .dat files in directory dir.dat.
# 
# Example:
# pre2.remove.genos.batch("/home/briollaislab/olga/curr/data/mach1in/hisformat", dir.out="/home/briollaislab/olga/curr/data/mach1in", dir.warning="/home/briollaislab/olga/curr/data/mach1in/warnings", empty="0/0")
#
# dir.dat - directory containing .ped and .dat files.
#    - If .ped files for some chromosomes are split into several files, these
#      files will be concatenated alphabetically.
#      Note: both CASE and CONTROL files must exist for each chromosome.
#    - .dat files must exist for all the chromosomes.  
#    Both files .ped and .dat must be tab separated. 
#        

if(missing(dir.dat)) stop("Name of input directory for .dat files must be provided.")
if(missing(dir.out)) stop("Name of output directory must be provided.")
if(missing(prefix.dat)) stop("Prefix of the .dat file name must be provided.")
if(missing(prefix.case)) stop("Prefix of the CASE .ped file name must be provided.")
if(missing(prefix.control)) stop("Prefix of the CONTROL .ped file name must be provided.")

# TODO: remove this line:
#source("pre2.remove.genos.R")
#source("get.file.name.R")
#source("get.chrom.num.R")


# *******************************************
# 1. Obtain all .dat, CASE, and CONTROL files
all.dat <- get.file.name(dir=dir.dat, prefix=prefix.dat, key=key.dat, ending=ending.dat)
all.case <- get.file.name(dir=dir.ped, prefix=prefix.case, key=key.case, ending=ending.case)
all.control <- get.file.name(dir=dir.ped, prefix=prefix.control, key=key.control, ending=ending.control)

if(length(all.dat) == 0 || length(all.case) == 0 || length(all.control) == 0)
	return()

# *******************************************
# 2. Combine multiple files for one chromosome into one:

new.case <- combine.same.chrom(dir.in=dir.ped, file.name=all.case, prefix=prefix.case, ending=ending.case)
new.control <- combine.same.chrom(dir.in=dir.ped, file.name=all.control, prefix=prefix.control, ending=ending.control)


# *******************************************
# 3. Match .ped and .dat and run the pre2.remove.genos()
 
chroms.case <- get.chrom.num(new.case, prefix=prefix.case)
chroms.control <- get.chrom.num(new.control, prefix=prefix.control)
chroms.dat <- get.chrom.num(all.dat, prefix=prefix.dat)

chroms.common <- intersect(chroms.case, chroms.control)
chroms.common <- intersect(chroms.common, chroms.dat)

i <- 1
while (i <= length(chroms.common)) {
	curr.chrom <- chroms.common[i]

	curr.dat <- all.dat[match(curr.chrom, chroms.dat)]
	curr.case <- new.case[match(curr.chrom, chroms.case)]
	curr.control <- new.control[match(curr.chrom, chroms.control)]

	pre2.remove.genos(file.dat=curr.dat, case.ped=curr.case, control.ped=curr.control, dir.dat=dir.dat, dir.out=dir.out, dir.warning=dir.warning, perc.snp=perc.snp, perc.patient=perc.patient, empty=empty, num.nonsnp.col=num.nonsnp.col)

        i <- i + 1
}

}



# Helper function to combine multiple files into one.
# Returns the list of resultant files, even those that haven't been combined.
combine.same.chrom <- function(dir.in, file.name, prefix, ending) {

	chroms <- get.chrom.num(file.name, prefix=prefix)
	uniq.shred <- unique(chroms)
	names <- rep("", length(uniq.shred))

	i <- 1
	while (i <= length(uniq.shred)) {
		shred.id <- which(chroms == uniq.shred[i])
		shred.name <- file.name[shred.id[1]]

		new.name <- paste(substr(shred.name, 1, nchar(shred.name) - nchar(ending)), ".all", ending, sep="")

		check.all <- grep(paste(".all", ending, "$", sep=""), file.name[shred.id])
		if(length(check.all) > 0)
			new.name <- file.name[check.all[1]]

		new.name.path <- paste(dir.in, "/", new.name, sep="")
		# Save the names of useful files to return later.
		if(length(shred.id) == 1)
			names[i] <- shred.name
		else
			names[i] <- new.name

		# Combine files for a Chromosome if it is saved in more than 1 file
		if(!file.exists(new.name.path) && length(shred.id) > 1) {
			shred.same <- file.name[shred.id]
			shred.comand <- paste(dir.in, "/", shred.same, sep="", collapse=" ")
			shred.comand <- paste("cat ", shred.comand, " > ", new.name.path, sep="")
			print(paste("Combining into: ", new.name, sep=""))
			system(shred.comand) # <----------------------------------------------------------------
		}
		i <- i + 1
	}

	return(names)


}



