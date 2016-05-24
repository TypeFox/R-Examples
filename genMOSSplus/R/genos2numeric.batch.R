genos2numeric.batch <- function(dir1, dir2, ending="", num.nonsnp.col=2, num.nonsnp.last.col=1, dir.dat, key.dat) {
# Provided with .ped files in dir1, whose names end with 'ending', 
# which have NO missing genotypes, 
# converts the genos into numeric values 1, 2, and 3. 
#
# num.nonsnp.col - how many first columns in the file format do not contain SNP information.
# 
# Example:
# genos2numeric.batch(dir1="/home/briollaislab/olga/curr/data/mach1in/result", dir2="/home/briollaislab/olga/curr/data/mach1out", ending="mlgeno")
#

# TODO: remove this line:
#source("genos2numeric.R")

# Obtain all files in the directory
all.files <- dir(dir1)

# 
file.id <- grep(ending, all.files)
num.file <- length(file.id)

if (num.file == 0) {
	print(paste("Error, no files with extension '", ending, "' were found in the directory: ", dir1, sep=""))
	return()
}

all.files <- all.files[file.id]

# now with .dat:
dat.files <- dir(dir.dat)
dat.id <- grep(key.dat, dat.files)
if(length(dat.id)== 0) {
	print(paste("Error, no files containing '", key.dat, "' were found in the directory: ", dir.dat, sep=""))
        return()
}
dat.files <- dat.files[dat.id]

# Check that there is a ".dat" in the filename
dat.id <- grep(".dat", dat.files)
if(length(dat.id)== 0) {
        print(paste("Error, no files containing '", key.dat, "' AND ending with .dat were found in the directory: ", dir.dat, sep=""))
        return()
}
dat.files <- dat.files[dat.id]

if(length(dat.files) != length(all.files)) {
	print(paste("Error: the number of data files is ", length(all.files), " but number of .dat files is: ", length(dat.files), sep=""))
	return()
}


i <- 1
while(i <= num.file) {
	print(paste("file: ", all.files[i], ", dat: ", dat.files[i], sep=""))
	genos2numeric(file.ped=all.files[i], dir.ped=dir1, dir.out=dir2, num.nonsnp.col=num.nonsnp.col, num.nonsnp.last.col=num.nonsnp.last.col, remove.bad.genos=TRUE, file.dat=dat.files[i], dir.dat=dir.dat)
	i <- i + 1
}


}
