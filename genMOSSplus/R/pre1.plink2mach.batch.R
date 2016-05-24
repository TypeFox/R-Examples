pre1.plink2mach.batch <- function(dir.in, dir.out, ending.ped=".ped", ending.map=".map", key.ped="", key.map="") {
# Provided with Plink-format .ped and .map files in dir.in, 
# re-format it into MACH pedirgree (file.ped) and data (file.dat) file fromats,
# and saves the reformatted files in dir.out.
# 
# Example:
# pre1.plink2mach.batch("/home/briollaislab/shared/Reshape", dir.out="/home/briollaislab/olga/curr/data/mach1in")
#

# TODO: remove this line:
#source("pre1.plink2mach.R")
#source("get.file.name.R")

if(missing(dir.in)) stop("Input directory (containing .ped and .map files) is required")
if(missing(dir.out)) stop("Output directory is required")

all.peds <- get.file.name(dir=dir.in, key=key.ped, ending=ending.ped)
all.maps <- get.file.name(dir=dir.in, key=key.map, ending=ending.map)

if(length(all.peds) == 0 || length(all.maps) == 0)
	return()

# Separately first do all the .ped files
i <- 1
while(i <= length(all.peds)) {
	print(paste("ped: ", all.peds[i], sep=""))
	pre1.plink2mach(all.peds[i], "", dir.in, dir.out)
	i <- i + 1
}

# Next do all the .map files
i <- 1
while(i <= length(all.maps)) {
        print(paste("map: ", all.maps[i], sep=""))
        pre1.plink2mach("", all.maps[i], dir.in, dir.out)
        i <- i + 1
}

}
