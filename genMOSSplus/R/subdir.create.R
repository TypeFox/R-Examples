subdir.create <- function(dir.out=".", out.name="subset1", prefix.dir="s") {
# Function to help create the recommended sub-subdirectory structure for the pre-processing
# of segments of the dataset.
# In dir.out a directory with name out.name will be created. 
# Inside of this out.name directory will be a 'working' directory that contains a
# set of subdirectories, whose names will begin with prefix.dir, followed by a number, 
# followed by short description of what the folder is designed to contain.
#
# Note: dir.out is intended to be 'd11_segment' folder of the work-dataset.
#
# Returns the list of all the full path names:
# out$s0 - directory of the prepared subset dataset will go
# out$s1 - subdirectory to which trimmed/small/relevantSNP data will go
# out$s2 - subdirectory to which re-computation of Mach1 with hapmap will go
# out$s3 - subdirectory into which CASE and CONTROL will be combined
# out$s4 - subdirectory into which genos will be converted to numeric values 

dir.main <- paste(dir.out, out.name, sep="/") 
d1  <- paste(dir.main, "/working/", prefix.dir, "01_trimmed", sep="")
d2  <- paste(dir.main, "/working/", prefix.dir, "02_machout", sep="")
d3  <- paste(dir.main, "/working/", prefix.dir, "03_combined", sep="")
d4  <- paste(dir.main, "/working/", prefix.dir, "04_numeric", sep="")
dw  <- paste(dir.main, "working", sep="/")


if(!file.exists(dir.main))
	dir.create(dir.main)
if(!file.exists(dw))
	dir.create(dw)
if(!file.exists(d1))
	dir.create(d1)
if(!file.exists(d2))
	dir.create(d2)
if(!file.exists(d3))
	dir.create(d3)
if(!file.exists(d4))
	dir.create(d4)

return(list(s0=dir.main, s1=d1, s2=d2, s3=d3, s4=d4))

}
