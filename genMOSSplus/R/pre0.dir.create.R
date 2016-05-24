pre0.dir.create <- function(dir.out=".", out.name="newdata", prefix.dir="d") {
# Function to help create the recommended subdirectory structure for the pre-processing.
# In dir.out a directory with name out.name will be created. 
# Inside of this out.name directory will be a set of subdirectories, 
# whose names will begin with prefix.dir, followed by a number, 
# followed by short description of what the folder is designed to contain.
#
#
# Returns the list of all the full path names:
# out$d0 - directory to which you should copy your original dataset in any format
# out$d1 - directory into which plink format of your original dataset should be stored;
#          note, there is no function in this package for this conversion,
#          write your own function, similar to ex2plink()
# out$d2 - 
# ...
# out$d11 -

dir.main <- paste(dir.out, out.name, sep="/") 
d0  <- paste(dir.main, "/", prefix.dir, "00_original", sep="")
d1  <- paste(dir.main, "/", prefix.dir, "01_plink", sep="")
d2  <- paste(dir.main, "/", prefix.dir, "02_machin", sep="")
d3  <- paste(dir.main, "/", prefix.dir, "03_removed", sep="")
d4  <- paste(dir.main, "/", prefix.dir, "04_ref", sep="")
d5  <- paste(dir.main, "/", prefix.dir, "05_machout", sep="")
d6  <- paste(dir.main, "/", prefix.dir, "06_combined", sep="")
d7  <- paste(dir.main, "/", prefix.dir, "07_numeric", sep="")
d8  <- paste(dir.main, "/", prefix.dir, "08_merged", sep="")
d9  <- paste(dir.main, "/", prefix.dir, "09_confound", sep="")
d10 <- paste(dir.main, "/", prefix.dir, "10_split", sep="")
d11 <- paste(dir.main, "/", prefix.dir, "11_subset", sep="")


dir.create(dir.main)
dir.create(d0)
dir.create(d1)
dir.create(d2)
dir.create(d3)
dir.create(d4)
dir.create(d5)
dir.create(d6)
dir.create(d7)
dir.create(d8)
dir.create(d9)
dir.create(d10)
dir.create(d11)
#dir.create(d12)

return(list(d0=d0, d1=d1, d2=d2, d3=d3, d4=d4, d5=d5, d6=d6, d7=d7, d8=d8, d9=d9, d10=d10, d11=d11)) #, d12=d12))

}
