call.mach.nohap.step1 <- function(file.dat, dir.dat, out.ped, num.iters=2, prefix1="ans.1", dir.out=".", mach.loc="/software/mach1") { 
# Helper function that calls the first step of MaCH1 algorithm.
# 
# file.dat: the data file as required for MaCH1, of the format:
#         M SNP1
#         M SNP2
# out.ped: the FULL path and file name of pedegree data file, 
#       should be randomly sampled set of 200-500 individuals.
# dir.dat: the directory where file.dat can be found
# prefix1: which prefix MaCH1 should be using for step1. 
# prefix: what prefix MaCH1 should use to output its files for step2.
# dir.out: the directory into which files with prefix1 should be placed.
# mach.loc: the location directory where "mach" executable can be found.


#prefix1 <- paste(prefix, ".1", sep="")

step1 <- paste("time ", mach.loc, "/mach1 -d ", dir.dat, "/", file.dat, " -p ", out.ped, " -r ", num.iters, " --compact --prefix ", dir.out, "/", prefix1, sep="")

print("*******************************************************")
print(step1)
print("*******************************************************")

system(step1)

}
