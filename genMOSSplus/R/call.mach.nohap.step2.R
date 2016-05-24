call.mach.nohap.step2 <- function(file.dat, dir.dat, file.ped, dir.file, prefix1="ans.1", prefix="ans", dir.out=".", mach.loc="/software/mach1") { 
# Helper function that calls the second step of MaCH1 algorithm.
#
# file.dat: the data file as required for MaCH1, of the format:
#         M SNP1
#         M SNP2
# file.ped: the pedegree data file name, should be reasonably small s.t.
#       MaCH1 does not run out of memory running it. 
# dir.dat: directory where file.dat is located.
# dir.file: the directory where file.ped can be found
# prefix1: which prefix MaCH1 has used during step1. 
# prefix: what prefix MaCH1 should use to output its files for step2.
# dir.out: the directory in which files with prefix1 are located, and to which
#       output of MaCH1 for step2 will go.
# mach.loc: the location directory where "mach" executable can be found.


step2 <- paste("time ", mach.loc, "/mach1 -d ", dir.dat, "/", file.dat, " -p ", dir.file, "/", file.ped, " --crossover ", dir.out, "/", prefix1, ".rec --errormap ", dir.out, "/", prefix1, ".erate --compact --mle --mldetails --prefix ", dir.out, "/", prefix, sep="")

print("*******************************************************")
print(step2)
print("*******************************************************")

system(step2)

# If the new version of MaCH is outputting .mlgeno.gz file names instead, unzip them to .mlgeno.
possible.mlgeno <- paste(dir.file, "/", prefix, ".mlgeno.gz", sep="")
if(file.exists(possible.mlgeno)) {
	try(system(paste("gzip -d ", possible.mlgeno,  sep="")))
}


}
