
library("sprint")
library("ff")
library("RUnit")

args <- commandArgs(trailingOnly = TRUE)
nreads <- as.numeric(args[1])
print("number of reads")
print(nreads)

filename_ <- args[2]
print("filename")
print(filename_)

# Fill first half of list of strings with 'hello' and second half with 'holla'
a <- rep(c("hello","holla"),2, len = nreads)

system(paste("rm ",filename_))

stime_sprint <- proc.time()["elapsed"]
actual_result <- pstringdistmatrix(a, a, method="h", filename=filename_)
etime_sprint <- proc.time()["elapsed"]

checkEquals(0, actual_result[1,1])
checkEquals(2, actual_result[1,2])
checkEquals(2, actual_result[2,1])
checkEquals(0, actual_result[nreads,nreads])
checkEquals(2, actual_result[nreads,nreads-1])
checkEquals(2, actual_result[nreads-1,nreads])


print(paste("SPRINT pstringDist time: ")); print(paste(etime_sprint-stime_sprint))

print(system(paste("ls -lh|grep ",filename_)))
system(paste("rm ",filename_))
pterminate()

