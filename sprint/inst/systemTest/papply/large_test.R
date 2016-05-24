
library("sprint")
library("ff")
library("RUnit")

#args <- commandArgs(trailingOnly = TRUE)
#total <- as.numeric(args[1])
#print("total")
#print(total)

#filename_ <- args[2]
#print("filename")
#print(filename_)

filename_ <- "papply.out"


ffobjectt = ff(sin(1:10000), vmode="double", dim=c(200,50))
matrixx = matrix(sin(1:10000), ncol=50)

dim(ffobjectt[,])
summary(matrixx)

ffobjectt[1,1]
matrixx[1,1]
ffobjectt[1,2]
matrixx[1,2]
ffobjectt[200,50]
matrixx[200,50]

#  expected_result = apply(ffobjectt, 1, mean)
       
	
#    checkEquals(as.vector(expected_result), unlist(papply_result), " papply mean function on matrix over rows")


system(paste("rm ",filename_))

expected_result <- apply(matrixx, 1, mean)
summary(expected_result)
expected_result

stime_sprint <- proc.time()["elapsed"]
actual_result <- papply(ffobjectt, mean, 1, out_filename=filename_)
etime_sprint <- proc.time()["elapsed"]
actual_result

papply_result = papply(matrixx, mean, 1)   
papply_result
class(papply_result)

checkEquals(as.vector(expected_result), unlist(actual_result[]), " papply mean function on matrix over rows")
checkEquals(unlist(papply_result), unlist(actual_result[]), " papply mean function on matrix over rows")

print(paste("SPRINT pstringDist time: ")); print(paste(etime_sprint-stime_sprint))

print(system(paste("ls -lh|grep ",filename_)))
system(paste("rm ",filename_))
pterminate()

