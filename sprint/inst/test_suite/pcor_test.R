library("sprint")
args <- commandArgs(trailingOnly = TRUE)
rows<- as.numeric(args[1])
print("number of rows")
print(rows)
cols<- as.numeric(args[2])
print("number of columns")
print(cols)
stime <- proc.time()["elapsed"]
my.matrix <- matrix(rnorm(10000,9,1.7), nrow=rows, ncol=cols)
etime <- proc.time()["elapsed"]
print(paste("Correlation time: ")); print(paste(etime-stime))
genecor <- pcor( t(my.matrix), filename_  ="pcor.out" )
print(system("ls -lh|grep pcor.out"))
#print(genecor)

