library("sprint")
library("ff")

args <- commandArgs(trailingOnly = TRUE)
rows<- as.numeric(args[1])
print("number of rows")
print(rows)
cols<- as.numeric(args[2])
print("number of columns")
print(cols)
filename<- args[3]
print("filename")
print(filename)


system(paste("rm ",filename))

my.matrix <- matrix(rnorm(500000,9,1.7), nrow=rows, ncol=cols)

genecor <- pcor( t(my.matrix), filename_  =filename )

print(system(paste("ls -lh|grep ",filename)))
system(paste("rm ",filename))