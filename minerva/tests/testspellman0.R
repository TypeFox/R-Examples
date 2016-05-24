library(minerva)

myalpha <- 0.6
myc <- 15
data(Spellman)
Spellman <- as.matrix(Spellman)

## Spellman[,2] <- 0
mydata <- Spellman[,1:10]

res <- mine(Spellman,master=1,n.cores=1,alpha=myalpha,C=myc)

res <- mine(mydata, n.cores=1, var.thr=0.0, use="pair")

## aval.cores <- detectCores()
## if (aval.cores > 1){
##   cat("Multicore: On this machinhe you have ",aval.cores," computational cores.\n")
## }
## if (aval.cores > 2){
##   cat("We suggest to exploit minerva parallel computing possibilities with ",aval.cores-1," cores.\n(where possible set 'n.cores = 3')\n")
## cat("Test ok!!!\n")
## }
