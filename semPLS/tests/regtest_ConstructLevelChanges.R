### Fixed Bug with version 1.0-3 (2011-10-06)
### as reply to e-mail by Pratyush Nidhi Sharma (2011-10-05)
library("semPLS")

clcData <- read.table("./ConstructLevelChanges/clcData.csv", header=TRUE)
sm <- as.matrix(read.table("./ConstructLevelChanges/clcStruc.csv", header=TRUE))
mm <- as.matrix(read.table("./ConstructLevelChanges/clcMeas.csv", header=TRUE))
CLC <- plsm(data=clcData, strucmod=sm, measuremod=mm)
clc <- sempls(model=CLC, data=clcData, maxit=100, wscheme="pathWeighting")

set.seed(1390544442) # Resulted in an error in bootstrap sample b=7
clcBoot <- bootsempls(clc, nboot=7, method="ConstructLevelChanges", verbose=FALSE)

