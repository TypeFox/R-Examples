#########################################################################
# File Name: buildGonnet.R
# Author: Ge Tan
# mail: gtan@me.com
# Created Time: Tue Jul 29 21:55:22 2014
#########################################################################
## This script is used to create the mutation matrix of GONNET at PAM 1.

transitionPAM1 <- read.table("daymatrix.txt", header=FALSE)
transitionPAM1 <- as.matrix(transitionPAM1)
AAOrder = c("A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V")
GONNETBF <- transitionPAM1[1,]
names(GONNETBF) <- AAOrder
GONNET <- transitionPAM1[2:21, ]
colnames(GONNET) <- AAOrder
rownames(GONNET) <- AAOrder

save(GONNET, GONNETBF, file="GONNET.rda", compress="xz")


