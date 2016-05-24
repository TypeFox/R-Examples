# ===========================================================================
# File: "runit_concordance.R"
#                        Created: 2012-11-06 20:02:30
#              Last modification: 2012-11-07 14:42:40
# Author: Bernard Desgraupes
# e-mail: <bdesgraupes@users.sourceforge.net>
# Unit test file for the R package clusterCrit.
# ===========================================================================
# Two artificial partitions with 3 and 5 clusters respectively
# set.seed(1)
# part1<-sample(1:3,150,replace=TRUE)
# part2<-sample(1:5,150,replace=TRUE)
# concordance(part1,part2)
# 			 y    n
# 		y  733 3026
# 		n 1504 5912


test.concordance <- function() {
	part1 <- c(1,2,2,3,1,3,3,2,2,1,1,1,3,2,3,2,3,3,2,3,3,1,2,1,1,2,1,2,3,2,2,2,2,1,3,3,3,1,3,2,3,2,3,2,2,3,1,2,3,3,2,3,2,1,1,1,1,2,2,2,3,1,2,1,2,1,2,3,1,3,2,3,2,2,2,3,3,2,3,3,2,3,2,1,3,1,3,1,1,1,1,1,2,3,3,3,2,2,3,2,2,2,1,3,2,1,1,2,3,2,3,3,2,2,1,1,3,1,2,2,3,2,2,1,3,2,2,1,1,2,2,1,1,2,3,2,2,2,3,2,3,2,1,1,3,2,1,3,1,3)
	part2 <- c(4,3,2,3,3,1,3,1,2,2,2,5,3,4,5,3,1,2,4,2,4,5,5,2,2,5,4,4,4,5,2,1,5,3,5,1,4,4,5,3,4,2,1,5,2,3,1,5,2,4,2,2,3,2,1,3,3,1,2,4,5,1,4,5,5,2,4,5,5,2,2,1,2,3,5,3,2,1,3,5,2,1,2,4,2,4,4,3,3,3,2,3,5,1,3,2,3,1,3,5,4,5,3,4,3,1,2,3,2,5,3,2,2,4,4,1,1,4,4,1,1,1,2,1,2,1,2,1,3,4,1,3,5,2,1,1,2,1,1,2,2,1,5,2,3,4,1,1,1,5)	

	expect <- matrix(c(733,1504,3026,5912), nrow=2)
	colnames(expect) <- c("y","n")
	rownames(expect) <- c("y","n")
	result <- concordance(as.integer(part1),as.integer(part2))

	checkEquals(result,expect)
}

