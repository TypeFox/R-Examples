# TODO: Add comment
# 
# Author: ahrnee-adm
###############################################################################

### INIT
if(!grepl("SafeQuant\\.Rcheck",getwd())){
	setwd(dirname(sys.frame(1)$ofile))
}
source("initTestSession.R")
### INIT END
### test functions

testGetImpuritiesMatrix <- function(){
	cat(" --- testGetImpuritiesMatrix --- \n")
	# 6-plex
	# old stopifnot(0.094 ==  getImpuritiesMatrix(6)[1,2])
	stopifnot(all.equal(0.004 , round(getImpuritiesMatrix(6)[1,2],3)))
	# 10-plex
	stopifnot(all.equal(0.004 , round(getImpuritiesMatrix(10)[1,2],3)))
	cat(" --- testGetImpuritiesMatrix: PASS ALL TEST --- \n")
	
	
	#getImpuritiesMatrix(test=T)
	# getImpuritiesMatrix(test=T)
	#       [,1]  [,2]  [,3]  [,4]  [,5]  [,6]
	# [1,] 0.939 0.061 0.000 0.000 0.000 0.000
	# [2,] 0.005 0.928 0.067 0.000 0.000 0.000
	# [3,] 0.000 0.011 0.947 0.042 0.000 0.000
	# [4,] 0.000 0.000 0.017 0.942 0.041 0.000
	# [5,] 0.000 0.000 0.000 0.016 0.963 0.021
	# [6,] 0.000 0.000 0.000 0.002 0.032 0.938
	# > cat("Synch1432219972384675000\n");

	
}

testPurityCorrectTMT <- function(){

	cat(" --- testPurityCorrectTMT --- \n")
	# 6-plex
	# old stopifnot( round(9.998839,4)  ==  round(purityCorrectTMT(tmtTestData6Plex,impurityMatrix=getImpuritiesMatrix(6))[1,1],4))
	stopifnot( all.equal(9.4965  , round(purityCorrectTMT(tmtTestData6Plex,impurityMatrix=getImpuritiesMatrix(6))[2,1],4)))
	
	# 10-plex
	stopifnot(all.equal(10.4493 ,  round(purityCorrectTMT(tmtTestData10Plex,impurityMatrix=getImpuritiesMatrix(10))[1,1],4)))
	cat(" --- testPurityCorrectTMT: PASS ALL TEST --- \n")
	
}

testCreateExpDesign <- function(){
	
	cat(" --- testCreateExpDesign --- \n")

	stopifnot(sum(createExpDesign("1,2,3:4,5,6",6)$isControl == c(T,T,T,F,F,F)) == 6 )
	stopifnot(sum(createExpDesign("10,2:3:4,5:6,7,8:9,1",10)$isControl == c(T,T,F,F,F,F,F,F,F,F)) == 10 )
	stopifnot(sum(createExpDesign("10,2:3:4,5:6,7,8:9,1",10)$condition[1:2] == c("Ctrl","Ctrl")) == 2 )
	stopifnot(length(unique(createExpDesign("10,2:3:4,5:6,7,8:9,1",10)$condition)) == 5 )
	
	#createExpDesign("1,4,10:2,5,8:3,6,9",10)
	cat(" --- testCreateExpDesign: PASS ALL TEST --- \n")
	
}

### compare our impurity correction to MSnbase
#comparePurityCorrectionToMsnbase <- function(){
#	
#	cat(" --- comparePurityCorrectionToMsnbase --- \n")
#	
#	library(MSnbase)
#	
#	impurities <- matrix(c(0.929, 0.059, 0.002, 0.000,
#					0.020, 0.923, 0.056, 0.001,
#					0.000, 0.030, 0.924, 0.045,
#					0.000, 0.001, 0.040, 0.923),
#			nrow = 4)
#	
#	qnt <- quantify(itraqdata[1:3,],reporters = iTRAQ4)
#	
#	exprs(qnt)[1,] <- rep(1,4)
#	exprs(qnt)[2,] <- 1:4
#	
#	exprs(qnt)[3,] <- rnorm(4,10)
#	
#	qnt.crct <- purityCorrect(qnt, impurities)
#	
#	testA <- solve(impurities,rep(1,4))
#	### our correction method
#	#testA <- as.vector(t(solve(impurities) %*% rep(1,4)))
#	#testB <- as.vector(t(solve(impurities) %*% 1:4))
#	testB <- solve(impurities,1:4)
#	#testC <- as.vector(t(solve(impurities) %*% exprs(qnt)[3,]))
#	testC <- solve(impurities,exprs(qnt)[3,])
#	
#	stopifnot( sum(round(testA,3) == round(exprs(qnt.crct)[1,],3)) == 4 )
#	stopifnot( sum(round(testB,3) == round(exprs(qnt.crct)[2,],3)) == 4 )
#	stopifnot( sum(round(testC,3) == round(exprs(qnt.crct)[3,],3)) == 4 )
#	
#	cat(" --- comparePurityCorrectionToMsnbase: PASS ALL TEST --- \n")
#	
#
#} 

### test functions end

# INIT

### CREATE TEST DATA

tmtTestData6Plex <- matrix(rep(10,24),ncol=6)
tmtTestData6Plex[2,1:3] <- c(9,9,9) 
tmtTestData6Plex[3,1:3] <- c(100,100,100) 
tmtTestData6Plex[4,c(1,3,5)] <- c(100,100,100) 

tmtTestData10Plex <- matrix(rep(10,100),ncol=10)

### CREATE TEST DATA END


### INIT END

### TESTS

testGetImpuritiesMatrix()
testPurityCorrectTMT()
testCreateExpDesign()
#comparePurityCorrectionToMsnbase()
### TESTS END


