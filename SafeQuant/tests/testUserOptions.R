# TODO: Add comment
# 
# Author: ahrnee-adm
###############################################################################


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


### TEST FUNCTIONS

testExpDesignTagToExpDesign <- function(){
	
	cat("--- testExpDesignTagToExpDesign: --- \n")
	
	stopifnot(all.equal(pData(eset),expDesign))
	stopifnot(all.equal(sampleNames(eset),colnames(m)))
	stopifnot(all.equal(nrow(exprs(eset)),nrow(m)))
	
	expDesignString1 <- "1,2,3:4,5,6"
	
	# 6-plex default: 1,2,3:4,5,6 
	#condition isControl
	#1 Condition 1      TRUE
	#2 Condition 1      TRUE
	#3 Condition 1     TRUE
	#4 Condition 2     FALSE
	#5 Condition 2     FALSE
	#6 Condition 2     FALSE
	
	
	expDesign <- data.frame(condition=paste("Condition",sort(rep(c(1,2),3))),isControl=sort(rep(c(T,F),3),decreasing=T) )
	
	expDesign1 <- expDesignTagToExpDesign(expDesignString1, expDesign)
	
	stopifnot(nrow(expDesign1) == 6 )
	stopifnot(length(unique(expDesign1$condition)) == 2 )
	stopifnot(sum(expDesign1$isControl) == 3 )
	
	expDesignString2 <- "1,4,7,10:2,5,8:3,6,9"
	
	# 10-plex default is "1,4,7,10:2,5,8:3,6,9"
	#condition isControl
	#1 Condition 1      TRUE
	#2 Condition 2      FALSE
	#3 Condition 3     FALSE
	#4 Condition 1     TRUE
	#5 Condition 2     FALSE
	#6 Condition 3     FALSE
	#7 Condition 1     TRUE
	#8 Condition 2     FALSE
	#9 Condition 3     FALSE
	#10 Condition 1     TRUE
	
	expDesign <- data.frame(condition=paste("Condition",c(1,2,3,1,2,3,1,2,3,1)),isControl=c(T,F,F,T,F,F,T,F,F,T) )
	expDesign2 <- expDesignTagToExpDesign(expDesignString2, expDesign)
	stopifnot(nrow(expDesign2) == 10 )
	stopifnot(length(unique(expDesign2$condition)) == 3 )
	stopifnot(sum(expDesign2$isControl) == 4 )
	
	### condition name assignment when mixing runs from different conditions
	expDesign <- data.frame(condition=paste("foo",c(1,1,1,2,2,3,3)),isControl=c(F,F,F,T,T,F,F) )
	stopifnot(all(grepl("foo" ,expDesignTagToExpDesign("1,2,3:4,5:6",expDesign)$condition)))
	stopifnot(all(grepl("Condition" ,expDesignTagToExpDesign("1:4,6:5",expDesign)$condition)))
	stopifnot(all(grepl("foo" ,expDesignTagToExpDesign("1,2,3:4,5",expDesign)$condition)))
	stopifnot(all(grepl("Condition" ,expDesignTagToExpDesign("1,2,4",expDesign)$condition)))
	stopifnot(all(grepl("foo" ,expDesignTagToExpDesign("2",expDesign)$condition)))
	stopifnot( length(unique(expDesignTagToExpDesign("1:2:3:4:5:6",expDesign)$condition)) == 6)
		
	cat("--- testExpDesignTagToExpDesign: PASS ALL TEST --- \n")
	
}


### TEST FUNCTIONS END

### TESTS
testExpDesignTagToExpDesign()













#names(expDesign) <- 1:ncol(expDesign)

