### R code from vignette source 'FeaLect_feature_scorer.Rnw'

###################################################
### code chunk number 1: FeaLect_example
###################################################
	library(FeaLect)
	data(mcl_sll)
	F <- as.matrix(mcl_sll[ ,-1])	# The Feature matrix
	L <- as.numeric(mcl_sll[ ,1])	# The labels
	names(L) <- rownames(F)
	message(dim(F)[1], " samples and ",dim(F)[2], " features.")

	FeaLect.result <-FeaLect(F=F,L=L,maximum.features.num=10,total.num.of.models=100,talk=TRUE)	


###################################################
### code chunk number 2: FeaLect_example
###################################################
	plot(FeaLect.result$log.scores, pch=19)


