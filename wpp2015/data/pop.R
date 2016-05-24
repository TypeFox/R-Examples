# This dataset is created on the fly as a sum of the age-specific population estimates popM and popF

pop <- local({
	library(plyr)
	source('popM.R')
	source('popF.R')

	sum.by.country <- function(dataset) {
		year.cols.idx <- grep('^[0-9]{4}', colnames(dataset))
		ddply(dataset[,c(which(colnames(dataset)=='country_code'), year.cols.idx)], "country_code", .fun=colwise(sum))
	}
	tpopM <- sum.by.country(popM)
	tpopF <- sum.by.country(popF)
	# The male and female dataset should be in the same format, 
	# i.e. the countries and years should be in the same order, but just to be sure
	# match columns and rows. It will fail if there are different sets of countries
	# in the two datasets.
	cols.to.sumM <- colnames(tpopM)[-match('country_code', colnames(tpopM))]
	cols.to.sumF <- colnames(tpopF)[-match('country_code', colnames(tpopF))]
	cols.to.sumF.idx <- match(cols.to.sumF, cols.to.sumM)
	rowsF.idx <- match(tpopF$country_code, tpopM$country_code)
	name.col <- grep('^name$|^country$', colnames(popM), value=TRUE)
	cbind(country_code=tpopM$country_code, name=popM[,name.col][match(tpopM$country_code, popM$country_code)],
			 tpopM[,cols.to.sumM] + tpopF[rowsF.idx, cols.to.sumF[cols.to.sumF.idx]])
})