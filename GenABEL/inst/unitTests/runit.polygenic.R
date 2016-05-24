### --- Test setup ---
#
# regression tests
#

if(FALSE) {
	## Not really needed, but can be handy when writing tests
	library(RUnit)
	library(GenABEL)
}

### do not run
#stop("SKIP THIS TEST")
###

### ---- common functions and data -----

#source(paste("../inst/unitTests/shared_functions.R"))
#source(paste(path,"/shared_functions.R",sep=""))

### --- Test functions ---

test.polygenic.Bug1322 <- function()
{
  require(GenABEL.data)
	data(ge03d2.clean)
	df <- ge03d2.clean[1:200,autosomal(ge03d2.clean)]
	gkin <- ibs(df,w="freq")
# get a SNP which is informative + no missings
	s <- summary(gtdata(df))
	maf <- pmin(s[,"Q.2"],1.-s[,"Q.2"])
	misSnp <- (rownames(s)[which(maf>0.2 & s[,"CallRate"]<1)])[1]
	snpX <- as.numeric(gtdata(df[, misSnp]))
	print(table(snpX,exclude=NULL))
	phdata(df)$snpX <- snpX
	pol1 <- polygenic(height~sex+age+snpX,df,kin=gkin,quiet=TRUE)
	print(pol1$esth2)
	nUsedIds1 <- sum(pol1$measuredIDs)
# introduce missing indicator
	snpInd <- rep(1,length(snpX))
	snpInd[is.na(snpX)] <- NA
	checkEquals(nUsedIds1,sum(!is.na(snpInd)))
# this does did not work!
	phdata(df)$snpInd <- snpInd
	pol0 <- polygenic(height~sex+age+snpInd,df,kin=gkin,quiet=TRUE)
# correct way to do comparison
	pol0 <- polygenic(height~sex+age,df[!is.na(snpX),],kin=gkin[!is.na(snpX),!is.na(snpX)],quiet=TRUE)
	nUsedIds0 <- sum(pol0$measuredIDs)
	checkEquals(nUsedIds1,nUsedIds0)
}

test.polygenic.eigenOfRel <- function()
{
  require(GenABEL.data)
	data(ge03d2.clean)
	completeIds <- complete.cases(phdata(ge03d2.clean)[,c("height","sex","age")])
	df <- ge03d2.clean[sample(which(completeIds),250),autosomal(ge03d2.clean)]
	gkin <- ibs(df,w="freq")
	pol1 <- polygenic(height~sex+age,df,kin=gkin,quiet=TRUE)
	gRel <- gkin
	gRel[upper.tri(gRel)] <- t(gkin)[upper.tri(gkin)]
	gRel <- gRel*2
	eigRes <- eigen(gRel)
	pol2 <- polygenic(height~sex+age,df,kin=gkin,quiet=TRUE,eigenOfRel=eigRes)
	checkEquals(pol1[which(names(pol1)!="call")],pol2[which(names(pol2)!="call")])
}
