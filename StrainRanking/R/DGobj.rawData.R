

### Constructor ###
DGobj.rawdata <- function(demographic.coord, demographic.measures, 
	genetic.coord, genetic.frequencies){
	## Construction of the genetic and demographic data sets
	geneticData=cbind(genetic.coord, genetic.frequencies)
	Z=log(1+as.numeric(demographic.measures[,2]))-log(1+as.numeric(demographic.measures[,1]))
	demographicData=cbind(demographic.coord,Z)[!is.na(Z),]
	nbStrains=ncol(genetic.frequencies)
	colnames(geneticData)=c("x1","x2",paste("strain",1:nbStrains,sep=""))
	colnames(demographicData)=c("x1","x2","growth")
	new("DGobj",demographic=as.matrix(demographicData),genetic=as.matrix(geneticData))
}
