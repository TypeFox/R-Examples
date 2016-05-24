`itemDescriptives` <-
	function(theta,rawData){
    itemMean <- colMeans(rawData, na.rm=TRUE)
    pBis      <- apply(rawData,2,function(XXX) cor(XXX,theta,use="pair"))
    bis       <- apply(rawData,2,function(XXX) polyserial(theta,XXX))	
	cbind(itemMean = itemMean, pBis = pBis, bis=bis)
}

	