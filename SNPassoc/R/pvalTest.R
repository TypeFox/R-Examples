`pvalTest` <-
function(dataX,Y,quantitative,type,genotypingRate)
 {
  pvalues<-t(data.frame(lapply(dataX,FUN=modelTest,Y=Y,
               quantitative=quantitative,type=type,
               genotypingRate = genotypingRate )))
  pvalues
 }

