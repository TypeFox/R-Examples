#do not edit, edit noweb/qmrparser.nw
pcAxisCubeToCSV <- function(prefix,pcAxisCube) {

  write.csv(pcAxisCube$pxCube              , file = paste(prefix,"pxCube.csv"              ,sep=""),row.names = FALSE)
  
  write.csv(pcAxisCube$pxCubeVariable      , file = paste(prefix,"pxCubeVariable.csv"      ,sep=""),row.names = FALSE)
  
  write.csv(pcAxisCube$pxCubeVariableDomain, file = paste(prefix,"pxCubeVariableDomain.csv",sep=""),row.names = FALSE)
  
  write.csv(pcAxisCube$pxCubeData          , file = paste(prefix,"pxCubeData.csv          ",sep=""),row.names = FALSE)

  for( name in names(pcAxisCube$pxCubeAttrN) )
    write.csv(pcAxisCube$pxCubeAttrN[[name]], file = paste(prefix,"pxCube",name,".csv ",sep=""),row.names = FALSE)      
}
