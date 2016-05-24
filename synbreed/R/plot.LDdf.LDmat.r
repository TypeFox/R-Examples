plot.LDdf <- function(x,gpData,plotType="dist",dense=FALSE,nMarker=TRUE,centr=NULL,
                      chr=NULL,type="p",breaks=NULL,n=NULL,file=NULL,fileFormat="pdf",onefile=TRUE,
                      colL=2,colD=1,...){
    if(plotType == "neighbour"){
      plotNeighbourLD(LD=x,gpData=gpData,dense=dense,nMarker=nMarker,centr=centr, file=file, fileFormat=fileFormat,...) 
    } else if(plotType == "dist"){ 
      LDDist(LDdf=x,chr=chr,type=type,breaks=breaks,n=n,file=file,fileFormat=fileFormat,onefile=onefile,colL=colL,colD=colD,...)
    } else stop("plotType has to be whether 'neighbour' or 'dist'!")
}

plot.LDmat <- function(x,gpData,plotType="map",dense=FALSE,nMarker=TRUE,centr=NULL,
                       chr=NULL,file=NULL,fileFormat="pdf",onefile=TRUE,...){
    if(plotType == "neighbour"){
      plotNeighbourLD(LD=x,gpData=gpData,dense=dense,nMarker=nMarker,centr=centr,file=file,fileFormat=fileFormat,...)
    } else if(plotType == "map"){
      LDMap(LDmat=x,gpData=gpData,chr=chr,file=file,fileFormat=fileFormat,onefile=onefile,...)
    } else stop("plotType has to be whether 'neighbour' or 'map'!")
}
