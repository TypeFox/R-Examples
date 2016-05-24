import.shape <-
function(filename, path=NULL, species.col, ID.col, area.col, dispersal)
  {
    pathfile <- paste(path, filename,sep="")
    sf <- readShapePoly(pathfile)
    df1 <- sf@data
    ctr <- gCentroid(sf, byid=TRUE)
    ctr<-as.data.frame(ctr)
    ID <- df1[, ID.col]
    area <- df1[, area.col]
    if(is.character(species.col)) species <- df1[, species.col]
    if(is.character(species.col)) df3 <- cbind(ID, ctr, area, species)
    else df3 <- cbind(ID,ctr,area)
    mapsize <- max(c(max(ctr$x)-min(ctr$x),max(ctr$y)-min(ctr$y)))
    return(convert.graph(dframe=df3,mapsize,dispersal))
  }
