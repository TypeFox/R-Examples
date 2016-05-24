fDrawfigure <-
function(x,iextract,idataxy=NULL,istudy=NULL)
{
  require(spatstat)

  if (class(iextract)  != "SpatialPolygonsDataFrame" )
    if(class(iextract)  != "SpatialGridDataFrame") 
      stop(paste(sQuote("iextract"), "must be either a", sQuote("SpatialPolygonsDataFrame"),"or a", 
                 sQuote("SpatialGridDataFrame"), "object", sep=" "))

  if (is.null(idataxy) != TRUE)
    if(class(idataxy)  != "SpatialPointsDataFrame") 
      stop(paste(sQuote("idataxy"), "must be either NULL or a", 
                 sQuote("SpatialPointsDataFrame"), "object", sep=" "))
  
  #Extract random points from all ids
  if(class(iextract) == "SpatialPolygonsDataFrame") #i.e. a 95% KUD
  {  
    ##select polygon only from 1 Roo
    HRpoly <- slot(iextract, "polygons")[[x]]
    HRid <- sapply(slot(iextract, "polygons"), function(i) slot(i, "ID"))
    HRpolySP <- SpatialPolygons(list(HRpoly))
    #Randomly sample 1 points within a polygon
    r.points <- spsample(HRpoly,n=10,type="random")[sample(1:10, 1)]
    ## Now plot the data
    if(is.null(istudy) == TRUE)
    {
      plot(iextract,col=NULL,border="white")
    }
    if(class(istudy)[1]=="RasterLayer")
    {
      image(istudy, col=colorRampPalette(c("black", "white"))(255))
    }
    if(class(istudy)[1]=="SpatialPolygons")
    {
      plot(istudy)
    }
    if(class(istudy)[1]=="SpatialPolygonsDataFrame")
    {
      plot(istudy)
    }
    if(is.null(idataxy) == FALSE) ###select actual locations from only these individuals and plot
    {
      xys <- idataxy[idataxy[[1]] %in% HRid[x],]
      plot(xys,col="green",add=TRUE,pch=16,cex=0.4)
    }
    plot(r.points,col="cyan",add=TRUE,pch=10,cex=1.5)
    plot(HRpolySP,border="red",lwd=1,lty=2,add=TRUE)                          
  }
  
  if(class(iextract) == "SpatialGridDataFrame")  #i.e. a utilisation distribution then draw figure of single individual xy, HR and probability weighted random points within ud
  { 
    #Randomly sample 1 points within a UD
    r.points <- fRanXY(x,iextract)[1:2]
    ## Now plot the data
    if(is.null(istudy) == TRUE)
    {
      image(iextract[x])
    }
    if(class(istudy)[1]=="RasterLayer")
    {
      image(istudy, col=colorRampPalette(c("black", "white"))(255))
    }
    if(class(istudy)[1]=="SpatialPolygons")
    {
      image(iextract[x])
      plot(istudy,add=T)
    }
    if(class(istudy)[1]=="SpatialPolygonsDataFrame")
    {
      image(iextract[x])
      plot(istudy,add=T)
    }  
    if(is.null(idataxy) != TRUE)
    {
      xys <- idataxy[idataxy[[1]] %in% names(iextract)[x],]
      plot(xys,col="green",add=TRUE,pch=16,cex=0.4)
    }
    
    points(r.points,col="cyan",pch=10,cex=1.5)
  }
}
