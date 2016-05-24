fAssocplot <-
function(dnndata,idataxy,iplotnames=FALSE)
{ 
  require(spdep)
  
  if(class(dnndata)  != "nb") 
    stop(paste(sQuote("dnndata"), "must be a", 
               sQuote("nb"), "object", sep=" "))
  
  if(class(idataxy)  != "SpatialPointsDataFrame") 
    stop(paste(sQuote("idataxy"), "must be a", 
               sQuote("SpatialPointsDataFrame"), "object", sep=" "))
  
  ## Display all random points and plot the associations
  plot(dnndata,coordinates(idataxy),col=2,lwd=2,
       xlim=c(min(coordinates(idataxy)[,1]-50),max(coordinates(idataxy)[,1]+50)),
       ylim=c(min(coordinates(idataxy)[,2]-50),max(coordinates(idataxy)[,2]+50))) 
  if(iplotnames==TRUE)
    text(idataxy,attr(dnndata, "region.id"),pos=rep(1,length(attr(dnndata, "region.id"))),cex=0.5,col="red")
  box("plot",col="blue")
}
