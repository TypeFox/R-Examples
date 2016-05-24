fRanXY <-
function(x,iextract)
{
  require(maptools)
  require(spatstat)

if (class(iextract)  != "SpatialPolygonsDataFrame" )
  if(class(iextract)  != "SpatialGridDataFrame") 
    stop(paste(sQuote("iextract"), "must be either a", sQuote("SpatialPolygonsDataFrame"),"or a", 
               sQuote("SpatialGridDataFrame"), "object", sep=" "))
   
  #Extract random points from all ids
  if(class(iextract) == "SpatialPolygonsDataFrame") #i.e. a 95% KUD
  {
    ranXY <- do.call(rbind,lapply(x,
                                  function(y)
                                  {
                                    ##select polygon only from 1 animal and randomly sample 1 times
                                    HRpoly <- slot(iextract, "polygons")[[y]]
                                    r.points <- spsample(HRpoly,n=10,type="random")[sample(1:10, 1)]
                                    return(data.frame(r.points,ID=iextract$id[y]))
                                  }
                                  ))
  }
  if(class(iextract) == "SpatialGridDataFrame")  #i.e. a utilisation distribution
  {
    ranXY <- do.call(rbind,lapply(x,
                                  function(y)
                                  {
                                    udsgdf1 <- iextract[y]
                                    udrange <- iextract[[y]]
                                    udrange[which(udrange==0)] <- NA            # Convert 0s to NAs
                                    udsgdf1$udrange <- as.vector(udrange)
                                    r.im <- as.im(udsgdf1[2])                   # convert to im object
                                    r.points <- rpoint(n=1,f=r.im)              # Sample points according to constant
                                    rpointdata <- data.frame(as.data.frame(r.points),ID=names(iextract)[y])
                                    return(rpointdata)
                                  }
                                  ))
  }
  return(ranXY)
}
