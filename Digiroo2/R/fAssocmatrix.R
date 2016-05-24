fAssocmatrix <-
function(sPerm,Gprox,iextract,iID)
{
  require(spdep)
  require(spatstat)
  
  if (class(iextract)  != "SpatialPolygonsDataFrame" )
    if(class(iextract)  != "SpatialGridDataFrame") 
      stop(paste(sQuote("iextract"), "must be either a", sQuote("SpatialPolygonsDataFrame"),"or a", 
                 sQuote("SpatialGridDataFrame"), "object", sep=" "))
  
  if (is.vector(iID) == FALSE)
    if(is.matrix(iID) == FALSE) 
      stop(paste(sQuote("iID"), " must be either a vector or a matrix"))
  
  fAssocmatrix1 <- function(iPerm,Gprox,iextract,iID)
  {
    if(is.vector(iID) == TRUE)
      ranXY <- fRanXY(iID,iextract)
    if(is.matrix(iID) == TRUE)
      ranXY <- fRanXY(iID[,iPerm],iextract)
    coordinates(ranXY) <- ~x+y
    dnn_digi <- dnearneigh(ranXY,0,Gprox,row.names=ranXY$ID) #Run NN function
    (tassoc <- data.frame(Permutation=iPerm,fAssoctable(dnn_digi))) #Run Association matrix convertor
    return(tassoc)
  }  
  lassoc <- lapply(sPerm,function(i) fAssocmatrix1(iPerm=i,
                                                   Gprox=Gprox,
                                                   iextract=iextract,
                                                   iID=iID))
  return(do.call(rbind,lassoc))
}
