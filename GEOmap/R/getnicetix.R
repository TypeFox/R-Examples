getnicetix<-function(lats, lons)
{
  GLON = niceLLtix(lons)

  GLAT = niceLLtix(lats)

  GLON$lab = paste(format(GLON$deg), format(GLON$min), format(GLON$sec)  , sep=":")

  GLAT$lab = paste(format(GLAT$deg), format(GLAT$min), format(GLAT$sec)  , sep=":")

  return(list( LAT=GLAT, LON=GLON ) ) 
}


