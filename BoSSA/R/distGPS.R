`distGPS` <-
function(input)
{
  dMat <- matrix(0,ncol=length(input$lon),nrow=length(input$lon))
  colnames(dMat) <- input$nom
  rownames(dMat) <- input$nom
  for(i in 1:(length(input$lon)-1))
  {
    for(j in (i+1):length(input$lon))
    {
      a <- geoDist(input$lat[i],input$lon[i],input$lat[j],input$lon[j])
      dMat[i,j] <- a
      dMat[j,i] <- a
    }
  }
dMat
}

