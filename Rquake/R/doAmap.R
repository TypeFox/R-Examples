doAmap <-
function(stas, doproj=TRUE)
  {
    proj = GEOmap::setPROJ(type=2, LAT0 =median(stas$lat) , LON0 = median(stas$lon) )
###  generic plotting routine for gMAP

    
    if(doproj)
      {
        XY = GEOmap::GLOB.XY(stas$lat, stas$lon, proj)

        BEX = GEOmap::expandbound(range(XY$x), 0.1)
        BEY  = GEOmap::expandbound(range(XY$y), 0.1)
        plot(BEX, BEY, type='n', xlab="km", ylab="km"  )
        points(XY, pch=25, bg=stas$col, cex=1.2)
        text(XY, labels=stas$name, pos=3, xpd=TRUE, cex=1.1)
      }
    else
      {
        plot(stas$lon, stas$lat, pch=25, bg=stas$col, cex=1.2)
        text(stas$lon, stas$lat, labels=stas$name, pos=3, xpd=TRUE, cex=1.1)
      }
    return(proj)
  }
