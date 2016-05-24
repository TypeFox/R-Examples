`SHOWTOMO` <-
function(MOD, colmap=topo.colors(100),  zlim=NULL,  MAP=NULL, I=1, J=2, bkgr="white" ,
         linelty=1, linelwd=1, ptpch=".", ptcex=1 )
{
  if(missing(colmap)) { colmap=tomo.colors(100) }
  if(missing(MAP)) { MAP=NULL }
  if(missing(I)) { I=1 }
  if(missing(J)) { J=length(MOD$MOD) }
  if(missing(zlim)) { zlim=NULL }
    if(missing(bkgr)) { bkgr="white" }
    

  
  opar = par(no.readonly = TRUE)
  
  NTOT = (J-I+1)
  KROW =  floor(NTOT/5)
  
  if(KROW>=1)
    {
    par(mfrow=c(KROW,5),mai=c(0, 0, 0,0))
  }
  else
    {
      par(mfrow=c(1, NTOT),mai=c(0, 0, 0,0))

    }

  KMAX = min(c(J, 25))
  for( i in I:KMAX)
    {
      
      pltomo(MOD$x,MOD$y,MOD$MOD,i, colmap=colmap, zlim=zlim,  bkgr= bkgr )
      
      ##  image(xo, yo, TM1$MOD[[i]], col=tomo.colors(100), axes=FALSE, ann=FALSE) 

      
      text(MOD$x[1],MOD$y[1], labels=paste(sep=" : ", paste(sep="", "LAY=", i) , paste(sep="", "Z=",MOD$D[i])), adj=c(0,0) ,xpd=TRUE, font=2)
      
###   HOZscale(MOD$MOD[[i]] , col= rainbow(100)  , units="1/q", SIDE=2)
###   image.SCALE(MOD$MOD[[i]] , col=rainbow(100), nlab=2)
###  PROJmap(JAPmap,  ADD=TRUE, COL=TRUE)

      box()
        if(!is.null(MAP)) GEOmap::plotGEOmapXY(MAP, PROJ=MAP$PROJ, add=TRUE, xpd=FALSE,   linelty=linelty, linelwd=linelwd, ptpch=ptpch, ptcex=ptcex  )
      
      
###     title(paste(sep=' ', "Depth", MOD$D[i], "-", MOD$D[i+1]))
      
    }

  newpar= par(no.readonly = TRUE)
  par(opar)

  
  invisible(newpar)
}

