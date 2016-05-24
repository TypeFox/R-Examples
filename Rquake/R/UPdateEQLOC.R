UPdateEQLOC <-
function(PF,  sol, vel, stas=NULL)
  {
######   update a hypocenter solution with a newmodel
    v = vel
    pstas  = PF$STAS

    if(!is.null(stas))
      {
        MASTA  = match(pstas$name, stas$name)
        
        stalats =  stas$lat[MASTA]
        stalons = stas$lon[MASTA]
        stanam  = stas$name[MASTA]
        staz = stas$z[MASTA]
        pstas$lat = stalats
        pstas$lon = stalons
        pstas$z = staz
      }

    proj = GEOmap::setPROJ(type=2, LAT0 =median(pstas$lat) , LON0 = median(pstas$lon) )
    
    XY = GEOmap::GLOB.XY(pstas$lat, pstas$lon, proj)
    elcor = rep(0, length(pstas$lat))
    
    DZ = pstas$z - mean(stas$z)
    
    elcor[pstas$phase=="P"] = DZ[pstas$phase=="P"]/vel$vp[1]
    elcor[pstas$phase=="S"] = DZ[pstas$phase=="S"]/vel$vs[1]
    
    XY$cor = elcor
    XY$phase = pstas$phase
    XY$sec = pstas$sec
    
    eqXY = GEOmap::GLOB.XY(sol[1], sol[2], proj)
    res =  EQXYresid(XY, vel=vel , h1=c(eqXY$x, eqXY$y, sol[2],sol[4]) , PLOT=FALSE)
    
    pstas$res = res
    
    PF$LOC$lat = sol[1]
    PF$LOC$lon = sol[2]
    PF$LOC$z = sol[3]
    PF$LOC$sec = PF$LOC$sec + sol[4]
    
    PF$H = PF$LOC
    
    PF$STAS = pstas
    return(PF)
    
######   eqsol = NLSlocate(GH, v=GH$vel,  PLOT=TRUE )
  }
