

data( fujitopo, package='geomapdata' )

data(japmap, package='geomapdata')
###data(ETOPO5)


######  set a location for the origin:
PLOC=list(LON=range(fujitopo$lon), x=range(fujitopo$lon),
  LAT=range(fujitopo$lat), y=range(fujitopo$lat))


####  set the projection
PROJ = GEOmap::setPROJ(type=2, LAT0=mean(PLOC$y) , LON0=mean(PLOC$x) )

########  all of Japan
plotGEOmapXY(japmap, PROJ=PROJ, add=FALSE)


readline("Hit Return>")

###  limited around Mt. Fuji

plotGEOmapXY(japmap, PROJ=PROJ, LIM=c(PLOC$LON[1], PLOC$LAT[1],PLOC$LON[2], PLOC$LAT[2] ) , add=FALSE)

readline("Hit Return>")

#########   set the color scheme for the perspective plot
calcol = settopocol()

###############  get a subset of the ETOPO5 data

data(ETOPO5, package='geomapdata')


### ZZ2 = subsetTOPO(ETOPO5, PLOC, PROJ)


ZZ = ETOPO::getetopo(ETOPO5, PLOC$LAT, PLOC$LON )

####  dim(ZZ2$z)
ZZ2 =list(x=attr(ZZ2, 'lon'), y=attr(ZZ, 'lat'), z=ZZ)

 image(ZZ2$x,ZZ2$y, ZZ2$z[ , rev(1:(dim(ZZ2$z)[2]))],  col=topo.colors(100), asp=TRUE , axes=FALSE, xlab="", ylab="" )


readline("Hit Return>")

####  create X-Y matrix for interpolation
G = setplotmat(ZZ2$x,ZZ2$y)



XY = GLOB.XY(rep(ZZ2$y[1], length(ZZ2$x)), ZZ2$x, PROJ)
XYs = GLOB.XY(ZZ2$y,rep(ZZ2$x[1], length(ZZ2$y)),  PROJ)

############  plot image, projected
image(XY$x, XYs$y, ZZ2$z[ , rev(1:(dim(ZZ2$z)[2]))],  col=topo.colors(100), asp=TRUE , axes=FALSE, xlab="", ylab="" )

########   add map data:
plotGEOmapXY(japmap, PROJ=PROJ, LIM=c(PLOC$LON[1], PLOC$LAT[1],PLOC$LON[2], PLOC$LAT[2] ) , add=TRUE)

readline("Hit Return>")
############   make a nice perspective plot of Mt. Fuji region....


##   PMAT = GEOpersp(japmap, jtop  ,calcol=calcol$calcol)

 CCOL = settopocol()
        calcol = CCOL$calcol
#####  this is an interpolation function from ETOPO package
interpETOPO <-function (b5, PROJ, nx = 500, ny = 500, nb = 4, mb = 4, hb = 8) 
{
    xlon = attr(b5, "lon")
    ylat = attr(b5, "lat")
    LLgrid = RPMG::meshgrid(xlon, ylat)
    GXY = GEOmap::GLOB.XY(LLgrid$y, LLgrid$x, PROJ)
    DF = cbind(as.vector(GXY$x), as.vector(GXY$y), as.vector(t(b5)))
    IZ = MBA::mba.surf(DF, nx, ny, n = nb, m = mb, h = hb, extend = TRUE)$xyz.est
    return(IZ)
}

IZ = interpETOPO(ZZ, PROJ, nx = 500, ny = 500, nb = 4, mb = 4, hb = 8)


Cmat = TOPOCOL(IZ$z, calcol)
Dcol = attr(Cmat, "Dcol")

PMAT = persp(IZ$x, IZ$y, IZ$z, theta = 0, phi = 90, 
            r = 4000, col = Cmat[1:(Dcol[1] - 1), 1:(Dcol[2] - 
                1)], scale = FALSE, ltheta = 120, lphi = 30, 
            shade = 0.75, border = NA, expand = 0.001, box = FALSE)
   

plotGEOmapXY(japmap, PROJ=PROJ, LIM=c(PLOC$LON[1], PLOC$LAT[1],PLOC$LON[2],
                                    PLOC$LAT[2] ) , PMAT=PMAT , add=TRUE)

