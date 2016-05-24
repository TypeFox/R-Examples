`subsetTOPO` <-
    function(TOPO, ALOC, PROJ, nx=500, ny=500, nb = 4, mb = 4, hb = 8)
{
    
#################   extract topographic information from ETOPO5 data base
####  this code works mainly with ETOPO5 and ETOPO2 from package ETOPO or geomapdata
    nn = names(ALOC)

    if(length(nn)<1)
        { return(NULL) }

    wlon = grep("lon", nn, ignore.case = TRUE)
    wlat = grep("lat", nn, ignore.case = TRUE)

    if(length(wlon)<1)   { print('subsetTOPO ERROR: no lon range input'); return(NULL) }
    if(length(wlat)<1)   { print('subsetTOPO ERROR: no lat range input');  return(NULL) }

    glon = ALOC[[wlon[1] ]]
    glat = ALOC[[wlat[1] ]]

    glon = RPMG::fmod(glon, 360)
    
    glon  =   range(glon)
    glat  =   range(glat)

    b5 = GEOmap::getETOPO(TOPO, glat =glat ,   glon =glon   )
   
####  uses the input of ETOPO::getetopo to convert to X-Y coordinates
        xlon=attr(b5, 'lon');
        ylat=attr(b5,'lat');
        ## zel=b5
         ##  zelda = t(b5 )

        LLgrid = RPMG::meshgrid(xlon, ylat)

####  now we have a grid of lat-lon - we will convert to X-Y

     ##    OLON = mean(xlon)
     ##    OLAT = mean(ylat)

        ##  PROJ = setPROJ(2, LAT0 = OLAT, LON0 = OLON)

        GXY = GEOmap::GLOB.XY(LLgrid$y, LLgrid$x, PROJ)


        DF = cbind( as.vector(GXY$x) , as.vector(GXY$y) ,  as.vector(t(b5))  )
        IZ =    MBA::mba.surf(DF, nx, ny, n = nb, m = mb, h = hb, extend=TRUE )$xyz.est

    return(IZ )
    
}

