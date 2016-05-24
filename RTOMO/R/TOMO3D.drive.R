`TOMO3D.drive` <-
function(MOD, COL=NULL, LIM=NULL, MAP=NULL, MAPLIM=NULL,  ZLIM=c(0, 30), STA=NULL, TOPO=NULL,  STDLAB=c("DONE", "QUIT") )
  {

    if(missing(COL)) { COL=tomo.colors(100) }
    if(missing(LIM)) { LIM=NULL }
    if(missing(MAP)) { MAP=NULL }
    if(missing(MAPLIM)) { MAPLIM=NULL }
    
    if(missing(STA)) { STA=NULL }
    if(missing(TOPO)) { TOPO=NULL }
    if(missing(ZLIM)) { ZLIM =c(0, 30) }
    
     
   
  if(missing(STDLAB))
    { STDLAB = c("DONE", "QUIT", "FRESH", "CLICKS", "AREA", "ASTATS", "LINE", "XSEC", "TOPO", "UP", "DOWN")}

    stdlab =  STDLAB

    BLABS = c(stdlab)
    NLABS = length(BLABS)
    NOLAB = NLABS +1000

    colabs = rep(1,length(BLABS))
    pchlabs = rep(4,length(BLABS))
    
    ilay = 1
    NLAY = length(MOD$MOD)
    MESHXY = meshgrid(MOD$x, MOD$y)
    PTSXY = cbind(as.vector(MESHXY$x), as.vector(MESHXY$y) )
    
### pltomo(MOD$x,MOD$y,MOD$MOD,ilay, tomocolors, zlim=c(-15, 15))


    proj = MAP$PROJ

    if(is.null(proj))
      {
        
        proj = MOD$proj

      }

    
    
    FANCY.TOMO(MOD, ilay, COL=COL, LIM=LIM, MAP=MAP, MAPLIM=MAPLIM, STA=STA, TIT = NULL)
    upar = par("usr")

    
    buttons = RPMG::rowBUTTONS(BLABS, col=colabs, pch=pchlabs)

    MAINdev = dev.cur()
    dev.set( MAINdev)
    u = par("usr")
    sloc = list(x=c(u[1],u[2]))

     iloc = locator(1, col=rgb(1,0.8, 0.8), type='p')
    zloc = iloc
  
    
    Nclick = length(zloc$x)
    
    if(is.null(iloc$x)) { return(0) }
    K = RPMG::whichbutt(iloc ,buttons)
  
  sloc = zloc

   
    DF = NULL
    MAINdev = dev.cur()

    
  PLOC = NULL

  while(TRUE)
    {
      
      if(K[Nclick] == match("DONE", BLABS, nomatch = NOLAB))
        {
          buttons = RPMG::rowBUTTONS(BLABS, col=rep(grey(.8), length(BLABS)), pch=rep("NULL", length(BLABS)))
          title("Return to Calling Program")
          
          
          break;
        }
      
      if(K[Nclick] == match("QUIT", BLABS, nomatch = NOLAB))
        {
          
          buttons = RPMG::rowBUTTONS(BLABS, col=rep(grey(.8), length(BLABS)), pch=rep("NULL", length(BLABS)))
          title("Return to Calling Program")
          
          return(NULL)
        }



        if(iloc$x<upar[1] & (iloc$y>upar[3] & iloc$y<upar[4]) )
          {
            FANCY.TOMO(MOD, ilay, COL=COL, LIM=LIM, MAP=MAP, MAPLIM=MAPLIM, STA=STA, TIT = NULL )
            upar = par("usr")

            buttons = RPMG::rowBUTTONS(BLABS, col=colabs, pch=pchlabs)
            zloc = list(x=NULL, y=NULL)
     

          }


      
     if(K[Nclick] == match("FRESH", BLABS, nomatch = NOLAB))
        {
          FANCY.TOMO(MOD, ilay, COL=COL, LIM=LIM, MAP=MAP, MAPLIM=MAPLIM, STA=STA, TIT = NULL)
          upar = par("usr")

          buttons = RPMG::rowBUTTONS(BLABS, col=colabs, pch=pchlabs)
          zloc = list(x=NULL, y=NULL)
        }

      
    if(K[Nclick] == match("UP", BLABS, nomatch = NOLAB))
        {
          ilay = ilay-1
          if(ilay<1) ilay = 1
          
          FANCY.TOMO(MOD, ilay, COL=COL, LIM=LIM, MAP=MAP, MAPLIM=MAPLIM, STA=STA)
          upar = par("usr")

          buttons = RPMG::rowBUTTONS(BLABS, col=colabs, pch=pchlabs)
          zloc = list(x=NULL, y=NULL)

        }

    if(K[Nclick] == match("DOWN", BLABS, nomatch = NOLAB))
        {
          ilay = ilay+1
          if(ilay>NLAY) ilay = NLAY
          

          FANCY.TOMO(MOD, ilay, COL=COL, LIM=LIM, MAP=MAP, MAPLIM=MAPLIM, STA=STA)
          upar = par("usr")

          buttons = RPMG::rowBUTTONS(BLABS, col=colabs, pch=pchlabs)
          zloc = list(x=NULL, y=NULL)

        }

    if(K[Nclick] == match("CLICKS", BLABS, nomatch = NOLAB))
        {
          NKLIX = length(zloc$x)-1
          ex = zloc$x[1:NKLIX]
          why = zloc$y[1:NKLIX]
          LL = GEOmap::XY.GLOB(ex, why, proj)

          kx = floor(ex/MOD$A$dx)
          ky = floor(why/MOD$A$dy)
          
          IM = MOD$MOD[[ilay]]
          print(paste(sep=" ", "Layer=", ilay))
          print(cbind(kx, ky))
          
          PLOC = cbind(ex, why, LL$lat, LL$lon, IM[cbind(kx, ky)])
          
          print( PLOC )
           zloc = list(x=NULL, y=NULL)
        }

    if(K[Nclick] == match("AREA", BLABS, nomatch = NOLAB))
        {
          NKLIX = length(zloc$x)-1
          ex = zloc$x[1:NKLIX]
          why = zloc$y[1:NKLIX]
          LL = GEOmap::XY.GLOB(ex, why, proj)
          polygon(ex, why, border='white')
          mypoly = cbind(ex, why)
          myarea = splancs::areapl(mypoly)
          print(paste(sep=' ', "AREA=",myarea))
          zloc = list(x=NULL, y=NULL)
        }


      if(K[Nclick] == match("ASTATS", BLABS, nomatch = NOLAB))
        {
          NKLIX = length(zloc$x)-1
          ex = zloc$x[1:NKLIX]
          why = zloc$y[1:NKLIX]
          LL = GEOmap::XY.GLOB(ex, why, proj)
          polygon(ex, why, border='white')

          
          mypoly = cbind(ex, why)

         
          
          INP =  splancs::inpip(PTSXY,mypoly )

          points(PTSXY[INP, 1] , PTSXY[INP, 2], pch='.', col='white')
          myarea = splancs::areapl(mypoly)

          IM = MOD$MOD[[ilay]]

          kx = floor(PTSXY[INP, 1]/MOD$A$dx)
          ky = floor(PTSXY[INP, 2]/MOD$A$dy)
         
          Ival = IM[cbind(kx, ky)]

          SVAL = sum(Ival)
          
          statval = jstats(Ival)
          
          print(paste(sep=' ', "AREA=",myarea))
          print(paste(sep=' ', "SUM=", SVAL))
          print("STATISTICS:")
          print(statval)
          zloc = list(x=NULL, y=NULL)
          
        }

      

      if(K[Nclick] == match("LINE", BLABS, nomatch = NOLAB))
        {
          NKLIX = length(zloc$x)-1
          ex = zloc$x[1:NKLIX]
          why = zloc$y[1:NKLIX]
          n = length(ex)
          LL = GEOmap::XY.GLOB(ex, why, proj)
          ##   check that n is even
           ##  print(paste(sep=' ', n, n%%2))
          if( identical(n%%2, 0))
             {
               x1 = ex[seq(from=1, to=n, by=2)]
               x2 = ex[seq(from=2, to=n, by=2)]
               y1 = why[seq(from=1, to=n, by=2)]
               y2 = why[seq(from=2, to=n, by=2)]
                
               segments(x1,y1, x2, y2, col='white')

               dis = sqrt((x2-x1)^2+(y2-y1)^2)
               print(dis)
             }
          
          zloc = list(x=NULL, y=NULL)
        }

      if(K[Nclick] == match("XSEC", BLABS, nomatch = NOLAB))
        {

          ##########  zloc = locator(col=rgb(1,0.8, 0.8), type='p' )
          ##########  Nclick = length(zloc$x)

         NKLIX = length(zloc$x)-1
          ex = zloc$x[1:NKLIX]
          why = zloc$y[1:NKLIX]
          n = length(ex)
          LL = GEOmap::XY.GLOB(ex, why, proj)
          ##   check that n is even
          ##   print(paste(sep=' ', n, n%%2))
          if( identical(n%%2, 0))
             {
               x1 = ex[seq(from=1, to=n, by=2)]
               x2 = ex[seq(from=2, to=n, by=2)]
               y1 = why[seq(from=1, to=n, by=2)]
               y2 = why[seq(from=2, to=n, by=2)]
                
               segments(x1,y1, x2, y2, col='white')
               text(x2, y2 , labels=1:length(x2), pos=3, col='white')

               dis = sqrt((x2-x1)^2+(y2-y1)^2)
               print(dis)

               zmax = ZLIM[2]
               for(jsec in 1:length(x1))
                 {
                   ############# X11(width = 10, height = 10*zmax/dis[jsec])
                   if(!is.null(TOPO))
                     {
                       topoZ = get2Drayblox(x1[jsec], y1[jsec], x2[jsec], y2[jsec], TOPO$x, TOPO$y , NODES=FALSE, PLOT=FALSE)
                       topo.trace =   TOPO$z[cbind(topoZ$ix, topoZ$iy ) ]/1000
                       topo.along  = sqrt( (topoZ$nodes$x-topoZ$nodes$x[1])^2+(topoZ$nodes$y-topoZ$nodes$y[1])^2)
                       TOPTRACE = list(x=topo.along[1:(length(topo.along)-1)], z = topo.trace)
                     }
                   print(paste(sep=" ", "x1=",x1,"; y1=",y1, "; x2=", x2, "; y2=",y2))

                 dev.new()
 
                   XSEC.drive(MOD, x1[jsec], y1[jsec], x2[jsec], y2[jsec] , zmax=zmax, COL=tomo.colors(100), LIM=LIM, STA=STA, TOP=TOPTRACE)
                 }

               
               
             }
          
          dev.set( MAINdev)
           zloc = list(x=NULL, y=NULL)
        }
      if(K[Nclick] == match("TOPO", BLABS, nomatch = NOLAB))
        {
          
          contour(TOPO, add=TRUE)
          zloc = list(x=NULL, y=NULL)
          
        }
      
      iloc = locator(1, col=rgb(1,0.8, 0.8), type='p')
      zloc  = list(x=c(zloc$x,iloc$x), y=c(zloc$y, iloc$y))
      Nclick = length(iloc$x)
      
      if(is.null(iloc$x)) { return(zloc) }
      K =  RPMG::whichbutt(iloc ,buttons)
        
    
      
      
    }
    return(PLOC)
  }

