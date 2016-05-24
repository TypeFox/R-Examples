`GXMA3DV` <-
function(name)
  {
    ##  GXMA3D function to get a tomographic inversion result from doinv
    ##  
    A = scan(file=name, nmax=1, list(lat=0, lon=0, nx=0, ny=0, nz=0, dx=0, dy=0))
    D = scan(file=name, skip=1, n=A$nz)
    V = scan(file=name, skip=2, n=A$nz)

    M = scan(file=name, skip=3)

    
    A$skip = 1000.0;
    
    M[M==A$skip] = NA

    x = A$dx*seq(from=0, length=A$nx)
    y = A$dy*seq(from=0, length=A$ny)
    


    MOD = as.list(1:A$nz)


    tot = A$nx*A$ny*A$nz
    
    toplen = A$nx*A$ny

    
    for(i in 1:A$nz)
      {
        k = (i-1)*toplen
        MOD[[i]] = t(matrix(M[(k+1):(k+toplen)], ncol=A$nx, nrow=A$ny, byrow=TRUE))

        ###  image(x,y,MOD[[i]], col=tomocolors)
        ###   locator()
        
      }


    ###  image(x,y,MOD[[3]], col=tomocolors)
    ###  image(x,y,MOD[[5]], col=tomocolors)

  ###     DAPROJ = PROJ.DATA
   ###    GEOmap::setPROJ(type=2, LAT0=A$lat, LON0=A$lon )

  ###     LL = GEOmap::XY.GLOB( x , rep(0,length(x)) )
  ###     lons = LL$lon
  ###     LL = GEOmap::XY.GLOB( rep(0,length(y)) , y )
  ###     lats = LL$lat

   ###    PROJ.DATA<<-DAPROJ

  ###     mx = GEOmap::GLOB.XY(A$lat, A$lon  )
  ###     kmx = mx$x+x
   ###    kmy = mx$y+y

    return(list(name=name, A=A, D=D, V=V, MOD=MOD, x=x, y=y) )

    
  }

