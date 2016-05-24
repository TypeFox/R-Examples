NLSlocate <-
function(GH, vel=list() , init=c(0,0,0,0) , PLOT=FALSE )
  {
    if(missing(PLOT))  { PLOT=TRUE }
    
    if(missing(vel))
      {
        ##   data(LITHOS.vel)
        LITHOS.vel=defaultVEL(2)
      
        vel= LITHOS.vel
      }

    ##########  get the station information from the pickfile
    STAS = GH$pickfile$STAS

    ###########  set up the projection - use the median station location

    proj = GEOmap::setPROJ(type=2, LAT0 =median(STAS$lat) , LON0 = median(STAS$lon) )

    #############    convert the LATLON of stations to X-Y
    XY = GEOmap::GLOB.XY(STAS$lat, STAS$lon, proj)

    ##########   set the initial guess: z= 6 km 
    initz = 6.0
    t0a = initz/vel$vp[1]
    ####################  use the station that has the earliest arrival time
    w1 = which.min(STAS$sec)
    if(missing(init))
      {
####  find station with earliest arrival time, set it below this station
###   at initial depth
        initz = 6.0
        t0a = initz/vel$vp[1]
        w1 = which.min(STAS$sec)
        h1 = c(XY$x[w1], XY$y[w1], initz,  STAS$sec[w1]-t0a )
      }
    else
      {
        h1 = init
      }
    
    ##############  set up the matrix 
    XY = XYSETUP(STAS, c(STAS$lat[w1],STAS$lon[w1], initz,  STAS$sec[w1]-t0a  ) , vel  )

    ############  set the function that gets the residuals
    residFun <- function(p, STAS, vel) EQXYresid(XY, vel=vel , h1=p , PLOT=FALSE)
    parStart <- h1
    ############   parstart is the initial guess

    #############   non-linear least square inversion
    nls.out <- minpack.lm::nls.lm(par=parStart, fn = residFun, STAS=XY, vel=vel, control = minpack.lm::nls.lm.control(nprint=1))

    #############  after done, get LAT-LON information 
    LL = GEOmap::XY.GLOB(nls.out$par[1], nls.out$par[2], proj)
    if(LL$lon>180)  LL$lon = LL$lon-360
    SOL = c(LL$lat, LL$lon, nls.out$par[3], nls.out$par[4] )

    OLDsol = c(GH$pickfile$LOC$lat, GH$pickfile$LOC$lon, GH$pickfile$LOC$z,  GH$pickfile$LOC$sec  )

    ###############  this is for comparison with older LQUAKE solution
    print( cbind(SOL, OLDsol    )) 

    #############   
    if(PLOT){

      plot(XY, pch=6, xlab="km", ylab="km")
      text(XY, labels=STAS$name, pos=3)
      points(h1[1], h1[2], pch=3, col='blue')

      gloc = GEOmap::GLOB.XY(GH$pickfile$LOC$lat, GH$pickfile$LOC$lon, proj)
      points(gloc$x, gloc$y, pch=8, col='green')

      points(nls.out$par[1], nls.out$par[2], pch=8, col='red')

       u = par('usr')
     LEG = RSEIS::jlegend( u[1], u[4], c("Init", "OLD", "NEW"), pch=c(3,8,8) , col=c('blue', 'green', 'red'), plot=FALSE  )



    }


    return(SOL)
  }
