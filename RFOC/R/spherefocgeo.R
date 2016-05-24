spherefocgeo<-function(CMTSOL, PROJ=NULL, icut=5, ndivs=10,  bbox=c(0,1, 0, 1), PLOT=TRUE, add=FALSE, RECT=FALSE, pal=terrain.colors(100))
  {

  if(missing(PLOT)) PLOT=TRUE
    if(missing(add)) add=FALSE
     if(missing(RECT)) RECT=FALSE

    if(missing(pal)) { pal=terrain.colors(100) }
    if(missing(PROJ))
      {

        PROJ = GEOmap::setPROJ(type=2, LAT0=mean(CMTSOL$lat) , LON0=mean(CMTSOL$lon) )   ##   utm


      }

    if (missing(icut)) {
        icut = 5
    }
    if (missing(ndivs)) {
        ndivs = 10
    }
   if (missing(bbox)) {
        bbox = par("usr")
    }
   




  XY = GEOmap::GLOB.XY(CMTSOL$lat, CMTSOL$lon, PROJ)
  PHAIT = prepFOCS(CMTSOL)
  

  u = par("usr")
#######  RI = RectDense( XY$x, XY$y, icut=5, u=u, ndivs=10)
  RI = RectDense( XY$x, XY$y, icut=icut, u=bbox, ndivs=ndivs)

  KPspat =  matrix(NA, nrow=length(RI$ipass), ncol=10)
  KTspat =  matrix(NA, nrow=length(RI$ipass), ncol=10)
  colnames(KTspat) = c("x", "y", "n", "Ir",  "Dr", "R", "K", "S", "Alph95" , "Kappa")
  colnames(KPspat) = c("x", "y", "n", "Ir",  "Dr", "R", "K", "S", "Alph95" , "Kappa")

  i=1
  sizy = RI$icorns[i,4]-RI$icorns[i,2]
  sizx = RI$icorns[i,3]-RI$icorns[i,1]
  siz = .5*min(c(sizy, sizx))
  

  
  if(RECT)   rect(RI$icorns[,1],RI$icorns[,2],RI$icorns[,3],RI$icorns[,4], col=NA, border='blue')
  
    for(i in 1:length(RI$ipass))
      {
        flag = XY$x>RI$icorns[i,1]& XY$y>RI$icorns[i,2] & XY$x<RI$icorns[i,3] & XY$y<RI$icorns[i,4]
        paz=PHAIT$Paz[flag]
        pdip=PHAIT$Pdip[flag]
         taz=PHAIT$Taz[flag]
        tdip=PHAIT$Tdip[flag]
        x=mean(RI$icorns[i,c(1,3)])
        y=mean(RI$icorns[i,c(2,4)])
        siz=(RI$icorns[1,3]-RI$icorns[1,1])/2.5

        PlotPTsmooth(paz, pdip, x=x, y=y, siz=siz, border=NA, bcol='white' , LABS=FALSE, add=FALSE, IMAGE=TRUE, CONT=FALSE, pal=pal)
          PlotPTsmooth(taz, tdip, x=x, y=y, siz=siz, border=NA, bcol='white' , LABS=FALSE, add=TRUE, IMAGE=FALSE, CONT=TRUE)

######dev.set(2)
######pnet(MN, add=FALSE)
  
   #### ALPH = alpha95(paz, pdip)
      ##    n = length( paz)
     ##  KPspat[i,] =   c(x, y, n, ALPH$Ir,  ALPH$Dr, ALPH$R, ALPH$K, ALPH$S, ALPH$Alph95 , ALPH$Kappa)

    ####    ALPH = alpha95(taz, tdip) 
     ##   KTspat[i,] =    c(x, y, n, ALPH$Ir,  ALPH$Dr, ALPH$R, ALPH$K, ALPH$S, ALPH$Alph95 , ALPH$Kappa)

         invisible(PHAIT)

      }













  }
