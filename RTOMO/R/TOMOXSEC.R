`TOMOXSEC` <-
function(MOD, x1, y1, x2, y2 , zmax=100, depth=c(-25, 0) , COL=rainbow(100), LIM=NULL, STA=NULL, PLOT=TRUE)
  {
    if(missing(COL)) { COL=tomo.colors(100) }
    if(missing(LIM)) { LIM=NULL }
    if(missing(STA)) { STA=NULL }
    if(missing(zmax)) { zmax = max(MOD$D) }
    if(missing(PLOT)) { PLOT=TRUE }
    if(missing(depth)) { depth=NULL }


    IYZ = get2Drayblox(x1, y1, x2, y2, MOD$x, MOD$y , NODES=FALSE, PLOT=FALSE)

    zed =  MOD$D[MOD$D<=zmax]

    LAUGH = matrix(NA, ncol=length(IYZ$ix), nrow=length(zed))
    
    Ixysel = cbind(IYZ$ix, IYZ$iy)

    for(ZED in 1:length(zed))
      {
        IM = MOD$MOD[[ZED]]
        LAUGH[ZED, ] = IM[Ixysel]

      }

    if(!is.null(LIM))
      {

        if(length(LIM)==2)
          {
            LAUGH[LAUGH<LIM[1]&LAUGH>LIM[2]] = NA
          }
      
      }

    along = sqrt( (IYZ$nodes$x-IYZ$nodes$x[1])^2+(IYZ$nodes$y-IYZ$nodes$y[1])^2)

    L1 = t(LAUGH)
    L1 = L1[,rev(1:dim(L1)[2])]

    xz = list(x=along, y=rev(-c(zed, zed[length(zed)]) ), z=L1)


    
    if(PLOT) {  PLOT.TOMOXSEC(xz, COL=COL, LIM) }
    invisible(xz)
  }

