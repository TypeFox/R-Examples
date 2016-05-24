ROTcovQUAKE<-function(XQ, proj, sel=1, domap=TRUE, cosomap=NA, add=FALSE)
  {

    if(missing(sel)) sel = 1:length(XQ)

    pct = 0.05
    pfact = 1-pct/2
    
    NEQ = length(sel)

    
    ZEDS = matrix(ncol=3, nrow=NEQ)
    POSq = matrix(ncol=3, nrow=NEQ)
    SIZq = matrix(ncol=3, nrow=NEQ)
    ANGq = matrix(ncol=3, nrow=NEQ)
    
    for(i in 1:NEQ)
      {
        j =  sel[i]
        ZEDS[i, ] = c(XQ[[j]]$EQ$lat,XQ[[j]]$EQ$lon, XQ[[j]]$EQ$z )
           }

    
    ZEXY = GEOmap::GLOB.XY(ZEDS[,1], ZEDS[,2], proj)
    
    for(i in 1:NEQ)
      {
        j = sel[i]
        KOV =  XQ[[j]]$ERR$cov[2:4, 2:4]

        
        ndf = XQ[[j]]$ERR$ndf

        U2 = eigen(solve(KOV))
        u = U2$vectors
        
        lam = (U2$values)
        
        delta= qt(pfact, ndf)

        zed =   (delta/  sqrt(lam))
        
        #####  EV = eigen(KOV)
        ##### dels = sqrt(diag(KOV))
       #####  sizes = dels[2:4]

       Eulerangles =getEulers(u)
        
        POSq[i,] = c(ZEXY$x[i] ,ZEXY$y[i] , -ZEDS[i, 3] )
        ANGq[i,] = Eulerangles[c(1,2,3) ]
        SIZq[i,] = zed
        
        
      }
    

if(add==FALSE)
  {
     rgl::rgl.open()
    rgl::par3d(windowRect=c(0, 0, 600, 600) )
   }

 cda::rgl.ellipsoids(POSq, SIZq
                       , ANGq, col="gold")

    

    if(domap)  {
      rglGEOmapXY(cosomap,   PROJ =proj, MAPstyle=2, add=TRUE )
    }
    else
      {

        cenx = mean(ZEXY$x)
        ceny = mean(ZEXY$y)
        cenz = mean( ZEDS[,3 ]) 
      }

        
## rgl.close()


    
  }
