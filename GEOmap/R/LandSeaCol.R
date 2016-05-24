LandSeaCol<-function(IZ, coastmap, PROJ, calcol=NULL)
  {
  setLANDcol<- function()
      {
##########  set up a color map good for topographic displays


        cmat = matrix( c(
          -20000,  70,      130 ,    180,  -3000,   141,     182,     205  ,
          -3000,   162,     181 ,    205,  -2000   ,  188  ,     210   ,    238     ,
          -2000,   202,     225 ,    255 , -100 ,   176 ,    196 ,    222  ,
          -100,    107,     142,     35   ,  0.1 ,    85 ,     107  ,   47 ,
          0.1,     85 ,     107 ,    47 ,   300 ,    143 ,    188  ,   143   ,
          300 ,    105  ,   139 ,    105 ,  600 ,    180,     238 ,    180 ,
          600 ,    193,     255 ,    193 ,  1000 ,   255 ,    211 ,    155     ,
          1000 ,   238,     197,     145,   2000 ,   255 ,    255 ,    255   ,
          2000 ,   255,     255,     255 ,  3500 ,   255 ,    255 ,    255  ),
          ncol = 8, byrow=TRUE)

        notes = c(
          '#SteelBlue,LightSkyBlue3',
          '#LightSteelBlue3,LightSteelBlue2',
          '#LightSteelBlue1,LightSteelBlue',
          '#OliveDrab,DarkOliveGreen',
          '#DarkOliveGreen,DarkSeaGreen',
          '#DarkSeaGreen4,DarkSeaGreen2',
          '#DarkSeaGreen1,burlywood1',
          '#burlywood2,White',
          '#White,White')


        calcol=list(z1=cmat[,1], r1=cmat[,2],g1=cmat[,3],b1=cmat[,4],
          z2=cmat[,5],  r2=cmat[,6],g2=cmat[,7],b2=cmat[,8], note=notes)


        coltab = cbind(calcol$r1, calcol$g1, calcol$b1,calcol$r2, calcol$g2, calcol$b2)

        coltab = rbind(coltab, coltab[length(calcol$r1),])

        return(list(calcol=calcol , coltab=coltab))


      }

 if(is.null(calcol)) {
   A = setLANDcol()

   calcol = A$calcol




 


   
 }

##   source("LandSeaCol.R")

 
    cat("Setting Colors....please wait....", file="", sep="\n")


    ##  image(x=xo, y=yo,   z=UZ, col=blues, asp=TRUE , axes=FALSE, xlab="", ylab="" )

    ##  image(x=xo, y=yo,   z=AZ, col=topo.colors(100), asp=TRUE , axes=FALSE, xlab="", ylab="", add=TRUE )

    ##  plotGEOmapXY(japmap, PROJ=PROJ,CZ LIM=c(A$LON[1], A$LAT[1],A$LON[2], A$LAT[2] ) , add=TRUE)


mxy = setplotmat(IZ$x, IZ$y)
 
FLAGland = matrix(FALSE, ncol=ncol(mxy$x) , nrow=nrow(mxy$x))


 
## plotGEOmapXY(coastmap, PROJ=PROJ)


##  plotGEOmapXY(coastmap, LIM = c(-180, -90, 180, 90), PROJ =PROJ)

##  coastmapXY(coastmap, LIM = c(-180, -90, 180, 90), PROJ =PROJ)


 rectoverlap<-function( x1,  y1,  x2,  y2, 
                       x3,  y3,  x4,  y4)
   {
#### /*  return true if two rectangles overlap  */
#### /*   rectangles are defined by lower left and upper right corners  */

     return((x2>=x3)&(x4>=x1)& (y2>=y3)&(y4>=y1))
   }



 MAPXY = GLOB.XY(coastmap$POINTS$lat ,  RPMG::fmod( coastmap$POINTS$lon, 360) , PROJ )

 mypoly  = splancs::as.points(as.vector(mxy$x) ,as.vector(mxy$y))
 
 if(TRUE)
   {
     ##  GG = GLOB.XY(c(-90, 90) , c(-180,  180) , PROJ)

     ## plot(GG, type='n', asp=1)
     rx = range(IZ$x)
     ry = range(IZ$y)
     
     RADmax = max( c( sqrt(rx[1]^2+ry[1]^2)  , sqrt(rx[2]^2+ry[2]^2) ))
     
     
     for(i in 1:length(coastmap$STROKES$index))
       {
       j1 = coastmap$STROKES$index[i]+1
       j2 = j1+coastmap$STROKES$num[i]-1
       
       
       JEC = j1:j2
       x = MAPXY$x[JEC]
       y = MAPXY$y[JEC]
       
       ##  lines(x,y, col=grey(.8) )
       
       hx = range(x)
       hy =range(y)
       
       rads = sqrt(x^2+y^2)


if(any(rads<RADmax) & rectoverlap(rx[1], ry[1], rx[2], ry[2], hx[1], hy[1], hx[2], hy[2])  )
{


      JJpoly = splancs::as.points(x, y)
      INTEMP = splancs::inout(mypoly,JJpoly, bound=TRUE )

     ##    points(mxy$x[INTEMP] ,mxy$y[INTEMP], pch=".", col='red') 
  FLAGland[INTEMP] = TRUE
     }

}


## plot(GG, type='n', asp=1)


## tem = trans3d(mxy$x[FLAGland==TRUE],   mxy$y[FLAGland==TRUE], rep(0, length(mxy$x[FLAGland==TRUE])), PMAT)

##  points(tem$x , tem$y, pch='.')


##   points(mxy$x[FLAGland==TRUE] ,mxy$y[FLAGland==TRUE], pch='.')


  

}


FLAGwat = t(FLAGland)

FLAGwat = FLAGwat[ , rev(seq(from=1, to=length(IZ$y))) ]

UZ = IZ$z
UZ[(FLAGwat)==TRUE ] = NA


## UZ =UZ[ rev(seq(from=1, to=length(IZ$y))), ]


## image(x=IZ$x, y=IZ$y, z=IZ$z, col=rainbow(100), add=FALSE)

## image(x=IZ$x, y=IZ$y, z=(UZ), col=blues, add=TRUE)

## image(x=IZ$x, y=IZ$y, z=(AZ), col=terrain.colors(100) , add=TRUE)





AZ = IZ$z

    AZ[(FLAGwat)!=TRUE] = NA

    blues = RPMG::shade.col(100, acol=as.vector(col2rgb("darkblue")/255)   , bcol= as.vector(col2rgb("paleturquoise")/255))





    CZ = AZ
########  TZ[TZ<0] = NA

    dz = dim(AZ)
    
    coltab = cbind(calcol$r1, calcol$g1, calcol$b1,calcol$r2, calcol$g2, calcol$b2)
    
    coltab = rbind(coltab, coltab[length(calcol$r1),])
    
    rng = range(AZ[!is.na(AZ)])
    ncol = 100
    levs = seq(from=rng[1], to=rng[2], length=100)
    
    
    brs = c(calcol$z1[1], calcol$z2,  calcol$z2[length(calcol$z2)]+10000)
    
    Cs = CZ

    naflag = which(is.na(Cs))
    
    fs = findInterval(Cs, brs, all.inside =TRUE)
    df = diff(brs)
    
    coldis = (Cs - brs[fs])/df[fs]
    
    newcol = list(r= round(coltab[fs,1]+coldis*(coltab[fs,4]-coltab[fs,1]  )),
      g = round(coltab[fs,2]+coldis*(coltab[fs,5]-coltab[fs,2]  )),
      b = round(coltab[fs,3]+coldis*(coltab[fs,6]-coltab[fs,3]  )))

    newcol$r[naflag] = 0
    newcol$g[naflag] = 0
    newcol$b[naflag] = 0

    
    Collist = rgb(newcol$r/255, newcol$g/255, newcol$b/255)

    
#### bluesrgb = col2rgb(blues)

####    Collist[IZ$z<0] = blues
    
####   Mollist = matrix(data=Collist, ncol=dz[2], nrow=dz[1])
    
####    PMAT = persp(jx, jy , TZ, theta = 0, phi = 90, r=4000, col=Mollist[1:479, 1:359] , scale = FALSE,
####      ltheta = 120, lphi=60, shade = 0.75, border = NA, expand=0.001, box = FALSE )


    if( all(is.na(UZ)) )
      {
        C2 = Collist
      }
    else
      {
        
        bluesrgb = col2rgb(blues)

        rngU = range(UZ[!is.na(UZ)])
        
        ncol = 100
        levs = seq(from=rngU[1], to=rngU[2], length=100)
        
        
        brs = levs
        Cs = UZ[!is.na(UZ)]
        fs = findInterval(Cs, brs, all.inside =TRUE)
        
        
####   df = diff(brs)
        
####  coldis = (Cs - brs[fs])/df[fs]
        
        
        Anewcol = list(r=bluesrgb[1,fs],
          g =bluesrgb[2,fs] ,
          b =bluesrgb[3,fs] )
        
        
        UCollist = rgb(Anewcol$r/255, Anewcol$g/255, Anewcol$b/255)
        
        C2 = Collist
        C2[!is.na(UZ)] = UCollist
        
      }
    
    CollorMatrix = matrix(data=C2, ncol=dz[2], nrow=dz[1])
    
    Dcol = dim(CollorMatrix)       
    attr(CollorMatrix, 'Dcol') <- Dcol
    
    invisible(list(Cmat=CollorMatrix, UZ=UZ, AZ=AZ) )

  }
