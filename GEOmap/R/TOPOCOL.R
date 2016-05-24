TOPOCOL<- function(IZ, calcol=NULL)
  {
    ######  source("/home/lees/Progs/R_PAX/GEOmap/R/TOPOCOL.R")
    if(missing(calcol)) {  CCOL = settopocol(); calcol = CCOL$calcol }
    
    if(is.null(calcol)) {  CCOL = settopocol(); calcol = CCOL$calcol  }
    
  
    if(!is.list(IZ)) { IZ = list(z=IZ) }
    
    UZ = IZ$z
    UZ[IZ$z>= .001 ] = NA

    AZ = IZ$z
    AZ[IZ$z<=-.001] = NA

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
    
    fs = findInterval(Cs, brs)
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
    
####   Cmat = matrix(data=Collist, ncol=dz[2], nrow=dz[1])
    
####    PMAT = persp(jx, jy , TZ, theta = 0, phi = 90, r=4000, col=Cmat[1:479, 1:359] , scale = FALSE,
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
       fs = findInterval(Cs, brs)
       
       
####   df = diff(brs)
    
####  coldis = (Cs - brs[fs])/df[fs]
       
       
       Anewcol = list(r=bluesrgb[1,fs],
         g =bluesrgb[2,fs] ,
         b =bluesrgb[3,fs] )
       
       
       UCollist = rgb(Anewcol$r/255, Anewcol$g/255, Anewcol$b/255)
       
       C2 = Collist
       C2[!is.na(UZ)] = UCollist
       
     }
    
    Cmat = matrix(data=C2, ncol=dz[2], nrow=dz[1])
    
    Dcol = dim(Cmat)       

    attr(Cmat, 'Dcol') <- Dcol
    invisible(Cmat)

  }
