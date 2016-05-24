# Author: Carlos A. Silva. Uses code by Remko Duursma from Maeswrap package, see Plotstand.
TreesModel<- function(crownshape=c("cone","ellipsoid","halfellipsoid","paraboloid","cylinder"),
                      nz=5, nalpha=5, CL=5, CW=5, HCB=10, x0=0, y0=0, z0=0, dbh = 0.3, crowncolor = "forestgreen", 
                      stemcolor = "chocolate4"
){
  
  crownshape <- match.arg(crownshape)
  
  z <- rep(seq(0,1,length=nz),each=nalpha)
  angs <- rep(seq(0,2*pi, length=nalpha),nz)
  
  if(crownshape == "cone")distfun <- (1-z)
  if(crownshape == "ellipsoid")distfun <- sqrt(1 - ((z-1/2)^2)/((1/2)^2))
  if(crownshape == "halfellipsoid")distfun <- sqrt(1 - z**2)
  if(crownshape == "paraboloid")distfun <- sqrt(1-z)
  if(crownshape == "cylinder")distfun <- 1
  H <- HCB + CL
  r <- CW/2
  x <- x0 + r*distfun*cos(angs)
  y <- y0 + r*distfun*sin(angs)
  z <- z0 + HCB + z*CL
  
  keep <- !duplicated(cbind(x,y,z))
  x <- x[keep]
  y <- y[keep]
  z <- z[keep]
  klj=matrix(cbind(x,y,z),ncol=3)
  
  mMatrix<-matrix(,ncol=3)[-1,]
  
  for ( i in 1:nrow(klj)){
    ln=i+nz
    
    if ( ln >= nrow(klj)) { ln2=nrow(klj) } else { ln2= ln}
    
    mMatrix<-rbind(mMatrix,rbind(klj[i,],klj[ln2,])) }
  
  
  kljzbase=subset(klj,klj[,3]==z[2])
  kljzbaseNew<-matrix(,ncol=3)[-1,]
  
  for ( i in 1:nrow(kljzbase)){
    kljzbaseNew<-rbind(kljzbaseNew,rbind(kljzbase[i,],c(x0,y0,HCB)))
    
  }
  
  newList<-rbind(kljzbaseNew,mMatrix,klj)
  plot3d(newList, type="l", col=crowncolor, add=T)
  m2 <- coord3dshape("cone", CW = dbh, CL = H, z0 = z0, x0 = x0, 
                     y0 = y0, nz = 50, nalpha = 50)
  interpol(m2, col = stemcolor)
  
}
