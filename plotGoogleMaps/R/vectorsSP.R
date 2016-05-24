vectorsSP <-
  function(SPDF,
           zcol=1:2,
           maxlength=30000,   # m max arrow length
           arrowSize=0.15, # 0.15*maxlength
           arrAng=30  # degrees
          ){
    SP=SPDF
    proj=SP@proj4string
    data=SP@data
    SP=spTransform(SP,CRS("+proj=merc +lon_0=0 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84"))  
    
    range01 <- function(x){(x-min(x))/(max(x)-min(x))}
    
    SP@data=SP@data[,zcol]
    SP@data[,1]=range01(abs(SP@data[,1])) +.2 # avoid 0
    
    SPVect=lapply(1:length(SP@data[,1]),function(j){
      nx= SP@coords[j,1]
      ny= SP@coords[j,2]
      
      r= SP@data[j,2] * pi/180   # direction 
      
      delta= SP@data[j,1]*maxlength
      
      mx = nx + sin( r ) * delta
      my = ny + cos( r ) * delta
      
      delta=arrowSize*delta
      
      bx = mx;
      by = my;
      dx = bx - sin( r ) * delta
      dy = by - cos( r ) * delta
      
      r1= r+arrAng/180*pi
      ax = dx-sin( r1 ) * delta
      ay = dy - cos( r1 ) * delta
      r2=r-arrAng/180*pi
      cx = dx - sin( r2 ) * delta
      cy = dy - cos( r2 ) * delta
      
      l1 = cbind(c(nx,bx,ax,bx,cx), c(ny,by,ay,by,cy))
      Sl1 = Line(l1)
      S1 = Lines(list(Sl1) ,ID = paste(j))
      S1 } )
    
    Sl = SpatialLines(SPVect)
    Sldf = SpatialLinesDataFrame(Sl, data = data, match.ID = F)
    Sldf@proj4string=SP@proj4string
    Sldf=spTransform(Sldf,proj)
    
  } # end of vectors
    
    
    