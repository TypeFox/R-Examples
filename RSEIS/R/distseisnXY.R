distseisnXY<-function(GH, sta=list(nam="", x=0 , y=0 , z=0) , LOC=list(x=0, y=0 , z=0)   )
{
  ###  given a RSEIS list get the distances to each station 
  m = match( GH$STNS,    sta$nam)
  X =  sta$x[m]
  Y =  sta$y[m]
  Z = sta$z[m]	
  dees =  sqrt((X-LOC$x)^2 + (Y-LOC$y)^2 +  (Z-LOC$z)^2 )

  return(dees)	
}
