areatriangles <-
function(line.slope,line.dist1) {
  area.tri1=c()
  ang=sin(-atan(line.slope[1])+pi/2)
  K=length(line.slope)
  for (i in 1:(K+1))
  {
      if (i==1)
      {
       area.tri1[i]=0.5*line.dist1[i]*ang
      }
      if (i==K+1)
      {
       area.tri1[i]=0.5*line.dist1[i-1]*ang
      }
      if (i>1 && i<K+1)
      {
       area.tri1[i]=0.5*line.dist1[i]*line.dist1[i-1]*ang
      }
  }
  area.auctri1=sum(area.tri1)
  answer=list(auctri=area.auctri1,areatri=area.tri1)
  return(answer)
}
