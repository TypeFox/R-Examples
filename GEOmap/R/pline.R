pline<-function(x1,y1,x2,y2,ex,ey)
{
  ##########   get the shortest distanvce to a line segment

dx1 = x2-x1
dy1 = y2-y1
gx = ex-x1
gy = ey-y1
hx = ex-x2
hy = ey-y2

d1 =sqrt(dx1^2 + dy1^2)


d2 = sqrt(gx^2 + gy^2)
d3  = sqrt(hx^2 + hy^2)

if(d1==0 | d2==0)
  {
     dis = d2
      px = x1
      py = y1
     zee=0
     dee = dis

     return(c(dis,dee,zee,px,py))

  }




cs = (dx1*gx +dy1*gy) / (d1*d2)
sn = (dx1*gy - gx*dy1)/(d1*d2)

zee = d2*cs
dee = abs(d2*sn)

if(zee<0 | zee>d1) {
  if(d2<d3)
    {
      dis = d2
      px = x1
      py = y1


    }
  else
    {
      dis=d3
      px = x2
      py= y2 


    }

return(c(dis,dee,zee,px,py))

  
}
else {
  dis = dee
  px = x1+zee*dx1/d1
  py = y1+zee*dy1/d1
}


return(c(dis,dee,zee,px,py))

}

##  source('/home/lees/Progs/R_PAX/GEOmap/R/pline.R')
