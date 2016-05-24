`bang` <-
function(x1,y1,x2,y2)
  {
    xp = x1*y2-y1*x2
    xp[xp>=0] =  1
    xp[xp<0]  = -1
    jang = x1*x2+y1*y2
    if(jang>1) jang=1
    ang = xp*acos(jang)
    return(ang)
  }

