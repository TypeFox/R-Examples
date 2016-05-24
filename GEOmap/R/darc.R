darc <-
function ( rad=1, ang1=0, ang2=360, x1=0, y1=0, n=1)
{

  if (missing(n)) n = 1
  if (missing(x1))  x1 = 0
  if (missing(y1))  y1 = 0
  if (missing(ang1))  ang1=0
  if (missing(ang2))  ang2=360

if(ang1>ang2 & n>0) n = -n

 i = pi * seq(from =ang1, to = ang2, by = n)/180
    cx = rad*cos(i)
    cy = rad*sin(i)
    C = list(x=x1+cx, y = y1+cy)
   return(C)
}

