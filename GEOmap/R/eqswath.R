eqswath<-function(x,y, z, L, width=1, PROJ=NULL)
{

  d1 =   sqrt( (L$x[1]-L$x[2])^2+ (L$y[1]-L$y[2])^2 )

  if(d1<=0) return(0)

  if(missing(width)) { width= 0.05*d1 }

  

 
cs1 =  (L$x[2]-L$x[1])/d1
sn1 =   (L$y[2]-L$y[1])/d1



  ##   translation
  x1 = x-L$x[1]
  y1 = y-L$y[1]

  ## first rotation


x2 = cs1* x1 + sn1*y1
y2 = -sn1* x1 + cs1*y1

  ##   dev.set(3)
  ##  plot(x2, y2, pch='.', asp=1)
  ##  segments(0, 0, 0, d1, lwd=2)

  gy1 = -width
  gy2 = width

  flag = (y2>=gy1) & (y2<=gy2) &  (x2>=0) & (x2<=d1)

  ##   segments(L$x[1], L$y[1], L$x[2], L$y[2], lwd=2)

 
boxpts = list(x=c(0,d1,d1,0), y = c(-width,-width,width,width))

x3 = cs1* boxpts$x - sn1* boxpts$y
y3 = sn1* boxpts$x + cs1* boxpts$y
InvBox = list(x=x3+L$x[1], y=y3+L$y[1])

 
  ## points(boxpts, col='blue')
  ## text(boxpts, labels=1:4, col='blue')
  ## points(rotpts, col='red')
  ## text(rotpts, labels=1:4, col='red')


  ## points(swathpts, col='blue')
  ## text(swathpts, labels=1:4, col='blue')

 ##  dev.set(3)




  ##  plot(r, -z2)

  invisible(list(r=x2[flag] , dh=y2[flag] , depth=z[flag] ,    flag=which(flag) , InvBox = InvBox, rdis=c(0, d1))   )  


}
