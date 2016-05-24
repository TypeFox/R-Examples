`fancyarrows` <-
function(x1,y1,x2,y2, thick=.08, headlength=.4, headthick=.2, col=grey(.5),  border='black')
  {
    if(missing(thick)) {  thick=.08}
    if(missing(headlength)) { headlength=.4 }
    if(missing(headthick)) {headthick=.2 }
    if(missing(col)) {col=grey(.5) }
    if(missing(border)) { border='black'}


    
    up = par("usr")
    ui = par("pin")

    ratx = (up[2]-up[1])/ui[1]
    raty=  (up[4]-up[3])/ui[2]
  ##   if(length(siz)<length(x)) { siz=rep(siz,length=length(x)) }
     if(length(col)<length(x1)) { col=rep(col,length=length(x1)) }

##     thick=.4; headlength=4; headthick=1
    
##  x1=-12;y1=-10;x2=12; y2=3
    
 for(i in 1:length(x1))
      {
    
    halfheadthickness = headthick/2
    halfthickness = thick/2
    dx = (x2[i]-x1[i])/ratx
    dy = (y2[i]-y1[i])/raty
    arrowlength = sqrt(dx^2 + dy^2)
    angle = atan2(dy, dx)
    base = arrowlength - headlength 

    iX = ratx*c(0, base, base, arrowlength, base, base, 0)
    iY=  raty*c(-halfthickness, -halfthickness, -halfheadthickness, 0, halfheadthickness, halfthickness, halfthickness)
    

    X = x1[i]+iX*cos(angle)-iY*sin(angle)
    Y =  y1[i]+iX*sin(angle)+iY*cos(angle)

    polygon(X,Y, col=col, border=border )
    
  }
    
  }

