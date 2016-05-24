`faultdip` <-
function(x,y, rot=0, h=1, lab='')
  {

    alpha = (90-rot)*pi/180
    arrows(x,y,x+h*cos(alpha), y+h*sin(alpha), length=h*.1)  
  }

