`teeth` <-
function(x,y,h, rot, col='black', border='black' )
  {

#########  plot(c(-1,1), c(-1,1), type='n', asp=1)
######   polygon(px,py)
######   polygon(ax,ay)
if(missing(col)) { col='black' }
if(missing(border)) {border =col }

###  if(missing(border)) { border='black' }

 ###     border=col 
#######  zim = cos(60*pi/180)
    zim =0.5
    px = c(0, -h*zim, h*zim)
    py=c(h, 0, 0 )
#####    polygon(px,py)


    cs = rot$cs
    sn = rot$sn

    for(i in 1:length(x))
      {
        
        ax = x[i]+(px*cs[i]-py*sn[i])
        ay = y[i]+(px*sn[i]+py*cs[i])
        polygon(ax,ay, col=col, border=border)
      }


  }

