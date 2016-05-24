`PPIX` <-
function(zloc, YN=NULL, col=1, lab='')
  {
    if(missing(YN)) { YN = 1 }
    if(missing(col)) { col = 1 }
    if(missing(lab)) { lab = NA }
    
    du = 1/(YN)
    j = floor((zloc$y)/du)
    y1 = j*du
    y2 = y1+du
    segments(zloc$x, y1, zloc$x, y2, col=col)
    if(!is.na(lab))
      {
        text(zloc$x, y2, labels=lab, pos=4)
      }
    
  }

