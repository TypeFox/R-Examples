`partmotnet` <-
function(temp , LINES=FALSE, STAMP=STAMP, COL= rainbow(100))
  {
    if(missing(STAMP)) { STAMP = " " }
    if(missing(LINES)) { LINES=FALSE }
    if(missing(COL)) { COL = rainbow(100)  }

    n = length(temp$x)
    ncol = length(COL)

    if( ncol >3 )
      {
        
        KR = n/(ncol-1)
        
        KX = floor(seq(from=1, length=n)/(KR))+1
        
        cols=COL[KX]
      }
    else
      {
        
        cols = 1

      }

    temp$x = temp$x-mean(temp$x)
    temp$y = temp$y-mean(temp$y)
    temp$z = temp$z-mean(temp$z)
    
    
    nrm = sqrt(temp$x^2+temp$y^2+temp$z^2)

    nrm[is.na(nrm)] = 1
    nrm[nrm==0] = 1
    
    bascd = list(x=temp$x/nrm, y=temp$y/nrm, z=temp$z/nrm)

   circ()
    segments(c(-1, 0), c(0, -1), c(1, 0), c(0, 1), col=rgb(0.85, 0.85, 0.85), lty=2)

  #####  if(LINES==TRUE)   lines(bascd$x, bascd$y, col=rgb(0.85, 0.85, 0.85))

    if(LINES==TRUE)
      { segments(bascd$x[1:(n-1)], bascd$y[1:(n-1)], bascd$x[2:n], bascd$y[2:(n)],  col=cols[2:n]) }
    
  
    zflag = bascd$z>=0

    points(bascd$x[zflag], bascd$y[zflag], type='p', pch=2, col=cols[zflag])
    
    points(bascd$x[!zflag], bascd$y[!zflag], type='p', pch=5, col=cols[!zflag])
    
   #####  if(LINES==TRUE)   lines(bascd$x, bascd$y, col=rgb(0.85, 0.85, 0.85))
    
    
    title(sub=STAMP)
    

  }

