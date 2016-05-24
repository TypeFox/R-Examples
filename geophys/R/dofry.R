dofry<-function(x,y, PLOT=FALSE)
  {
    if(missing(PLOT)) PLOT=FALSE
    rx = range(x)
    ry = range(y)

    EX = seq(from=rx[1]-diff(rx), to=rx[2]+diff(rx), length=100)
    WHY = seq(from=ry[1]-diff(ry), to=ry[2]+diff(ry), length=100)

    deltax = (EX[2]-EX[1])/2
    deltay = (WHY[2]-WHY[1])/2


###     MM = meshgrid(EX, WHY)
### desh(MM, PTS=FALSE, add=FALSE)

    
###    thefry = matrix(0, ncol=length(EX), nrow=length(WHY))


if(PLOT) plot(range(EX), range(WHY), type='n', asp=1)

    thefry = list()
    
    for(i in 1:length(x))
      {

        DX = x-x[i]+mean(EX)
        DY = y-y[i]+mean(WHY)
 ###       points(DX, DY, pch=3, cex=.6, col='green')
        
###   get rid of the center point
        DX =  DX[-i]
        DY =  DY[-i]

      if(PLOT)   points(DX, DY, pch=".", cex=2)
       
###       i1 = findInterval(DX, EX, rightmost.closed = FALSE, all.inside = FALSE)
###       i2 = findInterval(DY, WHY, rightmost.closed = FALSE, all.inside = FALSE)

     ###   rect(EX[i1]-deltax, WHY[i2]-deltay, EX[i1]+deltax, WHY[i2]+deltay, border='blue') 
###        thefry[cbind(i1, i2)] = thefry[cbind(i1, i2)]+1

        thefry[[i]] = list(x=DX, y=DY)
       
      }

    attr(thefry, "X") = EX
    attr(thefry, "Y") = WHY
    attr(thefry, "mx") = mean(EX)
    attr(thefry, "my") = mean(WHY)
    
    
    invisible(thefry)

  }
