`NSarrow` <-
function(x=NULL, y=NULL, R=1, col.arrow=1, col.N=1, col.circ=1, rot=0, PMAT=NULL)
  {
   ###   complex north-south arrow for maps

    if(is.list(x))
      {
        vane = x
      }
    else
      {
        vane = list(x=x, y=y)

      }

    if(is.null(vane$x)) {  vane = locator(2) }


    if(missing(y)) { y = vane$y }
    if(missing(col.arrow)) { col.arrow="black" }
    if(missing(col.N)) { col.N="black" }
    if(missing(col.circ)) { col.circ="white" }
    if(missing(rot)) { rot=0 }
    if(missing(PMAT)) { PMAT=NULL }


    if (!is.null(PMAT)) {
        tem1 = trans3d(vane$x, vane$y, rep(0, length(vane$y)), 
            PMAT)
    }
    else {
        tem1 = vane
    }


    
#######   rot is in degrees
    if(missing(R))
      {
        if(length(tem1$y)>1)
          {
            R = abs(tem1$y[2]-tem1$y[1])

          }
        else
          {
            R = 1
          }
      }



    
   C = RPMG::circle()
   #   rect(-.2, -1, .2, .5)

   ##  col.arrow="black"; col.N="black"; col.circ="white"
    #  plot(c(-1,1)  , c(-1,1)  , type='n', asp=TRUE)
   
   crad = 0.4 * R



   
    x1 = c(-.1, .1, .1, -.5, -.1)
    y1 = c(-1, -1, 1, .15,  .3)
   #   points(x1+x, y1+y)
   #  text(x1+x,y1+y,labels=1:length(y1), pos=4)

   N = list(x=R*c(-0.2, -0.2, -0.1, 0.1, 0.1, 0.2, 0.2, .1, -0.1, -0.1  ),
     y=R*(.05+c(-0.5, 0, 0, -0.3, 0, 0, -0.5, -0.5, -0.2, -0.5)))

       ROT = rotmat2D(rot*pi/180)
       arr = cbind(R*x1, R*y1)
       arr = arr %*%  ROT
       
   cir = cbind(crad*c(C$x, C$x[1]), crad*(c(C$y,C$y[1]))-R*0.2)
   cir = cir %*%  ROT

   NOR = cbind(N$x, N$y)
       NOR = NOR %*%  ROT
 
   polygon(arr[,1]+tem1$x[1], arr[,2]+tem1$y[1] , col =col.arrow)
   
   
   
   polygon(cir[,1]+tem1$x[1] , cir[,2]+tem1$y[1] , col=col.circ, lwd=1.2)


    polygon(NOR[,1]+tem1$x[1], NOR[,2]+tem1$y[1], col=col.N )

   # points(N)
   # text(N, labels = 1:length(N$x), pos=1 )


    #  plot(c(-1,1)  , c(-1,1)  , type='n', asp=TRUE)
   
   # text(0,0, labels="N", cex=10)

   return(vane)
  }

