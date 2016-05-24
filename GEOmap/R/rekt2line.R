rekt2line<-function(rekt, pnts)
  {
    ###  source('rekt2line.R')
     ###  source('/home/lees/XMDEMO_R/INDONESIA/rekt2line.R')

    
###  rekt is a rectangle with 4 points counterclockwise
    ###  pnts is a set of point pairs each representing a line
    ###  the purpose of the routine is
    ###  to find the points on the rectangle where the lines intersect.

    PEES = list(x=vector(), y=vector() )

    for(i in 1:length(pnts$x1))
      {
        ###   test each line of the rectangle


        DX = pnts$x2[i]  -pnts$x1[i]
        DY = pnts$y2[i]  -pnts$y1[i]
        if(DX!=0)
          {
            m = DY/DX
            b = pnts$y1[i]  - m * pnts$x1[i]

            y1 = rekt$x[1]*m + b
            y2 = rekt$x[2]*m + b
            if(m!=0)
              {
                x1 = (rekt$y[1]-b)/m
                x2 = (rekt$y[3]-b)/m
                ############     left right bottom top
                APNTS = list(x=c(rekt$x[1],rekt$x[2], x1, x2),
                  y=c(y1, y2, rekt$y[1], rekt$y[3]) )
              }
            else
              {
                ###########  horizontal line
                  x1 = rekt$x[1]
                  x2 = rekt$x[2]

                  APNTS = list(x=c(rekt$x[1],rekt$x[2], NA, NA),
                    y=c(b, b, NA , NA) )
              }
          }
        else
          {
            ###########  vertical line
             APNTS = list(x=c(NA,NA, pnts$x1[i], pnts$x1[i]),
                    y=c(NA, NA,rekt$y[1],  rekt$y[3]  ) )
                  
          }


        dis =  (APNTS$x-pnts$x2[i])^2 +  (APNTS$y-pnts$y2[i])^2
        w = which.min(dis)
        PEES$x[i] = APNTS$x[w]
        PEES$y[i] = APNTS$y[w]
        
      }


    return(PEES)



  }
