`get2Drayblox` <-
function(x1, y1, x2, y2, xo, yo, NODES=FALSE, PLOT=FALSE)
  {
#########    find the block indecies and lengths of rays in block model
    
    if(missing(PLOT) ) { PLOT = FALSE }
    if(missing(NODES) ) { NODES = FALSE }

    
####  xo, yo are nodes not block boundaries
    
    A = list(x=c(x1,x2), y=c(y1, y2))

    dx = diff(xo)[1]
    dy = diff(yo)[1]


  ####  M = meshgrid(xo, yo)

    if(NODES)
      {
        LX = c(xo-(dx/2), xo[length(xo)]+(dx/2))
        LY = c(yo-(dy/2), yo[length(yo)]+(dy/2))
      }
    else
      {
        LX = xo
        LY =yo
      }

    ####  get slope and intercept
    
    zm = lm(A$y ~ A$x)
    slope=as.numeric(zm$coefficients[2])
    intercept=as.numeric(zm$coefficients[1])

    if(is.na(slope))
      {
        WHY = c(min(A$y), c(LY[LY>=min(A$y) & LY<=max(A$y) ]), max(A$y))

        WHYX =  rep(A$x, length=length(WHY))
        EX = NULL
        EXY = NULL

      }
    else
      {
        EX = c(min(A$x),  LX[LX>=min(A$x) & LX<= max(A$x) ], max(A$x))
        EXY = intercept + EX* slope

        WHY = c(LY[LY>=min(A$y) & LY<=max(A$y) ])
        WHYX = (WHY-intercept)/slope

      }
    xp = c(EX, WHYX)
    yp = c(EXY, WHY)

    o = order(yp)

    kx = length(xp)
    ky = length(yp)

    

    nodes = list(x = xp[o], y = yp[o])

    if(identical(y2, y1))
      {
        o = order(xp)
        nodes = list(x = xp[o], y = yp[o])
      }
      
    if(y2<y1) { nodes = list(x = rev(xp[o]), y = rev(yp[o]))  }

    mids = list(x= (nodes$x[2:kx]+nodes$x[1:(kx-1)])/2,y= (nodes$y[2:ky]+nodes$y[1:(ky-1)])/2 )


    lengs = sqrt(  (nodes$x[2:kx]-nodes$x[1:(kx-1)])^2+(nodes$y[2:ky]-nodes$y[1:(ky-1)])^2)

    ix = findInterval(mids$x, LX)
    iy = findInterval(mids$y, LY)

    if(PLOT==TRUE)
      {
        blues = shade.col(100, acol=c(1,1,1) , bcol= as.vector(col2rgb("lightblue")/255))

        plot(xo, yo, type='n')
        segments(LX, rep(min(LY), length=length(LX)), LX, rep(max(LY), length=length(LX)))
        segments( rep(min( LX), length=length(LY) ) , LY,   rep(max( LX), length=length(LY) ) , LY )

        thecolors = 1+99*(lengs-min(lengs))/(max(lengs)-min(lengs))
        rect( LX[ix], LY[iy], LX[ix+1], LY[iy+1],col=blues[thecolors])

        segments(A$x[1], A$y[1], A$x[2], A$y[2])


      }

    return(list(ix=ix, iy=iy, lengs=lengs, mids=mids, nodes=nodes, LX=LX, LY=LY))
  }

