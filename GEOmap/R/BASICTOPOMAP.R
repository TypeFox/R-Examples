`BASICTOPOMAP` <-
function(xo, yo, DOIMG, DOCONT, UZ, AZ, IZ, perim, PLAT, PLON, PROJ=PROJ, pnts=NULL, GRIDcol=NULL)
  {
    if(missing(GRIDcol) ) { GRIDcol=NULL }
    if(missing(pnts) ) { pnts=NULL }


     blues = RPMG::shade.col(100, acol=as.vector(col2rgb("darkblue")/255)   , bcol= as.vector(col2rgb("paleturquoise")/255))

    plot(range(xo), range(yo),  type='n' , asp=TRUE , axes=FALSE, xlab="", ylab="")

    if(DOIMG==TRUE)
      {
        image(x=xo, y=yo,   z=UZ, col=(blues), add=TRUE )
        image(x=xo, y=yo,   z=AZ, col=terrain.colors(100), add=TRUE )
      }

    if(DOCONT==TRUE)
      {
        contour(x=xo, y=yo,   z=IZ$z, add=TRUE  )
        contour(x=xo, y=yo,   z=IZ$z, add=TRUE, levels=pretty(c(0, -1000) , n=5) )
      }

    lines(perim, col='red')

    antipolygon(perim$x, perim$y, col='white') 

    addLLXY(PLAT,  PLON,  LABS=TRUE, GRIDcol=GRIDcol  , TICS=c(.1,.1) , PROJ=PROJ )

    
   if(!is.null(pnts))
      {

        P = GLOB.XY(pnts$lat , pnts$lon, PROJ)
        points(P$x, P$y, col=pnts$PCOL, pch=pnts$PCH)
      }


  }

