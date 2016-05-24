"put.points.demo" <-
function( x=NULL, y=NULL, lsline=TRUE) {

  old.par <- par(no.readonly=T)
  on.exit(par(old.par))

  options(locatorBell=FALSE)

  mode='add'

  layout( matrix( c(2,1), nrow=1), widths=c(3,1) )


  repeat {

    ## right panel
    par(mar=c(0,0,0,0),usr=c(0,1,0,1))
    frame()
    box()
    abline(h=c(0.8,0.6))
    text( rep(0.5, 5), c(0.9, 0.725, 0.525, 0.325, 0.125),
         lab=c('End','LS Line','Add Point','Delete Point','Move Point') )
    lines( c(0.25,0.25,0.75,0.75,0.25), c(0.85,0.95,0.95,0.85,0.85) )

    points( rep(0.5,4), c(0.675,0.475,0.275,0.075),
           pch=c( ifelse(lsline,7,0),
             ifelse(mode=='add', 16, 1),
             ifelse(mode=='del', 16, 1),
             ifelse(mode=='mov', 16, 1)), cex=2.5 )

    ## left panel
    par(mar=c(5,4,4,1)+0.1)
    if(length(x) == 0) {
      plot(5,5,type='n', xlim=c(0,10), ylim=c(0,10),
           xlab='x', ylab='y')
    } else {
      plot(x,y, xlim=range(x,0,10), ylim=range(y,0,10),
           xlab='x', ylab='y')
      if( lsline && length(x) > 1 ){
        tmp.fit <- lm(y~x)
        abline(tmp.fit)
        title( paste( "r =", round(cor(x,y),2),
                      "r^2 =", round(cor(x,y)^2,2),
                     "\nSlope =",round(coef(tmp.fit)[2],2),
                      "Intercept =",round(coef(tmp.fit)[1],2)) )
      } else {
        title( paste( "r =", round(cor(x,y),4),
                      "r^2 =", round(cor(x,y)^2,4)))
      }
    }

    # get point
    pnt <- locator(1)

    if (pnt$x > par('usr')[2]) { ## clicked in left panel

#      pnt2 <- cnvrt.coords(pnt)$fig
      pnt2 <- list()
      pnt2$y <- grconvertY(pnt$y, to='nfc')

      if( pnt2$y > .8 ){
        break
      }
      if( pnt2$y > .6 ){
        lsline <- !lsline
        next
      }
      if( pnt2$y > .4 ){
        mode <- 'add'
        next
      }
      if( pnt2$y > .2 ){
        mode <- 'del'
        next
      }
      mode <- 'mov'
      next

    } else { ## clicked in right panel
      if( mode=='add' ) {
        x <- c(x,pnt$x)
        y <- c(y,pnt$y)
        next
      }
      if( mode=='del' ) {
        min.i <- which.min( (x-pnt$x)^2+(y-pnt$y)^2 )
        x <- x[-min.i]
        y <- y[-min.i]
        next
      }
      if( mode=='mov' ) {
        mov.i <- which.min( (x-pnt$x)^2+(y-pnt$y)^2 )
        points( x[mov.i], y[mov.i], pch=16 )
        pnt <- locator(1)
        x[mov.i] <- pnt$x
        y[mov.i] <- pnt$y
        next
      }
    }

  } ## end repeat

}

