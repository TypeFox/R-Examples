plot.MCmcmc <-
function( x,
        axlim = range( attr(x,"data")$y, na.rm = TRUE ),
        wh.cmp,
     lwd.line = c(3,1), col.line = rep("black",2), lty.line=rep(1,2),
          eqn = TRUE, digits = 2,
         grid = FALSE, col.grid=gray(0.8),
       points = FALSE,
      col.pts = "black", pch.pts = 16, cex.pts = 0.8,
          ... )
{
Got.coda <- require( coda )
if( !Got.coda )
  stop( "Using the plot.MCmcmc function requires that\n",
        "the package 'coda' is installed.\n",
        "All installed packages are shown if you type 'library()'." )

# Extract the conversion formulae
conv.array <- summary.MCmcmc(x)$conv.array

# A function to plot the single panels --- argument definition is
# limited, since the functions is defined inside plot.MCmcmc,
# and therefore uses the parameters from the calling environment,
# which is the frame of plot.MCmcmc
conv.pair <-
function( parm, Mnames, pts )
# Assumes that alpha and beta refer to the relation
#  Mnames[2] = alpha + beta*Mnames[1], which is printed
# in the upper left corner.
{
alpha <- parm[1]
 beta <- parm[2]
   sd <- parm[3]
if( is.logical(grid) ) if( grid ) grid <- pretty( axlim )
if( is.numeric(grid) ) abline( h=grid, v=grid, col=col.grid )
if( eqn )
  {
  y.x <- paste( Mnames[2], "=\n",
                formatC( alpha, format="f", digits=digits ), "+",
                formatC( beta, format="f", digits=digits ),
                Mnames[1], "\n  (",
                formatC( sd, format="f", digits=digits ), ")" )
  x.y <- paste( Mnames[1], "=\n",
                formatC( -alpha/beta, format="f", digits=digits ), "+",
                formatC( 1/beta, format="f", digits=digits ),
                Mnames[2], "\n  (",
                formatC( sd/beta, format="f", digits=digits ), ")" )
  # Heights and widths of the equations
  wul <- strwidth ( y.x )
  hul <- strheight( y.x )
  wlr <- strwidth ( x.y )
  hlr <- strheight( x.y )
  if( is.numeric(grid) )
    {
    rect( par("usr")[1], par("usr")[4],
          par("usr")[1]+wul+0.2*hul, par("usr")[4]-1.2*hul,
          border=col.grid, col="white" )
    rect( par("usr")[2], par("usr")[3],
          par("usr")[2]-wlr-0.2*hlr, par("usr")[3]+1.2*hlr,
          border=col.grid, col="white" )
    }
  text( par("usr")[1]+0.1*hul, par("usr")[4]-0.1*hul, y.x, adj=c(0,1) )
  text( par("usr")[2]-0.1*hlr, par("usr")[3]+0.1*hlr, x.y, adj=c(1,0) )
  }
if( !is.null(pts) ) points( pts[,1], pts[,2],
                            pch=pch.pts, cex=cex.pts, col=col.pts )
abline( alpha        , beta, lwd=lwd.line[1], col=col.line[1], lty=lty.line[1] )
if( lwd.line[2] > 0 )
{
abline( alpha-2.00*sd, beta, lwd=lwd.line[2], col=col.line[2], lty=lty.line[2] )
abline( alpha+2.00*sd, beta, lwd=lwd.line[2], col=col.line[2], lty=lty.line[2] )
}
box()
}

# Duplication to allow a single line specification
col.line <- rep( col.line, 2 )
lty.line <- rep( lty.line, 2 )
lwd.line <- rep( lwd.line, 2 )

# Here is the substantive part of the function
#-------------------------------------------------------------------------------
# Names for plotting
Mnames <- attr( x, "methods" )
     if( missing(      wh.cmp ) ) wh <- 1:length(Mnames)
else if( is.numeric(   wh.cmp ) ) wh <- wh.cmp
else if( is.character( wh.cmp ) ) wh <- match( wh.cmp, Mnames )
Nm <- length( wh )

# Points to be plotted
if( is.logical( points ) ) points <- if( points ) "repl" else ""
if( points=="repl" )
  {
  pts <- with( attr( x, "data" ), tapply( y, list( interaction(item,repl), meth ), mean ) )
  }
if( points=="mean" )
  {
  pts <- with( attr( x, "data" ), tapply( y, list( item, meth ), mean ) )
  }

# Set up the panel layout for plotting the relations
oldpar <- par( mfrow=c(Nm-1,Nm-1), mar=rep(0,4), oma=c(4,4,1,1) )
for( ir in 2:Nm )
for( ic in 1:(Nm-1) )
   {
   if( ir<=ic ) plot( 0, 0, type="n", axes=FALSE, xlab="", ylab="" )
   if( ir>ic )
     {
     plot( NA, xlim=axlim, ylim=axlim, xlab="", ylab="", axes=FALSE )
     conv.pair( conv.array[wh[ic],wh[ir],],
                Mnames[c(wh[ic],wh[ir])],
                pts = if( points %in% c("repl","mean") )
                        cbind( pts[,wh[ic]], pts[,wh[ir]] ) )
     if( ic==1 )
       {
       axis( side=2 )
       mtext( Mnames[wh[ir]], line=2.5, side=2, outer=FALSE )
       }
     if( ir==Nm )
       {
       axis( side=1 )
       mtext( Mnames[wh[ic]], line=2.5, side=1, outer=FALSE )
       }
     }
    }
# on.exit( par(oldpar) )
invisible( NULL )
}
