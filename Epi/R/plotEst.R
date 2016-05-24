get.ests <- function( ests, ... )
{
 # If a model object is supplied, extract the parameters and the
 # standard errors
 #
if( inherits( ests, c("glm","coxph","clogistic","gnlm","survreg") ) )
    ests <- ci.exp( ests, ... )
else
if( inherits( ests, c("lm","gls","lme","nls","polr",
                      "mer","MIresult","mipo") ) )
    ests <- ci.lin( ests, ... )[,-(2:4)]
ests
}

plotEst <-
function( ests,
             y = dim(ests)[1]:1,
           txt = rownames(ests),
        txtpos = y,
          ylim = range(y)-c(0.5,0),
          xlab = "",
          xtic = nice( ests[!is.na(ests)], log=xlog ),
          xlim = range( xtic ),
          xlog = FALSE,
           pch = 16,
           cex = 1,
           lwd = 2,
           col = "black",
       col.txt = "black",
      font.txt = 1,
     col.lines = col,
    col.points = col,
          vref = NULL,
          grid = FALSE,
      col.grid = gray(0.9),
   restore.par = TRUE,
           ...
         )
{
# Function to plot estimates from a model.
# Assumes that ests is a p by 3 matrix with estimate, lo and hi as
# columns OR a model object.

 # Extract the estimates if necessary
 #
ests <- get.ests( ests, ... )

 # Is it likley that we want a log-axis for the parameters?
 #
mult <- inherits( ests, c("glm","coxph","gnlm") )
if( missing(xlog) ) xlog <- mult

 # Create an empty plot in order to access the dimension so that
 # sufficient place can be made for the text in the margin
 #
plot.new()
mx <- max( strwidth( txt, units="in" ) )
oldpar <- par( mai=par("mai") + c(0,mx,0,0) )
if( restore.par ) on.exit( par( oldpar ) )

 # Set up the coordinate system witout advancing a frame
 #
plot.window( xlim = xlim, ylim = ylim, log = ifelse( xlog, "x", "") )

 # Draw a grid if requested
 #
if( !is.logical( grid ) )           abline( v = grid, col = col.grid )
if(  is.logical( grid ) & grid[1] ) abline( v = xtic, col = col.grid )

 # Draw a vertical reference line
 #
if( !missing( vref ) ) abline( v = vref )

 # Draw the estimates with c.i. s
 #
linesEst( ests, y, pch=pch, cex=cex, lwd=lwd,
          col.points=col.points, col.lines=col.lines )

 # Finally the x-axis and the annotation of the estimates
 #
axis( side=1, at=xtic )
mtext( side=1, xlab, line=par("mgp")[1], cex=par("cex")*par("cex.lab") )
axis( side=2, at=txtpos, labels=txt, las=1, lty=0, col.axis=col.txt )
invisible( oldpar )
  }

pointsEst <-
linesEst <-
function( ests,
             y = dim(ests)[1]:1,
           pch = 16,
           cex = 1,
           lwd = 2,
           col = "black",
     col.lines = col,
    col.points = col,
           ...
         )
{
# Function to add estimates from a model to a drawing.
# Assumes that ests is a p by 3 matrix with estimate, lo and hi as
# columns.

  # Extract the estimates if necessary
  #
ests <- get.ests( ests, ... )

  # Cut the confidence interval lines to fit inside the plot
  # before drawing them.
  #
xrng <- if( par("xlog") ) 10^par("usr")[1:2] else par("usr")[1:2]
segments( pmax(ests[, 2],xrng[1]), y, pmin(ests[, 3], xrng[2]), y,
          lwd = lwd, col=col.lines )

  # Then the point estimates on top of the lines.
  #
points( ests[, 1], y, pch = pch, cex = cex, col=col.points )
invisible()
}

nice <-
function( x,
        log = F,
       lpos = c(1,2,5), ... )
{
# Function to produce nice labels also for log-axes.
#
if( log )
  {
  fc <- floor( log10( min( x ) ) ):ceiling( log10( max( x ) ) )
  tick <- as.vector( outer( lpos, 10^fc, "*" ) )
  ft <- max( tick[tick<min(x)] )
  lt <- min( tick[tick>max(x)] )
  tick <- tick[tick>=ft & tick<=lt]
  if( length( tick ) < 4 & missing( lpos ) )
    tick <- nice( x, log=T, lpos=c(1:9) )
  if( length( tick ) > 10 & missing( lpos ) )
    tick <- nice( x, log=T, lpos=1 )
  tick
  }
else pretty( x, ... )
}
