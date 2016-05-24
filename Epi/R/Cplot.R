Cplot <-
function( rates,
          age = as.numeric( rownames( rates ) ),
          per = as.numeric( colnames( rates ) ),
         grid = FALSE,
       c.grid = grid,
        ygrid = grid,
     col.grid = gray( 0.9 ),
        c.lim = NULL,
         ylim = range( rates[rates>0], na.rm=TRUE ),
           at = NULL,
       labels = paste( at ),
        c.lab = names( dimnames( rates ) )[2],
         ylab = deparse( substitute( rates ) ),
         type = "l",
          lwd = 2,
          lty = 1,
          col = par( "fg" ),
       log.ax = "y",
          las = 1,
        xannx = 1/20,
          ann = FALSE,
      cex.ann = 0.8,
       a.thin = seq( 1, length( age ), 2 ),
          ...
          )
{
# First convert the age-period table to an age-cohort table
rt <- as.table( rates )
dimnames( rt ) <- list( age = age, per = per )
rtf <- data.frame( rt )
rtf$age <- as.numeric( as.character( rtf$age ) )
rtf$per <- as.numeric( as.character( rtf$per ) )
ac <- tapply( rtf$Freq, list( rtf$age, rtf$per-rtf$age ), mean )
coh <- as.numeric( colnames( ac ) )
if( is.null( c.lim ) )
  c.lim <- range( coh, na.rm=TRUE ) + c(0,diff( range( coh ) )/30) * ann

# Plot the frame
if( ann ) c.lim <- c.lim  - c(diff( range( coh ) ) * xannx,0)
matplot( coh, t(ac), type="n",
         xlim=c.lim, ylim=ylim, xlab=c.lab, ylab=ylab,
         log=log.ax, las=las, yaxt=if( !is.null( at ) ) "n" else "s" )
if( !is.null( at ) ) axis( side=2, at=at, labels=labels, yaxt="s", las=las ) 
# and the grid if required
if( !missing( c.grid ) | !missing( grid ) )
  {
  if( is.logical( c.grid ) & c.grid[1] ) c.grid <- nice( coh, log=par("xlog") )
  abline( v=c.grid, col=col.grid )
  }
if( !missing( ygrid ) | !missing( grid ) )
  {
  if( is.logical( ygrid ) & ygrid[1] ) ygrid <- nice( rates[!is.na(rates)], log=par("ylog") )
  abline( h=ygrid, col=col.grid )
  }
box()
# then the curves
matlines( coh, t(ac), lwd=lwd, lty=lty, col=col, type=type, ... )
# annotate them if required (every second by default )
if( ann )
  {
   nr <- nrow( ac )
   nc <- ncol( ac )
   # Find the cohorts for the last rates in each age-class
   c.end <- rev( per )[1] - age
   text( c.end[a.thin], rates[,ncol(rates)][a.thin],
         paste( "", age[a.thin] ), adj=c(0,0.5),
         cex=cex.ann, col=if( length(col)==1 ) col else col[a.thin] ) 
  }
}
