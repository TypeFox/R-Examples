Aplot <-
function( rates,
            age = as.numeric( dimnames( rates )[[1]] ),
            per = as.numeric( dimnames( rates )[[2]] ),
           grid = FALSE,
         a.grid = grid,
          ygrid = grid,
       col.grid = gray( 0.9 ),
          a.lim = range( age, na.rm=TRUE ),
           ylim = range( rates[rates>0], na.rm=TRUE ),
             at = NULL,
         labels = paste( at ),
          a.lab = names( dimnames( rates ) )[1],
           ylab = deparse( substitute( rates ) ),
           type = "l",
            lwd = 2,
            lty = 1,
            col = par( "fg" ),
         log.ax = "y",
            las = 1,
          c.col = col,
          p.col = col,
          c.ann = FALSE,
          p.ann = FALSE,
          xannx = 1/20,
        cex.ann = 0.8,
         c.thin = seq( 2, length( age ) + length( per ) - 1, 2 ),
         p.thin = seq( 1, length( per ), 2 ),
        p.lines = TRUE,
        c.lines = !p.lines,
            ...   # arguments passed on to matlines()
            )
{

# Plot the frame
if( p.ann ) a.lim <- a.lim + c(0,diff( range( age ) ) * xannx)
if( c.ann ) a.lim <- a.lim - c(  diff( range( age ) ) * xannx,0)
matplot( age, rates, type="n",
         xlim=a.lim, ylim=ylim, xlab=a.lab, ylab=ylab,
         log=log.ax, las=las, yaxt=if( !is.null( at ) ) "n" else "s" )
if( !is.null( at ) ) axis( side=2, at=at, labels=labels, yaxt="s", las=las )
# and the grid if required:
if( !missing( a.grid ) | !missing( grid ) )
  {
  if( is.logical( a.grid ) & a.grid[1] )
    a.grid <- nice( age, log=par("xlog") )
  abline( v=a.grid, col=col.grid )
  }
if( !missing( ygrid ) | !missing( ygrid ) )
    {
  if( is.logical( ygrid ) & ygrid[1] ) ygrid <- nice( rates[!is.na(rates)], log=par("ylog") )
  abline( h=ygrid, col=col.grid )
  }
box()
# What lines were required?
if( !missing( c.lines ) & missing( p.lines ) ) p.lines <- !c.lines

# Period curves:
if( p.lines ){
matlines( age, rates, type=type, lwd=lwd, lty=lty, col=p.col, ... )
# annotate them if required (every second by default):
if( p.ann )
  {
   nr <- nrow( rates )
   nc <- ncol( rates )
   text( rep( age[nr], nc )[p.thin], rates[nr,][p.thin],
         paste( "", per[p.thin] ), adj=c(0,0.5), cex=cex.ann,
         col=if( length(p.col)==1 ) p.col else p.col[p.thin] )  
  }
}

# Cohort curves:
if( c.lines ){
# First convert the age-period table to an age-cohort frame
rt <- as.table( rates )
dimnames( rt ) <- list( age = age, per = per )
rtf <- data.frame( rt )
rtf$age <- as.numeric( as.character( rtf$age ) )
rtf$per <- as.numeric( as.character( rtf$per ) )
ac <- tapply( rtf$Freq, list( rtf$age, rtf$per-rtf$age ), mean )
matlines( age, ac, type=type, lwd=lwd, lty=lty, col=c.col, ... )
# annotate them if required (every other by default):
if( c.ann )
  {
   nr <- nrow( rt )
   nc <- ncol( rt )
   # Find the ages, cohorts and rates where the cohort curves starts
   a.min <- c( rev( age ), rep( age[1], nc-1 ) )
   p.min <- c( rep( per[1], nr-1 ), per )
   c.min <- p.min - a.min           
   r.min <- c(rt[nr:1,1],rt[1,2:nc])
   text( a.min[c.thin], r.min[c.thin], paste( "", c.min[c.thin] ),
         adj=c(1,0.5), cex=cex.ann,
         col=if( length(c.col)==1 ) c.col else c.col[c.thin] )
  }
}

}
