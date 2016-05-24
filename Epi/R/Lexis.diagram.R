Lexis.diagram <- function( age = c( 0, 60),
                          alab = "Age",
                          date = c( 1940, 2000 ),
                          dlab = "Calendar time",
                           int = 5,
                       lab.int = 2*int,
                      col.life = "black",
                      lwd.life = 2,
                      age.grid = TRUE,
                     date.grid = TRUE,
                      coh.grid = FALSE,
                      col.grid = gray(0.7),
                      lwd.grid = 1,
                           las = 1,
                    entry.date = NA,
                     entry.age = NA,
                     exit.date = NA,
                      exit.age = NA,
                     risk.time = NA,
                    birth.date = NA,
                          fail = NA,
                      cex.fail = 1.1,
                      pch.fail = c(NA,16),
                      col.fail = rep( col.life, 2 ),
                          data = NULL,
                           ... )
{
# Function to plot a Lexis-diagram
#
# BxC, 2002, revsions in 2005

    ## Get variables from data argument, if supplied, or from parent
    ## frame if not.
    entry.date <- eval(substitute(entry.date), data)
    entry.age  <- eval(substitute(entry.age ), data)
    exit.date  <- eval(substitute(exit.date ), data)
    exit.age   <- eval(substitute(exit.age  ), data)
    risk.time  <- eval(substitute(birth.date), data)
    birth.date <- eval(substitute(birth.date), data)
    fail       <- eval(substitute(fail      ), data)

# First expand intervals to both dimensions
#
    int[1:2] <- c(    int,    int)[1:2]
lab.int[1:2] <- c(lab.int,lab.int)[1:2]

# Plot the diagram
#
plot( NA,
      xlim=date, xaxt="n", xaxs="i", xlab=dlab,
      ylim=age,  yaxt="n", yaxs="i", ylab=alab, ... )
axis( side=1, at=seq( date[1], date[2], lab.int[2] ), las=las )
axis( side=2, at=seq(  age[1],  age[2], lab.int[1] ), las=las )
box( col="white" ) # par("fg") )

# Then the required grids
#
if (  age.grid ) { abline( h = seq(  age[1],  age[2], int[1] ),
                           col=col.grid, lwd=lwd.grid ) }
if ( date.grid ) { abline( v = seq( date[1], date[2], int[2] ),
                           col=col.grid, lwd=lwd.grid ) }
ages  <- seq(  age[1],  age[2], min( int ) )
dates <- seq( date[1], date[2], min( int ) )
if ( coh.grid ) { segments( rep( date[1], length( ages ) ),
                            ages,
                            pmin( date[1] + ( age[2] - ages ), date[2] ),
                            pmin( ages + ( date[2] - date[1] ), age[2] ),
                            col=col.grid, lwd=lwd.grid )
                  segments( dates,
                            rep( age[1], length( dates ) ),
                            pmin( dates + ( age[2] - age[1] ), date[2] ),
                            pmin( age[1] + ( date[2] - dates ), age[2] ),
                            col=col.grid, lwd=lwd.grid )
                 }

# Check if data for lifelines is supplied and plot lifelines if so
#
if( sum( !is.na( list( entry.date,
  	                entry.age,
  	                 exit.date,
                         exit.age,
  	                birth.date,
                         risk.time ) ) ) > 2 )
  {
  LL <- Lexis.lines( entry.date = entry.date,
                      exit.date = exit.date,
                     birth.date = birth.date,
                      entry.age = entry.age,
                       exit.age = exit.age,
                      risk.time = risk.time,
                       col.life = col.life,
                       lwd.life = lwd.life,
                           fail = fail,
                       cex.fail = cex.fail,
                       pch.fail = pch.fail,
                       col.fail = col.fail,
                           data = data )
  invisible( LL )
  }
}









