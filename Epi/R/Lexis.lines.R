Lexis.lines <-
function( entry.date = NA,
           exit.date = NA,
          birth.date = NA,
           entry.age = NA,
            exit.age = NA,
           risk.time = NA,
            col.life = "black",
            lwd.life = 2,
                fail = NA,
            cex.fail = 1,
            pch.fail = c(NA,16),
            col.fail = col.life,
                data = NULL
          )
{
    ## Get variables from data argument, if supplied, or from parent
    ## frame if not.
    entry.date <- eval(substitute(entry.date), data)
    entry.age  <- eval(substitute(entry.age ), data)
    exit.date  <- eval(substitute(exit.date ), data)
    exit.age   <- eval(substitute(exit.age  ), data)
    risk.time  <- eval(substitute(birth.date), data)
    birth.date <- eval(substitute(birth.date), data)
    fail       <- eval(substitute(fail      ), data)

# If fail is numeric make it logical
if( is.numeric( fail ) ) fail <- ( fail > 0 )

# Complete the information on lifelines
XX <- Life.lines( entry.date = entry.date,
                   entry.age = entry.age,
                   exit.date = exit.date,
                    exit.age = exit.age,
                   risk.time = risk.time,
                  birth.date = birth.date )

# Expand lwd.life/col.life/pch.fail/col.fail/cex.fail
#
Np <- nrow( XX )

if( length( col.life )==1 ) col.life <- rep( col.life, Np  ) else
if( length( col.life )!=length(fail) ) stop("col.life must have length 1 or length(fail)" )

if( length( lwd.life )==1 ) lwd.life <- rep( lwd.life, Np ) else
if( length( lwd.life )!=length(fail) ) stop("lwd.life must have length 1 or length(fail)" )

if( length( col.fail )==1 ) col.fail <- rep( col.fail, Np ) else {
if( length( col.fail )==2 ) col.fail <- col.fail[fail+1] }
if( length( col.fail )!=length(fail) ) stop("col.fail must have length 1,2 or length(fail)" )

if( length( pch.fail )==1 ) pch.fail <- rep( pch.fail, Np ) else
if( length( pch.fail )==2 ) pch.fail <- pch.fail[fail+1]
if( length( pch.fail )!=length(fail) ) stop("pch.fail must have length 1,2 or length(fail)" )

if( length( cex.fail )==1 ) cex.fail <- rep( cex.fail, Np ) else
if( length( cex.fail )==2 ) cex.fail <- cex.fail[fail+1]
if( length( cex.fail )!=length(fail) ) stop("cex.fail must have length 1,2 or length(fail)" )

# Was XX returned as a Date-object?
# If so make a numerical version i LL, otherwise just a copy.
#
if( attr( XX, "Date" ) )
  {
  LL <- data.frame( lapply( XX, unclass ) )
  LL[,c(1,3,5)] <- LL[,c(1,3,5)] / 365.25 + 1970
  LL[,c(2,4,6)] <- LL[,c(2,4,6)] / 365.25
  } else LL <- XX

# Find age and date ranges in the current plot.
#
date <- par( "usr" )[1:2]
age  <- par( "usr" )[3:4]

# Plot the lifelines
  segments( LL[,1], LL[,2], LL[,3], LL[,4],
            lwd=lwd.life, col=col.life )
# If there are any non-NAs for pch.fail then blank out the space
# where they go before plotting the symbols
  if( any( !is.na(pch.fail) ) )
  points( LL[!is.na(pch.fail),3], LL[!is.na(pch.fail),4],
          pch=16,
          col="white", #par()$bg,
          cex=cex.fail[!is.na(pch.fail)] )
  points( LL[,3], LL[,4],
          pch=pch.fail,
          col=col.fail,
          cex=cex.fail )

# Return the untouched version of the completed dataframe
#
invisible( data.frame( XX, fail=fail ) )
}
