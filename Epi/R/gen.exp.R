# Acronyms:
# dop: date of purchase
# dpt: dose per time
# amt: amount
# dur: duration

use.amt.dpt <-
function( purchase,
          push.max = Inf,
            breaks,
              lags = NULL,
           lag.dec = 1 )
{
do.call( "rbind",
lapply( split( purchase, purchase$id ),
        function(set)
        {
        np <- nrow(set)
        if( np==1 ) return( NULL )
        set <- set[order(set$dop),]
        # Compute length of exposure periods
        drug.dur  <- set$amt / set$dpt
        # Put the exposed period head to foot
        new.start <- min( set$dop ) + c(0,cumsum(drug.dur[-np]))
        # Move them out so that the start of a period is never earlier than
        # the dop
        exp.start <- new.start + cummax( pmax(set$dop-new.start,0) )
        # Compute the pushes
        push.one <- exp.start - set$dop
        # Revise them to the maximally acceptable
        push.adj <- pmin( push.one, push.max )
        # Revise the starting dates of exposure
        exp.start <- exp.start - push.one + push.adj
        # Revise the durations to be at most equal to differences between the
        # revised starting dates
        drug.dur  <- pmin( drug.dur, c(diff(exp.start),Inf) )
        # Compute the end of the intervals
        exp.end   <- exp.start + drug.dur
        # Intervals in the middle not covered by the drug exposures - note
        # also that we make a record for the last follow-date
        followed.by.gap <- c( exp.start[-1]-exp.end[-length(exp.end)] > 0, TRUE )
        # To facilitate
        dfR <- rbind( data.frame( id = set$id[1],
                                 dof = exp.start,
                                 dpt = set$dpt ),
                      data.frame( id = set$id[1],
                                 dof = exp.end[followed.by.gap],
                                 dpt = 0 ) )
        dfR <- dfR[order(dfR$dof),]
        # We now compute the cumulative dose at the end of the interval using
        # interval length and dpt:
        dfR$cum.amt <- with( dfR, cumsum( c(0, diff(dof)*dpt[-length(dpt)]) ) )
        return( dfR )
        } ) )
}

use.only.amt <-
function( purchase,
          pred.win = Inf,
            breaks,
              lags = NULL,
           lag.dec = 1 )
{
# Compute the cumulative dose at all purcase dates and at the last
# (unknown) future expiry date, computed based on previous
# consumption.  The resulting data frame has one more line per person
# than no. of purchases.
do.call( "rbind",
lapply( split( purchase, purchase$id ),
        function(set)
        {
        np <- nrow(set)
        if( np==1 ) return( NULL )
        set <- set[order(set$dop),]
        # The points to include in the calculation:
        # All dates after pred.win before last purchase,
        # but at least the last two purchase dates,
        wp <- ( set$dop > pmin( max(set$dop)-pred.win,
                               sort(set$dop,decreasing=TRUE)[2] ) )
        # Cumulative amount consumed at each dop
        cum.amt <- cumsum(c(0,set$amt))
        # Average slope to use to project the duration last purchase
        avg.slp <- diff(range(cum.amt[c(wp,FALSE)]))/
                   diff(range(set$dop[wp]))
        # Purchase dates and the date of last consumption
        dof <- c( set$dop, set$dop[np]+set$amt[np]/avg.slp )
        return( data.frame( id = set$id[1],
                           dof = dof,
                       cum.amt = cum.amt ) )
        } ) )
}

gen.exp <-
function( purchase, id="id", dop="dop", amt="amt", dpt="dpt",
                fu, doe="doe", dox="dox",
            breaks,
           use.dpt = ( dpt %in% names(purchase) ),
              lags = NULL,
          push.max = Inf,
          pred.win = Inf,
           lag.dec = 1 )
{
# Make sure that the data fames have the right column names
wh <- match( c(id,dop,amt), names(purchase) )
if( any( is.na(wh) ) ) stop("Wrong column names for the purchase data frame")
names( purchase )[wh] <- c("id","dop","amt")
wh <- match( c(id,doe,dox), names(fu) )
if( any( is.na(wh) ) ) stop("Wrong column names for the follow-up data frame")
names( fu )[wh] <- c("id","doe","dox")

if( use.dpt )
  {
  # This is to allow dpt to be entered as numerical scalar common for all
  if( is.numeric(dpt) )
    {
    if( length(dpt) > 1 ) stop( "If dpt is numeric it must have length 1" )
    purchase$dpt <- dpt
    }
  else
  names( purchase )[match(c(dpt),names(purchase))] <- "dpt"
  tmp.dfr <- use.amt.dpt( purchase,
                              lags = lags,
                          push.max = push.max,
                           lag.dec = lag.dec )
  }
else
  tmp.dfr <- use.only.amt( purchase,
                               lags = lags,
                           pred.win = pred.win,
                            lag.dec = lag.dec )


# Merge in the follow-up period for the persons
tmp.dfr <- merge( tmp.dfr, fu, all=T )

# Interpolate to find the cumulative doses at the dates in breaks
do.call( "rbind",
lapply( split( tmp.dfr, tmp.dfr$id ),
        function(set)
        {
        # All values of these are identical within each set (=person)
        doe <- set$doe[1]
        dox <- set$dox[1]
        # The first and last date of exposure according to the assumption
        doi <- min(set$dof)
        doc <- max(set$dof)
        # Get the breakpoints and the entry end exit dates
        breaks <- sort( unique( c(breaks,doe,dox) ) )
        xval   <- breaks[breaks>=doe & breaks<=dox]
        dfr    <- data.frame( id = set$id[1],
                             dof = xval )
        dfr$tfi  <- pmax(0,xval-doi)
        dfr$tfc  <- pmax(0,xval-doc)
        dfr$cdos <- approx( set$dof, set$cum.amt, xout=xval, rule=2 )$y
        for( lg in lags )
           dfr[,paste( "ldos",
                       formatC(lg,format="f",digits=lag.dec),
                       sep="." )] <-
           approx( set$dof, set$cum.amt, xout=xval-lg, rule=2 )$y
        dfr
        } ) )
}

