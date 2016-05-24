##
## Code originally from Frank Harrell's 'Hmisc' library: 
##   http://biostat.mc.vanderbilt.edu/twiki/bin/view/Main/Hmisc
## Copied with permission on 2007-08-04
##

## Determine if variable is a date, time, or date/time variable in R
## or S-Plus.  The following 2 functions are used by describe.vector
## timeUsed assumes is date/time combination variable and has no NAs

#' @include AFirst_lib.R

testDateTime <- function(x, what=c('either','both','timeVaries'))
{
  what <- match.arg(what)
  cl <- class(x)  # was oldClass 22jun03
  if(!length(cl))
    return(FALSE)

  dc <- if(.R.)
          c('Date', 'POSIXt','POSIXct','dates','times','chron')
        else
          c('timeDate','date','dates','times','chron')
  
  dtc <- if(.R.)
           c('POSIXt','POSIXct','chron')
         else
           c('timeDate','chron')
  
  switch(what,
         either = any(cl %in% dc),
         both   = any(cl %in% dtc),
         timeVaries = {
           if('chron' %in% cl || 'Date' %in% cl || !.R.) { 
             ## chron or S+ timeDate
             y <- as.numeric(x)
             length(unique(round(y - floor(y),13))) > 1
           }
           else if(.R.)
             length(unique(format(x,'%H%M%S'))) > 1
           else
             FALSE
         })
}
