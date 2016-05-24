##
## Code originally from Frank Harrell's 'Hmisc' library: 
##   http://biostat.mc.vanderbilt.edu/twiki/bin/view/Main/Hmisc
## Copied with permission on 2007-08-04
##

## $Id: AFirst.lib.s,v 1.6 2005/09/26 15:44:17 dupontct Exp $
under.unix <- !(version$os=='Microsoft Windows' ||
                version$os=='Win32' || version$os=='mingw32')

.R.   <- TRUE
.SV4. <- FALSE

.noGenenerics <- TRUE  # faster loading as new methods not used

if(!exists('existsFunction')) {
  existsFunction <- function(...) exists(..., mode='function')
}

if(.R.) {  
  ## create some function definitions just to avoid R CMD CHECK warnings
  timeDate <- function(...) stop("Not Implemented")
  dates <- function(...) stop("Not Implemented")
}

