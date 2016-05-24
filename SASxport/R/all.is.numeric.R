##
## Code originally from Frank Harrell's 'Hmisc' library: 
##   http://biostat.mc.vanderbilt.edu/twiki/bin/view/Main/Hmisc
## Copied with permission on 2007-08-04
##

all.is.numeric <- function(x, what=c('test','vector'), extras=c('.','NA'))
{
  what <- match.arg(what)
  old <- options(warn=-1)
  on.exit(options(old))
  ##.Options$warn <- -1  6Aug00
  x <- sub('[[:space:]]+$', '', x)
  x <- sub('^[[:space:]]+', '', x)
  xs <- x[x %nin% c('',extras)]
  isnum <- !any(is.na(as.numeric(xs)))
  if(what=='test')
    isnum
  else if(isnum)
    as.numeric(x)
  else x
}
