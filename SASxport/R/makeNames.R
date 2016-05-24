##
## Code originally from Frank Harrell's 'Hmisc' library: 
##   http://biostat.mc.vanderbilt.edu/twiki/bin/view/Main/Hmisc
## Copied with permission on 2007-08-04
##

#' @include AFirst_lib.R

makeNames <- function(names, unique=FALSE, allow=NULL)
{
  ## Runs make.names with exceptions in vector allow
  ## By default, R 1.9 make.names is overridden to convert _ to . as
  ## with S-Plus and previous versions of R.  Specify allow='_' otherwise.
  if(!.R. & length(allow))
    stop('does not apply for S-Plus')
  n <- make.names(names, unique)
  if(!length(allow))
    n <- gsub('_', '.', n)
  n
}
