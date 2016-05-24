# Function from package 'fda' (c) 2014

rangechk <- function(rangeval)
{
#  check a range vector argument

#  last modified 29 September 2008 by Jim Ramsay

  nrangeval = length(rangeval)
  OK <- TRUE
  if (!is.numeric(rangeval))          OK <- FALSE
  if (!is.vector(rangeval))           OK <- FALSE
  if (nrangeval < 1 || nrangeval > 2) OK <- FALSE
  if (rangeval[1] >= rangeval[2])     OK <- FALSE
  return(OK)
}
