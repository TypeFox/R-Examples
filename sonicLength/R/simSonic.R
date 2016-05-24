simSonic <- function(theta , phi )
{
  ## Purpose: simulate sonicLength data
  ## ----------------------------------------------------------------------
  ## Arguments: theta - a vector poisson parameters
  ##            phi - a vector (or a list with a phi component) whose
  ##            names are like 'rep len', where 'rep' defines the
  ##            replicate and 'len' defines the sonicant length
  ## ----------------------------------------------------------------------
  ## Author: Charles Berry, Date:  9 Jun 2011, 12:54
  if (is.list(phi)) phi <- phi$phi
  locnames <-
    if( (length(nm <- names(theta)) && !any(duplicated(nm)))) nm else seq( along = theta )
  loc.cnt <- rpois(length( theta ), theta )
  nloc <- sum(loc.cnt)
  phi.pos <- suppressWarnings(sample( names(phi), nloc, prob=phi, replace=TRUE))
  phi.len <- as.numeric(sub("^.*[ ]","", phi.pos ))
  phi.rep <- as.numeric(sub("[ ].*$","", phi.pos))
  unique(data.frame(locations=rep(locnames,loc.cnt),lengths=phi.len,replicates=phi.rep))
}
