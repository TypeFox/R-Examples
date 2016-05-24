adjacentreaches <- function(x, ...) UseMethod("adjacentreaches")


adjacentreaches.rivernet <- function(x,crit.reach,thresh.length=0,...)
{
  extend.region <- function(ind,crit.reach,thresh.length,dist,
                            region,regions,attrib.reach)
  {
    if ( !is.na(regions[ind]) ) return(regions)
    if ( crit.reach[ind] ) d <- 0
    else                   d <- dist + attrib.reach$length[ind]
    if ( crit.reach[ind] | d <= thresh.length )
    {
      regions[ind] <- region
      r <- c(which(attrib.reach$node_start==attrib.reach$node_start[ind] |
                   attrib.reach$node_end  ==attrib.reach$node_start[ind] |
                   attrib.reach$node_start==attrib.reach$node_end[ind] |
                   attrib.reach$node_end  ==attrib.reach$node_end[ind] ))
      r <- r[-match(ind,r)]
      if ( length(r) > 0 )
      {
        while ( sum(is.na(regions[r])) > 0 )
        {
          i <- match(NA,regions[r])
          regions <- extend.region(r[i],crit.reach,thresh.length,d,
                                   region,regions,attrib.reach)
        }
      }
    }
    else
    {
      regions[ind] <- 0
    }
    
    return(regions)
  }
  
  rivernet <- x
  n.reach <- length(rivernet$reaches)
  if ( length(crit.reach) != n.reach ) return(NA)
  regions <- rep(NA,n.reach)
  region <- 0
  while ( sum(is.na(regions)) > 0 )
  {
    region <- region + 1
    ind <- match(NA,regions)
    regions <- extend.region(ind,crit.reach,thresh.length,0,
                             region,regions,rivernet$attrib.reach)
  }
  return(regions)
}



