upstreamconnectivity <- function(x, ...) UseMethod("upstreamconnectivity")


upstreamconnectivity.rivernet <- function(x,crit.reach,crit.node,thresh.length=0,...)
{
  rivernet <- x
  if ( length(rivernet$paths) == 0 )
  {
    rivernet <- analyze(rivernet)
    if ( length(rivernet$paths) == 0 ) return(NA)  # unable to analyze network structure
  }
  paths.reachable <- list()
  firstorder.reachable <- rep(FALSE,length(rivernet$paths))
  for ( i in 1:length(rivernet$paths) )
  {
    path <- rivernet$paths[[i]]
    path.reachable <- numeric(0)
    dist.bad <- 0
    for ( j in length(path):1 )
    {
      reach <- path[j]
      if ( crit.reach[reach] ) 
      {
        path.reachable <- c(reach,path.reachable)
        dist.bad <- 0
        if ( rivernet$attrib.reach$streamorder[reach]==1 ) firstorder.reachable[i] <- TRUE
        if ( rivernet$attrib.reach$downstream[reach] )
        {
          downstream.node <- rivernet$attrib.reach$node_start[reach]
        }
        else
        {
          downstream.node <- rivernet$attrib.reach$node_end[reach]
        }
        if ( !crit.node[downstream.node] ) break
      }
      else
      {
        dist.bad <- dist.bad + rivernet$attrib.reach$length[reach]
        if ( dist.bad > thresh.length ) break
      }
    }
    paths.reachable[[i]] <- path.reachable
  }
  return(list(paths.reachable      = paths.reachable,
              firstorder.reachable = firstorder.reachable,
              fract.firstorder.reachable = sum(firstorder.reachable)/length(paths.reachable)))
}
