#' scc --- find the strongly connected component
#' INPUT = an adjacency matrix
#' OUTPUT = list of membership and values
#' S. Borrett | July 2011
#' ------------------------------------

scc <- function(A="adjacency"){
                                        #Check for network class
  if (class(A) != 'matrix'){warning('A is not a matrix class object')}
  n <- dim(A)[1]
  c <- component.dist(A) # finds strong components in A (from sna package)
  no.scc <- length(c$csize)  # numer of scc
  j <- which(c$csize>1)  # finds scc > 1
  no.scc.big <- length(j)    # number of scc > 1
  pscc <- sum(c$csize[j])/n  # percent of nodes participating in a scc
  sp <- c("no.scc"=no.scc,"no.scc.big"=no.scc.big,"pscc"=pscc)
  y <- list("sp"=sp,"membership"=c$membership,"scc.id"=j-1)
  return(y)
}
