#' environCentrality --- calculates the centrality of 
#' flow network environs
#' INPUT = environ matrix
#' OUTPUT = in-going, out-going and average centralities
#' 
#' M. Lau | July 2011
#' ---------------------------------------------------

environCentrality <- function(x='matrix'){
  if (class(x) != 'matrix'){warning('x is not a matrix class object')}
  ECin <- rowSums(x)/sum(rowSums(x))
  ECout <- colSums(x)/sum(rowSums(x))
  AEC <- (ECin + ECout)/2
  names(ECin) <- rownames(x)
  names(ECout) <- rownames(x)
  names(AEC) <- rownames(x)
  return(list('ECin'=ECin,'ECout'=ECout,'AEC'=AEC))
}
