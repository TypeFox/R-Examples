#' enaStructure --- performes strucutral analysis of the
#' network graph (see Borrett et al. 2007)
#' INPUT = network object
#' OUTPUT = list of structure statistics
#'
#' S. Borrett and M. Lau | March 2011
#' ---------------------------------------------------

enaStructure <- function(x = 'network object'){
                                        #Check for network class
  if (class(x) != 'network'){warning('x is not a network class object')}
  Flow <- t(as.matrix(x,attrname = 'flow')) #get flows
  A <- sign(Flow)   # get adjacency matrix
  sp <- structure.statistics(A)    # calls structure.statistics helper function
                                          #Output orientation
  orient <- get.orient()
  if (orient=='rc'){A <- t(A)}else{}
  return(list('A'=A,'ns'=sp))  # "A" is the adjacency matrix oriented
                                        # column to row and "sp" is a list of
                                        # structural network staistics
}
