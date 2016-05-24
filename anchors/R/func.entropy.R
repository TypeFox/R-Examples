#######################################################################
##
## Function: entropy helpers
## Author  : Jonathan Wand <wand(at)stanford.edu>
##           http://wand.stanford.edu
##
#######################################################################

  entropy.max <- function(n.cat) {
    p <- rep(1/n.cat,n.cat)
    return( -sum( p*log(p)))
  }
    
  entropy.calc <- function(p) {
    p <- p[ p > 0]
    -sum(p*log(p))
  }
