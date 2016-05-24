#######################################################################
##
## Function: minimum.entropy()
## Author  : Jonathan Wand <wand@stanford.edu>
##
#######################################################################
minimum.entropy <- function( obj , debug=0 ) {
  sa <- summary.anchors.rank.type( obj )
  me <- minimum.entropy.calc(sa, debug)
  
  tt  <- sa$interval
  idx <- tt$from != tt$to 

  ## find cases in 'me' that are applicable to 'tt'
  if (!is.null(me)) {

    if (debug > 0) {
      print(tt)
      print(me)
    }
    fidx<- names(me) %in% row.names(tt)[idx]
    ## and do the replacement
    tt$from[idx] <- me[fidx]
  }
  
  ## don't need 'to' column
  tt <- tt[,-ncol(tt)]

  cc <- colnames(tt)
  cc[length(cc)] <- "MinEnt"
  colnames(tt) <- cc

  tt <- as.data.frame(tt)
#  class(tt) <- "minium.entropy"
  return(tt)
}


