################################################################################
#  
#  Check to see if basis object requested has been pre-computed and is available
#  to load from disk.
#
#  Currently only works for several small examples included with CRAN package;
#  will be extended to work with a number of larger examples which are too large
#  to include on CRAN.
#  
#  (C) wburr, July 2014
#
#  Licensed under GPL-2
#
################################################################################

checkSaved <- function(N, W, K) {

  # load list of saved objects
  slpSavedObjects <- NULL
  data("slpSavedObjects", envir = environment())

  check <- unlist(lapply(slpSavedObjects, FUN = function(x) { 
                    (x[1] == N & x[2] == W & x[3] == K)  }))
  
  if(any(check)) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}
