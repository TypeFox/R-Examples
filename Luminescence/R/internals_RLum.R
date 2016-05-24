####################################################################################################
##                     INTERNAL HELPER FUNCTIONS                                                  ##
####################################################################################################

#' Create unique object identifier
#'
#' This function creates a unique object identifier. This function has no further parameters, as
#' an md5 hash is calcualted from the system timestamp adding an abitrary value using \code{runif}
#'
#' @return
#' Returns a UID
#'
#' @examples
#' .create_UID()
#'
#' @noRd
.create_UID <- function(){
  digest::digest(object = paste(as.character(Sys.time()), runif(1)), algo = "md5")
}
