# utilities to reset par() and options() to 'factory defaults'
PARDEFAULTS <- list()

.onLoad <- function(libname, pkgname){
  tmpdev <- tempfile()
  on.exit(unlink(tmpdev))
  pdf(file = tmpdev)
  PARDEFAULTS <<- par(no.readonly=TRUE)
  dev.off()
  PARDEFAULTS$new <<- FALSE
  PARDEFAULTS <<- PARDEFAULTS[!names(PARDEFAULTS) %in% c("bg", "fin", "mai", "new", "pin", "plt", "ps")]
}



# set defaults at build time.
OPTIONSDEFAULTS <- options()

#' Reset graphical options in 'par' to factory defaults.
#' 
#' Reset the \code{\link[graphics]{par}} to R's defaults. 
#' 
#' Some of \code{par}'s settings are readonly. These are obviously not reset.
#' 
#' Settings stored in \code{\link[graphics]{par}} are device-dependent. In practice,
#' most settings in \code{par} are initially the same accross devices. Exceptions
#' we noted are: 
#' \itemize{
#' \item{\code{bg}: background color}
#' \item{\code{fin}: figure region dimensions}
#' \item{\code{mai}: margin size (inches)}
#' \item{\code{pin}: current plot dimensions (inches)}
#' \item{\code{plt}: coordinates of the plot region as fractions of the current figure region}
#' \item{\code{ps}: point size of text (but not symbos)}
#' }
#' Consequently, these options are currently not reset by calling \code{reset_par()}  
#' 
#' @seealso \code{\link{reset_options}}, \code{\link[graphics]{par}}
#' @export 
reset_par <- function(){
  par(PARDEFAULTS)
}

#' Reset general options in 'options' to factory defaults.
#' 
#' @seealso \code{\link{reset_par}}
#' @export 
reset_options <- function(){
  options(OPTIONSDEFAULTS)
}


