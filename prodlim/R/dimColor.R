##' This function calls first \code{\link{col2rgb}} on a color name and then 
##' uses \code{\link{rgb}} to adjust the intensity of the result.
##' 
##' @title Dim a given color to a specified density
##' @param col Color name or number passed to \code{\link{col2rgb}}.
##' @param density Integer value passed as alpha coefficient to
##' \code{\link{rgb}} between 0 and 255
##' @return A character vector with the color code. See \code{rgb} for details.
##' @seealso rgb col2rgb
##' @examples 
##' dimColor(2,33)
##' dimColor("green",133)
##' @export 
##' @author Thomas A. Gerds <tag@@biostat.ku.dk>
dimColor <- function(col,density=55){
    ccrgb=as.list(grDevices::col2rgb(col,alpha=TRUE))
    names(ccrgb) <- c("red","green","blue","alpha")
    ccrgb$alpha=density
    do.call(grDevices::rgb,c(ccrgb,list(max=255)))
}
