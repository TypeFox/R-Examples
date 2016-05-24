#' Color fade out
#' 
#' Make existing colors fade away to white
#' 
#' @param col Vector of color names (\code{\link{colors}}), hexadecimal or integer that can be interpreted by \code{\link{col2rgb}}
#' @param fade Level of fading towards target. between 0 (target) and 1 (col). Can also be a vector. DEFAULT: 0.3
#' @param target Target color that should be faded into. DEFAULT: "white"
#' @param \dots Further arguments passed to \code{\link{colorRamp}}

#' @return character matirx with hexadecimal color codes.
#' @author Berry Boessenkool, \email{berry-b@@gmx.de}, Feb 2016
#' @seealso \code{\link{addAlpha}}, \code{\link{colorRamp}}, \code{\link{colors}}
#' @keywords dplot color
#' @export
#' @examples
#' 
#' plot(1:11, pch=16, cex=3, col=addFade(2, 10:0/10))
#' plot(1:11, pch=16, cex=3, col=addFade(2, 10:0/10, target="blue"))
#' plot(1:11, pch=16, cex=3, col=addFade(2, 10:0/10, target=3:4))
#' plot(1:21, pch=16, cex=3, col=addFade(2:3, 10:0/10))
#' plot(1:21, pch=16, cex=3, col=addFade(2:3, 10:0/10, target=4:5))
#' NewColors <- addFade(c("red","blue","yellow","green", "purple"), 0:200/200)
#' plot(runif(1000), col=NewColors, pch=16, cex=2)
#' 
addFade <- function(
  col, 
  fade=0.3,
  target="white", 
  ... 
  )
{
if(any(fade<0 | fade>1)) stop("fade must be between 0 and 1, not ",
                                 pastec(fade[fade<0|fade>1]))
cR <-  function(fade, col, target)
  {
  x <- colorRamp(c(col, target))(1-fade)   # , ...
  x2 <- if (ncol(x) == 4L)
             rgb(x[, 1L], x[, 2L], x[, 3L], x[, 4L], maxColorValue = 255)
        else rgb(x[, 1L], x[, 2L], x[, 3L], maxColorValue = 255)
  x2
  }
# output <- sapply(fade, FUN=cR)
output <- sapply(col, FUN=function(co) sapply(fade, FUN=cR, col=co, target=target))
return(output)
}
