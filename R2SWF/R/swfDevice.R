#' SWF graphics device
#' 
#' This function opens a SWF device that produces Flash animation
#' in SWF format. Every time you call a high level plotting function
#' like \code{\link[graphics]{plot}()}, the movie will create a new
#' frame and draw following shapes on it.
#' 
#' @param file a character string giving the output SWF file
#' @param width the width of the device in inches
#' @param height the height of the device in inches
#' @param bg the background color of the SWF file
#' @param fg initial foreground color
#' @param frameRate how many frames to be played in 1 second
#' 
#' @export
#' 
#' @author Yixuan Qiu <\url{http://yixuan.cos.name/}>
#' 
#' @examples \dontrun{
#' ## A demonstration of K-means clustering, using animation package
#' if(require(animation)) {
#'     swf("kmeans.swf", frameRate = 1)
#'     kmeans.ani()
#'     dev.off()
#' }
#' 
#' ## Test built-in fonts in sysfonts package
#' swf("fonts.swf", 8, 8)
#' plot(1, type = "n")
#' 
#' par(family = "sans", cex = 2)
#' text(0.7, 1.3, "Sans-R", font = 1)
#' text(0.7, 1.1, "Sans-B", font = 2)
#' text(0.7, 0.9, "Sans-I", font = 3)
#' text(0.7, 0.7, "Sans-BI", font = 4)
#' 
#' par(family = "serif")
#' text(1.0, 1.3, "Serif-R", font = 1)
#' text(1.0, 1.1, "Serif-B", font = 2)
#' text(1.0, 0.9, "Serif-I", font = 3)
#' text(1.0, 0.7, "Serif-BI", font = 4)
#' 
#' par(family = "mono")
#' text(1.3, 1.3, "Mono-R", font = 1)
#' text(1.3, 1.1, "Mono-B", font = 2)
#' text(1.3, 0.9, "Mono-I", font = 3)
#' text(1.3, 0.7, "Mono-BI", font = 4)
#' 
#' dev.off()
#' }
swf <- function(file = "Rplots.swf", width = 7, height = 7, bg = "white",
                fg = "black", frameRate = 12)
{
    .Call("swfDevice", as.character(file),
          as.double(width), as.double(height),
          c(col2rgb(bg)), c(col2rgb(fg)),
          as.double(frameRate), PACKAGE='R2SWF');
    invisible(NULL);
}
