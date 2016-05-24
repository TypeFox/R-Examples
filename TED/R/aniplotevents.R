#' Generate a gif to visualise the event detection process
#' 
#' This function generates a gif file demonstrating how the event detection process is implemented.

#'
#' @param x a vector or a time series.
#' @param w a scalar specifying the size of the sliding window.
#' @param noiseType background noise type assumed for x. There are two options: white noise or red noise.
#' @param alpha the significance level. When the noise test p value of the subsequence is smaller than this significance level,
#' it is defined as a potential event.
#' @param main title of the animiation plot; default is `Animation plot of event detection'.
#' @param xlab x label of the animation plot; default is `t'.
#' @param ylab y label of the animation plot; default is `x'.
#' @param movie.name name of the output gif file; default is `animation.gif'.
#' @param interval a positive number to set the time interval of the animation (unit in seconds); default is 0.05.
#' @param ani.width width of the gif file (unit in px), default is 1000.
#' @param ani.height height of the gif file (unit in px); default is 400.
#' @param outdir character: specify the output directory when exporting the animations; default to be the 
#' current working directory.
#' @return ...
#' @seealso \code{\link{noiseTests}}, \code{\link{eventExtraction}}, \code{\link{plotevents}}
#' @references Yihui Xie (2013). Animation: An R Package for Creating Animations and Demonstrating Statistical Methods. 
#' \emph{Journal of Statistical Software}, \bold{53}(1), 1-27. \url{http://www.jstatsoft.org/v53/ i01/}.
#' @export
#' @examples
#' set.seed(123)
#' # generate an artificial time series
#' x=c(rnorm(128),cbfs(type='box'),rnorm(128),cbfs(type='rc'),rnorm(128))
#' # generate a gif file to show the event detection process
#' \dontrun{
#' aniplotevents(x,w=128,noiseType='white',outdir=getwd())
#' }

aniplotevents <- function(x, w, noiseType = c("white", "red"), alpha = 0.05, main = "Animation plot of events", xlab = "t", ylab = "x", 
    movie.name = "animation.gif", interval = 0.05, ani.width = 1000, ani.height = 400, outdir = getwd()) {
    noiseType <- match.arg(noiseType)
    tests = c()
    eventsFound = c()
    pbar <- tkProgressBar("test progress bar", "Some information in %", 0, 100, 50)
    saveGIF(for (i in 1:(length(x) - w)) {
        info <- sprintf("%d%% done", round(i/(length(x) - w) * 100))
        setTkProgressBar(pbar, round(i/(length(x) - w) * 100), sprintf("test (%s)", info), info)
        xsub = x[(i + 1):(w + i)]
        testx = noiseTests(xsub, w, noiseType)
        tests = c(tests, testx)
        if (testx <= alpha) {
            if (is.null(eventsFound$start)) {
                plot(c(1:(w + i)), x[1:(w + i)], ty = "l", xlab = xlab, ylab = ylab, col = "#9FC8DC", main = main)
            } else {
                plot(c(1:(w + i)), x[1:(w + i)], ty = "l", xlab = xlab, ylab = ylab, col = "#9FC8DC", main = main)
                for (j in 1:length(a)) {
                  lines(c(a[j]:b[j]), x[a[j]:b[j]], xlab = xlab, ylab = ylab, col = "#ff53a9")
                  xline(a[j], lty = 2, col = "#ff53a9")
                  xline(b[j], lty = 2, col = "#ff53a9")
                }
            }
            lines(c((i + 1):(w + i)), x[(i + 1):(w + i)], xlab = xlab, ylab = ylab, col = "#ff53a9")
        } else if (testx > alpha) {
            eventsFound = eventExtraction(tests, w, alpha)
            if (is.null(eventsFound$start)) {
                plot(c(1:(w + i)), x[1:(w + i)], ty = "l", xlab = xlab, ylab = ylab, col = "#9FC8DC", main = main)
            } else {
                a = ceiling((eventsFound$start + eventsFound$end)/2)
                b = a + w - 1
                plot(c(1:(w + i)), x[1:(w + i)], ty = "l", xlab = xlab, ylab = ylab, col = "#9FC8DC", main = main)
                for (j in 1:length(a)) {
                  lines(c(a[j]:b[j]), x[a[j]:b[j]], xlab = xlab, ylab = ylab, col = "#ff53a9")
                  xline(a[j], lty = 2, col = "#ff53a9")
                  xline(b[j], lty = 2, col = "#ff53a9")
                }
            }
        }
    }, movie.name = movie.name, interval = interval, ani.width = ani.width, ani.height = ani.height, out.dir = outdir)
} 
