#' Generate an artificial event with white noise
#' 
#' This function generates a box, cliff-ramp, ramp-cliff or a sine function with different levels of white noise as the background noise. Length
#' of the generated event is 128. Generation of events are similar to that of Cylinder-Bell-Funnel dataset in the reference below (Keogh and Lin 2005).
#' 
#' @param type type of the event to be generated. There are four options: `box', `rc',`cr',`sine' representing 
#' a box, cliff-ramp, ramp-cliff or a sine function.
#' @param A amplitude of the event; default is 10.
#' @param sigma a scalar specifying the level of white noise. Default is 1, which means the standard deviation of noise is 1.
#' @return an artificial event with white noise.
#' @references Eamonn Keogh and Jessica Lin (2005). Clustering of time-series subsequences is meaningless: implications for previous and future research. 
#' \emph{Knowl. Inf. Syst.}, \bold{8}(2), 154-177. \url{http://dblp.uni- trier.de/db/journals/kais/kais8.html#KeoghL05}.
#' @references Yanfei Kang, Kate Smith-Miles, Danijel Belusic (2013). How to extract meaningful shapes from noisy time-series subsequences? \emph{2013 IEEE Symposium on 
#' Computational Intelligence and Data Mining}, Singapore, 65-72. \url{http://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=6597219&isnumber=6597208}. 
#' @references Yanfei Kang, Danijel Belusic, Kate Smith-Miles (2014). Detecting and Classifying Events in Noisy Time Series. \emph{J. Atmos. Sci.}, \bold{71}, 1090-1104.
#' \url{http://dx.doi.org/10.1175/JAS-D-13-0182.1}.
#' @export
#' @examples
#' # generate a box function with white noise
#' set.seed(123)
#' x1 = cbfs(type = 'box', sigma = 1)
#' # generate a box function with higher level noise
#' set.seed(123)
#' x2 = cbfs(type = 'box', sigma = 3)
#' # plot them
#' par(mfrow=c(1,2))
#' plot(x1,type='l',xlab='t',ylab=expression(x[1]))
#' plot(x2,type='l',xlab='t',ylab=expression(x[2]))

cbfs <- function(type = c("box", "rc", "cr", "sine"), A = 10, sigma = 1) {
    t = 0:127
    # simulate backgroud white noise
    noise = rnorm(128) * sigma
    type <- match.arg(type)
    if (type == "box") {
        # simulate an event
        y = rep(0, length(t))
        a = ceiling(runif(1) * 16) + 16
        b = ceiling(runif(1) * 16) + 112
        y[a:b] = rep(A, b - a + 1)
        # add noise to the event
        finalevent = noise + y
    }
    if (type == "rc") {
        # simulate an event
        y = rep(0, length(t))
        a = ceiling(runif(1) * 16) + 16
        b = ceiling(runif(1) * 16) + 112
        y[a:b] = A * (t[a:b] - a)/(b - a)
        # add noise to the event
        finalevent = noise + y
    }
    if (type == "cr") {
        # simulate an event
        y = rep(0, length(t))
        a = ceiling(runif(1) * 16) + 16
        b = ceiling(runif(1) * 16) + 112
        y[a:b] = A * (b - t[a:b])/(b - a)
        # add noise to the event
        finalevent = noise + y
    }
    if (type == "sine") {
        # simulate an event
        y = A/2 * sin(t/127 * 2 * pi - pi/2) + A/2
        # add noise to the event
        finalevent = noise + y
    }
    return(finalevent)
} 
