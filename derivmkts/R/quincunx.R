#' @title Quincunx simulation
#'
#' @name quincunx
#'
#' @description \code{quincunx} simulates balls dropping down a
#' pegboard with a 50\% chance of bouncing right or left at each
#' level. The balls accumulate in bins. If enough balls are dropped,
#' the distribution approaches normality. This device is called a
#' quincunx. See \url{http://www.mathsisfun.com/data/quincunx.html}
#'
#'
#' @param n  Integer The number of peg levels, default is 3
#' @param numballs Integer The number of balls dropped, default is 20
#' @param delay Numeric Number of seconds between ball drops. Set
#' delay > 0 to see animation with \code{delay} seconds between
#' dropped balls. If \code{delay < 0}, the simulation will run to
#' completion without delays. If \code{delay == 0}, the user must
#' hit <return> for the next ball to drop. The default is 0.1 second
#' and can be set with the \code{delay} parameter.
#' @param probright Numeric The probability the ball bounces to the
#' right; default is 0.5
#' @param plottrue Boolean If \code{TRUE}, the display will indicate
#' bin levels if the distribution were normal. Default is TRUE
#'
#' @importFrom graphics lines barplot plot par points
#' 
#' @examples
#' 
#' ## These examples will not display correctly within RStudio unless
#' ## the plot window is large
#' quincunx(delay=0)
#' quincunx(n=10, numballs=200, delay=0)
#' quincunx(n=20, numballs=200, delay=0, probright=0.7)
#'
#' @importFrom stats rbinom
#' 

#' @export
quincunx <- function(n=3, numballs=20, delay=0.1,
                     probright=0.5, plottrue=TRUE) {
    ## inputs: n: number of binomial steps, number of balls, 
    ## time delay for plotting
    nlev <- n + 1 ## number of levels
    tbl <- array(0, dim=nlev)
    names(tbl) <- 0:(n)
    x <- array(dim=nlev*(nlev + 1)/2)
    y <- array(dim=nlev*(nlev + 1)/2)
    center <- round((nlev + 0.01)/2)
    ## enumerate the pegs
    x[1] <- center
    y[1] <- nlev
    for (i in 2:nlev) {
        w <- ((i*(i-1)/2)+1):(i*(i+1)/2)
        x[w] <- center + (1:i) - ((i+1)/2)
        y[w] <- rep(nlev+1-i, i)
    }
    ## graph the outcome
    par(mfcol=c(2,1))
    for (j in 1:numballs) {
        s <- cumsum(rbinom(nlev-1, size=1, prob=probright))
        tbl[s[nlev-1]+1] <- tbl[s[nlev-1]+1] + 1
        if (delay != 0 | j==numballs) {
            peglook <- function(a, b)
                c(rep(a, n*(n+1)/2), rep(b, n+1))
            ## Need last row of plot to look different, so vectorize
            ## characteristics
            plot(x, y, axes=FALSE, xlab='', ylab='',
                 pch=peglook(1, 0),
                 cex=peglook(1, 1.25),
                 col=peglook('black', 'red'))
            path <- c(1, cumsum(1:(nlev))[-1] - 1:(nlev-1) + s)
            lines(x[path], y[path], lty=2)
            h <- barplot(tbl, ylim=c(0, round(numballs/2.8)),
                    xlab='Number of Rightward Bounces')
            if (delay > 0) Sys.sleep(delay)
            if (delay < 0) readline()
        }
    }
    ## theoretical probability distribution
    if (plottrue) {
        k <- 0:n
        num <- choose(n, k)*probright^k*(1-probright)^(n-k)*numballs
        points(h, num, cex=3, xlim=c(0, n), pch=1)
    }
}

