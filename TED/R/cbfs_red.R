#' Generate an artificial event with red noise
#' 
#' This function generates a box, cliff-ramp, ramp-cliff or a sine function with red noise (AR(1)) as the background noise.  Length
#' of the generated event is 128.
#' 
#' @param type type of the event to be generated. There are four options: ``box', ``rc',``cr',``sine' representing 
#' a box, cliff-ramp, ramp-cliff or a sine function.
#' @param A amplitude of the event; default is 10.
#' @param s standard deviation of the AR(1) model innovations.  Default is 1.
#' @param coeff coefficient of the AR(1) process, which is used to control the level of red noise. Default is 0.5.
#' @return an artificial event with red noise.
#' @export
#' @examples
#' # generate a box function with red noise
#' set.seed(123)
#' x = cbfs_red(type = 'box', coeff=0.5, s=1, A=10)
#' # plot it
#' plot(x,type='l',xlab='t',ylab='x')
cbfs_red <- function(type = c("box", "rc", "cr", "sine"), A = 10, s = 1, coeff = 0.5) {
    n = 128
    t = 1:n
    # simulate backgroud red noise
    noise = arima.sim(list(order = c(1, 0, 0), ar = coeff), n = n, sd = s)
    
    if (type == "box") {
        
        # simulate an event
        y = rep(0, length(t))
        a = ceiling(runif(1) * 16) + 16
        b = ceiling(runif(1) * 16) + 112
        y[a:b] = rep(A, b - a + 1)
        # add noise to the event
        finalevent = y + noise
        
    }
    if (type == "rc") {
        # simulate an event
        y = rep(0, length(t))
        a = ceiling(runif(1) * 16) + 16
        b = ceiling(runif(1) * 16) + 112
        y[a:b] = A * (t[a:b] - a)/(b - a)
        # add noise to the event
        finalevent = y + noise
    }
    
    if (type == "cr") {
        # simulate an event
        y = rep(0, length(t))
        a = ceiling(runif(1) * 16) + 16
        b = ceiling(runif(1) * 16) + 112
        y[a:b] = A * (b - t[a:b])/(b - a)
        # add noise to the event
        finalevent = y + noise
        
    }
    
    if (type == "sine") {
        # simulate an event
        y = A/2 * sin(t/127 * 2 * pi - pi/2) + A/2
        # add noise to the event
        finalevent = y + noise
        
    }
    
    return(finalevent)
    
} 
