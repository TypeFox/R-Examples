#' Generate auto-regressive order 1 sequence
#' @param length size of the sequence
#' @param mean mean value of the sequence \deqn{\mu}
#' @param sd standard deviation of the sequence \deqn{\sigma}
#' @param rho auto-correlation factor
#' @param rand function to generate one random step
#' @return Auto Regressive ("AR") sequence \deqn{X=\mu * \left\{ x_n \right\}_{1 \leq n \leq N}, x_1 = \sigma * u_1 / \sqrt{1-\rho^2}, x_n=\rho * x_{n-1} + \sigma * u_n}
#' @examples
#' x = AR1()
#' plot(x,type='l',col=rgb(.5,0,0,.5))
AR1 <- function(length=1000, mean=0, sd=1, rho=.1, rand = function(){rnorm(1)}) {
    X = array(NaN, length)
    X[1] = sd * rand() / sqrt(1-rho^2)
    for (i in 2:length) {
        X[i] = X[i-1] * rho + sd * rand()
    }
    return(X + mean)
}
