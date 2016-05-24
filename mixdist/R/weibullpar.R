## last modified June 2002

weibullpar <- function(mu, sigma, loc = 0) 
{
    weibullpar1 <- function(weibpar) {
        mu <- weibpar[1]
        sigma <- weibpar[2]
        loc <- weibpar[3]
        if ((mu - loc) <= 0 | sigma < 0) {
            shape <- NA
            scale <- NA
        }
        else {
            cv <- sigma/(mu - loc)
            if (cv < 1e-06) {
                nu <- cv/(sqrt(trigamma(1)) - cv * digamma(1))
                shape <- 1/nu
                scale <- (mu - loc)/(1 + nu * digamma(1))
            }
            else {
                aa <- log(cv^2 + 1)
                nu <- 2 * cv/(1 + cv)
                repeat {
                  gb <- (lgamma(1 + 2 * nu) - 2 * lgamma(1 + 
                    nu) - aa)/(2 * (digamma(1 + 2 * nu) - digamma(1 + 
                    nu)))
                  nu <- nu - gb
                  if (abs(gb) < 1e-12) 
                    break
                }
                shape <- 1/nu
                scale <- exp(log(mu - loc) - lgamma(1 + nu))
            }
        }
        c(shape, scale, loc)
    }
    wpar <- data.frame(t(apply(cbind(mu, sigma, loc), 1, weibullpar1)))
    names(wpar) <- c("shape", "scale", "loc")
    wpar
}
