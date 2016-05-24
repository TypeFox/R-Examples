z.test <- function(x, mu=0, stdev,
                   alternative = c("two.sided", "less", "greater"),
                   sd=stdev, n=length(x),
                   conf.level = 0.95, ... ) {

  if(missing(stdev) && missing(sd)) stop("You must specify a Standard Deviation of the population")

  alternative <- match.arg(alternative)

  z <- (mean(x)-mu)/(sd/sqrt(n))

  out <- list(statistic=c(z=z))
  class(out) <- 'htest'

  out$parameter <- c(n=n,"Std. Dev." = sd,
                     "Std. Dev. of the sample mean" = sd/sqrt(n))

  out$p.value <- switch(alternative,
                    two.sided = 2*pnorm(abs(z),lower.tail=FALSE),
                    less = pnorm(z),
                    greater = pnorm(z, lower.tail=FALSE) )

  out$conf.int <- switch(alternative,
                         two.sided = mean(x) +
                           c(-1,1)*qnorm(1-(1-conf.level)/2)*sd/sqrt(n),
                         less = c(-Inf, mean(x)+qnorm(conf.level)*sd/sqrt(n)),
                         greater = c(mean(x)-qnorm(conf.level)*sd/sqrt(n), Inf)
                         )
  attr(out$conf.int, "conf.level") <- conf.level

  out$estimate <- c("mean of x" = mean(x))
  out$null.value <- c("mean" = mu)
  out$alternative <- alternative
  out$method <- "One Sample z-test"
  out$data.name <- deparse(substitute(x))
  names(out$estimate) <- paste("mean of", out$data.name)

  return(out)
}
