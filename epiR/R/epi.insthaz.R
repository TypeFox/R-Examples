"epi.insthaz" <- function(survfit.obj, conf.level = 0.95){

   N <- 1 - ((1 - conf.level) / 2)
   z <- qnorm(N, mean = 0, sd = 1)
        
   time <- survfit.obj$time
   time0 <- c(0, time[-length(time)])
   interval <- (time - time0)
   
   a <- survfit.obj$n.event
   n <- survfit.obj$n.risk
   p <- a/n
   a. <- n/(n + z^2)
   b. <- a/n
   c. <- z^2/(2 * n)
   d. <- (a * (n - a)) / n^3
   e. <- z^2 / (4 * n^2)
   
   est <- p / interval     
   low <- (a. * (b. + c. - (z * sqrt(d. + e.)))) / interval
   up <- (a. * (b. + c. + (z * sqrt(d. + e.)))) / interval
   
   rval <- as.data.frame(cbind(time, est, low, up))
   names(rval) <- c("time", "est", "lower", "upper")
   return(rval)
}