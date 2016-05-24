epi.edr <- function(dat, n = 4, conf.level = 0.95, nsim = 99, na.zero = TRUE){
   
   N. <- 1 - ((1 - conf.level) / 2)
   alpha <- 1 - conf.level
   z <- qnorm(N., mean = 0, sd = 1)
   
   num.sum <- 0; num.sd <- 0; num.n <- 0
   den.sum <- 0; den.sd <- 0; den.n <- 0   
   
   start <- 2 * n

   for (i in start:length(dat)){
      top.start <- (i - (n - 1))
      top.finish <- i
      bot.start <- (i - (2 * n)) + 1
      bot.finish <- i - n

      # Vector of outbreak counts for numerator and denominator:
      num.tmp <- dat[top.start:top.finish]
      den.tmp <- dat[bot.start:bot.finish]
      
      num.sum <- c(num.sum, sum(num.tmp))
      num.sd <- c(num.sd, sd(num.tmp))
      num.n <- c(num.n, length(num.tmp))

      den.sum <- c(den.sum, sum(den.tmp))
      den.sd <- c(den.sd, sd(den.tmp))
      den.n <- c(den.n, length(den.tmp))      
      }
      
   # Remove the initiating zero and add a vector of zeroes to the start:
   num.sum <- c(rep(0, times = (start - 1)), num.sum[-1])
   num.sd <- c(rep(0, times = (start - 1)), num.sd[-1])
   num.n <- c(rep(0, times = (start - 1)), num.n[-1])
   
   den.sum <- c(rep(0, times = (start - 1)), den.sum[-1])
   den.sd <- c(rep(0, times = (start - 1)), den.sd[-1])
   den.n <- c(rep(0, times = (start - 1)), den.n[-1])

   # Work out the standard error of numerator and denominator:
   # SE_total = (n * SE_mean):
   # num.se <- num.n * (num.sd / sqrt(num.n))
   # den.se <- den.n * (den.sd / sqrt(den.n))
   
   num.mat <- matrix(rep(0, times = length(num.sum) * nsim), nrow = length(num.sum))
   den.mat <- matrix(rep(0, times = length(num.sum) * nsim), nrow = length(num.sum))
   
   for(i in 1:nsim){
      num.mat[,i] <- rpois(n = length(num.sum), lambda = num.sum) 
      den.mat[,i] <- rpois(n = length(den.sum), lambda = den.sum) 
   }

   edr.p <- num.sum / den.sum
   edr.mat <- num.mat / den.mat   
   edr.mat[is.na(edr.mat)] <- 0
   
   quant <- function(x, probs) quantile(x, probs, na.rm = TRUE)
   edr.l <- apply(edr.mat, MARGIN = 1, FUN = quant, probs = alpha/2)
   edr.u <- apply(edr.mat, MARGIN = 1, FUN = quant, probs = 1 - alpha/2)

   # Work out EDR and confidence intervals of EDR: 
   # Source: http://www.agron.missouri.edu/mnl/55/34kowles.html
   # edr.sed <- sqrt(num.se^2 + den.se^2)
   # edr.var <- (num.se^2 / den.sum^2) + (num.sum^2 / den.sum^4) * den.se^2

   # Method 1 - use of extremes:
   # edr.l <- (num.sum - (z * num.se)) / (den.sum + (z * den.se))
   # edr.u <- (num.sum + (z * num.se)) / (den.sum - (z * den.se))
   
   # Method 2 - standard error of the difference between means:
   # edr.l <- 1 + ((num.sum - den.sum) - edr.sed) / (den.sum - (z * den.se))
   # edr.u <- 1 + ((num.sum - den.sum) - edr.sed) / (den.sum + (z * den.se))
   
   # Method 3 - approximate variance of the error of the ratios:
   # edr.l <- edr.p - (z * sqrt(edr.var))
   # edr.l[edr.l < 0] <- 0
   # edr.u <- edr.p + (z * sqrt(edr.var))
   
   if(na.zero == FALSE) {
      rval <- as.data.frame(cbind(edr.p, edr.l, edr.u))
      names(rval) <- c("est", "lower", "upper")
      }
   
   else if(na.zero == TRUE) {
      id <- is.na(edr.p)
      edr.p[id] <- 0
      edr.l[id] <- 0
      edr.u[id] <- 0
      
      id <- is.infinite(edr.p)
      edr.p[id] <- 0
      edr.l[id] <- 0
      edr.u[id] <- 0
 
      rval <- as.data.frame(cbind(edr.p, edr.l, edr.u))
      names(rval) <- c("est", "lower", "upper")
      }
rval
}

