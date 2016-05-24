epi.betabuster <- function(mode, conf, greaterthan, x, conf.level = 0.95, max.shape1 = 100, step = 0.001){

   shape1 <- seq(from = 1, to = max.shape1, by = step)
   shape2 <- 2 - shape1 + (shape1 - 1) / mode
   p.vec <- pbeta(q = x, shape1 = shape1, shape2 = shape2)
  
   # What value of a has the lowest (abs(p.vec-(1 - q)))?
   if(greaterthan){
      index <- which((abs(p.vec - (1 - conf))) == min(abs(p.vec - (1 - conf))))
    }
   else{
      index <- which((abs(p.vec - conf)) == min(abs(p.vec - conf)))
    }
  
   shape1 <- shape1[index]
   shape2 <- shape2[index]
  
   #  In general, if an experiment resulted in 's' successes (e.g. no. test-positive animals) 
   #  recorded in 'n' trials (e.g. number of truly infected animals), 
   #  use of a beta (a, b) distribution with a = s+1 and b = n-s+1 is an appropriate choice to model the uncertainty in that parameter.
   s <- shape1 - 1
   n <- shape1 + shape2 - 2
   .mode <- (shape1 - 1) / (shape1 + shape2 - 2)
   .mean <- shape1 / (shape1 + shape2)
   .var <- shape1 * shape2 / (((shape1 + shape2)^2) * (shape1 + shape2 + 1))
   .median <- qbeta(p = 0.5, shape1 = shape1, shape2 = shape2)
  
   lower <- qbeta(p = (1 - conf.level) / 2, shape1 = shape1, shape2 = shape2)
   upper <- qbeta(p = 1 - ((1 - conf.level) / 2), shape1 = shape1, shape2 = shape2)      
  
  # dens <- dbeta(x = seq(from = 0, to = 1,by = 0.001), shape1 = a, shape2 = b)
  # beta.plot <- plot(x = seq(from = 0, to = 1, by = 0.001), y = dens, type = 'l', xlab = "Proportion", ylab = "Density")
  rval <- list(shape1 = shape1, shape2 = shape2, mode = .mode, mean = .mean, median = .median, lower = lower, upper = upper, variance = .var)
  rval
  
  # Example:
  # fred <- epi.betabuster(mode = 0.25, conf.level = 0.95, greaterthan = FALSE, x = 0.30, max.a = 100, step = 0.001); fred$a; fred$b;
  # plot(seq(from = 0, to = 1,by = 0.001), dbeta(x = seq(from = 0, to = 1,by = 0.001), shape1 = fred$a, shape2 = fred$b), type = 'l', xlab = "Proportion", ylab = "Density")
}
