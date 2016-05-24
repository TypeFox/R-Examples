######################################################################
# Qini 
######################################################################

qini <- function(x, ...)  UseMethod("qini")

qini.default <- function(x, ...)
  stop("uplift: No method implemented for this class of object")

qini.performance <- function(x, direction = 1, plotit = TRUE, ...) {
  
  if (!inherits(x, "performance"))
    stop("uplift: x is not of class performance")
  
  ### check valid arguments
  if (!direction %in% c(1, 2))
    stop("uplift: direction must be either 1 or 2")
  
  perf <- x
  groups <- nrow(perf)
  
  if (direction == 1) {
  
    ### Model Incremental gains 
    inc.gains <- cumsum(perf[, 4] - perf[, 5] * sum(perf[, 2]) / sum(perf[, 3])) / sum(perf[, 2])
  
    ### Overall incremental gains
    overall.inc.gains <- sum(perf[, 4]) / sum(perf[, 2]) - sum(perf[, 5]) / sum(perf[, 3])
  
  
  } else {
    
    ### Model Incremental gains 
    inc.gains <- cumsum(-1 * (perf[, 4] - perf[, 5] * sum(perf[, 2]) / sum(perf[, 3]))) / sum(perf[, 2])
      
    ### Overall incremental gains
    overall.inc.gains <- sum(perf[, 5]) / sum(perf[, 3]) - sum(perf[, 4]) / sum(perf[, 2]) 
    
 }
  
  ### Random incremental gains
  random.inc.gains <- cumsum(rep(overall.inc.gains / groups, groups))
  
  ### Compute area under the model incremental gains (uplift) curve 
  x <- seq(1 / groups, 1, 1 / groups)
  y <- inc.gains
  
  auuc <- 0
  for (i in 2:length(x)) {
    auuc <- auuc + 0.5 * (x[i] - x[i-1]) * (y[i] + y[i-1])
  }
  
  ### Compute area under the random incremental gains curve
  y.rand <- random.inc.gains
  auuc.rand <- 0
  for (i in 2:length(x)) {
    auuc.rand <- auuc.rand + 0.5 * (x[i] - x[i-1]) * (y.rand[i] + y.rand[i-1])
  }
  
  ### Compute the difference between the areas (Qini coefficient)
  Qini <- auuc - auuc.rand
  miny <- 100 * min(c(random.inc.gains, inc.gains))
  maxy <- 100 * max(c(random.inc.gains, inc.gains))
  
  if (plotit) {
    plot(inc.gains * 100 ~ seq(100 / groups, 100, 100 / groups), type ="b",
         col = "blue", lty = 2, xlab = "Proportion of population targeted (%)", 
         ylab = "Cumulative incremental gains (pc pt)", ylim = c(miny, maxy), ...)
    lines(random.inc.gains * 100 ~ seq(100 / groups, 100, 100 / groups), type = "l", col = "red", lty = 1)
    legend("topright", c("Model", "Random"), 
           col=c("blue", "red"), lty=c(2,1))
  }  
    
  res <- list(Qini = Qini,
              inc.gains = inc.gains,
              random.inc.gains = random.inc.gains)
  
  return(res)
  
}


### END FUN
  
    
  