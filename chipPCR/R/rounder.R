rounder <- function(object, cyc = 1) {
  if (class(object) != "der") 
    stop("Object must have class 'der'.")
  
  #number of cycles
  n.cyc <- round(range(object[, cyc]), 0)
  #outline values appropriate for different cycle
  cuts <- cut(object[, cyc], 
              breaks = c(n.cyc[1], seq(n.cyc[1] + 0.5, n.cyc[2] - 0.5, by = 1), n.cyc[2]), 
              include.lowest = TRUE)
  res <- cbind(cyc = seq(n.cyc[1], n.cyc[2]), t(sapply(levels(cuts), function(i) colMeans(object[cuts == i, -cyc]))))
  rownames(res) <- NULL
  new("der", '.Data' = res, 'method' = slot(object, "method"))
}
