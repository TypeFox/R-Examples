`tss` <-
function(x) {
  n <- length(as.numeric(x))
  sum(as.numeric(x)^2)-n*mean(as.numeric(x))^2
}

