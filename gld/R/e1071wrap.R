# use e1071 to calculate moments for estimation
# allow passing na.rm
# don't export

moments4data <- function(x, ...){
  mean <- mean(x)
  var <- var(x)
  skewr <- e1071::moment(x,order = 3,center = TRUE, ...)/(var^(1.5))
  kurtr <- e1071::moment(x,order = 4,center = TRUE, ...)/(var^(2))
  c(mean,var,skewr,kurtr)
}