###
### R routines for the R package mvmeta (c) Antonio Gasparrini 2012-2014
#
sumlist <-
function(list) {
#
################################################################################
#
  res <- 0
  for(i in seq(list)) res <- res + list[[i]]
#
  res
}
