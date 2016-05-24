retrieve.skeletal <-
function(profile_data, a, col) {
  skeletal_retrieved <- profile_data[[a]][[3]]
  small_shapes <- skeletal_retrieved[[1]]
  for (i in 1:length(small_shapes)) {
    clast_temp <- small_shapes[[i]]
    polygon(clast_temp, col=col)
    ## names(small_shapes) <- 'extra_shapes'
  }
  if (length(skeletal_retrieved)!=1) {
    big_shapes <- skeletal_retrieved[[2]]
    for (i in 1:length(big_shapes)) {
      jj <- big_shapes[[i]]
      polygon(jj, col=col)
    }
  }
}
