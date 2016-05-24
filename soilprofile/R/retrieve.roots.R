retrieve.roots <-
function(profile_data, a) {
  ## a list containing the 4 basic root shapes
  ## data(root_unit, package='soilprofile')
  roots_retrieved <- profile_data[[a]][[2]]
  if (length(roots_retrieved)!=0) {
    for (i in 1:length(roots_retrieved)) {
      ##final lwd for the given root
      lwd_final <- roots_retrieved[[i]][2]$lwd
      ##the choosen root
      the_root <- roots_retrieved[[i]][1]$root_coordinates
      ##a large root for plotting cases when root is larger than 5 cm
      ## plotting the root
      polygon(the_root, lwd=lwd_final, border='#24110A', col='#24110A')
    }
  }
}
