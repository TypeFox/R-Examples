roots <-
function(horizon, root.abundance, root.dimension, orientation='h') {
  if (root.abundance=='') {root.abundance <- 'absent'}
  the_roots <- list()
  if (root.abundance!='absent') {
    ## a list containing the 4 basic root shapes
    ## data(root_unit, package='soilprofile')
    area <- areapl(horizon)
    x_shift <- mean(horizon[,1])
    y_shift <- mean(horizon[,2])
    for (jj in 1:4) {
      data <- root_unit[[jj]]
      ## increase the size of single shapes
      data <- data*4
      ##shifts the roots to the center of the horizon 
      root_shifted <- cbind(data[,1]+x_shift, data[,2]+y_shift)
      root_unit[[jj]] <- root_shifted
    }
    ##assessing root abundance
    if (root.abundance=='absent') {n_roots <- 0}
    if (root.abundance=='few') {n_roots <- 20}
    if (root.abundance=='common') {n_roots <- 50}
    if (root.abundance=='many') {n_roots <- 100}
    if (root.abundance=='abundant') {n_roots <- 200}
    n_roots <- trunc(n_roots*area/300)
    ##assessing root dimension (cm)
    if (root.dimension<0.1) {lwd=rep(1,6)}
    if (root.dimension>=0.1 |root.dimension<1) {lwd=c(1,1,1.5,1.5,1.5,1.5)/2}
    if (root.dimension>=1 |root.dimension<5) {lwd=c(1,1.5,2,2,2,2)/2}
    if (root.dimension>=5) {lwd=c(1,1.5,2,4,4,4)/2}
    ## computes random coordinates within the horizon surface
    rdm_coordinates <-cbind(runif(10000, min=min(horizon[,1])+1, max=max(horizon[,1])-1), runif(10000, min=min(horizon[,2]), max=max(horizon[,2])))
    ## set the root orientation
    if (orientation=='v') {to_sample <- c(3:4)} else {to_sample <- c(1:2)} 
    sample <- sample(to_sample, size=10000, replace=T)
    sample2 <- sample(1:6, size=10000, replace=T)
    ## the loop for plotting roots
    n_final <- 0
    for (i in 1:10000) {
      ##final lwd for the given root
      lwd_final <- lwd[sample2[i]]
      ##the choosen root
      the_root <- root_unit[[sample[i]]]
      ##a large root for plotting cases when root is larger than 5 cm
      the_big_root <- cbind((the_root[,1]-x_shift)*2+x_shift, (the_root[,2]-y_shift)*2+y_shift)
      if (lwd_final==2) {the_root <- the_big_root}
      ##computation for final shift of the root
      rdm_coordinate <- rdm_coordinates[i,]
      distance <- c(x_shift-rdm_coordinate[1],y_shift-rdm_coordinate[2])
      root_temp <- cbind(the_root[,1]-distance[1], the_root[,2]-distance[2])
      points_in <- point.in.polygon(na.omit(root_temp[,1]), na.omit(root_temp[,2]), horizon[,1], horizon[,2])
      out <- which(points_in!=1)
      ##if the coordinate leads to a root tht overlaps the border it is discarded
      if (length(out)!=0) {next()}
      n_final <- n_final+1
      ## plotting the root
      polygon(root_temp, lwd=lwd[sample2[i]], border='#24110A', col='#24110A')
      the_root_tmp <- list(root_temp, lwd[sample2[i]])
      names(the_root_tmp) <- c('root_coordinates', 'lwd')
      the_roots[[n_final]] <- the_root_tmp
      names(the_roots) <- 'roots'
      ## when final number is reached, break the loop
      if (n_final==n_roots) {break()}
    }
  } else {print('no roots')}
  roots_dataframe <- the_roots
invisible(roots_dataframe)
}
