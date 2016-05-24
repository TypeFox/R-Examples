skeletal <-
function(horizon, clast_dimension, type, abundance, col) {
  if (is.na(clast_dimension)==FALSE | is.na(abundance)==FALSE) {
    ## data(subangular, package='soilprofile')
    ## data(subcircle, package='soilprofile')
    ## data(channer, package='soilprofile')
    horizon_reduced <- horizon*0.95
    ##the coordinates to shift the shape in the center of the horizon
    x_shift <- mean(horizon[,1])
    y_shift <- mean(horizon[,2])
    x1_shift <- mean(horizon_reduced[,1])
    y1_shift <- mean(horizon_reduced[,2])
    diff_x <- x1_shift-x_shift
    diff_y <- y1_shift-y_shift
    hor_red <- cbind(horizon_reduced[,1]-diff_x, horizon_reduced[,2]-diff_y)
    ##structures : subangular, subcircle, channer
    ##computes area of polygon
    area <- areapl(horizon)
    ## factor for skeletal dimension (three different dimensions)
    red_factor <- clast_dimension/3
    ##select shape type
    if (type=='subangular' | type=='subcircle' | type=='channer') {
      shape <- get(type)
      ##build units for smaller clasts
      unit <- shape/2
      unit2 <- shape/3
      if (type!='channer') {
        unit_tmp <- cbind(unit[,1]-0.17, unit[,2]-0.17)
        unit_tmp2 <- cbind(unit2[,1]+0.3, unit2[,2]+0.17)
        shape_2 <- rbind(unit_tmp,c(NA,NA), unit_tmp2,c(NA,NA))
        unit_tmp <- cbind(unit[,1]-0.17, unit[,2]+0.17)
        unit_tmp2 <- cbind(unit2[,1]+0.2, unit2[,2]-0.2)
        shape_3 <- rbind(unit_tmp,c(NA,NA), unit_tmp2,c(NA,NA))
      } else {
        unit_tmp <- cbind(unit[,1]-0.17, unit[,2]-0.07)
        unit_tmp2 <- cbind(unit2[,1]+0.3, unit2[,2]+0.13)
        shape_2 <- rbind(unit_tmp,c(NA,NA), unit_tmp2,c(NA,NA))
        unit_tmp <- cbind(unit[,1]-0.17, unit[,2]+0.1)
        unit_tmp2 <- cbind(unit2[,1]+0.2, unit2[,2]-0.1)
        shape_3 <- rbind(unit_tmp,c(NA,NA), unit_tmp2,c(NA,NA))
      }
      ##scale according to clast dimension (cm)
      shape_scaled <- shape*red_factor
      shape_scaled2 <- shape_2*red_factor
      shape_scaled3 <- shape_3*red_factor
      ##compute tolerance to coordinates for not having overlapped polygons
      x_tolerance <- shape_scaled[which.max(0-shape_scaled[,1]),1]*1.1
      y_tolerance <- shape_scaled[which.max(0-shape_scaled[,2]),2]*1.1
      ##computes area
      area_shape <- areapl(shape_scaled)
      ##computes expected area to cover
      skeletal_pc <- abundance
      area_exp <- area*skeletal_pc
      ##clast dimension higher than 40
      if (clast_dimension>=20 | skeletal_pc>0.3) {
        if (skeletal_pc>0.7) {surface <- area*2*skeletal_pc}
        else {surface <- area}
        extra_coordinates <-cbind(runif(trunc(100*surface/800), min=min(hor_red[,1])+1, max=max(hor_red[,1])-1), runif(trunc(100*surface/800), min=min(hor_red[,2]), max=max(hor_red[,2])))
        small_shp <- list()
        for (bb in 1:length(extra_coordinates[,1])) {
          extra_coordinate <- extra_coordinates[bb,]
          xdistance <- c(x_shift-extra_coordinate[1],y_shift-extra_coordinate[2])
          small_shape <- shape_scaled/4
          clast_centered <- cbind(small_shape[,1]+x_shift, small_shape[,2]+y_shift)
          clast_temp <- cbind(clast_centered[,1]-xdistance[1], clast_centered[,2]-xdistance[2])
          points_in <- point.in.polygon(clast_temp[,1], clast_temp[,2], hor_red[,1], hor_red[,2])
          out <- which(points_in!=1)
          if (length(out)==0) {
            small_shapes <- clast_temp
            polygon(clast_temp, col=col)
            ## names(small_shapes) <- 'extra_shapes'
          }
          else {small_shapes <- NA}
          small_shp[[bb]] <- small_shapes
        }
      } else (small_shp <- NA)
      ##computes n of polygons to plot
      n_pol <- trunc(area_exp/area_shape)
      ##the coordinates to shift the shape in the center of the horizon
      x_shift <- mean(horizon[,1])
      y_shift <- mean(horizon[,2])
      shape_shifted <- cbind(shape_scaled[,1]+x_shift, shape_scaled[,2]+y_shift)
      shape_shifted2 <- cbind(shape_scaled2[,1]+x_shift, shape_scaled2[,2]+y_shift)
      shape_shifted3 <- cbind(shape_scaled3[,1]+x_shift, shape_scaled3[,2]+y_shift)
      ##choose fist coordinate within the horizon
      rdm_coordinate <-cbind(runif(1, min=min(hor_red[,1]), max=max(hor_red[,1])), runif(1, min=min(hor_red[,2]), max=max(hor_red[,2])))
      coord_ok <- rdm_coordinate
      ##the final number that stops iteractions
      final_n <- 0
      ##get into the loop
      for (xx in 1:5000) {
        ##second coordinate
        rdm_coordinate <-c(runif(1, min=min(hor_red[,1]), max=max(hor_red[,1])), runif(1, min=min(hor_red[,2]), max=max(hor_red[,2])))
        ##check if the coordinate gives a polygon that is within the given horizon
        x_shift <- mean(horizon[,1])
        y_shift <- mean(horizon[,2])
        distance <- c(x_shift-rdm_coordinate[1],y_shift-rdm_coordinate[2])
        shape_temp <- cbind(shape_shifted[,1]-distance[1], shape_shifted[,2]-distance[2])
        points_in <- point.in.polygon(shape_temp[,1], shape_temp[,2], hor_red[,1], hor_red[,2])
        out <- which(points_in!=1)
        ##if the coordinate leads to a polygon tht overlaps the border it is discarded
        if (length(out)!=0) {next()}
        for (p in 1:length(coord_ok[,1])) {
          ##check the distance between the already existing and the new coordinates
          x_difference <- abs(rdm_coordinate[1]- coord_ok[p,1])
          y_difference <- abs(rdm_coordinate[2]- coord_ok[p,2])
          if (x_difference< 2*abs(x_tolerance) & y_difference< 2*abs(y_tolerance)) {break()}
          if (p==length(coord_ok[,1])) {
            coord_ok <- rbind(rdm_coordinate, coord_ok)
            final_n <- final_n+1
          }
          ##print(p)
        }
        if (final_n==n_pol) {break()}
        ##print(xx)
      }
      ##condition for plotting different skeletal dimensions randomly
      if (clast_dimension<5) {sample <- rep(1, length(coord_ok[,1]))}
      else {
        sample <- sample(c(1:3), size=length(coord_ok[,1]), replace=T)}
      the_polygons <- list()
      for (i in 1:length(coord_ok[,1])) {
        the_coord <- coord_ok[i,]
        the_distance <-c(x_shift-the_coord[1],y_shift-the_coord[2])
        the_shape <- list(shape_shifted, shape_shifted2, shape_shifted3)
        jj <- cbind(the_shape[[sample[i]]][,1]-the_distance[1], the_shape[[sample[i]]][,2]-the_distance[2])
        ##plot the polygons
        polygon(jj, col=col)
        the_polygons[[i]] <- jj
      }
      names(the_polygons) <- 'skeletal_polygons'
      skeletal_list <- list(small_shp, the_polygons)
    } else {skeletal_list <- NA}
  } else {skeletal_list <- NA}
  invisible(skeletal_list)
}
