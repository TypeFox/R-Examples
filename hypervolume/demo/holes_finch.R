if (exists('doHypervolumeHolesFinchDemo')==TRUE)
{
  data(morphSnodgrassHeller)
  
  # select data for only Isabela Island
  finch_isabela <- morphSnodgrassHeller[morphSnodgrassHeller$IslandID=="Isa_Alb",]
  
  # select trait axes
  trait_axes <- c("BodyL","WingL","TailL","BeakW")
  trait_data <- finch_isabela[,trait_axes]
  # keep complete cases only
  trait_data <- na.omit(trait_data)
  # convert length units to comparable scales
  trait_data_scaled <- scale(trait_data)
  
  getholes <- function(bw)
  {
    # get overall community hypervolume
    hv_finch <- hypervolume(trait_data_scaled,bandwidth=bw,name="Isabela Island finches")
    # compute convex expectation
    ec_finch <- expectation_convex(hv_finch, check_memory=F)
    # find holes
    holes_finch <- hypervolume_holes(hv_finch, ec_finch, set_check_memory=F)
    
    # return combined result 
    return(list(hv=hv_finch, ec=ec_finch, holes=holes_finch))
  }
  
  # the below line is commented out because of approximate two-hour runtime
  # bw_plugin <- estimate_bandwidth(trait_data_scaled, method="plug-in")
  # instead, use pre-computed output for that function
  bw_plugin <- c(0.5278828, 0.4883812, 0.5951435, 0.4480163)
  
  # do holes calculation
  result <- getholes(bw_plugin)
  
  # extract volume statistics
  vol_hv <- result$hv@Volume
  vol_ec <- result$ec@Volume
  vol_holes <- result$holes@Volume
  
  # calculate fraction of volume that is holey
  print(hole_volume_ratio <- vol_holes / vol_ec)
  # calculate approximate length of axis occupied
  print(length_ratio <- hole_volume_ratio ^ (1/4))
  
  # plot holes
  plot(hypervolume_join(result$hv,result$holes),
       showcentroid=F,darkfactor=0,col=c('purple','green'), npmax_random=2000,
       names=c("Body length","Wing length","Tail length", "Beak width"),
       legend=F,cex.names=1.5,contour.lwd=3)
  
  
  # calculate (in transformed coordinates) the centroid of the holes	
  print(holepos <- apply(result$holes@RandomUniformPointsThresholded,2,mean))
  print(holepos)
  
  # calculate (in untransformed coordinates) the centroid of the holes
  cpos <- attr(trait_data_scaled, "scaled:center")
  csca <- attr(trait_data_scaled, "scaled:scale")
  hole_origcoords <- holepos * csca + cpos
  print(hole_origcoords)
  
  # determine which other species in the dataset is most similar to the hole
  # calculate species mean trait values
  speciesmeans <- as.data.frame(do.call("rbind",by(morphSnodgrassHeller[, trait_axes], morphSnodgrassHeller $TaxonOrig, colMeans,na.rm=T)))
  # calculate rescaled distances
  scaled_diffs <- scale(speciesmeans - hole_origcoords)
  speciesmeans$dist <- apply(scaled_diffs, 1, function(x) { sqrt(sum(x^2))})
  # reorder species list by distance
  speciesmeans <- speciesmeans[order(speciesmeans $dist),]
  # identify 'top candidates' for filling the hole
  print(head(speciesmeans))
  
} else
{
  message('Demo does not run by default to meet CRAN runtime requirements.')
  message('This demo requires approximately 1 minutes to run.')  
  message('To run the demo, type')
  message('\tdoHypervolumeHolesFinchDemo=TRUE')
  message('\tdemo(holes_finch)')
  message('at the R command line prompt.')
}