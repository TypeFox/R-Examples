ParEx <-
function(img.file, out.file, section.title){
  #Read image
  my.img <- readImage(img.file)
  #Object fitting operation
  my.obj <- bwlabel(my.img)
  #Extract parameters into data frame
  my.par <- data.frame(computeFeatures.moment(my.obj[,,1], my.img[,,1]))
  #Determine minor axis length from major axis and eccentricity
  my.par$m.minoraxis <- sqrt(my.par$m.majoraxis^2 * (1 - my.par$m.eccentricity^2))
  #Determine ellipse area
  my.par$m.area <- pi * (my.par$m.majoraxis / 2) * (my.par$m.minoraxis / 2)
  #Convert from image coordinates to Cartesian
  my.par$m.cy <- my.img@dim[2] - my.par$m.cy
  #Calculate axial ratio
  my.par$m.ratio <- my.par$m.majoraxis / my.par$m.minoraxis
  #Determine rake in radians based on ellipsoid 2003 convention
  my.par$m.rake.rad <- sapply(my.par$m.theta, function(x){ifelse(x < 0, x + pi, x)})
  #Determine rake in degrees based on ellisoid 2003 convention
  my.par$m.rake.deg <- sapply(my.par$m.rake.rad, function(x){x * (180 / pi)})
  #Open PDF device for ellipse fit output
  pdf(file = out.file, width = 6, height = 6, useDingbats = FALSE)
  #Set PDF and plotting parameters
  par(family = 'serif', xaxs = 'i', yaxs = 'i', omi = c(0,0,0,0), mai = c(.5,.5,.5,.5))
  plot(x = 0, y = 0, pch = '', xlim = c(0 , my.img@dim[1]), ylim = c(0, my.img@dim[2]), asp = 1, ann = FALSE, axes = FALSE)
  title(paste("Fitted ellipses for section: ", section.title))
  polygon(c(0, my.img@dim[1], my.img@dim[1], 0, 0 ),
    c(0, 0, my.img@dim[2], my.img@dim[2], 0), border = 'gray', lwd = 1)
  par(ps = 6)
  #Run loop to add ellipses to plot
  i <- 0
  while(i < length(my.par$m.cx)){
    i <- 1 + i
    #Determine ellipse node coordinates
    tt <- seq(0, 2 * pi, length = 72)
    sa <- sin(-1 * my.par$m.theta[i])
    ca <- cos(-1 * my.par$m.theta[i])
    ct <- cos(tt)
    st <- sin(tt)
    x <- my.par$m.cx[i] + my.par$m.majoraxis[i]/2 * ct * ca - my.par$m.minoraxis[i]/2 * st * sa
    y <- my.par$m.cy[i] + my.par$m.majoraxis[i]/2 * ct * sa + my.par$m.minoraxis[i]/2 * st * ca
    #Plot ellipse
    polygon(x, y, col = "gray", border = "white", lwd = .01)
    #Add numerical lables to objects
    text(x = my.par$m.cx[i], y = my.par$m.cy[i], i)
  }
  #Destroy device connection
  dev.off()
  #Remove image
  remove(my.img)
  #Return parameters corresponding to labels
  return(my.par)
}
