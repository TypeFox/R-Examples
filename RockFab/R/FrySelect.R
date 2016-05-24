FrySelect <-
function(fry.data, max.dim, out.file, normalized = FALSE, select = TRUE){
  #Remove Fry coordinates outside of max dimension
  point.data <- fry.data@fryParams[which(fry.data@fryParams$dist <= max.dim),]
  
  #Populate new objects wiith abridged Fry coords
  x.coords <- point.data$x.coords
  y.coords <- point.data$y.coords
  
  #Determine coordinates to plot outer ring
  p.ring <- seq(from = 0, to = 2 * pi, length = 72)
  x.ring <- max.dim * cos(p.ring)
  y.ring <- max.dim * sin(p.ring)
  if(select){
    #Open interactive Fry plot for apogee determination
    dev.new()
    par(mai = c(0, 0, .25, 0), omi = c(0, 0, .5, 0), family = "serif")
    plot(x.coords, y.coords, ann = FALSE, axes = FALSE, asp = 1, pch = 19, , xlim = c(-1 * max.dim - (max.dim * .1), max.dim + (max.dim * .1)))
    lines(x.ring, y.ring, col = "gray", lwd = 1.5)
    mtext(text = "Select central void apogee", side = 3, outer = TRUE, cex = 1.25)
    points(0, 0, pch = 3, col = "red")
    #Select central void apogee
    apogee <- locator(n = 1, type="n")
    
    #Open interactive Fry plot for perigee determination
    plot(x.coords, y.coords, ann = FALSE, axes = FALSE, asp = 1, pch = 19, xlim = c(-1 * max.dim, max.dim))
    lines(x.ring, y.ring, col = "gray", lwd = 1.5)
    mtext(text = "Select central void perigee", side = 3, outer = TRUE, cex = 1.25)
    points(0, 0, pch = 3, col = "red")
    #Add apogee selection
    lines(x = c(-1 * as.numeric(apogee[1]), as.numeric(apogee[1])),
      y = c(-1 * as.numeric(apogee[2]), as.numeric(apogee[2])),
      col = "red", lwd = 1.5)
    #Select central void perigee
    perigee <- locator(n = 1, type="n")
    #Add perigee section
    lines(x = c(-1 * as.numeric(perigee[1]), as.numeric(perigee[1])),
      y = c(-1 * as.numeric(perigee[2]), as.numeric(perigee[2])),
      col = "blue", lwd = 1.5)
  
    #Combine coordinates into matrix
    elli.coords <- matrix(rbind(apogee, perigee, deparse.level = 0), ncol = 2)
    
    #Calculate semi-major and minor axial lengths
    a.half.axis <- sqrt(as.numeric(elli.coords[1,1])^2 + as.numeric(elli.coords[1,2])^2)
    b.half.axis <- sqrt(as.numeric(elli.coords[2,1])^2 + as.numeric(elli.coords[2,2])^2)
    
    #Run algorith to determine rake azimuth
    if(as.numeric(elli.coords[1,1]) >= 0 & as.numeric(elli.coords[1,2]) > 0){
      rake.azimuth <- abs(atan(as.numeric(elli.coords[1,1]) / as.numeric(elli.coords[1,2])))
    }
    if(as.numeric(elli.coords[1,1]) >= 0 & as.numeric(elli.coords[1,2]) <= 0){
      rake.azimuth <- (pi / 2) + abs(atan(as.numeric(elli.coords[1,2]) / as.numeric(elli.coords[1,1])))
    }
    if(as.numeric(elli.coords[1,1]) < 0 & as.numeric(elli.coords[1,2]) <= 0){
      rake.azimuth <- (pi) + abs(atan(as.numeric(elli.coords[1,1]) / as.numeric(elli.coords[1,2])))
    }
    if(as.numeric(elli.coords[1,1]) < 0 & as.numeric(elli.coords[1,2]) > 0){
      rake.azimuth<-((3 / 2) * pi) + abs(atan(as.numeric(elli.coords[1,2]) / as.numeric(elli.coords[1,1])))
    }
    #Convert to ellipsoid azimuth for plotting
    elli.azimuth <- pi / 2 + rake.azimuth
    
    #Determine ellipse node coords and add to plot
    tt <- seq(0, 2 * pi, length = 72)
    sa <- sin(-1 * elli.azimuth)
    ca <- cos(-1 * elli.azimuth)
    ct <- cos(tt)
    st <- sin(tt)
    x <- a.half.axis * ct * ca - b.half.axis * st * sa
    y <- a.half.axis * ct * sa + b.half.axis * st * ca
    polygon(x, y, border = "black", lwd = 1.5)
    
    #populate Fry data object with parameters
    fry.data@rsAxes <- c(2 * a.half.axis, 2 *  b.half.axis)
    fry.data@strainRatio <- a.half.axis / b.half.axis
    fry.data@meanObjectArea <- pi * a.half.axis * b.half.axis
    fry.data@vectorMean <- (elli.azimuth%%pi) * (180 / pi)
  }else{
    #Determine ellipse node coords and add to plot
    tt <- seq(0, 2 * pi, length = 72)
    sa <- sin(-1 * fry.data@vectorMean * (pi / 180))
    ca <- cos(-1 * fry.data@vectorMean * (pi / 180))
    ct <- cos(tt)
    st <- sin(tt)
    x <- (fry.data@rsAxes[1] / 2) * ct * ca - (fry.data@rsAxes[2] / 2) * st * sa
    y <- (fry.data@rsAxes[1] / 2) * ct * sa + (fry.data@rsAxes[2] / 2) * st * ca
  }
  
  #Create standardized Fry plot output to PDF
  pdf(file = out.file, width = 3, height = 4.5, useDingbats = FALSE, family = 'serif')
  par(omi = c(0,0,.5,0), mai = c(0, 0, .25, 0), xaxs = 'r', yaxs = 'r')
  layout(matrix(c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 2, 2, 2, 0), ncol = 5, byrow = TRUE))
  plot(x.coords, y.coords, axes = FALSE, asp = 1, pch = 19, , xlim = c(-1 * max.dim, max.dim), main = expression("Scaled Fry plot"),
    xlab = '', ylab = '', cex = .5, col = "gray")
  points(0, 0, pch = 3, col = "black")
  lines(x.ring, y.ring, col = "black", lwd = 1.5)
  polygon(x, y, border = "black", lwd = 1.5)
  plot(0,0,
    xlim = c(0,1), 
    ylim = c(0, 1), 
    main = expression("Parameters"), 
    axes = FALSE, pch = '', ylab = '')
  text(0, .965, "Sample size:", adj = 0)
  text(1, .965, fry.data@sampleSize, adj = 1)
  text(0, .813,  "Fabric ratio:", adj = 0)
  text(1, .813, round(fry.data@strainRatio, 2), adj = 1)
  text(0, .66, "Rake:", adj = 0)
  text(1, .66, round(fry.data@vectorMean, 0), adj = 1)
  text(0, .508, "Void area*:", adj = 0)
  text(1, .508, round(fry.data@meanObjectArea, 0), adj = 1)
  text(.5, .25, "*Central void area in pixles")
  if(normalized){
    mtext(paste("Normalized Fry results for plane: ", fry.data@sectionName, sep = ""), side = 3, line = 1, outer = TRUE)
  }else{
    mtext(paste("Standard Fry results for plane: ", fry.data@sectionName, sep = ""), side = 3, line = 1, outer = TRUE)
  }
  dev.off()
  return(fry.data)
}
