RfPhi <-
function(my.par, out.file, section.title, weight.vec = TRUE, unit.area = "px"){
  #Caclculate the vector mean components
  if(weight.vec){
    mean.sin <- mean(my.par$m.eccentricity * sin(2 * my.par$m.rake.rad))
    mean.cos <- mean(my.par$m.eccentricity * cos(2 * my.par$m.rake.rad))
  } else{
    mean.sin <- mean(sin(2 * my.par$m.rake.rad))
    mean.cos <- mean(cos(2 * my.par$m.rake.rad))
  }
  #Calculate vector mean
  if(mean.sin > 0 & mean.cos > 0){
    v.mean.rad <- atan(mean.sin / mean.cos)
  }
  if(mean.cos < 0){
    v.mean.rad <- atan(mean.sin / mean.cos) + pi
  }
  if(mean.sin < 0 & mean.cos > 0){
    v.mean.rad <- atan(mean.sin / mean.cos) + 2 * pi
  }
  v.mean.rad <- v.mean.rad / 2
  v.mean.deg <-  v.mean.rad * (180 / pi)
  
  my.par$centered.rad <- my.par$m.rake.rad - v.mean.rad
  my.par$centered.rad[which(my.par$centered.rad < (-1 * pi / 2))] <- my.par$centered.rad[which(my.par$centered.rad < (-1 * pi / 2))] + pi
  my.par$centered.rad[which(my.par$centered.rad > (pi / 2))] <- my.par$centered.rad[which(my.par$centered.rad > (pi / 2))] - pi
  my.par$centered.deg <- my.par$centered.rad * (180 / pi)
  
  #Calculate harmonic mean
  h.mean <- (length(my.par$m.ratio)) / (sum(my.par$m.ratio^(-1)))
  
  #Calculate index of symmetry
  quad.a <- length(which(my.par$centered.deg <= 0 & my.par$m.ratio >= h.mean))
  quad.b <- length(which(my.par$centered.deg > 0 & my.par$m.ratio >= h.mean))
  quad.c <- length(which(my.par$centered.deg <= 0 & my.par$m.ratio < h.mean))
  quad.d <- length(which(my.par$centered.deg > 0 & my.par$m.ratio < h.mean))
  i.sym <- 1 - ((abs(quad.a - quad.b) + abs(quad.c - quad.d)) / length(my.par$m.rake.rad))
  
  #Create Rs destraining sequence
  rs.seq <- seq(from = 1.1, to = 16, by = .1)
  #Define chi squared object
  chi.sq <- vector()
  
  #Loop through each rs value
  for(j in seq(rs.seq)){
    phi <- vector()

r.s <- 1 / rs.seq[j]

#Destrain each object based in rs value
for(i in 1:length(my.par$m.rake.rad)){
  r.i <- my.par$m.ratio[i]
  theta <- my.par$centered.rad[i]
  
  phi[i] <- .5 * atan((2 * r.s * (r.i^2 - 1) * sin(2 * theta)) / ((r.i^2 + 1) * (r.s^2 - 1) + (r.i^2 - 1) * (r.s^2 + 1) * cos(2 * theta)))
  if(theta > 0 & phi[i] < 0){
    phi[i] <- phi[i] - pi / 2
  }
  if(theta < 0 & phi[i] > 0){
    phi[i] <- phi[i] + pi / 2
  }
}
phi[which(phi < 0)] <- phi[which(phi < 0)] + pi
phi[which(phi > pi)] <- phi[which(phi > pi)] - pi

n.bins <- floor(length(phi) / 5)
bin.breaks <- seq(from = 0, to = pi, length = n.bins + 1)
ex.freq <- length(phi) / n.bins

chi.stat <- vector()
for(k in 1:n.bins){
  chi.stat[k] <- (length(which(phi > bin.breaks[k] & phi <= bin.breaks[k + 1])) - ex.freq)^2 / ex.freq
}
chi.sq[j] <- sum(chi.stat)
  }
  my.r.s <- rs.seq[which(chi.sq == min(chi.sq))[1]]
  
  #Populate new class with data
  rfphi.data <- new("RFPHI", 
    sectionName = section.title,
    vectorMean = v.mean.deg,
    harmonicMean = h.mean,
    strainRatio = my.r.s,
    indexSymmetry = i.sym,
    sampleSize = length(my.par$m.rake.rad),
    meanObjectArea = sum(my.par$m.area) / length(my.par$m.rake.rad),
    rsAxes = c(sqrt((4 * my.r.s * (sum(my.par$m.area) / length(my.par$m.rake.rad))) / (pi)), 
      sqrt((4 * (sum(my.par$m.area) / length(my.par$m.rake.rad))) / (my.r.s * pi))),
    chiSquare = data.frame(rs.seq, chi.sq)
  )
  
  #Create standardized plot of sectional RfPhi data
  pdf(file = out.file, width = 6, height = 4.5, useDingbats = FALSE, family = 'serif')
  par(omi = c(0,0,.5,0))
  layout(matrix(c(1, 2, 2, 1, 3, 4), ncol = 3, byrow = TRUE))
  #Plot Rf/Phi graph
  plot(1, 1,
    main = expression(paste(R[f], " / ", phi, sep = "")), axes = FALSE,
    xlim = c(-90, 90),
    ylim = c(1, 20),
    pch = '',
    cex = .5,
    log='y',
    xlab = expression(paste("Recentered ", phi, sep = '')),
    ylab = expression(R[f]))
  #Add center and harmonic mean lines, points, and axes
  lines(c(1, 1), c(1, 20), col = 'gray')
  lines(c(-90, 90), c(h.mean, h.mean), col = 'gray')
  points(my.par$centered.deg,
    my.par$m.ratio,
    pch = 19,
    cex = .5)
  axis(1, at = c(-90, 0, 90), col = 'gray')
  axis(2, at = c(1, 2, 5, 10, 20), col = 'gray')
  #Plot Chi squared graph
  plot(rfphi.data@chiSquare, type = 'l', lwd = 1,
    xlab = expression("R"["s"]),
    ylab = expression(Chi^2),
    main = expression(paste("R"["s"], " vs. ", Chi^2, sep = '')), 
    axes = FALSE,
    ylim = c(0, ceiling(max(chi.sq))))
  axis(1, at = seq(from = 2, to = 16, by = 2), col = "gray")
  axis(2, at = c(0, round(max(chi.sq), 0)), col = 'gray')
 #Calculate nodes for Rs ellipse
  tt <- seq(0, 2 * pi, length = 72)
  sa <- sin(-1 * v.mean.rad)
  ca <- cos(-1 * v.mean.rad)
  ct <- cos(tt)
  st <- sin(tt)
  x <- sqrt(my.r.s)/2 * ct * ca - sqrt(1 / my.r.s) / 2 * st * sa
  y <- sqrt(my.r.s)/2 * ct * sa + sqrt(1 / my.r.s) / 2 * st * ca
  #Plot calculated Rs ellipse
  plot(0, 0, 
    pch = '', xlim = c(-2, 2), ylim = c(-2, 2),
    asp = 1, axes = FALSE,
    xlab = "Normalized ellipse", ylab = '',
    main = expression(paste(, R[s], " ellipse", sep = '')))
  #Add ellipse, reference line and center point
  polygon(x, y, border = 'black', lwd = .5)
  lines(c(-2, 2), c(0, 0), lwd = .5)
  points(0, 0, pch = 19, cex = .5)
  box(col = 'gray')
  #Plot list of parameters
  plot(0,0,
    xlim = c(0,1), 
    ylim = c(0, 1), 
    main = expression("Parameters"), 
    axes = FALSE, pch = '', ylab = '',
    xlab = paste("*Average area in squared", unit.area, sep = " "))
  text(0, .965, "Sample size:", adj = 0)
  text(1, .965, length(my.par$m.ratio), adj = 1)
  text(0, .813,  "Fabric ratio:", adj = 0)
  text(1, .813, my.r.s, adj = 1)
  text(0, .66, "Vector mean:", adj = 0)
  text(1, .66, round(v.mean.deg, 0), adj = 1)
  text(0, .508, "Harmonic mean:", adj = 0)
  text(1, .508, round(h.mean, 1), adj = 1)
  text(0, .355, "Symmetry:", adj = 0)
  text(1, .355, round(i.sym, 2), adj = 1)
  text(0, .203, "Min. chi sq:", adj = 0)
  text(1, .203, round(min(chi.sq)[1], 2), adj = 1)
  text(0, .05, "Mean area*:", adj = 0)
  text(1, .05, round(sum(my.par$m.area) / length(my.par$m.rake.rad), 1), adj = 1)
  #Add figure title
  mtext(substitute(paste("Sectional ", R[f], " / ", phi, " results for plane: ", z, sep = ''), list(z = section.title)), side = 3, cex = 1.25, line = 1, outer = TRUE)
  #Destroy device and return data
  dev.off()
  return(rfphi.data)
}
