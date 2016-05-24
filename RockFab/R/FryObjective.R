FryObjective <-
function(object.data, n.pass = 15, pie.step = 12, expansion = 1.5, pie.pts = 1, section.name, ave.piepts = FALSE, norm = TRUE){ 
  #Function to fit 2D ellipse to points modified from Michael Bedward
  EllipFit <- function(dat, npts = 180){
    #Determine least squares ellipse coeff.
    D1 <- cbind(dat$x * dat$x, dat$x * dat$y, dat$y * dat$y)
    D2 <- cbind(dat$x, dat$y, 1)
    S1 <- t(D1) %*% D1
    S2 <- t(D1) %*% D2
    S3 <- t(D2) %*% D2
    T <- -solve(S3) %*% t(S2)
    M <- S1 + S2 %*% T
    M <- rbind(M[3,] / 2, -M[2,], M[1,] / 2)
    evec <- eigen(M)$vec
    cond <- 4 * evec[1,] * evec[3,] - evec[2,]^2
    a1 <- evec[, which(cond > 0)]
    f <- c(a1, T %*% a1)
    
    #Determine semi-axes. Note a is not always > b
    a <- sqrt((2 * (f[1] * (.5 *f[5])^2 + f[3] * (.5 * f[4])^2 + f[6] * (.5 * f[2])^2 - 2 * (.5 * f[2]) * (.5 * f[4]) * (.5 * f[5]) - f[1] * f[3] * f[6])) / 
      (((.5 * f[2])^2 - f[1] * f[3]) * (sqrt((f[1] - f[3])^2 + 4 * (.5 * f[2])^2) - (f[1] + f[3]))))
    b <- sqrt((2 * (f[1] * (.5 *f[5])^2 + f[3] * (.5 * f[4])^2 + f[6] * (.5 * f[2])^2 - 2 * (.5 * f[2]) * (.5 * f[4]) * (.5 * f[5]) - f[1] * f[3] * f[6])) / 
      (((.5 * f[2])^2 - f[1] * f[3]) * (-1 * sqrt((f[1] - f[3])^2 + 4 * (.5 * f[2])^2) - (f[1] + f[3]))))
    
    #Find angle based on the inequalities of f[1] and f[3]; and axes a and b
    if(f[1] < f[3] & a > b){
      angle <- .5 * atan(f[2] / (f[1] - f[3]))
    }
    if(f[1] > f[3] & a > b){
      angle <- pi / 2 + .5 * atan(f[2] / (f[1] - f[3]))
    }
    if(f[1] < f[3] & a < b){
      angle <- pi / 2 + .5 * atan(f[2] / (f[1] - f[3]))
    }
    if(f[1] > f[3] & a < b){
      angle <- .5 * atan(f[2] / (f[1] - f[3]))
    }
    
    #Define true major and minor axes
    semi.maj <- ifelse(a > b, a, b)
    semi.min <- ifelse(b < a, b, a)
    
    #Determine XY plotting coords
    tt <- seq(from = 0, to = 2 * pi, length = npts)
    sa <- sin(angle)
    ca <- cos(angle)
    ct <- cos(tt)
    st <- sin(tt)
    x <- semi.maj * ct * ca - semi.min * st * sa
    y <- semi.maj * ct * sa + semi.min * st * ca
    
    #Calculate angle of semi.maj from 0 to pi clockwise from +x-axis
    angle <- ifelse((angle * -1) < 0 , (angle * -1) + pi, (angle * -1))
    
    #Create list of parameters
    my.ellipse <- list()
    my.ellipse[[1]] <- c(semi.maj, semi.min)
    my.ellipse[[2]] <- angle
    my.ellipse[[3]] <- data.frame(x, y)
    return(my.ellipse)
  }

  #Define empty objects for coordinates
  x.coords <- NULL
  y.coords <- NULL
  radii <- NULL
  
  #Determine object radius based on circle with equal area
  object.data$radius <- sqrt(object.data$m.area / pi)
  
  #Loop through each point to determine Fry cooordinates
  for(j in 1:length(object.data$m.cx)){
    x.coords <- c(x.coords, object.data$m.cx[j] - object.data$m.cx)
y.coords <- c(y.coords, object.data$m.cy[j] - object.data$m.cy)
radii <- c(radii, object.data$radius[j] + object.data$radius)
  }
  
  #Construct data frame of Fry coordinates and remove origin points
  my.data <- data.frame(x.coords, y.coords, radii)
  my.data <- my.data[which(my.data$x.coords != 0 & my.data$y.coords != 0),]
  
  #Determine center to point distances, normalize if defined
  my.data$dist <- sqrt(my.data$x.coords^2 + my.data$y.coords^2)
  if(norm){
    norm.dist <- my.data$dist / my.data$radii
    my.data$x.coords <- my.data$x.coords * (norm.dist / my.data$dist) * sqrt(mean(object.data$m.area)/pi)
    my.data$y.coords <- my.data$y.coords * (norm.dist / my.data$dist) * sqrt(mean(object.data$m.area)/pi)
    my.data$dist <- norm.dist * sqrt(mean(object.data$m.area)/pi)
  }
  
  #Sort by distance from origin
  my.data <- my.data[order(my.data$dist),]
  
  #Determine center to point angle
  my.data$angle <- rep(NA, times = length(my.data$x.coords))
  idx <- which(my.data$x.coords > 0 & my.data$y.coords <= 0)
  my.data[idx,]$angle <- abs(atan(my.data[idx,]$y.coords / my.data[idx,]$x.coords)) * (180 / pi)
  idx <- which(my.data$x.coords <= 0 & my.data$y.coords < 0)
  my.data[idx,]$angle <- abs(atan(my.data[idx,]$x.coords / my.data[idx,]$y.coords)) * (180 / pi) + 90
  idx <- which(my.data$x.coords < 0 & my.data$y.coords >= 0)
  my.data[idx,]$angle <- abs(atan(my.data[idx,]$y.coords / my.data[idx,]$x.coords)) * (180 / pi) + 180
  idx <- which(my.data$x.coords >= 0 & my.data$y.coords > 0)
  my.data[idx,]$angle <- abs(atan(my.data[idx,]$x.coords / my.data[idx,]$y.coords)) * (180 / pi) + 270
  
  #Loop through n.pass iterations of ellipse fitting
  for(k in 1:n.pass){
    #Create empty objects
x.coord <- vector()
y.coord <- vector()
min.point <- vector()

#On first iteration create equally spaced wedges
if(k == 1){
  wedges <- seq(from = 0, to = 360, by = pie.step)
  wedges <- wedges[-length(wedges)]
}
    
#Loop through each wedge increment and select Fry points
for(j in 1:length(wedges)){
  #On first iteration select points greater than last bin break and less then first bin.
  if(j == 1){
    temp.data <- my.data[which(my.data$angle <= wedges[j] | my.data$angle > wedges[length(wedges)]),]
  } else{
    temp.data <- my.data[which(my.data$angle <= wedges[j] & my.data$angle > wedges[j - 1]),]
  }
  #Conditional if mulitple pts are to be averaged into a moment
  if(ave.piepts){
    x.coord <- c(x.coord, mean(temp.data$x.coords[1:(pie.pts)], na.rm = TRUE))
    y.coord <- c(y.coord, mean(temp.data$y.coords[1:(pie.pts)], na.rm = TRUE))
  } else{
    x.coord <- c(x.coord, temp.data$x.coords[1:(pie.pts)])
    y.coord <- c(y.coord, temp.data$y.coords[1:(pie.pts)])
  }
  #Define minimum Fry point distance
  min.point <- c(min.point, temp.data$dist[1])
}
    
    #Create dataframe object for selected point coordinates and remove NA rows
dat <- data.frame(x = x.coord, y = y.coord)
dat <- dat[complete.cases(dat),]

#Fit ellipse to selected points and determine node locatons
dat <- EllipFit(dat = dat, npts = length(wedges) + 1)
node.data <- dat[[3]]

#Determine angle intervals between nodes to degine wedge sequence
    node.data$angle <- rep(NA, times = length(node.data$x))

    idx <- which(node.data$x > 0 & node.data$y <= 0)
    node.data[idx,]$angle <- abs(atan(node.data[idx,]$y / node.data[idx,]$x)) * (180 / pi)

    idx <- which(node.data$x <= 0 & node.data$y < 0)
    node.data[idx,]$angle <- abs(atan(node.data[idx,]$x / node.data[idx,]$y)) * (180 / pi) + 90

    idx <- which(node.data$x < 0 & node.data$y >= 0)
    node.data[idx,]$angle <- abs(atan(node.data[idx,]$y / node.data[idx,]$x)) * (180 / pi) + 180

    idx <- which(node.data$x >= 0 & node.data$y > 0)
    node.data[idx,]$angle <- abs(atan(node.data[idx,]$x / node.data[idx,]$y)) * (180 / pi) + 270
    
#Sort new wedge intervals
    wedges <- sort(node.data$angle)
wedges <- wedges[-length(wedges)]
  }
  
  #Populate FRY class with fitted ellipse parameters
  fry.data <- new("FRY")
  fry.data@rsAxes <- c(2 * dat[[1]][1], 2 * dat[[1]][2])
  fry.data@strainRatio <- dat[[1]][1] / dat[[1]][2]
  fry.data@meanObjectArea <- dat[[1]][1] * dat[[1]][2] * pi
  fry.data@vectorMean <-  dat[[2]] * (180 / pi)
  fry.data@sectionName <- section.name
  fry.data@sampleSize <- length(object.data$m.cx)
  fry.data@fryParams <- my.data
  fry.data@voidScale <- expansion * max(min.point, na.rm = TRUE)
  
  return(fry.data)
}
