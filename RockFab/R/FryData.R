FryData <-
function(object.data, pie.step = 5, expansion = 2, section.name){
  #Define empty objects for coordinates
  x.coords <- NULL
  y.coords <- NULL
  
  #Loop through each point to determine Fry cooordinates
  for(j in 1:length(object.data$m.cx)){
    x.coords <- c(x.coords, object.data$m.cx[j] - object.data$m.cx)
y.coords <- c(y.coords, object.data$m.cy[j] - object.data$m.cy)
  }
  #Construct data frame of Fry coordinates and remove origin points
  my.data <- data.frame(x.coords, y.coords)
  my.data <- my.data[which(my.data$x.coords != 0 & my.data$y.coords != 0),]
  
  #Determine center to point distances and sort by distance
  my.data$dist <- sqrt(my.data$x.coords^2 + my.data$y.coords^2)
  my.data <- my.data[order(my.data$dist),]
  
  #Determine center to point angle
  my.data$angle <- atan(my.data$y.coords / my.data$x.coords) * (180 / pi)
  my.data$angle[which(my.data$angle < 0)] <- my.data$angle[which(my.data$angle < 0)] + 180
  
  #Create sequence of wedges to approixmate central void distance
  wedges <- seq(from = 0, to = 180, by = pie.step)
  #Search through wedges for largest distance given by the closest point
  min.point <- NULL
  for(j in 1:(length(wedges) - 1)){
    min.temp <- min(my.data$dist[which(my.data$angle >= wedges[j] & my.data$angle < wedges[j + 1])])
    min.point <- c(min.point, min.temp)
  }

  #Populate new object with available data
  fry.data <- new("FRY", 
    sectionName = section.name,
    sampleSize = length(object.data$m.cx),
fryParams = my.data,
    voidScale = expansion * max(min.point)
  )
  return(fry.data)
}
