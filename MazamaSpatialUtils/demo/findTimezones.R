# findTimezones demo

library(MazamaSpatialUtils)
library(sp)

# Vector of lons and lats
lons <- seq(-120,-60,5)
lats <- seq(20,80,5)

# Get Olson timezone names
timezones <- getTimezone(lons, lats)
print(timezones)

# Get all information in the dataset
timezoneDF <- getTimezone(lons, lats, allData=TRUE)
print(timezoneDF)

# Subset the SpatialPolygonsDataFrame to only include our timezones
timezoneMask <- SimpleTimezones@data$timezone %in% timezones

# Plot the timezone polygons
plot(SimpleTimezones[timezoneMask,],col='gray90',border='gray70')
# Add countries from the 'maps' package
map::map('world',add=TRUE,col='gray80')
# Add our points in red
points(lons,lats,pch=16,col='red')
# Add text to the right
timezoneText <- ifelse(is.na(timezones),'water',timezones)
text(lons,lats,timezoneText,pos=4)
# Add a title
title('Timezones in North America')

