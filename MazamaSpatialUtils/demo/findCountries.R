# findCountries demo

library(MazamaSpatialUtils)
library(sp)

# Vector of lons and lats
lons <- seq(5,25,2)
lats <- seq(30,50,2)

# Get country names
countryNames <- getCountryName(lons, lats)
print(countryNames)

# Get all information in the dataset
countryDF <- getCountryName(lons, lats, allData=TRUE)
print(countryDF)

# Subset the SpatialPolygonsDataFrame to only include our countries
countryMask <- SimpleCountries@data$countryName %in% countryNames

# Plot the country polygons
plot(SimpleCountries[countryMask,],col='gray90',border='gray70')
# Add countries from the 'maps' package
maps::map('world',col='gray80',add=TRUE)
# Add our points in red
points(lons,lats,pch=16,col='red')
# Add text to the right
countryText <- ifelse(is.na(countryNames),'water',paste0(countryDF$countryCode,' = ',countryDF$countryName))
text(lons,lats,countryText,pos=4)
# Add a title
title('Country Codes and Names')

