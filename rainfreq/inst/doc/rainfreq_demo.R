## ----, echo = FALSE, message = FALSE-------------------------------------
knitr::opts_chunk$set(
  comment = "#>",
  error = FALSE,
  tidy = FALSE)

## ----, message = FALSE---------------------------------------------------
require(rainfreq)
require(RCurl)
require(SDMTools)
require(raster)
require(maps)

## ------------------------------------------------------------------------
x_se <- extract_freq()
print(x_se)

## ------------------------------------------------------------------------
x_mw <- extract_freq(region_name = "mw", storm_RP = 1000, storm_duration = "48h")
print(x_mw)

## ------------------------------------------------------------------------
x_hi <- extract_freq(region_name = "hi", storm_RP = 10, storm_duration = "6h")
print(x_hi)

## ------------------------------------------------------------------------
data(rain_max_usa)
head(rain_max_usa)
data(rain_max_world)
head(rain_max_world)

## ------------------------------------------------------------------------
x_se <- x_se * 0.001
x_mw <- x_mw * 0.001
x_hi <- x_hi * 0.001

## ------------------------------------------------------------------------
# southeast
plot(x_se, breaks = c(6, 9, 12, 15, 18), 
     col = c("red", "yellow", "green", "blue"), 
     main = "100-year 24-hour Rainfall for Southeast USA (inches)")
map('state', region = c('florida', 'arkansas', 'louisiana', 'mississippi', 
                        'alabama', 'georgia'), add = TRUE)

## ------------------------------------------------------------------------
# midwest
plot(x_mw, breaks = c(2, 5, 10, 15, 20), 
     col = c("red", "yellow", "green", "blue"),
     main = "1000-year 48-hour Rainfall for Midwest USA (inches)")
map('state', region = c('colorado', 'north dakota', 'south dakota', 'nebraska', 
                        'oklahoma', 'minnesota', 'iowa', 'missouri', 
                        'wisconsin', 'michigan'), add = TRUE)

## ------------------------------------------------------------------------
# hawaii
plot(x_hi, breaks = c(1, 3, 6, 9, 12), 
     col = c("red", "yellow", "green", "blue"),
     main = "10-year 6-hour Rainfall for Hawaii (inches)")

