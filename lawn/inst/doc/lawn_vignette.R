## ----echo=FALSE----------------------------------------------------------
library("knitr")
hook_output <- knitr::knit_hooks$get("output")
knitr::knit_hooks$set(output = function(x, options) {
   lines <- options$output.lines
   if (is.null(lines)) {
     return(hook_output(x, options))  # pass to default hook
   }
   x <- unlist(strsplit(x, "\n"))
   more <- "..."
   if (length(lines)==1) {        # first n lines
     if (length(x) > lines) {
       # truncate the output, but add ....
       x <- c(head(x, lines), more)
     }
   } else {
     x <- c(if (abs(lines[1])>1) more else NULL,
            x[lines],
            if (length(x)>lines[abs(length(lines))]) more else NULL
           )
   }
   # paste these lines together
   x <- paste(c(x, ""), collapse = "\n")
   hook_output(x, options)
 })

knitr::opts_chunk$set(
  comment = "#>",
  collapse = TRUE,
  warning = FALSE,
  message = FALSE
)

## ----eval=FALSE----------------------------------------------------------
#  devtools::install_github("rstudio/leaflet")

## ----eval=FALSE----------------------------------------------------------
#  install.packages("lawn")

## ----eval=FALSE----------------------------------------------------------
#  devtools::install_github("ropensci/lawn")

## ------------------------------------------------------------------------
library("lawn")

## ------------------------------------------------------------------------
lawn_point(c(-74.5, 40))

## ------------------------------------------------------------------------
rings <- list(list(
  c(-2.275543, 53.464547),
  c(-2.275543, 53.489271),
  c(-2.215118, 53.489271),
  c(-2.215118, 53.464547),
  c(-2.275543, 53.464547)
))
lawn_polygon(rings)

## ------------------------------------------------------------------------
lawn_count(polygons = lawn_data$polygons_count, points = lawn_data$points_count)

## ------------------------------------------------------------------------
lawn_average(polygons = lawn_data$polygons_average,
             points = lawn_data$points_average,
             field = 'population')

## ------------------------------------------------------------------------
from <- '{
 "type": "Feature",
 "properties": {},
 "geometry": {
   "type": "Point",
   "coordinates": [-75.343, 39.984]
 }
}'
to <- '{
  "type": "Feature",
  "properties": {},
  "geometry": {
    "type": "Point",
    "coordinates": [-75.534, 39.123]
  }
}'

## ------------------------------------------------------------------------
lawn_distance(from, to)

## ------------------------------------------------------------------------
lawn_random(n = 2)

## ------------------------------------------------------------------------
lawn_random(n = 5)

## ------------------------------------------------------------------------
gr_position()

## ------------------------------------------------------------------------
gr_point(2)

## ------------------------------------------------------------------------
gr_polygon(n = 1, vertices = 5, max_radial_length = 5)

## ------------------------------------------------------------------------
dat <- lawn_data$points_average
lawn_sample(dat, 1)

## ------------------------------------------------------------------------
lawn_sample(dat, 2)

## ------------------------------------------------------------------------
lawn_sample(dat, 3)

## ------------------------------------------------------------------------
lawn_extent(lawn_data$points_average)

## ------------------------------------------------------------------------
lawn_within(lawn_data$points_within, lawn_data$polygons_within)

## ------------------------------------------------------------------------
dat <- '{
 "type": "Feature",
 "properties": {},
 "geometry": {
     "type": "Polygon",
     "coordinates": [[
       [-112.072391,46.586591],
       [-112.072391,46.61761],
       [-112.028102,46.61761],
       [-112.028102,46.586591],
       [-112.072391,46.586591]
     ]]
   }
}'
lawn_buffer(dat, 1, "miles")

## ------------------------------------------------------------------------
dat <- '{
  "type": "FeatureCollection",
  "features": [
    {
      "type": "Feature",
      "properties": {
        "population": 200
      },
      "geometry": {
        "type": "Point",
        "coordinates": [10.724029, 59.926807]
      }
    },
      {
      "type": "Feature",
      "properties": {
        "population": 600
      },
      "geometry": {
        "type": "Point",
        "coordinates": [10.715789, 59.904778]
      }
    }
  ]
}'
lawn_extent(dat)

## ----eval=FALSE----------------------------------------------------------
#  dat <- '{
#    "type": "FeatureCollection",
#    "features": [
#      {
#        "type": "Feature",
#        "properties": {
#          "population": 200
#        },
#        "geometry": {
#          "type": "Point"
#        }
#      },
#        {
#        "type": "Feature",
#        "properties": {
#          "population": 600
#        },
#        "geometry": {
#          "type": "Point",
#          "coordinates": [10.715789, 59.904778]
#        }
#      }
#    ]
#  }'
#  lawn_extent(dat, lint = TRUE)
#  
#  #> Error: Line 1 - "coordinates" property required

## ----eval=FALSE----------------------------------------------------------
#  view(lawn_data$points_average)

## ----eval=FALSE----------------------------------------------------------
#  lawn_sample(lawn_data$points_average, 2) %>% view()

