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
#  install.packages("wellknown")

## ----eval=FALSE----------------------------------------------------------
#  install.packages("devtools")
#  devtools::install_github("ropensci/wellknown")

## ------------------------------------------------------------------------
library("wellknown")

## ------------------------------------------------------------------------
point <- list('type' = 'Point', 'coordinates' = c(116.4, 45.2, 11.1))
geojson2wkt(point)

## ------------------------------------------------------------------------
mp <- list(type = 'MultiPoint',
           coordinates = list( c(100.0, 3.101), c(101.0, 2.1), c(3.14, 2.18)
))
geojson2wkt(mp)

## ------------------------------------------------------------------------
st <- list(type = 'LineString',
            coordinates = list(c(0.0, 0.0, 10.0), c(2.0, 1.0, 20.0),
                              c(4.0, 2.0, 30.0), c(5.0, 4.0, 40.0)))
geojson2wkt(st, fmt = 0)

## ------------------------------------------------------------------------
multist <- list(type = 'MultiLineString',
      coordinates = list(
        list(c(0.0, -1.0), c(-2.0, -3.0), c(-4.0, -5.0)),
        list(c(1.66, -31023.5), c(10000.9999, 3.0), c(100.9, 1.1), c(0.0, 0.0))
      ))
geojson2wkt(multist)

## ------------------------------------------------------------------------
poly <- list(type = 'Polygon',
      coordinates = list(
        list(c(100.001, 0.001), c(101.12345, 0.001), c(101.001, 1.001), c(100.001, 0.001)),
        list(c(100.201, 0.201), c(100.801, 0.201), c(100.801, 0.801), c(100.201, 0.201))
))
geojson2wkt(poly)

## ------------------------------------------------------------------------
mpoly <-
  list(type = "MultiPolygon",
    coordinates = list(list(list(c(30, 20), c(45, 40), c(10, 40), c(30, 20))),
        list(list(c(15, 5), c(40, 10), c(10, 20), c(5 ,10), c(15, 5))))
  )
geojson2wkt(mpoly, fmt = 1)

## ------------------------------------------------------------------------
gmcoll <- list(type = 'GeometryCollection',
   geometries = list(
     list(type = "Point", coordinates = list(0.0, 1.0)),
     list(type = 'LineString', coordinates = list(c(-100.0, 0.0), c(-101.0, -1.0))),
     list(type = 'MultiPoint',
          'coordinates' = list(c(100.0, 3.101), c(101.0, 2.1), c(3.14, 2.18)))
  )
)
geojson2wkt(gmcoll, fmt = 0)

## ------------------------------------------------------------------------
library("jsonlite")
(json <- toJSON(list(type = "Point", coordinates = c(-105, 39))))

## ------------------------------------------------------------------------
geojson2wkt(json)

## ------------------------------------------------------------------------
str <- '{"type":["LineString"],"coordinates":[[0,0,10],[2,1,20],[4,2,30],[5,4,40]]}'
geojson2wkt(str)

## ----output.lines=1:10---------------------------------------------------
str <- "POINT (-116.4000000000000057 45.2000000000000028)"
wkt2geojson(str)

## ------------------------------------------------------------------------
wkt2geojson(str, feature = FALSE)

## ----output.lines=1:10---------------------------------------------------
str <- 'MULTIPOINT ((100.000 3.101), (101.000 2.100), (3.140 2.180))'
wkt2geojson(str, feature = FALSE)

## ----output.lines=1:10---------------------------------------------------
str <- "POLYGON ((100 0.1, 101.1 0.3, 101 0.5, 100 0.1), (103.2 0.2, 104.8 0.2, 100.8 0.8, 103.2 0.2))"
wkt2geojson(str, feature = FALSE)

## ----output.lines=1:10---------------------------------------------------
str <- "MULTIPOLYGON (((40 40, 20 45, 45 30, 40 40)),
    ((20 35, 45 20, 30 5, 10 10, 10 30, 20 35), (30 20, 20 25, 20 15, 30 20)))"
wkt2geojson(str, feature = FALSE)

## ----output.lines=1:10---------------------------------------------------
wkt2geojson("LINESTRING (0 -1, -2 -3, -4 5)", feature = FALSE)

## ----eval=FALSE----------------------------------------------------------
#  lint("POINT (1 2)")
#  #> [1] TRUE
#  lint("LINESTRING EMPTY")
#  #> [1] TRUE
#  lint("MULTIPOINT ((1 2), (3 4), (-10 100))")
#  #> [1] TRUE
#  lint("POLYGON((20.3 28.6, 20.3 19.6, 8.5 19.6, 8.5 28.6, 20.3 28.6))")
#  #> [1] TRUE
#  lint("MULTIPOLYGON (((30 20, 45 40, 10 40, 30 20)), ((15 5, 40 10, 10 20, 5 10, 15 5)))")
#  #> [1] TRUE
#  lint("POINT (1 2 3 4 5)")
#  #> [1] FALSE
#  lint("LINESTRING (100)")
#  #> [1] FALSE
#  lint("MULTIPOLYGON (((30 20, 45 40, 10 40, 30 20)), ((15 5, a b, 10 20, 5 10, 15 5)))")
#  #> [1] FALSE

