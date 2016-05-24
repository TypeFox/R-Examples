# list to SpatialPolygons
as_SpatialPolygons <- function(x) {
  if (has_features(x)) {
    if (length(x$features) == 1) {
      polys <- list(make_poly(x$features[[1]]$geometry$coordinates[[1]]))
    } else {
      polys <- Filter(Negate(is.null), lapply(x$features, function(z) {
        switch(z$geometry$type,
               Polygon = {
                 make_poly(z$geometry$coordinates[[1]])
               },
               Point = {
                 message("only polygons supported")
                 NULL
               }
        )
      }))
    }
  } else {
    polys <- list(make_poly(x$coordinates[[1]]))
  }

  SpatialPolygons(list(Polygons(polys, length(polys))))
}

has_features <- function(x) "features" %in% names(x)

make_poly <- function(x) {
  Polygon(cbind(vapply(x, function(w) w[[1]], 1),
                vapply(x, function(w) w[[2]], 1)))
}
