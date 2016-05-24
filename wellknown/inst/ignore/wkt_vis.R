#' Visualize well-known text area's on a map.
#'
#' @import rgeos whisker
#' @export
#'
#' @param x Input well-known text area (character)
#' @param zoom Zoom level, defaults to 6 (numeric)
#'
#' @examples \dontrun{
#' poly <- 'POLYGON((-111.06 38.84, -110.80 39.37, -110.20 39.17, -110.20 38.90,
#'      -110.63 38.67, -111.06 38.84))'
#' wkt_vis(poly)
#'
#' poly2 <- 'POLYGON((-125 38.4,-125 40.9,-121.8 40.9,-121.8 38.4,-125 38.4))'
#' wkt_vis(poly2)
#' }

wkt_vis <- function(x, zoom = 6)
{
  long = lat = group = NULL
  stopifnot(!is.null(x))
  stopifnot(is.character(x))

  poly_wkt <- readWKT(x)
  df <- fortify(poly_wkt)

  pts <- apply(df, 1, function(x) as.list(x[c('long','lat')]))
  centroid <- poly_wkt@polygons[[1]]@labpt
  rend <- whisker.render(map)
  foot <- sprintf(footer, centroid[2], centroid[1])
  res <- paste(rend, foot)
  tmpfile <- tempfile(pattern = 'spocc', fileext = ".html")
  write(res, file = tmpfile)
  browseURL(tmpfile)
}

map <- '
<!DOCTYPE html>
<html>
<head>
<meta charset=utf-8 />
<title>spocc WKT Viewer</title>
<meta name="viewport" content="initial-scale=1,maximum-scale=1,user-scalable=no" />
<script src="https://api.tiles.mapbox.com/mapbox.js/v1.6.4/mapbox.js"></script>
<link href="https://api.tiles.mapbox.com/mapbox.js/v1.6.4/mapbox.css" rel="stylesheet" />
<style>
  body { margin:0; padding:0; }
  #map { position:absolute; top:0; bottom:0; width:100%; }
</style>
</head>
<body>

<div id="map"></div>

<script>
var geojson = [
{
    "type": "Feature",
    "geometry": {
        "type": "Polygon",
        "coordinates": [
        [
            {{#pts}}
            [ {{long}}, {{lat}} ],
            {{/pts}}
        ]
    ]
    },
    "properties": {
        "title": "Polygon"
    }
}
];
'

footer <- '
L.mapbox.map("map", "examples.map-i86nkdio")
  .setView([ %s , %s ], 6)
  .featureLayer.setGeoJSON(geojson);
</script>

  </body>
  </html>
'
