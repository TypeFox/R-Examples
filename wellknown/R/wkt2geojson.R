#' Convert WKT to GeoJSON-like objects.
#'
#' @export
#'
#' @param str A GeoJSON-like object representing a Point, LineString, Polygon, MultiPolygon, etc.
#' @param fmt Number of digits to display after the decimal point when formatting
#' coordinates.
#' @param feature (logical) Make a feature geojson object. Default: TRUE
#' @details Should be robust against a variety of typing errors, including extra spaces between
#' coordinates, no space between WKT type and coordinates. However, some things won't pass,
#' includingn lowercase WKT types, no spaces between coordinates.
#' @seealso \code{\link{geojson2wkt}}
#' @examples
#' # point
#' str <- "POINT (-116.4000000000000057 45.2000000000000028)"
#' wkt2geojson(str)
#' wkt2geojson(str, feature=FALSE)
#'
#' # multipoint
#' str <- 'MULTIPOINT ((100.000 3.101), (101.000 2.100), (3.140 2.180))'
#' wkt2geojson(str, fmt = 2)
#' wkt2geojson(str, fmt = 2, feature=FALSE)
#'
#' # polygon
#' str <- "POLYGON ((100 0.1, 101.1 0.3, 101 0.5, 100 0.1),
#'    (103.2 0.2, 104.8 0.2, 100.8 0.8, 103.2 0.2))"
#' wkt2geojson(str)
#' wkt2geojson(str, feature=FALSE)
#'
#' # multipolygon
#' str <- "MULTIPOLYGON (((40 40, 20 45, 45 30, 40 40)),
#'    ((20 35, 45 20, 30 5, 10 10, 10 30, 20 35), (30 20, 20 25, 20 15, 30 20)))"
#' wkt2geojson(str)
#' wkt2geojson(str, feature=FALSE)
#'
#' # linestring
#' str <- "LINESTRING (100.000 0.000, 101.000 1.000)"
#' wkt2geojson(str)
#' wkt2geojson(str, feature=FALSE)
#' wkt2geojson("LINESTRING (0 -1, -2 -3, -4 5)")
#' wkt2geojson("LINESTRING (0 1 2 3, 4 5 6 7)")
#'
#' # multilinestring
#' str <- "MULTILINESTRING ((30 1, 40 30, 50 20)(10 0, 20 1))"
#' wkt2geojson(str)
#'
#' str <- "MULTILINESTRING (
#'    (-105.0 39.5, -105.0 39.5, -105.0 39.5, -105.0 39.5, -105.0 39.5, -105.0 39.5),
#'    (-105.0 39.5, -105.0 39.5, -105.0 39.5),
#'    (-105.0 39.5, -105.0 39.5, -105.0 39.5, -105.0 39.5, -105.0 39.5),
#'    (-105.0 39.5, -105.0 39.5, -105.0 39.5, -105.0 39.5))"
#' wkt2geojson(str)
#'
#' # Geometrycollection
#' str <- "GEOMETRYCOLLECTION (
#'    POINT (0 1),
#'    LINESTRING (-100 0, -101 -1),
#'    POLYGON ((100.001 0.001, 101.1235 0.0010, 101.001 1.001, 100.001 0.001),
#'            (100.201 0.201, 100.801 0.201, 100.801 0.801, 100.201 0.201)),
#'    MULTIPOINT ((100.000 3.101), (101.0 2.1), (3.14 2.18)),
#'    MULTILINESTRING ((0 -1, -2 -3, -4 -5),
#'                    (1.66 -31023.50 1.10, 10001.0 3.0 2.2, 100.9 1.1 3.3, 0.0 0.0 4.4)),
#'    MULTIPOLYGON (((100.001 0.001, 101.001 0.001, 101.001 1.001, 100.001 0.001),
#'                  (100.201 0.201, 100.801 0.201, 100.801 0.801, 100.201 0.201)),
#'                  ((1 2 3 4, 5 6 7 8, 9 10 11 12, 1 2 3 4))))"
#' wkt2geojson(str)

wkt2geojson <- function(str, fmt = 16, feature = TRUE){
  type <- get_type(str)
  res <- switch(type,
         Point = load_point(str, fmt, feature),
         Multipoint = load_multipoint(str, fmt, feature),
         Polygon = load_polygon(str, fmt, feature),
         Multipolygon = load_multipolygon(str, fmt, feature),
         Linestring = load_linestring(str, fmt, feature),
         Multilinestring = load_multilinestring(str, fmt, feature),
         Geometrycollection = load_geometrycollection(str, fmt, feature)
  )
  structure(res, class = "geojson")
}

types <- c("POINT",'MULTIPOINT',"POLYGON","MULTIPOLYGON",
           "LINESTRING","MULTILINESTRING","GEOMETRYCOLLECTION",
           "TRIANGLE","CIRCULARSTRING","COMPOUNDCURVE")

get_type <- function(x, ignore_case=FALSE){
  type <- cw(types[sapply(types, grepl, x = x, ignore.case = ignore_case)], onlyfirst = TRUE)
  if (length(type) > 1) {
    grep(tolower(strextract(x, "[A-Za-z]+")), type, ignore.case = TRUE, value = TRUE)
  } else {
    type
  }
}

load_point <- function(str, fmt = 16, feature = TRUE){
  str_coord <- gsub("POINT|\\(|\\)", "", str)
  coords <- strsplit(gsub("[[:punct:]]$", "", str_trim_(str_coord)), "\\s")[[1]]
  coords <- nozero(coords)
  # iffeat('Point', as.numeric(coords, fmt), feature)
  iffeat('Point', format_num(coords, fmt), feature)
}

format_num <- function(x, fmt) {
  sprintf(paste0("%.", fmt, "f"), as.numeric(x))
}

load_multipoint <- function(str, fmt = 16, feature = TRUE){
  str_coord <- str_trim_(gsub("MULTIPOINT\\s?", "", str))
  str_coord <- gsub("^\\(|\\)$", "", str_coord)
  str_coord <- strsplit(str_coord, "\\),")[[1]]
  coords <- unname(sapply(str_coord, function(z){
    pairs <- strsplit(strsplit(gsub("\\(|\\)", "", str_trim_(z)), ",|,\\s")[[1]], "\\s")
    lapply(pairs, function(x) {
      format_num(nozero(x), fmt)
    })
  }))
  iffeat('MultiPoint', coords, feature)
}

load_polygon <- function(str, fmt = 16, feature = TRUE){
  str_coord <- str_trim_(gsub("POLYGON\\s?", "", str))
  str_coord <- gsub("^\\(|\\)$", "", str_coord)
  str_coord <- strsplit(str_coord, "\\),")[[1]]
  coords <- lapply(str_coord, function(z){
    pairs <- strsplit(strsplit(gsub("\\(|\\)", "", str_trim_(z)), ",|,\\s")[[1]], "\\s")
    lapply(pairs, function(x) {
      format_num(nozero(x), fmt)
    })
  })
  iffeat('Polygon', coords, feature)
}

load_multipolygon <- function(str, fmt = 16, feature = TRUE){
  str <- gsub("\n", "", str)
  str_coord <- str_trim_(gsub("MULTIPOLYGON\\s?", "", str))
  str_coord <- gsub("^\\(|\\)$", "", str_coord)
  str_coord <- strsplit(str_coord, "\\)),")[[1]]
  coords <- lapply(str_coord, function(z){
    pairs <- strsplit( gsub("\\(|\\)", "", strsplit(str_trim_(z), "\\),")[[1]]), ",|,\\s")
    lapply(pairs, function(zz){
      unname(lapply(sapply(str_trim_(zz), strsplit, split = "\\s"), function(x) {
        format_num(nozero(x), fmt)
      }))
    })
  })
  iffeat('MultiPolygon', coords, feature)
}

load_linestring <- function(str, fmt = 16, feature = TRUE){
  str_coord <- str_trim_(gsub("LINESTRING\\s?", "", str))
  str_coord <- gsub("^\\(|\\)$", "", str_coord)
  str_coord <- strsplit(str_coord, "\\),")[[1]]
  coords <- lapply(str_coord, function(z){
    pairs <- strsplit(strsplit(gsub("\\(|\\)", "", str_trim_(z)), ",|,\\s")[[1]], "\\s")
    lapply(pairs, function(x) {
      format_num(nozero(x), fmt)
    })
  })[[1]]
  iffeat('LineString', coords, feature)
}

load_multilinestring <- function(str, fmt = 16, feature = TRUE){
  str_coord <- str_trim_(gsub("MULTILINESTRING\\s?", "", str))
  str_coord <- gsub("^\\(|\\)$", "", str_coord)
  str_coord <- strsplit(str_coord, "\\),|\\)\\(")[[1]]
  coords <- lapply(str_coord, function(z){
    pairs <- strsplit(strsplit(str_trim_(gsub("\\(|\\)", "", str_trim_(z))), ",|,\\s")[[1]], "\\s")
    lapply(pairs, function(x) {
      format_num(nozero(x), fmt)
    })
  })
  iffeat('MultiLineString', coords, feature)
}

load_geometrycollection <- function(str, fmt = 16, feature = TRUE){
  str_coord <- str_trim_(gsub("GEOMETRYCOLLECTION\\s?", "", gsub("\n", "", str)))
  str_coord <- gsub("^\\(|\\)$", "", str_coord)
  matches <- noneg(sort(sapply(types, regexpr, text = str_coord)))
  out <- list()
  for (i in seq_along(matches)) {
    end <- if (i == length(matches)) nchar(str_coord) else matches[[i + 1]] - 1
    strg <- substr(str_coord, matches[[i]], end)
    out[[ i ]] <- get_load_fxn(tolower(names(matches[i])))(strg, fmt)
  }
  list(type = 'GeometryCollection', geometries = out)
}

get_load_fxn <- function(type){
  switch(type,
         point = load_point,
         multipoint = load_multipoint,
         linestring = load_linestring,
         multilinestring = load_multilinestring,
         polygon = load_polygon,
         multipolygon = load_multipolygon,
         geometrycollection = load_geometrycollection)
}

noneg <- function(x) x[!x < 0]

iffeat <- function(type, cd, feature){
  tmp <- list(type = type, coordinates = cd)
  if (feature) {
    list(type = "Feature", geometry = tmp)
  } else {
    tmp
  }
}
