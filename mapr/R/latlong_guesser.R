###### code adapted from the leaflet package - source at github.com/rstudio/leaflet
guess_latlon <- function(x, lat = NULL, lon = NULL) {
  nms <- names(x)

  if (is.null(lat)) {
    lats <- nms[grep(sprintf("^(%s)$", paste0(lat_options, collapse = "|")), nms, ignore.case = TRUE)]

    if (length(lats) == 1) {
      if (length(nms) > 2) {
        message("Assuming '", lats, "' is latitude")
      }
      names(x)[names(x) %in% lats] <- lat_var <-  "latitude"
    } else {
      stop("Couldn't infer latitude column, please specify with the 'lat' parameter",
           call. = FALSE)
    }
  } else {
    if (!any(names(x) %in% lat)) stop("'", lat, "' not found in your data", call. = FALSE)
    names(x)[names(x) %in% lat] <- lat_var <- "latitude"
  }

  if (is.null(lon)) {
    lngs <- nms[grep(sprintf("^(%s)$", paste0(lon_options, collapse = "|")), nms, ignore.case = TRUE)]

    if (length(lngs) == 1) {
      if (length(nms) > 2) {
        message("Assuming '", lngs, "' is longitude")
      }
      names(x)[names(x) %in% lngs] <- lon_var <- "longitude"
    } else {
      stop("Couldn't infer longitude column, please specify with 'lon' parameter",
           call. = FALSE)
    }
  } else {
    if (!any(names(x) %in% lon)) stop("'", lon, "' not found in your data", call. = FALSE)
    names(x)[names(x) %in% lon] <- lon_var <- "longitude"
  }

  structure(x, lat_var = lat_var, lon_var = lon_var)
}

lat_options <- c("lat", "latitude", "decimallatitude", "y")
lon_options <- c("lon", "lng", "long", "longitude", "decimallongitude", "x")
