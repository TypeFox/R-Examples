#' Create a spatial points dataframe from a spocc search
#'
#' @export
#'
#' @param x The resuslts of a spocc search called by occ()
#' @param coord_string A valid EPGS cooridate string from the sp package, the default is WSGS 84
#' @param just_coords Return data frame with specios names and provenance or just a spatial points
#' object, which is the default.
#'
#' @details This function will return either a spatial points dataframe or spatial points object.
#' Conversion to spatial points objects allows spocc searches to interact with other spatial
#' data sources. More coordinate system codes can be found at the EPGS registry:
#' \url{http://www.epsg-registry.org/}
#'
#' @examples \dontrun{
#' ### See points on a map
#' library("maptools")
#' library("spocc")
#' data(wrld_simpl)
#' plot(wrld_simpl[wrld_simpl$NAME == "United States", ], xlim = c(-70, -60))
#' out <- occ(query = "Accipiter striatus", from = c("vertnet", "gbif"), limit = 50)
#' xx <- occ2sp(out, just_coords = TRUE)
#' points(xx, col = 2)
#' }
occ2sp <- function(x, coord_string = "+proj=longlat +datum=WGS84", just_coords = FALSE) {
  x <- spocc::occ2df(x)

  # numerics
  x$longitude <- as.numeric(x$longitude)
  x$latitude <- as.numeric(x$latitude)

  # remove NA rows
  x <- x[complete.cases(x$latitude, x$longitude), ]

  # check valid coords
  index <- 1:dim(x)[1]
  index <- index[(x$longitude < 180) & (x$longitude > -180) & !is.na(x$longitude)]
  index <- index[(x$latitude[index] < 90) & (x$latitude[index] > -90) & !is.na(x$latitude[index])]

  spobj <- sp::SpatialPoints(as.matrix(x[index,c('longitude', 'latitude')]), proj4string = sp::CRS(coord_string))

  sp_df <- sp::SpatialPointsDataFrame(spobj, data = data.frame(x[index, c('name', "prov")]))
  if (just_coords) spobj else sp_df
}
