#' Generate reliable elevation data.
#'
#' Using the latitude and longitude columns of the supplied \emph{formatted}
#' data, a vector of elevation values is returned of the same length. If no
#' elevation data files exist within the working directory, files are first
#' downloaded. Note that NAs in the data will return corresponding NAs in the
#' corrected elevation.
#'
#' @param data a dataset with longitude ("lng") and lattitude ("lon") columns.
#' @param country character string; the country to which the data pertain, given
#'   as an ISO3 code (see \code{raster::getData("ISO3")})
#'
#' @return a vector of elevation values. If there is an error at any stage, a
#'   vector of NAs is returned.
#'
#' @examples \dontrun{
#' data(ridedata)
#'
#' ## When run the first time, geographical data will need to be downloaded.
#' ridedata$elevation.corrected <- elevation_correct(ridedata, "GBR")
#'
#' ## A Bland-Altman-type plot.
#' difference <- ridedata$elevation.m - ridedata$elevation.corrected
#' plot(difference ~ ridedata$timer.min, cex = 0.2, ylab = "raw minus corrected")
#' m <- mean(difference, na.rm = TRUE); stdev <- sd(difference, na.rm = TRUE)
#' abline(h = c(m + c(-stdev, 0, stdev)), lty = c(1, 2, 1), col = "red")
#' }
#'
#' @seealso \code{\link{download_elev_data}}.
#'
#' @export
elevation_correct <- function(data, country) {
  if (!requireNamespace("raster", quietly = TRUE))
    stop("Raster package not installed.", call. = FALSE)
  stopifnot("lat" %in% names(data), # Nice error.
            "lon" %in% names(data))

  position_data <- data.frame(lng = data$lon, lat = data$lat)

  exit <- 0
  if (!country %grep_in% list.files()) exit <- download_elev_data(country)

  if (as.logical(exit)) {
    message("Error, returning NAs.")
    return(rep_len(NA, nrow(data)))
  }

  e <- raster::getData("alt", country = country, download = FALSE)
  # Extract interpolated elevation data from raster object.
  raster::extract(e, position_data, method = "bilinear")
}

#' Download geographical elevation data.
#'
#' Downloads elevation data files to the working directory for use with
#' \code{\link{elevation_correct}}. Requires package \code{raster} to be
#' installed.
#'
#' @param country character string; the ISO3 country code (see
#'   \code{raster::getData("ISO3")}) for which to download the data. If "all",
#'   then all available data is downloaded - this may take some time.
#'
#' @return nothing, files are downloaded to the working directory.
#'
#' @seealso \code{\link{elevation_correct}}.
#'
#' @export
download_elev_data <- function(country = "all") {
  if (country == "all") {
    response <- readline(paste(
      "This will download ALL available elevation data, and will take some time. Continue [Y/N]? "
    ))
    if (!grepl(pattern = "y", response, ignore.case = TRUE))
      stop("Download aborted.", call. = FALSE)

    mx <- raster::getData("ISO3")
    for (i in mx[, 1]) {
      message(rep("-", times = 50))
      message(paste(mx[match(i, mx), 2]))
      message(rep("=", times = 50))
      tryCatch(raster::getData("alt", country = i, download = TRUE, mask = TRUE),
               error   = function(e) message(paste(i, "not available")),
               finally = message("done"))
    }
  } else {
    if (country %in% raster::getData("ISO3")[, 1]) {
      message(paste("Downloading elevation data for", country, "to", getwd()))
      test <- try(suppressWarnings(
        raster::getData("alt", country = country, download = TRUE, mask = TRUE)
      ))
    } else {
      message("Invalid country argument or no data available.")
      return(1)
    }
    if (inherits(test, "try-error")) {
      message("No internet connection.")
      return(1)
    } else
      return(0)
  }
}
