#' Fetches image from Earth Imagery API
#'
#' Calls NASA's Earth Imagery API and returns list with identification information and image.
#'
#' @param key Key for API authentication.
#' @param lon Longitud of coordinate position.
#' @param lat Latitud of coordinate position.
#' @param date In YYYY-MM-DD format. The API wil return the image that is closest to this date.
#' @param cloud_score Gives a score of percentage of cloud cover, via algorithm (see official documentation). Defaults to TRUE.
#' @param plot If TRUE will plot the image via generic plot function.
#' @param meta_only if TRUE will only download the meta data for the image.
#'
#' @importFrom png readPNG
#' @importFrom jsonlite fromJSON
#' @importFrom utils download.file
#' @importFrom graphics plot.new
#' @importFrom graphics rasterImage
#' @examples
#'\dontrun{
#' key <- "123key"
#' img <- earth_image(key, -100.31008, 25.66779, "2016-01-01")
#'}
#' @export
earth_image <- function(key, lon, lat, date, cloud_score = TRUE, plot = FALSE, meta_only = FALSE)
  {
  tryCatch({
    date <- as.Date(date)
  },
    error = function(e){
      stop("Date parameter must be YYYY-MM-DD")
    })

  # Validate a few things
  if(!is.numeric(lon)){
    stop("Lon parameter must be numeric")
  }
  if(!is.numeric(lat)){
    stop("Lat parameter must be numeric")
  }

  # Change true to NASA-True
  tryCatch({
    if(cloud_score){
      cloud_score <- "True"
    }else{
      cloud_score <- "False"
    }
  }, error = function(e){
    stop("Cloud_score parameter must be TRUE or FALSE")
  })

  # ojo con el "true"
  h <- "https://api.nasa.gov/planetary/earth/imagery?"
  query <- paste0(h,
                  "lon=", lon, "&",
                  "lat=", lat, "&",
                  "date=", date, "&",
                  "cloud_score=", cloud_score, "&",
                  "api_key=", key)

  # Download json object...
  s <- jsonlite::fromJSON(query)

  if("error" %in% names(s)){
    stop(cat(paste0("NASA API Error \n",
                    "The following is the output: ", s$error , "\n",
                    "You can use earth_asset to review availability for location"
                    )))
  }

  if(is.null(s$cloud_score)){
    s$cloud_score <- "NA"
  }

  image_data <- data.frame("Date" = s$date,
                           "URL" = s$url,
                           "CloudScore" = s$cloud_score,
                           "ID" = s$id)
  if(!meta_only){
    z <- tempfile()
    suppressMessages(
      download.file(s$url, z, mode = "wb", quiet = TRUE)
    )
    image <- png::readPNG(z)
    file.remove(z)
  }

  l <- list("image_data" = image_data,
            "image_png" = image)
  if(plot){
    plot.new()
    rasterImage(image,0,0,1,1)
  }
  return(l)
}
#' Plots the image to device
#'
#' To avoid S4 Classes and methods, this small wrapper simply plots an image from NASA.
#' If the purpose is to this interactively on one image, set the parameter plot = TRUE in \code{earth_image}
#'
#' @param image_png image downloaded using earth_image.
#' @seealso earth_image
#' @importFrom graphics plot.new
#' @importFrom graphics rasterImage
#' @examples
#'\dontrun{
#' key <- "123key"
#' img <- earth_image(key, -100.31008, 25.66779, "2016-01-01")
#' plot_earth_image(img$image_png)
#'}
#' @export
plot_earth_image <- function(image_png){
  plot.new()
  rasterImage(image_png, 0,0,1,1)
}
