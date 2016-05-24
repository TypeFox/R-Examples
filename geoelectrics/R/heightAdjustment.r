#' Adjust Profile Height
#' 
#' Adjusts the height of a single profile.
#' GPS measurement heights might differ otherwise.
#' 
#' @param Profile a single Profile.
#' @param deltaMeter positive or negative value.
#' @export
#' @seealso \code{\link{GpsCoordinates-class}}, \code{\link{Profile-class}}
#' @examples 
#' # p3 <- new("Profile",
#' #            title = "Profile 3",
#' #            xyzData = 
#' #              new("XyzData", 
#' #                  address = "../example/xyzFiles/p3_DipolDipol_S-N.xyz"),
#' #            rawData = 
#' #              new("RawData",
#' #                  address = "../example/rawdata/p3_DipolDipol_S-N.dat"),
#' #            measurementType = "DipolDipol",
#' #            gpsCoordinates = 
#' #              new("GpsCoordinates",
#' #                  address = "../example/gps/p3.txt")) 
#'
#' # p3 <- heightAdjustment(p3, -10)
heightAdjustment <- function(Profile, deltaMeter) {
  Profile@xyzData@heightAdaption$depth <- 
    Profile@xyzData@heightAdaption$depth + deltaMeter
  return(Profile)
}