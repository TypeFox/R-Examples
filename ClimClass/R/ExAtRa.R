NULL
#' @description Calculates Extra-Atmospheric Radiation. Called by function \code{\link{arid}} for Thornthwaite's index.
#' 
#' @param DOY day of the year.
#' @param latitude latitude in degrees (negative for S emishpere).
#' @param Gsc solar constant in MJ m-2 min-1 (default: 0.0820).
#' @param unit unit for solar radiation. Accepted values are "mm" and "MJ".
#' @param T temperature in degrees C. Default is 12.
#'
#' @title Extra-Atmospheric Radiation
#' 
#' @author Emanuele Eccel
#'
#' @return The daily extra-atmospheric solar radiation energy, espressed either in MJ or in mm of evaporated water.
#'
#' @details 
#' If \code{unit} = "mm", the calculated value represents the water height evaporated by solar radiation, calculated by the latent heat for vaporization. Otherwise (\code{unit} = "MJ") output is the solar radiation energy in MJ.
#' Temperature \code{T} is used only for the assessment of latent heat of vaporization, when \code{unit} = "mm". 
#' 
#' @export
#' 
#' @examples
#' 
#' data(Trent_climate)
#' # creates a vector with middle days for every month in a year
#' quinci   <- paste(15,"/",1:12,"/",2014,sep="")
#' posixlt  <- strptime(quinci, format="%d/%m/%Y")
#' yDay     <- posixlt$yday+1  # field yday starts from 0
#' latitude<-46  
#' 
#' # generates 12 values, one for each month
#' coeff_rad<- ExAtRa(DOY=yDay,latitude=latitude, unit="mm")
#' 
#'   @seealso \code{\link{arid}}



ExAtRa     <- function(DOY,latitude,Gsc= 0.0820, unit="mm", T=12)
  
{
  phi<-latitude*pi/180
  delta <- 0.4093 * sin(2*pi*(284+DOY)/365) # declination
  dr    <- 1 + 0.033*cos(2*pi*DOY/365)      # relative distance Earth - Sun
  omegaS  <- acos(-tan(phi)*tan(delta))     # Sun angle at sunset  
  R_MJ<- (24*(60)/pi)*Gsc*dr*((omegaS)*sin(phi)*sin(delta) + cos(phi)*cos(delta)*sin(omegaS))
  if(unit=="mm")
  {
    lambda  <- 2.501 - 2.361E-3 * T  # latent heat of vaporization
    R_H2O <- R_MJ / lambda
    return(R_H2O)
  } else if(unit== "MJ") return(R_MJ) else print("Incorrect unit!", quote=FALSE)
  
}