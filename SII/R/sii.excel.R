sii.excel <- function(
                      THDI,
                      freq=c(250, 500, 1000, 2000, 3000, 4000, 6000, 8000),
                      matlab.spline=TRUE
                )
{
  if(!require("splines"))
    stop("'splines' package must be installed to run this function.")
  
  if(length(THDI) != length(freq))
    stop("THDI must have the same length as freq")

  # Interpolate unobserved frequencies 
  sii.freqs <- sii.constants[,1]
  
  if(matlab.spline)
    {
      ispl <- splines::interpSpline( THDI ~ freq )
      THDI <- predict(ispl, sii.freqs)$y
      THDI[sii.freqs<freq[1]] <- THDI[1]
    }
  else
    {
      approx.l <- function(x,y,xout,...)
        { 
          retval <- approx(log(x), y, log(xout), ...)
          retval$x <- xout
          retval
        }

      THDI <- approx.l(x=freq, y=THDI, method="linear",  
                       xout=sii.freqs, rule=2)$y
    }

  # Set up working table
  sii.tab <- sii.constants
  
  # Remove example data
  sii.tab[,c("Ti'.THDN", "Xi'.Int.Noise", "Di.Equiv.Dist", "ki.Temp.Var", "Ai.BandAud","Si.SII.band.")] <- NA
  
  # Perform the calculations
  
  sii.tab[,"Ti'.THDN"]      <- THDI
  sii.tab[,"Xi'.Int.Noise"] <- sii.tab[,"Ti'.THDN"] + sii.tab[,"Xi.Ref.Int.Nz"]
  sii.tab[,"Di.Equiv.Dist"] <- sii.tab[,"Xi'.Int.Noise"]
  sii.tab[,"ki.Temp.Var"]   <- ifelse((sii.tab[,"Ei.Raised"] - sii.tab[,"Di.Equiv.Dist"] + 15)/30 >1, 1, 0)
  sii.tab[,"Ai.BandAud"]    <- sii.tab[,"ki.Temp.Var"]
  sii.tab[,"Si.SII.band."]  <- sii.tab[,"I.BandImp"] * sii.tab[,"Ai.BandAud"]

  sii.val <- sum(sii.tab[,"Si.SII.band."])

  retval <- list()

  retval$speech    <- THDI
  #retval$noise     <- noise
  #retval$threshold <- threshold
  #retval$loss      <- loss
  retval$freq      <- freq
  retval$sii       <- sii.val
  retval$table     <- sii.tab

  class(retval) <- "SII"
  retval
}
