NightDay <- function(time, timezone){

   is.wholenumber <- function(timezone, tol = .Machine$double.eps^0.5){
      abs(timezone - round(timezone)) < tol
   }


   if(!is.wholenumber(timezone)) {
     stop("The variable timezone has to be an integer, e.g. a number between -11 and +11.")
   }
   if(timezone < -11 || timezone > 11){
     stop("The variable timezone needs to be a number between -11 and +11.")
   }



   time <- as.POSIXlt(time)
   Y <- time$year + 1900
   M <- time$mon + 1
   D <- time$mday
   H <- (time$hour - timezone) + time$min/60  + time$sec/3600

   K <- pi/180.0

   ## Deklination berechnen
   Dec <- function(D, M, Y, H) {

      N <- 365 * Y + D + 31 * M - 46
      if (M < 3){
        N <- N + as.integer((Y - 1) / 4)
      }else{
        N <- N - as.integer(0.4 * M + 2.3) + as.integer(Y / 4)
      }

      X <- (N - 693960) / 1461.0
      X <- (X - as.integer(X)) * 1440.02509 + as.integer(X) * 0.0307572
      X <- X + H/24 * 0.9856645 + 356.6498973
      X <- X + 1.91233 * sin(0.9999825 * X * K)
      X <- (X + sin(1.999965 * X * K) / 50 + 282.55462)/360

      X <- (X - as.integer(X)) * 360

      Y2000 <- (Y - 2000)/100
      Ekliptik <- 23.43929111 - (46.8150 + (0.00059 - 0.001813 * Y2000) * Y2000) * Y2000 / 3600

      X <- sin(X * K) *  sin(K * Ekliptik)

      return(atan(X / sqrt(1 - X^2)) / K + 0.00075)
    }

   ## Greenwhich hour Angle berechnen
   GHA <- function(D, M, Y, H) {

      N <- 365 * Y + D + 31 * M - 46
      if (M < 3){
        N <- N + as.integer((Y - 1) / 4)
      }else{
        N <- N - as.integer(0.4 * M + 2.3) + as.integer(Y / 4)
      }

      P <- H/24
      X <- (P + N - 7.22449E5) * 0.98564734 + 279.306
      X <- X * K
      XX <- -104.55 * sin(X) - 429.266 * cos(X) + 595.63 * sin(2.0 * X) - 2.283 * cos(2.0 * X)
      XX <- XX + 4.6 * sin(3.0 * X) + 18.7333 * cos(3.0 * X)
      XX <- XX - 13.2 * sin(4.0 * X) - cos(5.0 * X) - sin(5.0 * X) / 3.0 + 0.5 * sin(6.0 * X) + 0.231
      XX <- XX / 240 + 360 * (P + 0.5)
      if (XX > 360){
        XX = XX - 360
      }
      return(XX)
    }

   ## Breitengrade berechnen
   Lat <- function(longitude, dec = Dec(D, M, Y, H)){

     tang <- -cos(longitude*K)/tan(dec*K)
     itan <- atan(tang)
     itan = itan/K

     return(round(itan))
   }

   x0 <- 360

   x <- round(GHA(D, M, Y, H))

   if (x > 0){
     von <- x*(-1)
   }else{
     von <- x
   }

   yy <- vector()

   for ( i in von:x0-1 ){
     yy[i] <- Lat(longitude = i, Dec(D, M, Y, H))
   }

   l <- list(time, timezone, yy, Dec(D, M, Y, H), GHA(D, M, Y, H))
   names(l) <- c('Time', 'tz', 'Latitude', 'Declination', 'GHA')

   class(l) <- "NightDay"
   return(l)
}

