Dat2Jul <- function( yr, mo, dy, hr=12 ) { #
  if ( missing(mo)) { # yr = c(yr, mo, dy, hr)
    hr <-  if(length(yr)==4) yr[4] else 12;  dy <- yr[3];  mo <- yr[2];  yr <- yr[1]  
  }
  if ( mo <= 2 )  { yr <- yr-1; mo <-  mo+12 } # Jan, Feb
  H <-  hr/24
    catn(" H =", H)
  if ( (yr*100 + mo )*100 + dy >= 15821015 ) {  # i.e. 1582-10-15
  # Gregorian calendar:
    A <- yr %/% 100  # centuries
    B <- 2 - A + A %/% 4 # correction for yr=k*400
  } else if (  (yr*100 + mo )*100 + dy <= 15821004 )  {
  # Julian calendar:
#    if (yr < 0) yr <- yr + 1  # no yr == 0
    B <- 0
  } else  {  # "no legal date, hiatus of calendars
    JD <- NA
    return(JD)
  }
  JD <- floor(365.25*(yr+4716))  # days of julian centuries since 4715 BC
  JD <- JD + floor(30.6001*(mo+1))  # days since March 1, except for Jan, Feb
  JD <- JD + dy + H + B
  JD <- JD - 1524.5
  return( JD )
} #  Dat2Jul

Jul2Dat <- function( JD ) { #
  stopifnot(is.numeric(JD)) 
  CC16 <- 2305447.5 # 1.1.1600 0h
  CC82 <- 2299160.5 # 1582   10    15    0
  BLp   <- c(366, 365, 365, 365) # for normal leap blocks, sum = 1461
  if ( JD < CC16 ) { #  1600
    if (JD >= CC82)  JD <-  JD + 10   #  make up for 10d gap
    yr   <- -4712
    rest <- JD + 0.5 # noon
    BL   <- BLp
  } else { # Gregorian calendar
    JD <- JD - CC16  # y1600
    yr   <- 1600
    rest <- JD
# block 400  ||||YnnnYnnn...|NnnnYnnnYnnn...|Nnnnnnn...|Nnnnynnn...
    BL <- c(36525, 36524, 36524, 36524)
    M <- dsm( rest, BL )
    yr <- yr +M[1]*400 + M[2]*100
    rest <- M[3]
# block 100   # |.nnnYnnnYnnn...
    BL <- BLp
    if ( M[2] == 0 ) {  # all quartets have c(366, 365, 365, 365)
      M <- dsm( rest, BL )
      yr <- yr + M[1]*4 + M[2]
      rest <- M[3]
    } else {  # |NnnnYnnnYnnnYnnnYnnn...
        if ( rest >= sum(BL) )  { # blocks with c(366, 365, 365, 365) 
      } else {         # first block with c(365, 365, 365, 365)
          M <- dsm( rest, BL )
        yr <-  yr + M[1]
        rest <- M[3]
      } # END blocks with c(366
    }  # all quartets have c(366
# end block 100
  } # Gregorian calendar
# block 4
  if ( yr < 1600 || yr %% 400 == 0 ) {
    BL <- BLp          # isLeap  Ynnn
  } else {
    BL <- rep(365,4)   # NO leap nnnn
  } # END leap
  M <- dsm( rest, BL )
  yr <- yr + M[1]*4 + M[2]
  rest   <- M[3]
  isLeap <- (yr %% 4 == 0) && ( (yr < 1600) || (yr %% 100 != 0) || ( yr %% 400 == 0 ) )
  md <- mdiny( rest, isLeap )
  mo <- md[1]
  if (mo > 12) {
    mo <- mo-12; yr <-  yr + 1
  }
  dy <- floor(md[2]) + 1 ### 201503021453
  hr <- (md[2] - dy)*24.0
  res <- c(yr, mo, dy, hr)    
  names(res) <- c("yr", "mo", "dy", "hr")
    return ( res )
} # Jul2Dat
monthsN <- function(leap)  {
    res <- c( 31, 28+leap, 31, diff(floor( ((1:10)*153+2)/5 )))
    return( res)
} # monthsN
Mnames <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
Dnames <- c("Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday", "Sunday") 
mdiny <- function( dk, leap ) submod( dk, cumsum( monthsN(leap ) )) + c(1,0) # 
Wday <- function( JD ) { d <- JD %% 7; if (d == 0) d <- 7; return( Dnames[d] ) }
Yday <- function(mo, dy, leap ) if (mo==1) 0 else cumsum( monthsN( leap ) )[mo-1] + dy
