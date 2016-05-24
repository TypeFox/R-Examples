##  IAU1976
##  http://www.kayelaby.npl.co.uk/general_physics/2_7/2_7_2.html

cJDJ2000      <- 2451545.0  ## Julian day number of the epoch J2000.0
cDAYPJULCENT  <- 36525.0             ##  days per julian century
cDAYPYEARTROP <- 365.242193   ## - 0.0000061 T,  days per tropical year  )
cDAYPYEARSID  <- 365.25636042    ##  days per sidereal year  )
cDAYPMONSYN   <- 29.530859      ##  days per synodical month )
cDAYPMONSID   <- 27.321661          ##  days per sidereal month    )

cK       <- 0.01720209895  ##  Gravitational constant, GAUSSian definition
cC       <- 299792458.0      ##  [m/s] defined speed of light  )
cRE      <- 6378140.0          ##  [m]    radius of earth.s equator )
cMY      <- 0.01230002        ##  ratio mass of moon/mass of earth  )
cPRECESS <- 5029.0966          ##  [arc sec] precession per year at 2000.0 )
cEPSOBL  <- 23.43929111       ##  [deg]  inclination of ecliptic at 2000.0   )
cAE      <- 1.49597870E11       		##  [m]  distance Earth to Sun
cSBYE    <- 332946.0            ##  ratio mass of Sun/mass of Earth    )
cSBYEMY  <- 328900.5             ##  ratio mass of Sun/mass of ## Earth+Moon) )
cSBYME   <- 6023600.0          ##  ratio mass of Sun/mass of Mercury    )
cSBYVE   <- 408523.5            ##  ratio mass of Sun/mass of Venus     )
cSBYMA   <- 3098710.0          ##  ratio mass of Sun/mass of Mars      )
cSBYJU   <- 1047.355            ##  ratio mass of Sun/mass of Jupiter   )
cSBYSA   <- 3498.5                ##  ratio mass of Sun/mass of Saturn    )
cSBYUR   <- 22869.0              ##  ratio mass of Sun/mass of Uranus    )
cSBYNE   <- 19314.0              ## * ratio mass of Sun/mass of Neptun    *)
cSBYPL   <- 130000000.0      ## ** ratio mass of Sun/mass of Pluto     *)
cSOLBYSID <- 1.00273790934  ## ** ratio solar/sidereal day       *)
cSIDBYSOL <- 0.99726956634  ## ** ratio sidereal/solar day        *)

DPY <-  cDAYPJULCENT/100.0; ## Tage/jul.Jahr      *)
DAYINMONTH <- c(31,28,31,30,31,30,31,31,30,31,30,31,31)
