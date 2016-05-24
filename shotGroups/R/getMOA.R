## input range checking for angular diameter <-> absolute size conversion
checkSAngDst <-
function(s, ang, dst, type=c("deg", "rad", "MOA", "SMOA", "mrad", "mil")) {
    type <- match.arg(toupper(type),
                      choices=c("DEG", "RAD", "MOA", "SMOA", "MRAD", "MIL"))
    if(!missing(s)) {
        stopifnot(is.numeric(s))
        is.na(s) <- s < 0
    } else { s <- NULL }

    if(!missing(ang)) {
        ## pi/21600
        cnstPi21600    <- 0.00014544410433286079807697423070738439278690599071181
        ## 1 / atan(1/7200)
        cnst1atan17200 <- 7200.0000462962960581466263965930916898334026996030
        stopifnot(is.numeric(ang))
        is.na(ang) <- ang < 0
        is.na(ang) <- switch(type,
             DEG=(ang >= 180),
             MOA=(ang >= (60*180)),
            SMOA=(ang >= (cnstPi21600*cnst1atan17200*60*180)),
             RAD=(ang >= pi),
            MRAD=(ang >= (1000*pi)),
             MIL=(ang >= 3200))
    } else { ang <- NULL }

    if(!missing(dst)) {
        stopifnot(is.numeric(dst))
        is.na(dst) <- dst <= 0
    } else { dst <- NULL }

    list(s=s, ang=ang, dst=dst, type=type)
}

## convert size to angular measure
## angle in degree
## angle in radian
## MOA (= arcmin)
## Shooter's MOA (SMOA = IPHY, inches per hundred yards)
## milliradian mrad (1/1000 of a rad)
## mil -> NATO convention
getMOA <-
function(x, dst, conversion="m2cm", type=c("deg", "rad", "MOA", "SMOA", "mrad", "mil")) {
    arg  <- checkSAngDst(s=x, dst=dst, type=type)
    x    <- arg$s
    dst  <- arg$dst
    type <- arg$type

    ## convert distance measure to the unit of the (x,y)-coordinates
    dstCommon <- getConvFac(conversion) * dst

    ## MOA
    ## calculate angle in rad, convert to deg, convert to MOA
    ## one degree = 1/360 of a circle's arc
    ## one degree = 60 MOA -> 1 MOA = 1/60 degree
    ## rad <- 2*atan(x/(2*dstCommon))  # size as angle in rad
    ## deg <- (360/(2*pi)) * rad       # size as angle in degree
    ## MOA <- 60*deg                   # degree in MOA

    ## mrad
    ## calculate angle in rad, convert to milliradian
    ## 1 rad = 1000 mils -> 1 mrad = 1/1000 rad
    ## rad  <- atan(x/dstCommon)      # size as angle in rad
    ## mrad <- 1000*rad

    ## constants
    cnst360pi      <-  114.59155902616464175359630962821034066481094493313
    cnst21600pi    <- 6875.4935415698785052157785776926204398886566959877
    cnst1atan17200 <- 7200.0000462962960581466263965930916898334026996030
    cnstmil        <- 2037.1832715762602978417121711681838340410834654778

#     angle <- switch(type,
#          MOA =      (21600/pi)*atan(x/(2*dstCommon)),  # size in arcmin
#          SMOA=(1/atan(1/7200))*atan(x/(2*dstCommon)),  # size in SMOA
#          mrad=            2000*atan(x/(2*dstCommon)),  # size in milliradian
#           mil=       (6400/pi)*atan(x/(2*dstCommon)))  # size in NATO mil
    atanArg <- x/(2*dstCommon)
    angle   <- switch(type,
            DEG=     cnst360pi*atan(atanArg), # size in degree
            RAD=             2*atan(atanArg), # size in radian
            MOA=   cnst21600pi*atan(atanArg), # size in arcmin
           SMOA=cnst1atan17200*atan(atanArg), # size in SMOA
           MRAD=          2000*atan(atanArg), # size in milliradian
            MIL=       cnstmil*atan(atanArg)) # size in NATO mil

    return(angle)
}

## convert angular measure to size
## MOA (= arcmin)
## Shooter's MOA (SMOA = IPHY, inches per hundred yards)
## milliradian mrad (1/1000 of a rad)
fromMOA <-
function(x, dst, conversion="m2cm", type=c("deg", "rad", "MOA", "SMOA", "mrad", "mil")) {
    arg  <- checkSAngDst(ang=x, dst=dst, type=type)
    x    <- arg$ang
    dst  <- arg$dst
    type <- arg$type

    ## convert distance measure to the unit of the (x,y)-coordinates
    dstCommon <- getConvFac(conversion) * dst

    ## MOA
    ## convert MOA to degree, convert to rad, calculate size
    ## one degree = 1/360 of a circle's arc
    ## one degree = 60 MOA -> 1 MOA = 1/60 degree
    ## deg <- MOA / 60                   # MOA as angle in degree
    ## rad <- (2*pi/360) * deg           # deg as angle in rad
    ## x   <- 2*dstCommon * tan(rad/2)   # angle in rad as size

    ## SMOA
    ## get SMOA for 1 MOA = conversion factor for MOA2SMOA
    ## MOA2SMOA <- (21600/pi)*atan(1/7200)
    MOA2SMOA    <- 0.95492965241113508367872462693595042591621952511092
    cnstPi360   <- 0.0087266462599716478846184538424430635672143594427086
    cnstPi21600 <- 0.00014544410433286079807697423070738439278690599071181
    cnstmilI    <- 0.00049087385212340519350978802863742232565580771865236
    
    ## mrad
    ## convert milliradian to rad, calculate size
    ## rad <- (1/1000)*mrad              # milliradian as rad
    ## x   <- 2*dstCommon * tan(rad/2)   # angle in rad as size

#     size <- switch(type,
#            MOA =           2*dstCommon*tan(x*pi/21600), # size from arcmin
#            SMOA=MOA2SMOA * 2*dstCommon*tan(x*pi/21600), # size from SMOA
#            mrad=           2*dstCommon*tan(0.0005*x),   # size from milliradian
#             mil=           2*dstCommon*tan(x*pi/6400)) # size from NATO mil
    size <- switch(type,
             DEG=         2*dstCommon*tan(x*cnstPi360),    # size from degree
             RAD=         2*dstCommon*tan(x/2),            # size from radian
             MOA=         2*dstCommon*tan(x*cnstPi21600),  # size from arcmin
            SMOA=MOA2SMOA*2*dstCommon*tan(x*cnstPi21600),  # size from SMOA
            MRAD=         2*dstCommon*tan(x/2000),         # size from milliradian
             MIL=         2*dstCommon*tan(x*cnstmilI))     # size from NATO mil

    return(size)
}

#####---------------------------------------------------------------------------
## get distance from absolute and angular size
getDistance <-
function(x, angular, conversion="m2cm",
         type=c("deg", "rad", "MOA", "SMOA", "mrad", "mil")) {
    arg     <- checkSAngDst(s=x, ang=angular, type=type)
    x       <- arg$s
    angular <- arg$ang
    dst     <- arg$dst
    type    <- arg$type

    ## SMOA2MOA <- (21600/pi)*atan(1/7200)
    SMOA2MOA      <- 0.95492965241113508367872462693595042591621952511092
    ## pi/360
    cnstPi360     <- 0.0087266462599716478846184538424430635672143594427086
    ## pi/(2*60*180)
    cnstPi21600   <- 0.00014544410433286079807697423070738439278690599071181
    ## atan(1/7200)
    cnstAtan7200  <- 0.00013888888799582762807755515064054772884717306648502
    cnstmilI      <- 0.00049087385212340519350978802863742232565580771865236

    ## distance in unit of (x,y)-coordinates
    dstCommon <- switch(type,
           DEG=(x/2)*(1/tan(angular*cnstPi360)),
           RAD=(x/2)*(1/tan(angular/2)),
           MOA=(x/2)*(1/tan(angular*cnstPi21600)),
          SMOA=(x/2)*(1/tan(angular*cnstAtan7200)),
          MRAD=(x/2)*(1/tan(angular/2000)),
           MIL=(x/2)*(1/tan(angular*cnstmilI)))
    
    ## convert distance measure from the unit of the (x,y)-coordinates
    dst <- (1/getConvFac(conversion)) * dstCommon
    return(dst)
}

#####---------------------------------------------------------------------------
## get all angular measures for sizes x
makeMOA <-
function(x, dst, conversion) {
    if(length(x) < 2L) {
        y <- c(x,
               getMOA(x, dst=dst, conversion=conversion, type="MOA"),
               getMOA(x, dst=dst, conversion=conversion, type="SMOA"),
               getMOA(x, dst=dst, conversion=conversion, type="mrad"))
        setNames(y, c("unit", "MOA", "SMOA", "mrad"))
    } else {
        rbind(unit=x,
               MOA=getMOA(x, dst=dst, conversion=conversion, type="MOA"),
              SMOA=getMOA(x, dst=dst, conversion=conversion, type="SMOA"),
              mrad=getMOA(x, dst=dst, conversion=conversion, type="mrad"))
    }
}
