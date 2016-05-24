##
##  d e g 2 r a d . R  Angle Conversion
##


deg2rad <- function(deg) {
    stopifnot(is.numeric(deg))
    ( rad <- (pi/180)*deg )
}


rad2deg <- function(rad) {
    stopifnot(is.numeric(rad))
    ( deg <- rad/(pi/180) )
}
