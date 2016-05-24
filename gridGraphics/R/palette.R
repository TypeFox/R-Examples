
# NOTE that there is only one palette for all graphics devices in R
# so only need to call these functions to replicate their effect
# (does not matter which graphics device is current)

C_palette <- function(x) {
    do.call("palette", x[-1])    
}

C_palette2 <- function(x) {
    hex <- sprintf("%08X", x[[2]])
    alpha <- substring(hex, 1, 2)
    blue <- substring(hex, 3, 4)
    green <- substring(hex, 5, 6)
    red <- substring(hex, 7, 8)
    palette(paste0("#", red, green, blue, alpha))
}

