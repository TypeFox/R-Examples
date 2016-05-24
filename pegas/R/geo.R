## geo.R (2014-11-17)

##   Tools for Geographic Data

## Copyright 2014 Emmanuel Paradis

## This file is part of the R-package `pegas'.
## See the file ../DESCRIPTION for licensing issues.

geoTrans <- function(x, degsym = NULL, minsym = "'", secsym = "\"")
{
    if (is.null(degsym)) degsym <- "\u00b0"
    x <- as.character(x)
    sow <- grep("[SOW]", x)
    x <- gsub("[SNOWE]", "", x)
    x <- strsplit(x, paste0("[", degsym, minsym, secsym, "]"))
    foo <- function(x) {
        x <- as.numeric(x)
        sum(x[1:3] / c(1, 60, 3600), na.rm = TRUE)
    }
    res <- sapply(x, foo)
    if (length(sow)) res[sow] <- -res[sow]
    res
}

geod <- function(lon, lat = NULL, R = 6371)
{
    deg2rad <- function(x) x * pi / 180
    if (is.null(lat)) {
        lat <- lon[, 2L]
        lon <- lon[, 1L]
    }
    n <- length(lat)
    lat <- deg2rad(lat)
    lon <- deg2rad(lon)
    absdiff <- function(x, y) abs(x - y)
    Delta_lon <- outer(lon, lon, absdiff)
    Delta_lat <- outer(lat, lat, absdiff)
    tmp <- cos(lat) # store the cosinus(lat) for all individuals
    A <- (sin(Delta_lat / 2))^2 + rep(tmp, n) * rep(tmp, each = n) * (sin(Delta_lon / 2))^2
    R * 2 * asin(sqrt(A))
### alternative:
### tmp <- sin(lat); tmp2 <- cos(lat)
### A <- acos(outer(tmp, tmp) + outer(tmp2, tmp2) * cos(Delta_lon))
### R * A
}

