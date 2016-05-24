#' Calculate the grid cella area for a centered lat/lon grid
#'
#' Calculate the grid cella area for a centered lat/lon grid.
#'
#' @param lat latitude coordinates of the grid centers
#' @param lon longitude coordinates of the grid centers
#' @param verbose logical. Print info as we go?
#' @return The grid cell area in m^2 (meter*meter)
#' @details Currently the lon must be uniform but the lat does not need to be.
#' @keywords internal
#' @note This is an internal RCMIP5 function and not exported.
calcGridArea<- function(lon, lat, verbose=FALSE) {

    # Sanity checks - parameter classes and lengths
    stopifnot(is.numeric(lat))
    stopifnot(is.numeric(lon))
    stopifnot(length(verbose)==1 & is.logical(verbose))

    if(verbose) cat('Calculating grid cell areas...\n')
    numLat <- length(lat)
    numLon <- length(lon)

    # Check that the latitudes are centered in the grids
    if(any(abs(lat) == 90)) {
        stop('No grid cell can be centered at pole.')
    }

    # If for some reason we have a -180:180 lon base, reset to span 0:360
    if(any(lon < 0)) {
        lon <- lon + 180
    }

    # Calculate the longitude degrees spanned by a grid cell
    # ... modulo 360 to deal with wrapping boundries
    deltaLon <- (lon[c(2:length(lon),1)] - lon[1:length(lon)]) %% 360

    # Calculate the min/max latitude for each grid cell
    edgeLat <- (lat[2:length(lat)]+lat[2:length(lat)-1])/2
    minLat <- c(-90, edgeLat)
    maxLat <- c(edgeLat, 90)

    # Convert from degree to radius
    deltaLon <- matrix(deltaLon/180*pi,  nrow=length(lon), ncol=length(lat))
    minLat <- matrix(minLat/180*pi,  nrow=length(lon), ncol=length(lat), byrow=TRUE)
    maxLat <- matrix(maxLat/180*pi,  nrow=length(lon), ncol=length(lat), byrow=TRUE)
    lat <- matrix(lat/180*pi,   nrow=length(lon), ncol=length(lat), byrow=TRUE)

    # Assume the radius of the earth: 6371e3 meter
    R <- 6371e3 # meters

    # Calculate the east/west edges by assuming the earth is spherical and
    # ...east/west edges are defined by latitude arc lengths
    # ... => R*(maxLat-minLat)
    # Calculate the north/south edges by assuming the arc length of longitude
    # ...is the lattitude corrected radius (R*cos(lat)) times the change in lon
    # ... => (R*cos(lat))*deltaLon
    return( R*(maxLat-minLat) * (R*cos(lat))*deltaLon)

    #old formulation for reference (updated 29 September 2014)
    # ...no significant difference but harder to explain
    #R^2*(sin(maxLat)-sin(minLat))*deltaLon

} # calcGridArea
