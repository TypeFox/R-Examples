##
##  c a r t 2 s p h . R  Coordinate Transformations
##


cart2sph <- function(xyz) {
    stopifnot(is.numeric(xyz))

    # Transform cartesian to spherical coordinates
    if (is.vector(xyz) && length(xyz) == 3) {
        x <- xyz[1];  y <- xyz[2];  z <- xyz[3]
        m <- 1
    } else if (is.matrix(xyz) && ncol(xyz) == 3) {
        x <- xyz[, 1];  y <- xyz[, 2];  z <- xyz[, 3]
        m <- nrow(xyz)
    } else
        stop("Input must be a vector of length 3 or a matrix with 3 columns.")
    
    hypotxy <- hypot(x, y)
    r       <- hypot(hypotxy, z)
    phi     <- atan2(z, hypotxy)
    theta   <- atan2(y, x)

    if (m == 1) tpr <- c(theta, phi, r)
    else        tpr <- cbind(theta, phi, r)
    return(tpr)
}


sph2cart <- function(tpr) {
    stopifnot(is.numeric(tpr))

    # Transform spherical to cartesian coordinates
    if (is.vector(tpr) && length(tpr) == 3) {
        theta <- tpr[1];  phi <- tpr[2];  r <- tpr[3]
        m <- 1
    } else if (is.matrix(tpr) && ncol(tpr) == 3) {
        theta <- tpr[, 1];  phi <- tpr[, 2];  r <- tpr[, 3]
        m <- nrow(tpr)
    } else
        stop("Input must be a vector of length 3 or a matrix with 3 columns.")

    z   <- r * sin(phi)
    tmp <- r * cos(phi)
    x   <- tmp * cos(theta)
    y   <- tmp * sin(theta)

    if (m == 1) xyz <- c(x, y, z)
    else        xyz <- cbind(x, y, z)
    return(xyz)
}


cart2pol <- function(xyz) {
    stopifnot(is.numeric(xyz))

    # Transform cartesian to cylindrical or polar coordinates
    if (is.vector(xyz) && (length(xyz) == 2 || length(xyz) == 3)) {
        x <- xyz[1];  y <- xyz[2]
        m <- 1; n <- length(xyz)
    } else if (is.matrix(xyz) && (ncol(xyz) == 2 || ncol(xyz) == 3)) {
        x <- xyz[, 1];  y <- xyz[, 2]
        m <- nrow(xyz); n <- ncol(xyz)
    } else
        stop("Input must be a vector of length 3 or a matrix with 3 columns.")

    phi <- atan2(y, x)
    r   <- hypot(x, y)

    if (n == 2) {
        if (m == 1) prz <- c(phi, r)
        else        prz <- cbind(phi, r)
    } else {
        if (m == 1) {
            z <- xyz[3]
            prz <- c(phi, r, z)
        } else {
            z <- xyz[, 3]
            prz <- cbind(phi, r, z)
        }
    }
    return(prz)
}


pol2cart <- function(prz) {
    stopifnot(is.numeric(prz))

    # Transform polar or cylindrical to cartesian coordinates
    if (is.vector(prz) && (length(prz) == 2 || length(prz) == 3)) {
        phi <- prz[1];  r <- prz[2]
        m <- 1; n <- length(prz)
    } else if (is.matrix(prz) && (ncol(prz) == 2 || ncol(prz) == 3)) {
        phi <- prz[, 1];  r <- prz[, 2]
        m <- nrow(prz); n <- ncol(prz)
    } else
        stop("Input must be a vector of length 3 or a matrix with 3 columns.")

    x <- r * cos(phi)
    y <- r * sin(phi)

    if (n == 2) {
        if (m == 1) xyz <- c(x, y)
        else        xyz <- cbind(x, y)
    } else {
        if (m == 1) {
            z <- prz[3]
            xyz <- c(x, y, z)
        } else {
            z <- prz[, 3]
            xyz <- cbind(x, y, z)
        }
    }
    return(xyz)
}
