SphericalStandardError3D <- function (coord) 
{
    x <- coord[, 1]
    y <- coord[, 2]
    z <- coord[, 3]
    if (length(x) >= 25) {
        n_elements = length(x)
        R = sqrt((sum(x) * sum(x)) + (sum(y) * sum(y)) + (sum(z) * 
            sum(z)))
        meanX <- sum(x)/R
        meanY <- sum(y)/R
        meanZ <- sum(z)/R
        x <- x * meanX
        y <- y * meanY
        z <- z * meanZ
        suma <- x + y + z
        suma2 <- suma * suma
        d <- 1 - (1/n_elements) * sum(suma2)
        Mm <- MeanModule3D(coord)
        sigma <- sqrt(Mod(d/(n_elements * Mm * Mm)))
        ea <- -log(0.05)
        Q <- asin(sigma * sqrt(ea))
        Q <- ToSexagesimal3D(Q)
        return(Q)
    }
    else {
        print("Number of Data < 25")
        return("---")
    }
}
