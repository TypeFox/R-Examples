`ordiequilibriumcircle` <-
function(pca,ordiplot,...) {
    `drawcircle` <- 
    function(x0=0, y0=0, radius=1, npoints=100,...) {
        a <- seq(0, 2*pi, len=npoints)    
        c <- array(dim=c(2,npoints))    
        c[1,] <- x0+cos(a)*radius
        c[2,] <- y0+sin(a)*radius
        for (i in 1:(npoints-1)) {graphics::segments(c[1, i], c[2, i], c[1, 1+i], c[2, 1+i],...)}
        graphics::segments(c[1, i], c[2, i], c[1, npoints], c[2, npoints],...)
    }
    eigen <- pca$CA$eig
    p <- length(eigen)   
    n <- nrow(pca$CA$u)
    tot <- sum(eigen)
    const <- ((n-1)*tot)^0.25
    radius <- (2/p)^0.5
    radius <- radius * const
    result <- list(radius=radius, constant=const)
    drawcircle(radius=radius,...)
    speciescoord <- scores(ordiplot, display="species")
    for (i in 1:nrow(speciescoord)) {
        length <- (speciescoord[i,1]^2+speciescoord[i,2]^2)^0.5
        if (length > radius) {graphics::arrows(0, 0, speciescoord[i,1], speciescoord[i,2],...)}
    }
    return(result)
}

