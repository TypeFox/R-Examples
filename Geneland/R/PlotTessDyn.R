PlotTessDyn <-
function (coordinates = NULL, path.mcmc, nxgrid, nygrid) 
{
    param <- as.matrix(read.table(paste(path.mcmc, "parameters.txt", 
        sep = "")))
    delta.coord <- as.numeric(param[param[, 1] == "delta.coord", 
        3])
    npopmax <- as.numeric(param[param[, 1] == "npopmax", 3])
    nppmax <- as.numeric(param[param[, 1] == "nb.nuclei.max", 
        3])
    nit <- as.numeric(param[param[, 1] == "nit", 3])
    thinning <- as.numeric(param[param[, 1] == "thinning", 3])
    nindiv <- as.numeric(param[param[, 1] == "nindiv", 3])
    spatial <- as.logical(param[param[, 1] == "spatial", 3])
    npp <- scan(paste(path.mcmc, "nuclei.numbers.txt", sep = ""))
    c <- matrix(ncol = nit/thinning, byrow = FALSE, scan(paste(path.mcmc, 
        "color.nuclei.txt", sep = "")))
    u <- array(dim = c(2, nppmax, nit/thinning), data = scan(paste(path.mcmc, 
        "coord.nuclei.txt", sep = "")))
    indgrid <- rep(-999, nxgrid * nygrid)
    distgrid <- rep(-999, nxgrid * nygrid)
    coorddom <- matrix(nrow = 2, ncol = nxgrid * nygrid, -999)
    dom <- matrix(nrow = nxgrid * nygrid, ncol = npopmax, -999)
    ccur <- rep(-999, nppmax)
    ucur <- matrix(nrow = 2, ncol = nppmax, data = -999)
    if (is.null(coordinates)) {
        if (spatial) {
            stop("Please give spatial coordinates of individuals or set argument spatial to FALSE")
        }
        else {
            n.int <- ceiling(sqrt(nindiv))
            x <- rep(seq(from = 0, to = 1, length = n.int), n.int)
            y <- rep(seq(from = 0, to = 1, length = n.int), n.int)
            y <- as.vector(t(matrix(nrow = n.int, ncol = n.int, 
                y, byrow = FALSE)))
            coordinates <- cbind(x, y)[1:nindiv, ]
        }
    }
    else {
        if (ncol(coordinates) != 2) 
            stop("matrix of coordinates does not have 2 columns")
    }
    coordinates <- as.matrix(coordinates)
    for (where in seq(0/thinning, nit/thinning - 1, 1)) {
        burnin = where
        ninrub = nit/thinning - (where + 1)
        res.tess <- .Fortran(name = "tessdyn", PACKAGE = "Geneland", 
            as.integer(npopmax), as.integer(nppmax), as.integer(nxgrid), 
            as.integer(nygrid), as.integer(indgrid), as.integer(npp), 
            as.integer(c), as.integer(ccur), as.integer(nit/thinning), 
            as.integer(burnin), as.integer(ninrub), as.double(dom), 
            as.double(distgrid), as.double(coorddom), as.double(u), 
            as.double(ucur), as.integer(nindiv), as.double(t(coordinates)), 
            as.double(rep(0, 2)), as.double(rep(0, 2)), as.double(delta.coord))
        dom <- matrix(nrow = (nygrid * nxgrid), ncol = npopmax, 
            data = res.tess[[12]], byrow = FALSE)
        map.dom <- t(apply(dom, 1, order))[, npopmax]
        add <- ifelse(where > 1, TRUE, FALSE)
        image(seq(min(coordinates[, 1]), max(coordinates[, 1]), 
            length = nxgrid), seq(min(coordinates[, 2]), max(coordinates[, 
            2]), length = nygrid), matrix(map.dom, nrow = nxgrid, 
            ncol = nygrid, byrow = TRUE), xlab = "Easting", ylab = "Northing", 
            cex.main = 1.8, col = terrain.colors(npopmax), add = add)
        points(coordinates)
    }
}
