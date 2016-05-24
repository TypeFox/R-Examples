PostProcessChain <-
function (coordinates = NULL, path.mcmc, nxdom, nydom, burnin) 
{
    print("Reading MCMC parameter file")
    param <- as.matrix(read.table(paste(path.mcmc, "parameters.txt", 
        sep = "")))
    delta.coord <- as.numeric(param[param[, 1] == "delta.coord", 
        3])
    npopmax <- as.numeric(param[param[, 1] == "npopmax", 3])
    nb.nuclei.max <- as.numeric(param[param[, 1] == "nb.nuclei.max", 
        3])
    nit <- as.numeric(param[param[, 1] == "nit", 3])
    thinning <- as.numeric(param[param[, 1] == "thinning", 3])
    nindiv <- as.numeric(param[param[, 1] == "nindiv", 3])
    use.geno2 <- as.logical(param[param[, 1] == "use.geno2", 
        3])
    use.geno1 <- as.logical(param[param[, 1] == "use.geno1", 
        3])
    use.qtc <- as.logical(param[param[, 1] == "use.qtc", 3])
    use.ql <- as.logical(param[param[, 1] == "use.ql", 3])
    spatial <- as.logical(param[param[, 1] == "spatial", 3])
    if (substring(path.mcmc, first = nchar(path.mcmc), last = nchar(path.mcmc)) != 
        "/") {
        path.mcmc <- paste(path.mcmc, "/", sep = "")
    }
    short.path <- substring(path.mcmc, first = 1, last = nchar(path.mcmc) - 
        1)
    if (!file_test("-d", short.path)) 
        stop(paste("Directory ", path.mcmc, "does not exist."))
    if (nit/thinning <= burnin) {
        print(paste("nit=", nit))
        print(paste("thinning=", thinning))
        print(paste("Number of saved iterations (nit/thinning) =", 
            nit/thinning))
        print(paste("burnin=", burnin))
        stop("*** burnin has to be smaller than the number of saved iterations ***")
    }
    if (use.geno2) {
        nloc.geno2 <- as.numeric(param[param[, 1] == "nloc.geno2", 
            3])
        allele.numbers.geno2 <- scan(paste(path.mcmc, "allele.numbers.geno2.txt", 
            sep = ""))
        dominance <- param[param[, 1] == "dominance", 3]
        filter.null.alleles <- as.logical(param[param[, 1] == 
            "filter.null.alleles", 3])
    }
    else {
        nloc.geno2 <- 1
        allele.numbers.geno2 <- rep(-999, nloc.geno2)
    }
    if (use.geno1) {
        nloc.geno1 <- as.numeric(param[param[, 1] == "nloc.geno1", 
            3])
        allele.numbers.geno1 <- scan(paste(path.mcmc, "allele.numbers.geno1.txt", 
            sep = ""))
    }
    else {
        nloc.geno1 <- 1
        allele.numbers.geno1 <- rep(-999, nloc.geno1)
    }
    if (use.ql) {
        nql <- as.numeric(param[param[, 1] == "nql", 3])
        allele.numbers.ql <- scan(paste(path.mcmc, "number.levels.ql.txt", 
            sep = ""))
    }
    else {
        nql <- 1
        allele.numbers.ql <- rep(-999, nql)
    }
    if (use.qtc) {
        nqtc <- as.numeric(param[param[, 1] == "nqtc", 3])
    }
    else {
        nqtc <- 1
        qtc <- matrix(nrow = nindiv, ncol = nqtc, data = -999)
    }
    nalt <- c(allele.numbers.geno2, allele.numbers.geno1, allele.numbers.ql)
    nalmax <- max(1, max(nalt))
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
    filenpop <- paste(path.mcmc, "populations.numbers.txt", sep = "")
    filenpp <- paste(path.mcmc, "nuclei.numbers.txt", sep = "")
    fileu <- paste(path.mcmc, "coord.nuclei.txt", sep = "")
    filec <- paste(path.mcmc, "color.nuclei.txt", sep = "")
    filef <- paste(path.mcmc, "frequencies.txt", sep = "")
    fileperm <- paste(path.mcmc, "perm.txt", sep = "")
    filedom <- paste(path.mcmc, "proba.pop.membership.txt", sep = "")
    filedomperm <- paste(path.mcmc, "proba.pop.membership.perm.txt", 
        sep = "")
    filemeanqtc <- paste(path.mcmc, "mean.qtc.txt", sep = "")
    filemeanf <- paste(path.mcmc, "mean.freq.txt", sep = "")
    print("Estimating number of populations")
    npop <- scan(paste(path.mcmc, "populations.numbers.txt", 
        sep = ""))
    npop.est <- order(hist(npop[-(1:burnin)], breaks = seq(0.5, 
        npopmax + 0.5, 1), plot = FALSE)$counts, decreasing = TRUE)[1]
    lpd <- scan(paste(path.mcmc, "log.posterior.density.txt", 
        sep = ""))
    subK <- (npop == npop.est) & ((1:(nit/thinning)) > burnin)
    pivot <- order(lpd[subK], decreasing = TRUE)[1]
    print(paste("Iteration with highest posterior density:", 
        pivot))
    ncolt <- nloc.geno2 + nloc.geno1 + nql
    dom <- matrix(nrow = nxdom * nydom, ncol = npopmax, data = 0)
    domperm <- matrix(nrow = nxdom * nydom, ncol = npopmax, data = 0)
    coorddom <- matrix(nrow = 2, ncol = nxdom * nydom, data = -999)
    indvois <- numeric(nxdom * nydom)
    distvois <- numeric(nxdom * nydom)
    orderf <- orderftmp <- 1:npopmax
    u <- matrix(nrow = 2, ncol = nb.nuclei.max, data = -999)
    c <- rep(times = nb.nuclei.max, -999)
    xlim <- ylim <- rep(-999, 2)
    f <- fpiv <- fmean <- array(dim = c(npopmax, ncolt, nalmax), 
        -999)
    fmean[1:npop.est, , ] <- 0
    meanqv <- meanqvpiv <- matrix(nrow = npopmax, ncol = nqtc, 
        -999)
    ninrub = 0
    nitsaved <- nit/thinning
    out.orderf <- matrix(nrow = nit/thinning, ncol = npopmax, 
        data = -999)
    print("Calling Fortran function postprocesschain2")
    out.res <- .Fortran("postprocesschain2", PACKAGE = "Geneland", 
        as.integer(nxdom), as.integer(nydom), as.integer(burnin), 
        as.integer(ninrub), as.integer(npopmax), as.integer(nb.nuclei.max), 
        as.integer(nindiv), as.integer(nloc.geno2), as.integer(nloc.geno1), 
        as.integer(nql), as.integer(ncolt), as.integer(nalt), 
        as.integer(nalmax), as.double(xlim), as.double(ylim), 
        as.double(delta.coord), as.integer(nit), as.integer(thinning), 
        as.character(filenpop), as.character(filenpp), as.character(fileu), 
        as.character(filec), as.character(filef), as.character(fileperm), 
        as.character(filedom), as.character(filemeanqtc), as.character(filemeanf), 
        as.double(t(coordinates)), as.double(u), as.integer(c), 
        as.double(f), as.integer(pivot), as.double(fpiv), as.double(fmean), 
        as.double(dom), as.double(coorddom), as.integer(indvois), 
        as.double(distvois), as.integer(orderf), as.integer(orderftmp), 
        as.integer(npop.est), as.integer(use.geno2), as.integer(use.geno1), 
        as.integer(use.ql), as.integer(use.qtc), as.integer(nqtc), 
        as.double(meanqv), as.double(meanqvpiv), as.integer(nitsaved), 
        as.integer(out.orderf))
    print("End of Fortran function postprocesschain2")
    out.orderf <- matrix(nrow = nit/thinning, ncol = npopmax, 
        data = out.res[[50]])
    write.table(out.orderf, file = paste(path.mcmc, "perm.txt", 
        sep = ""), row.names = FALSE, col.names = FALSE, append = FALSE)
    fmean <- array(dim = c(npopmax, ncolt, nalmax), data = out.res[[34]])
    for (iloc in 1:ncolt) {
        append <- ifelse(iloc == 1, FALSE, TRUE)
        write.table(t(fmean[, iloc, ]), paste(path.mcmc, "mean.freq.txt", 
            sep = ""), row.names = FALSE, col.names = FALSE, 
            append = append)
    }
    dom <- matrix(nrow = nxdom * nydom, ncol = npopmax, data = out.res[[35]])
    coorddom <- matrix(nrow = 2, ncol = nxdom * nydom, data = out.res[[36]])
    write.table(cbind(t(coorddom), dom), paste(path.mcmc, "proba.pop.membership.txt", 
        sep = ""), append = FALSE)
    coord.grid <- read.table(filedom)[, 1:2]
    pmbr <- as.matrix(read.table(filedom)[, -(1:2)])
    pmp <- rep(NA, dim(pmbr)[1])
    for (ipix in 1:(dim(pmbr)[1])) {
        pmp[ipix] <- order(pmbr[ipix, ], decreasing = TRUE)[1]
    }
    write.table(cbind(coord.grid, pmp), file = paste(path.mcmc, 
        "modal.pop.txt", sep = ""), quote = FALSE, row.names = FALSE, 
        col.names = FALSE)
    indvois <- numeric(nindiv)
    distvois <- numeric(nindiv)
    u <- matrix(nrow = 2, ncol = nb.nuclei.max, data = -999)
    c <- rep(times = nb.nuclei.max, -999)
    pmp <- matrix(nrow = nindiv, ncol = npopmax, data = 0)
    out.res <- .Fortran("pppmindiv2", PACKAGE = "Geneland", as.integer(nindiv), 
        as.double(t(coordinates)), as.integer(npopmax), as.integer(nb.nuclei.max), 
        as.integer(indvois), as.double(distvois), as.double(u), 
        as.integer(c), as.double(pmp), as.character(filenpop), 
        as.character(filenpp), as.character(fileu), as.character(filec), 
        as.character(fileperm), as.integer(nit), as.integer(thinning), 
        as.integer(burnin), as.integer(orderf), as.integer(npop.est), 
        as.integer(pivot))
    pmp <- matrix(nrow = nindiv, ncol = npopmax, data = out.res[[9]])
    mod.pop.indiv <- t(apply(pmp, 1, order))[, npopmax]
    write.table(cbind(coordinates, pmp), file = paste(path.mcmc, 
        "proba.pop.membership.indiv.txt", sep = ""), quote = FALSE, 
        row.names = FALSE, col.names = FALSE)
    write.table(cbind(coordinates, mod.pop.indiv), file = paste(path.mcmc, 
        "modal.pop.indiv.txt", sep = ""), quote = FALSE, row.names = FALSE, 
        col.names = FALSE)
    param <- c(paste("nxdom :", nxdom), paste("nydom :", nydom), 
        paste("burnin :", burnin))
    write.table(param, file = paste(path.mcmc, "postprocess.parameters.txt", 
        sep = ""), quote = FALSE, row.names = FALSE, col.names = FALSE)
}
