MCMC <-
function (coordinates = NULL, geno.dip.codom = NULL, geno.dip.dom = NULL, 
    geno.hap = NULL, qtc = NULL, qtd = NULL, ql = NULL, path.mcmc, 
    rate.max, delta.coord = 0, shape1 = 2, shape2 = 20, npopmin = 1, 
    npopinit, npopmax, nb.nuclei.max, nit, thinning = 1, freq.model = "Uncorrelated", 
    varnpop = TRUE, spatial = TRUE, jcf = TRUE, filter.null.alleles = TRUE, 
    prop.update.cell = 0.1, write.rate.Poisson.process = FALSE, 
    write.number.nuclei = TRUE, write.number.pop = TRUE, write.coord.nuclei = TRUE, 
    write.color.nuclei = TRUE, write.freq = TRUE, write.ancestral.freq = TRUE, 
    write.drifts = TRUE, write.logposterior = TRUE, write.loglikelihood = TRUE, 
    write.true.coord = TRUE, write.size.pop = FALSE, write.mean.quanti = TRUE, 
    write.sd.quanti = TRUE, write.betaqtc = FALSE, miss.loc = NULL) 
{
    ploidy <- 2
    if (!is.null(geno.dip.dom) & !is.null(geno.hap)) {
        stop("It is currently not possible to analyze jointly diploid dominant genotypes and haploid genotypes")
    }
    if (is.null(geno.dip.dom) & is.null(geno.hap)) 
        ploidy <- 2
    if (!is.null(geno.dip.dom)) 
        ploidy <- 2
    if (!is.null(geno.hap)) 
        ploidy <- 1
    geno2 <- geno.dip.codom
    if (ploidy == 2) {
        geno1 <- geno.dip.dom
    }
    if (ploidy == 1) {
        geno1 <- geno.hap
    }
    if (substring(path.mcmc, first = nchar(path.mcmc), last = nchar(path.mcmc)) != 
        "/") {
        path.mcmc <- paste(path.mcmc, "/", sep = "")
    }
    short.path <- substring(path.mcmc, first = 1, last = nchar(path.mcmc) - 
        1)
    if (!file_test("-d", short.path)) 
        stop(paste("Directory ", path.mcmc, "does not exist."))
    if ((nit%%thinning) != 0) 
        stop("nit/thinning is not an integer")
    if (missing(npopmax)) 
        stop("Argument npopmax is missing with no default")
    if (missing(npopinit)) 
        npopinit <- npopmax
    if (npopinit > npopmax) 
        stop("npopinit > npopmax")
    if (npopinit < npopmin) 
        stop("npopinit < npopmin")
    if ((freq.model != "Correlated") & (freq.model != "Uncorrelated")) {
        stop(paste("Error:", freq.model, "is not a frequency model. Check spelling (case sensitive) "))
    }
    if ((ploidy != 1) & (ploidy != 2)) {
        stop(paste("ploidy = ", ploidy, " is not a valid value."))
    }
    if (!is.null(geno1) & is.null(geno2)) {
        if (filter.null.alleles) {
            if (ploidy == 1) {
                stop("Algorithm for filtering null alleles not compatible with haploid data")
            }
            if (ploidy == 2) {
                stop("Algorithm for filtering null alleles not compatible with dominant markers")
            }
        }
    }
    if (!is.null(geno1) & (ploidy == 2)) {
        geno1 <- geno1 + 1
    }
    if (is.null(geno1)) {
        nloc.geno1 <- 1
        nindiv.geno1 <- 0
    }
    else {
        nloc.geno1 <- ncol(geno1)
        nindiv.geno1 <- nrow(geno1)
        print(c("in MCMC.R nindiv.geno1=", nindiv.geno1))
    }
    if (is.null(geno2)) {
        nloc.geno2 <- 1
        nindiv.geno2 <- 0
        print(c("in MCMC.R nindiv.geno2=", nindiv.geno2))
    }
    else {
        nloc.geno2 <- ncol(geno2)/2
        nindiv.geno2 <- nrow(geno2)
    }
    if (is.null(qtc)) {
        nqtc <- 1
        nindivqtc <- 0
    }
    else {
        nqtc <- ncol(qtc)
        nindivqtc <- nrow(qtc)
    }
    if (is.null(qtd)) {
        nqtd <- 1
        nindivqtd <- 0
    }
    else {
        nqtd <- ncol(qtd)
        nindivqtd <- nrow(qtd)
    }
    if (is.null(ql)) {
        nql <- 1
        nindivql <- 0
    }
    else {
        nql <- ncol(ql)
        nindivql <- nrow(ql)
    }
    nnn <- c(nindiv.geno1, nindiv.geno2, nindivqtc, nindivqtd, 
        nindivql)
    sub <- nnn > 0
    if (length(unique(nnn[sub])) > 1) {
        print(paste("nindiv.geno1 = ", nindiv.geno1))
        print(paste("nindiv.geno2 = ", nindiv.geno2))
        print(paste("nindivqtc = ", nindivqtc))
        print(paste("nindivqtd = ", nindivqtd))
        print(paste("nindivql = ", nindivql))
        stop("Number of rows of data matrices do not match")
    }
    else {
        nindiv <- nnn[sub][1]
    }
    use.geno1 <- nindiv.geno1 > 0
    use.geno2 <- nindiv.geno2 > 0
    use.qtc <- nindivqtc > 0
    use.qtd <- nindivqtd > 0
    use.ql <- nindivql > 0
    print("defining dummy data arrays")
    if (nindiv.geno1 == 0) {
        geno1 <- matrix(nrow = nindiv, ncol = 1, data = NA)
    }
    if (nindiv.geno2 == 0) {
        geno2 <- matrix(nrow = nindiv, ncol = 2, data = NA)
    }
    if (nindivqtc == 0) {
        qtc <- matrix(nrow = nindiv, ncol = nqtc, data = NA)
    }
    if (nindivqtd == 0) {
        qtd <- matrix(nrow = nindiv, ncol = nqtd, data = NA)
    }
    if (nindivql == 0) {
        ql <- matrix(nrow = nindiv, ncol = nql, data = NA)
    }
    print("defining dummy coordinates if coord are missing")
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
        if (nrow(coordinates) != nindiv) {
            print(paste("number of individuals in data matrix =", 
                nindiv))
            print(paste("number of individuals in coordinate matrix =", 
                nrow(coordinates)))
            stop("Number of rows in coordinate matrix and data matrices do not match ")
        }
    }
    if (!is.matrix(coordinates)) 
        coordinates <- as.matrix(coordinates)
    print("defining dummy matrix indicating genuinely missing data in geno2")
    if (is.null(miss.loc)) {
        miss.loc <- matrix(nrow = nindiv, ncol = ncol(geno2)/2, 
            data = 0)
    }
    if (missing(rate.max)) 
        rate.max <- nindiv
    if (missing(nb.nuclei.max)) {
        nb.nuclei.max <- ifelse(spatial, 2 * nindiv, nindiv)
    }
    if (nb.nuclei.max < nindiv) {
        stop("nb.nuclei.max should be at least equal to the number of individuals")
    }
    if (spatial & (nb.nuclei.max < 2 * rate.max)) {
        stop("nb.nuclei.max is too small as compared to rate.max")
    }
    if (nindiv.geno1 > 0) {
        res <- FormatGenotypes(as.matrix(geno1), ploidy = 1)
        geno1 <- res$genotypes
        allele.numbers.geno1 <- res$allele.numbers
        print(c("In R function MCMC, allele.numbers.geno1=", 
            allele.numbers.geno1))
        if (sum(allele.numbers.geno1 == 1) > 0) {
            stop("Some of the markers do not display any polymorphism")
        }
        print(paste("Number of missing data in matrix geno1: ", 
            sum(is.na(geno1))))
    }
    else {
        allele.numbers.geno1 <- -999
    }
    if (nindiv.geno2 > 0) {
        res <- FormatGenotypes(as.matrix(geno2), ploidy = 2)
        geno2 <- res$genotypes
        allele.numbers.geno2 <- res$allele.numbers
        print(c("In R function MCMC, allele.numbers.geno2=", 
            allele.numbers.geno2))
        if (sum(allele.numbers.geno2 == 1) > 0) {
            stop("Some of the markers do not display any polymorphism")
        }
        print(paste("Number of missing data in matrix geno2: ", 
            sum(is.na(geno2))))
    }
    else {
        allele.numbers.geno2 <- -999
    }
    if (nindivql > 0) {
        res <- FormatGenotypes(as.matrix(ql), ploidy = 1)
        ql <- res$genotypes
        allele.numbers.ql <- res$allele.numbers
        print(paste("Number of missing obs. for qualitative variables: ", 
            sum(is.na(ql))))
    }
    else {
        allele.numbers.ql <- -999
    }
    if (nindivqtc > 0) {
        print(paste("Number of missing obs. for quantitative continuous variables: ", 
            sum(is.na(qtc))))
    }
    if (nindivqtd > 0) {
        print(paste("Number of missing obs. for quantitative discrete variables: ", 
            sum(is.na(qtd))))
    }
    nchar.path <- nchar(path.mcmc)
    if (filter.null.alleles) {
        allele.numbers.geno2 <- allele.numbers.geno2 + 1
    }
    nalt <- c(allele.numbers.geno2, allele.numbers.geno1, allele.numbers.ql)
    nalmax <- max(1, max(nalt))
    print("Defining working arrays for parameters of spatial model...")
    u <- matrix(nrow = 2, ncol = nb.nuclei.max, data = -999)
    utemp <- matrix(nrow = 2, ncol = nb.nuclei.max, data = -999)
    c <- rep(times = nb.nuclei.max, -999)
    ctemp <- rep(times = nb.nuclei.max, -999)
    t <- matrix(nrow = 2, ncol = nindiv, data = -999)
    ttemp <- matrix(nrow = 2, ncol = nindiv, data = -999)
    indcell <- rep(times = nindiv, -999)
    indcelltemp <- rep(times = nindiv, -999)
    distcell <- rep(times = nindiv, -999)
    distcelltemp <- rep(times = nindiv, -999)
    xlim <- ylim <- rep(-999, times = 2)
    print("Defining working arrays for parameters of genetic model...")
    ncolt <- nloc.geno2 + nloc.geno1 + nql
    f <- array(dim = c(npopmax, ncolt, nalmax), data = -999)
    ftemp <- array(dim = c(npopmax, ncolt, nalmax), data = -999)
    fa <- array(dim = c(ncolt, nalmax), data = -999)
    drift <- rep(-999, npopmax)
    drifttemp <- rep(-999, npopmax)
    n <- array(dim = c(npopmax, ncolt, nalmax), data = -999)
    ntemp <- array(dim = c(npopmax, ncolt, nalmax), data = -999)
    a <- rep(times = nalmax, -999)
    ptemp <- rep(times = nalmax, -999)
    cellclass <- rep(times = nb.nuclei.max, -999)
    listcell <- rep(times = nb.nuclei.max, -999)
    fmodel <- ifelse(freq.model == "Correlated", 1, 0)
    kfix <- 1 - as.integer(varnpop)
    full.cond.y <- matrix(nrow = nalmax, ncol = 2, 0)
    print("Defining working arrays for parameters of quantitative variables...")
    meanqtc <- sdqtc <- meanqtctmp <- sdqtctmp <- matrix(nrow = npopmax, 
        ncol = nqtc, data = -999)
    nnqtc <- sqtc <- ssqtc <- matrix(nrow = npopmax, ncol = nqtc, 
        data = -999)
    ksiqtc <- kappaqtc <- alphaqtc <- betaqtc <- gbeta <- hbeta <- rep(-999, 
        nqtc)
    if (use.qtc) {
        ksiqtc <- apply(qtc, 2, mean, na.rm = TRUE)
        kappaqtc <- hbeta <- 2/(apply(qtc, 2, max, na.rm = TRUE) - 
            apply(qtc, 2, min, na.rm = TRUE))^2
        alphaqtc <- rep(2, nqtc)
        gbeta <- rep(0.5, nqtc)
        betaqtc <- rgamma(n = nqtc, shape = gbeta, rate = hbeta)
    }
    print(paste("ksiqtc", ksiqtc))
    print(paste("kappaqtc", kappaqtc))
    print(paste("alphaqtc", alphaqtc))
    print(paste("hbeta", hbeta))
    print(paste("gbeta", gbeta))
    print(paste("betaqtc", betaqtc))
    geno2.999 <- geno2
    geno2.999[is.na(geno2)] <- -999
    geno1.999 <- geno1
    geno1.999[is.na(geno1.999)] <- -999
    qtc.999 <- qtc
    qtc.999[is.na(qtc.999)] <- -999
    qtd.999 <- qtd
    qtd.999[is.na(qtd.999)] <- -999
    ql.999 <- ql
    ql.999[is.na(ql.999)] <- -999
    true.geno <- cbind(geno2.999, matrix(nrow = nindiv, ncol = nloc.geno1 * 
        2, data = -999))
    sub <- seq(nloc.geno2 * 2 + 1, nloc.geno2 * 2 + nloc.geno1 * 
        2 - 1, 2)
    true.geno[, sub] <- geno1.999
    sub <- seq(nloc.geno2 * 2 + 2, nloc.geno2 * 2 + nloc.geno1 * 
        2, 2)
    true.geno[, sub] <- geno1.999
    integer.par <- c(write.rate.Poisson.process, write.number.nuclei, 
        write.number.pop, write.coord.nuclei, write.color.nuclei, 
        write.freq, write.ancestral.freq, write.drifts, write.logposterior, 
        write.loglikelihood, write.true.coord, write.size.pop, 
        write.mean.quanti, write.sd.quanti, write.betaqtc, fmodel, 
        kfix, spatial, jcf, filter.null.alleles, ploidy, nchar.path, 
        nit, thinning, use.geno1, use.geno2, use.qtc, use.qtd, 
        use.ql)
    integer.par <- c(integer.par, rep(-999, 100 - length(integer.par)))
    double.par <- c(rate.max, delta.coord, shape1, shape2, prop.update.cell)
    double.par <- c(double.par, rep(-999, 100 - length(double.par)))
    nitsaved <- nit/thinning
    out1 <- array(dim = c(nitsaved, 5), data = -999)
    outspace <- array(dim = c(nitsaved, 5, nb.nuclei.max), data = -999)
    outfreq <- array(dim = c(nitsaved, 3, npopmax, ncolt, nalmax), 
        data = -999)
    outqtc <- array(dim = c(nitsaved, 2, npopmax, nqtc), data = -999)
    out.res <- .Fortran("mcmcgld", PACKAGE = "Geneland", as.double(t(coordinates)), 
        as.integer(geno2.999), as.integer(miss.loc), as.integer(geno1.999), 
        as.integer(ql.999), as.integer(nql), as.double(qtc.999), 
        as.integer(nqtc), as.character(path.mcmc), as.integer(integer.par), 
        as.double(double.par), as.integer(nindiv), as.integer(nloc.geno2), 
        as.integer(nloc.geno1), as.integer(ncolt), as.integer(nalt), 
        as.integer(nalmax), as.integer(nb.nuclei.max), as.integer(npopinit), 
        as.integer(npopmin), as.integer(npopmax), as.double(xlim), 
        as.double(ylim), as.integer(indcell), as.integer(indcelltemp), 
        as.double(distcell), as.double(distcelltemp), as.double(t), 
        as.double(ttemp), as.double(u), as.double(utemp), as.integer(c), 
        as.integer(ctemp), as.double(f), as.double(ftemp), as.double(fa), 
        as.double(drift), as.double(drifttemp), as.integer(n), 
        as.integer(ntemp), as.double(a), as.double(ptemp), as.double(meanqtc), 
        as.double(sdqtc), as.double(meanqtctmp), as.double(sdqtctmp), 
        as.integer(nnqtc), as.double(sqtc), as.double(ssqtc), 
        as.double(ksiqtc), as.double(kappaqtc), as.double(alphaqtc), 
        as.double(betaqtc), as.double(gbeta), as.double(hbeta), 
        as.integer(cellclass), as.integer(listcell), as.integer(true.geno), 
        as.double(full.cond.y), as.integer(nitsaved), as.double(out1), 
        as.double(outspace), as.double(outfreq), as.double(outqtc))
    param <- c(paste("nindiv :", nindiv), paste("rate.max :", 
        rate.max), paste("nb.nuclei.max :", nb.nuclei.max), paste("nit :", 
        nit), paste("thinning :", thinning), paste("varnpop :", 
        varnpop), paste("npopmin :", npopmin), paste("npopinit :", 
        npopinit), paste("npopmax :", npopmax), paste("spatial :", 
        spatial), paste("delta.coord :", delta.coord), paste("use.geno1 :", 
        use.geno1), paste("use.geno2 :", use.geno2), paste("use.qtc :", 
        use.qtc), paste("use.qtd :", use.qtd), paste("use.ql :", 
        use.ql))
    if (use.geno1) {
        param <- c(param, paste("nloc.geno1 :", nloc.geno1), 
            paste("filter.null.alleles :", filter.null.alleles))
    }
    if (use.geno2) {
        param <- c(param, paste("nloc.geno2 :", nloc.geno2), 
            paste("filter.null.alleles :", filter.null.alleles))
    }
    if ((use.geno1 | use.geno2) | use.ql) {
        param <- c(param, paste("nalmax :", nalmax), paste("freq.model :", 
            freq.model))
    }
    if (use.geno1) {
        param <- c(param, paste("ploidy :", ploidy))
    }
    if (use.qtc) {
        param <- c(param, paste("nqtc :", nqtc))
    }
    if (use.qtd) {
        param <- c(param, paste("nqtd :", nqtd))
    }
    if (use.ql) {
        param <- c(param, paste("nql :", nql))
    }
    write.table(param, file = paste(path.mcmc, "parameters.txt", 
        sep = ""), quote = FALSE, row.names = FALSE, col.names = FALSE)
    write.table(allele.numbers.geno1, file = paste(path.mcmc, 
        "allele.numbers.geno1.txt", sep = ""), quote = FALSE, 
        row.names = FALSE, col.names = FALSE)
    write.table(allele.numbers.geno2, file = paste(path.mcmc, 
        "allele.numbers.geno2.txt", sep = ""), quote = FALSE, 
        row.names = FALSE, col.names = FALSE)
    write.table(allele.numbers.ql, file = paste(path.mcmc, "number.levels.ql.txt", 
        sep = ""), quote = FALSE, row.names = FALSE, col.names = FALSE)
    print("Writing MCMC outputs in external text files")
    out1 <- array(dim = c(nitsaved, 5), data = out.res[[61]])
    write.table(file = paste(path.mcmc, "Poisson.process.rate.txt", 
        sep = ""), out1[, 1], row.names = FALSE, col.names = FALSE)
    write.table(file = paste(path.mcmc, "nuclei.numbers.txt", 
        sep = ""), out1[, 2], row.names = FALSE, col.names = FALSE)
    write.table(file = paste(path.mcmc, "populations.numbers.txt", 
        sep = ""), out1[, 3], row.names = FALSE, col.names = FALSE)
    write.table(file = paste(path.mcmc, "log.likelihood.txt", 
        sep = ""), out1[, 4], row.names = FALSE, col.names = FALSE)
    write.table(file = paste(path.mcmc, "log.posterior.density.txt", 
        sep = ""), out1[, 5], row.names = FALSE, col.names = FALSE)
    outspace <- array(dim = c(nitsaved, 5, nb.nuclei.max), data = out.res[[62]])
    for (iitstor in 1:(nit/thinning)) {
        append <- ifelse(iitstor > 1, TRUE, FALSE)
        write.table(file = paste(path.mcmc, "coord.nuclei.txt", 
            sep = ""), t(outspace[iitstor, 1:2, ]), row.names = FALSE, 
            col.names = FALSE, append = append)
        write.table(file = paste(path.mcmc, "color.nuclei.txt", 
            sep = ""), outspace[iitstor, 3, ], row.names = FALSE, 
            col.names = FALSE, append = append)
        write.table(file = paste(path.mcmc, "hidden.coord.txt", 
            sep = ""), t(outspace[iitstor, 4:5, ]), row.names = FALSE, 
            col.names = FALSE)
    }
    outfreq <- array(dim = c(nitsaved, 3, npopmax, ncolt, nalmax), 
        data = out.res[[63]])
    for (iitstor in 1:(nit/thinning)) {
        append <- ifelse(iitstor > 1, TRUE, FALSE)
        write.table(file = paste(path.mcmc, "ancestral.frequencies.txt", 
            sep = ""), outfreq[iitstor, 2, 1, , ], row.names = FALSE, 
            col.names = FALSE, append = append)
        write.table(file = paste(path.mcmc, "drifts.txt", sep = ""), 
            t(outfreq[iitstor, 3, , 1, 1]), row.names = FALSE, 
            col.names = FALSE, append = append)
        for (iloc in 1:ncolt) {
            append <- ifelse(iitstor > 1 | iloc > 1, TRUE, FALSE)
            write.table(file = paste(path.mcmc, "frequencies.txt", 
                sep = ""), t(outfreq[iitstor, 1, , iloc, ]), 
                row.names = FALSE, col.names = FALSE, append = append)
        }
    }
    outqtc <- array(dim = c(nitsaved, 2, npopmax, nqtc), data = out.res[[64]])
    for (iitstor in 1:(nit/thinning)) {
        append <- ifelse(iitstor > 1, TRUE, FALSE)
        write.table(file = paste(path.mcmc, "mean.qtc.txt", sep = ""), 
            outqtc[iitstor, 1, , ], row.names = FALSE, col.names = FALSE, 
            append = append)
        write.table(file = paste(path.mcmc, "sd.qtc.txt", sep = ""), 
            outqtc[iitstor, 2, , ], row.names = FALSE, col.names = FALSE, 
            append = append)
    }
}
