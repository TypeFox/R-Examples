HZ <-
function (coordinates, geno.dip.codom = NULL, geno.dip.dom = NULL, 
    geno.hap = NULL, dist.IC = NULL, allele.freq = NULL, ncluster = NULL, 
    cluster.indiv = NULL, path.mcmc.noadm = NULL, a.init = NULL, 
    b.init = NULL, c.init = 1, a.max = 10, b.max = NULL, c.max = 1, 
    estimate.a = TRUE, estimate.b = TRUE, estimate.c = FALSE, 
    common.param = TRUE, nit, thinning, path.mcmc.adm = NULL) 
{
    rdirichlet <- function(param, n.real = 1) {
        param = t(param)
        nb = length(param)
        Y <- matrix(nrow = nb, ncol = n.real, data = NA)
        X <- matrix(nrow = nb, ncol = n.real, data = NA)
        for (icol in 1:n.real) {
            Y[, icol] <- rgamma(n = nb, shape = param)
        }
        for (icol in 1:n.real) {
            X[, icol] <- Y[, icol]/sum(Y[, icol])
        }
        return(t(X))
    }
    EstimateFreq <- function(path.mcmc.noadm) {
        fileparam <- paste(path.mcmc.noadm, "parameters.txt", 
            sep = "")
        param <- as.matrix(read.table(fileparam))
        nit <- as.numeric(param[param[, 1] == "nit", 3])
        thinning <- as.numeric(param[param[, 1] == "thinning", 
            3])
        filter.null.alleles <- as.logical(param[param[, 1] == 
            "filter.null.alleles", 3])
        use.geno1 <- as.logical(param[param[, 1] == "use.geno1", 
            3])
        use.geno2 <- as.logical(param[param[, 1] == "use.geno2", 
            3])
        allele.numbers.dip <- scan(paste(path.mcmc.noadm, "allele.numbers.geno2.txt", 
            sep = ""))
        allele.numbers.hap <- scan(paste(path.mcmc.noadm, "allele.numbers.geno1.txt", 
            sep = ""))
        nb.levels.ql <- scan(paste(path.mcmc.noadm, "number.levels.ql.txt", 
            sep = ""))
        nall <- c(allele.numbers.dip, allele.numbers.hap, nb.levels.ql)
        filemeanf <- paste(path.mcmc.noadm, "mean.freq.txt", 
            sep = "")
        fmean <- as.matrix(read.table(filemeanf))
        npopmax <- ncol(fmean)
        nlocd <- length(allele.numbers.dip)
        nloch <- length(allele.numbers.hap)
        nql <- length(nb.levels.ql)
        ncolt <- nlocd + nloch + nql
        nallmax <- max(nall[nall > 0])
        res <- array(dim = c(npopmax, ncolt, nallmax))
        irow <- 1
        for (iloc in 1:ncolt) {
            for (iall in 1:nallmax) {
                res[, iloc, iall] <- fmean[irow, ]
                irow <- irow + 1
            }
        }
        if (use.geno2) {
            res <- res[, 1:nlocd, ]
        }
        if (use.geno1) {
            res <- res[, nlocd + (1:nloch), ]
        }
        list(Estimated.freq = res)
    }
    if (!is.null(geno.dip.codom) | !is.null(geno.dip.dom)) {
        ploidy <- 2
    }
    else {
        ploidy <- 1
    }
    nindiv <- nrow(coordinates)
    if (!is.null(geno.dip.codom)) {
        data.fmt <- FormatGenotypes(geno.dip.codom, ploidy = ploidy)
        if (is.null(allele.freq)) {
            geno.dip.codom.fmt <- data.fmt$genotypes
        }
        else {
            geno.dip.codom.fmt <- geno.dip.codom
        }
        nalmax.dip.codom <- max(data.fmt$allele.numbers)
        nlocd <- ncol(geno.dip.codom)/2
        use.codom <- TRUE
    }
    else {
        nalmax.dip.codom <- 1
        use.codom <- FALSE
        nalmax.dip.codom <- -999
        geno.dip.codom.fmt <- matrix(nrow = nindiv, ncol = 2, 
            -999)
    }
    if (!is.null(geno.dip.dom)) {
        nlocd <- ncol(geno.dip.dom)
        use.dom <- TRUE
        nalmax.dip.dom <- 2
        geno.dip.dom.fmt <- geno.dip.dom + 1
    }
    else {
        use.dom <- FALSE
        nalmax.dip.dom <- -999
        geno.dip.dom.fmt <- matrix(nrow = nindiv, ncol = 1, -999)
    }
    if (!is.null(geno.hap)) {
        nlocd <- 1
    }
    if (!is.null(geno.hap)) {
        data.fmt <- FormatGenotypes(geno.hap, ploidy = ploidy)
        if (is.null(allele.freq)) {
            geno.hap.fmt <- data.fmt$genotypes
        }
        else {
            geno.hap.fmt <- geno.hap
        }
        nalmax.hap <- max(data.fmt$allele.numbers)
        nloch <- ncol(geno.hap)
        use.hap <- TRUE
    }
    else {
        nloch <- 1
        geno.hap.fmt <- matrix(nrow = nindiv, ncol = nloch, -999)
        nalmax.hap <- 1
        use.hap <- FALSE
        nalmax.hap <- -999
    }
    nalmax <- max(nalmax.hap, nalmax.dip.codom, nalmax.dip.dom)
    if (!is.null(path.mcmc.noadm)) {
        param <- as.matrix(read.table(paste(path.mcmc.noadm, 
            "parameters.txt", sep = "")))
        freq.est <- EstimateFreq(path.mcmc.noadm = path.mcmc.noadm)$Estimated.freq
        cluster.indiv <- read.table(paste(path.mcmc.noadm, "modal.pop.indiv.txt", 
            sep = ""))[, 3]
        npopmax <- as.numeric(param[param[, 1] == "npopmax", 
            3])
        file.npop <- paste(path.mcmc.noadm, "populations.numbers.txt", 
            sep = "")
        npop.mcmc <- scan(file.npop)
        burnin <- floor(length(npop.mcmc) * 0.3)
        npop.est <- order(hist(npop.mcmc[-(1:burnin)], breaks = seq(0.5, 
            npopmax + 0.5, 1), plot = FALSE)$counts, decreasing = TRUE)[1]
        filter.null.alleles <- as.logical(param[param[, 1] == 
            "filter.null.alleles", 3])
        if (filter.null.alleles == TRUE) {
            nalmax <- nalmax + 1
        }
    }
    if (!is.null(allele.freq)) {
        freq.est <- allele.freq
    }
    if (is.null(path.mcmc.noadm)) {
        npop.est <- length(unique(cluster.indiv))
    }
    freq.est <- freq.est[1:npop.est, , ]
    dist.indiv <- as.matrix(dist(coordinates))
    if (is.null(dist.IC)) {
        dist.IC <- matrix(nrow = nindiv, ncol = npop.est)
        for (iindiv in 1:nindiv) {
            for (ipop in 1:npop.est) {
                dd <- 1e+300
                for (j in which(cluster.indiv == ipop)) {
                  dd <- min(dd, dist.indiv[iindiv, j])
                }
                dist.IC[iindiv, ipop] <- dd
            }
        }
    }
    if (estimate.a) {
        if (is.null(a.max)) 
            stop("Argument a.max missing ")
    }
    if (!estimate.a) {
        if (is.null(a.init)) 
            stop("Argument a.init missing ")
    }
    if (is.null(a.init)) {
        if (common.param) {
            a.init <- rep(runif(n = 1, 0, a.max), npop.est)
        }
        else {
            a.init <- runif(n = npop.est, 0, a.max)
        }
    }
    else {
        a.init <- rep(a.init, npop.est)
    }
    if (!estimate.b) {
        if (is.null(b.init)) 
            stop("Argument b.init missing ")
    }
    if (estimate.b) {
        if (is.null(b.max)) {
            b.max <- 10 * max(range(coordinates[, 1]), range(coordinates[, 
                2]))
        }
    }
    if (is.null(b.init)) {
        if (common.param) {
            b.init <- rep(runif(n = 1, 0, b.max), npop.est)
        }
        else {
            b.init <- runif(n = npop.est, 0, b.max)
        }
    }
    else {
        b.init <- rep(b.init, npop.est)
    }
    if (estimate.c) {
        if (is.null(c.max)) 
            stop("Argument c.max missing ")
    }
    if (!estimate.c) {
        if (is.null(c.init)) 
            stop("Argument c.init missing ")
    }
    if (is.null(c.init)) {
        if (common.param) {
            c.init <- rep(runif(n = 1, 0, c.max), npop.est)
        }
        else {
            c.init <- runif(n = npop.est, 0, c.max)
        }
    }
    else {
        c.init <- rep(c.init, npop.est)
    }
    alphadmix <- alphadmixtmp <- matrix(nrow = nindiv, ncol = npop.est, 
        data = a.init * exp(-(dist.IC/b.init)^c.init))
    q.init <- q.tmp <- matrix(nrow = nindiv, ncol = npop.est, 
        data = NA)
    for (iindiv in 1:nindiv) {
        q.init[iindiv, ] <- rdirichlet(param = alphadmix[iindiv, 
            ])
        q.tmp[iindiv, ] <- rdirichlet(param = alphadmixtmp[iindiv, 
            ])
    }
    r.coord <- range(coordinates)
    delta.b <- 0.05 * (r.coord[2] - r.coord[1])
    estimate.q <- TRUE
    a.dum <- b.dum <- c.dum <- rep(-999, npop.est)
    nchar.path.adm <- nchar(path.mcmc.adm)
    geno.dip.codom.fmt[is.na(geno.dip.codom.fmt)] <- -999
    geno.dip.dom.fmt[is.na(geno.dip.dom.fmt)] <- -999
    geno.hap.fmt[is.na(geno.hap.fmt)] <- -999
    nitstor <- nit/thinning
    qout <- array(dim = c(nitstor, nindiv, npop.est), data = -999)
    aout <- matrix(nrow = nitstor, ncol = npop.est, data = -999)
    bout <- matrix(nrow = nitstor, ncol = npop.est, data = -999)
    cout <- matrix(nrow = nitstor, ncol = npop.est, data = -999)
    print("coucou avant .Fortran")
    res <- .Fortran("mcmchz", PACKAGE = "Geneland", as.double(q.init), 
        as.double(q.tmp), as.integer(geno.dip.codom.fmt), as.integer(geno.dip.dom.fmt), 
        as.integer(geno.hap.fmt), as.integer(use.codom), as.integer(use.dom), 
        as.integer(use.hap), as.integer(npop.est), as.integer(npop.est), 
        as.integer(nindiv), as.double(freq.est), as.integer(nlocd), 
        as.integer(nloch), as.integer(nalmax), as.double(alphadmix), 
        as.double(alphadmixtmp), as.double(a.init), as.double(b.init), 
        as.double(c.init), as.double(a.dum), as.double(b.dum), 
        as.double(c.dum), as.double(a.max), as.double(b.max), 
        as.double(c.max), as.double(dist.IC), as.integer(nchar.path.adm), 
        as.character(path.mcmc.adm), as.integer(nit), as.integer(thinning), 
        as.integer(estimate.a), as.integer(estimate.b), as.integer(estimate.c), 
        as.integer(estimate.q), as.double(delta.b), as.integer(common.param), 
        as.integer(nitstor), as.double(qout), as.double(aout), 
        as.double(bout), as.double(cout))
    qout <- array(dim = c(nitstor, nindiv, npop.est), data = res[[39]])
    aout <- matrix(nrow = nitstor, ncol = npop.est, data = res[[40]])
    bout <- matrix(nrow = nitstor, ncol = npop.est, data = res[[41]])
    cout <- matrix(nrow = nitstor, ncol = npop.est, data = res[[42]])
    for (iitstor in 1:(nit/thinning)) {
        append <- ifelse(iitstor > 1, TRUE, FALSE)
        write.table(qout[iitstor, , ], paste(path.mcmc.adm, "q.txt", 
            sep = ""), row.names = FALSE, col.names = FALSE, 
            append = append)
    }
    write.table(aout, paste(path.mcmc.adm, "a.txt", sep = ""), 
        row.names = FALSE, col.names = FALSE)
    write.table(bout, paste(path.mcmc.adm, "b.txt", sep = ""), 
        row.names = FALSE, col.names = FALSE)
    write.table(cout, paste(path.mcmc.adm, "c.txt", sep = ""), 
        row.names = FALSE, col.names = FALSE)
    param.adm <- c(paste("nit :", nit), paste("thinning :", thinning), 
        paste("npop :", npop.est))
    write.table(param.adm, file = paste(path.mcmc.adm, "parameters.hz.txt", 
        sep = ""), quote = FALSE, row.names = FALSE, col.names = FALSE)
}
