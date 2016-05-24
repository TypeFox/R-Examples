# allows dist="dice"
# gives out skipping indicator
# (and does all the geco and abundances stuff too) 
prabinit <- function (file = NULL, prabmatrix = NULL, rows.are.species = TRUE, 
    neighborhood = "none", nbbetweenregions=TRUE,
                      geodist = NULL, gtf = 0.1, distance = "kulczynski", 
    toprab = FALSE, toprabp = 0.05, outc = 5.2) 
{
    if (is.null(prabmatrix)) 
        m1 <- as.matrix(read.table(file))
    else m1 <- as.matrix(prabmatrix)
    if (rows.are.species) 
        m1 <- t(m1)
    regperspec <- apply(m1, 2, sum)
    specperreg <- apply(m1, 1, sum)
    nonempty.species <- regperspec > 0
    m1 <- as.matrix(m1[, regperspec > 0])
    regperspec <- regperspec[regperspec > 0]
    n.species <- ncol(m1)
    n.regions <- nrow(m1)
    nb <- list()
    if (nbbetweenregions)
          nreg <- n.regions
    else
          nreg <- n.species
    if (is.list(neighborhood)) 
        nb <- neighborhood
    else {
        if (identical(neighborhood,"none")) 
            for (i in 1:nreg) nb[[i]] <- numeric(0)
        else for (i in 1:nreg) nb <- c(nb, list(scan(file = neighborhood, 
            skip = i - 1, nlines = 1, quiet = TRUE)))
    }
    nbtest(nb, nreg)
    if (toprab) {
        nscut <- toprabp * specperreg
#        if (toprabp == 0) 
#            nscut <- min(m1[m1 > 0])/10
        minno <- c()
        for (i in 1:n.species) {
            g0 <- log(m1[, i][m1[, i] > 0])
            minno[i] <- median(g0) - mad(g0, constant = outc)
        }
        for (i in 1:n.regions) for (j in 1:n.species) m1[i, j] <- as.integer(m1[i, 
            j] > nscut[i] | log(m1[i, j]) > minno[j])
        regperspec <- apply(m1, 2, sum)
        specperreg <- apply(m1, 1, sum)
    }
    distmat <- switch(distance, kulczynski = kulczynski(m1), 
        jaccard = jaccard(m1), geco = geco(m1, geodist, tf = gtf), 
        qkulczynski = qkulczynski(m1), logkulczynski = qkulczynski(m1, 
            log.distance = TRUE), dice=dicedist(m1),
                      none = NULL)
    if (identical(levels(as.factor(m1)),c("0","1"))){
      prabm <- as.matrix(as.logical(m1))
      dim(prabm) <- dim(m1)
      m1 <- prabm
    }
    out <- list(distmat = distmat, prab = m1, nb = nb, regperspec = regperspec, 
        specperreg = specperreg, n.species = n.species, n.regions = n.regions, 
        distance = distance, geodist = geodist, gtf = gtf, spatial = (!identical(neighborhood, 
            "none")), nonempty.species=nonempty.species,
                nbbetweenregions=nbbetweenregions)
    class(out) <- "prab"
    out
}

# 
# prabinit <- function (file = NULL, prabmatrix = NULL, rows.are.species = TRUE, 
#     neighborhood = "none", geodist = NULL, gtf = 0.1, distance = "kulczynski", 
#     toprab = FALSE, toprabp = 0.05, outc = 5.2) 
# {
#     if (is.null(prabmatrix)) 
#         m1 <- as.matrix(read.table(file))
#     else m1 <- as.matrix(prabmatrix)
#     if (rows.are.species) 
#         m1 <- t(m1)
#     regperspec <- apply(m1, 2, sum)
#     specperreg <- apply(m1, 1, sum)
#     nonzero <- regperspec>0
#     m1 <- as.matrix(m1[, regperspec > 0])
#     regperspec <- regperspec[regperspec > 0]
#     n.species <- ncol(m1)
#     n.regions <- nrow(m1)
#     nb <- list()
#     if (is.list(neighborhood)) 
#         nb <- neighborhood
#     else {
#         if (neighborhood == "none") 
#             for (i in 1:n.regions) nb[[i]] <- numeric(0)
#         else for (i in 1:n.regions) nb <- c(nb, list(scan(file = neighborhood, 
#             skip = i - 1, nlines = 1, quiet = TRUE)))
#     }
#     nbtest(nb, n.regions)
#     if (toprab) {
#         nscut <- toprabp * specperreg
#         if (toprabp == 0) 
#             nscut <- min(m1[m1 > 0])/10
#         minno <- c()
#         for (i in 1:n.species) {
#             g0 <- log(m1[, i][m1[, i] > 0])
#             minno[i] <- median(g0) - mad(g0, constant = outc)
#         }
#         for (i in 1:n.regions)
#           for (j in 1:n.species)
#             m1[i, j] <- as.integer(m1[i,j] >= nscut[i] | log(m1[i, j]) >= minno[j])
#         regperspec <- apply(m1, 2, sum)
#         specperreg <- apply(m1, 1, sum)
#     }
#     distmat <- switch(distance, kulczynski = kulczynski(m1), 
#         jaccard = jaccard(m1), geco = geco(m1, geodist, tf = gtf), 
#         qkulczynski = qkulczynski(m1), logkulczynski=qkulczynski(m1,
#                                          log.distance=TRUE), none = NULL)
#     out <- list(distmat = distmat, prab = m1, nb = nb, regperspec = regperspec, 
#         specperreg = specperreg, n.species = n.species, n.regions = n.regions, 
#         distance = distance, geodist = geodist, gtf = gtf,
#                 spatial = (!identical(neighborhood, "none")),
#                 nonzero=nonzero)
#     class(out) <- "prab"
#     out
# }
