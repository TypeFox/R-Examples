condit <-
    function(y,
             x,
             y.names = colnames(y),
             x.names = colnames(x), 
             id) {

        x.cent <- .Call("center", x, id, PACKAGE = "drgee")
        colnames(x) <- x.names
        colnames(x.cent) <- x.names

        ## Assuming that the observations are sorted by id
        
        ## Find the total number of observations
        n.obs <- length(id)

        ## Find the number of parameters
        n.params <- ncol(x.cent) 

        ## Indices in the original matrix
        idx <- seq_len(n.obs)

        y.dt <- data.table(idx = idx,
                           id = id,
                           y = y)

        names(y.dt) <- c("idx", "id", y.names)
        
        setkeyv(y.dt, "id")
        ## with(y.dt, setkeyv(y.dt, id) )

        ## Find cluster sizes
        ## clust.size <- y.dt[, .N, by = "id"][[2]]
        clusters <- data.table( with(y.dt,
                                     y.dt[, list(clust.size = .N,
                                                 ysum = sum(get(y.names)),
                                                 min.idx = min(idx)), 
                                          by = id]) )

        ## Find cluster sums
        ## ysum <- y.dt[, sum(y), by = "id"][[2]]
        ## with(y.dt, y.dt[, ysum := sum(y), by = id])

        ## Find cluster sums
        ## min.idx <- y.dt[, min(idx), by = "id"][[2]]
        ## with(y.dt, y.dt[, min.idx := min(idx), by = id])

        with(clusters, clusters[, disc := TRUE])
        with(clusters, clusters[ysum == 0 | ysum == clust.size, disc := FALSE])

        with(clusters, clusters[, inv := 0L])
        with(clusters, clusters[ysum >  0.5 * clust.size, inv := 1L])

        disc.clusters <- with(clusters, clusters[which(disc), ])
        
        ## disc.idx <- with(clusters, clusters[which(disc), min.idx])
        ## disc.ids <- with(clusters, clusters[which(disc), id])
        ## Find number of clusters
        ## n.clust <- length(ysum)
        ## n.clust <- nrow(clusters)
        
        ## Find concordant clusters
        ## disc <- rep(TRUE, n.clust)
        ## disc[clusters$ysum == 0 | clusters$ysum == clusters$clust.size] <- FALSE

        ## Find cluster with more than half being cases
        ## inv <- integer(n.clust)
        ## inv[clusters$ysum > 0.5 * clusters$clust.size] <- 1L

        ## ysum.disc <- ysum[disc]
        ## clust.size.disc <- clust.size[disc]
        ## min.idx.disc <- min.idx[disc]
        ## inv.disc <- inv[disc]

        ## yd.data <- cbind(y, x, id)
        ## colnames(yd.data) <- c(y.names, x.names, "id")

        yd.data <- data.frame(y, x, id)
        names(yd.data) <- c(y.names, x.names, "id")

        oformula <- formula( paste(y.names,
                                   " ~ ",
                                   paste(x.names,
                                         collapse = " + "),
                                   " + strata(id)", 
                                   sep = "") )

        ## ########################################################
        ## Fit using clogit
        ## ########################################################

        fit.clogit <- try(survival::clogit(oformula,
                                            data = yd.data, 
                                            method = "exact",
                                            subset = id %in% disc.clusters$id ) )

        if (class(fit.clogit)[1] == "try-error") {

            beta.hat <- rep(NA, n.params)

            eq.res <- rep(NA, n.obs)

            d.eq.res <- matrix( rep(NA, n.obs*n.params), nrow = n.obs )

            naive.var = NULL

        } else {

            beta.hat <- coef(fit.clogit)

            resids <- .Call("conditRes",
                            beta.hat,
                            disc.clusters$ysum,
                            disc.clusters$clust.size,
                            disc.clusters$min.idx,
                            disc.clusters$inv,
                            as.vector(y),
                            x.cent,
                            PACKAGE = "drgee")

            eq.res <- as.vector(resids$res)

            d.eq.res <- resids$dres

            naive.var <- vcov(fit.clogit)

        }

        return( list(coefficients = beta.hat,
                     res = eq.res,
                     d.res = d.eq.res,
                     eq.x = x.cent,
                     optim.object = NULL,
                     naive.var = naive.var) )

    }
