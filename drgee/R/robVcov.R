robVcov <-
    function(u, d.u, id = NULL) {

        n.obs <- nrow(u)

        n.par <- ncol(d.u)

        if(!is.null(id)) {
            n.clust <- length(unique(id))
        } else {
            n.clust <- 1
        }

        inv.d.u <- try(solve(d.u))

        if (class(inv.d.u) == 'try-error') {

            vcov <- matrix(rep( NA, n.par^2 ), ncol = n.par )

        } else {
            
            ## If data is not clustered we divide by
            ## the number of observations
            if (is.null(id) | n.clust == n.obs) {
                
                correction.term <- 1 / n.obs
                
            } else {
                
                ## If data is clustered we use the clustersums of u
                ## and divide by the number of clusters
                u.dt <- as.data.table(cbind(id, u))
                setkey(u.dt, id)
                u <- as.matrix( u.dt[, lapply(.SD, sum), by = id] )[, -1]

                correction.term <- n.clust / n.obs^2
            }

            vcov <- inv.d.u %*% cov(u, u) %*% t(inv.d.u) * correction.term
            
        }

        return( vcov )
    }

