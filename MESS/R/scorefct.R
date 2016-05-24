scorefct <- function(o, beta=NULL, testidx=NULL, sas=FALSE) {

    # Check that ids are correctly ordered
    if (!ordered.clusters(o$id)) {
        stop("clusters in the gee model are not ordered and contiguous. They really should be since otherwise geepack will consider non-contiguous clusters with same id as different. Reorder your dataset or redefine the cluster id variable and run the gee fit again.")
    }

    if (any(o$weights != 1)) {
        stop("Haven't thought about if there is a problem with weights so will not do any computations")
    }
    
    clusters <- unique(o$id)
    nclusters <- length(clusters)
    if (is.null(beta)) {
        beta <- coef(o)
    }
    
    # Offsets handled correctly?
    y <- o$y
    x <- model.matrix(o)
    linear.predictors <- x%*%beta
    if (!is.null(o$offset))
        linear.predictors <- linear.predictors + o$offset
    mui <- o$family$linkinv(linear.predictors)

    invert <- if ("MASS" %in% loadedNamespaces()) {
        MASS::ginv
    } else { solve }
   
    myres <- lapply(clusters, function(cluster) {
        # Individiuals in cluster
        idx <- (o$id == cluster)  

        # Cluster size
        csize <- sum(idx)

        # Di is r*k
        Di <- t(x[idx,,drop=FALSE]) %*% diag(o$family$mu.eta(linear.predictors[idx]), nrow=csize)
        A  <- diag(sqrt(o$family$variance(mui[idx])), nrow=csize)
        Rmat <- diag(csize)
        Ralpha <- switch(o$corstr,
                         independence = Rmat,
                         exchangeable = matrix(rep(o$geese$alpha, csize^2), csize),
                         ar1 = o$geese$alpha^abs(row(Rmat) - col(Rmat)))
        if (o$corstr=="exchangeable")
            diag(Ralpha) <- 1

        V <- outer(diag(A), diag(A))*Ralpha*o$geese$gamma  # Ok

        # V inverse
        Vinv <- invert(V)
        DiVinv <- Di %*% Vinv   
        list(score = DiVinv %*% (y[idx] - mui[idx]),
             DUD   = DiVinv %*% t(Di))
    })

    ## Speed improvement? 
    # S <- apply(sapply(myres, function(oo) oo[[1]]), 1, sum)
    S <- rowSums(sapply(myres, function(oo) oo[[1]]))

    Vsand <- Reduce("+", lapply(myres, function(oo) { oo[[1]] %*% t(oo[[1]])})) # I_1
    VDUD  <- Reduce("+", lapply(myres, function(oo) { oo[[2]] })) #  I_0
    iVDUD <- invert(VDUD)
    
    if(is.null(testidx)) {
        cat("Should not really be here")
        as.numeric(S %*% invert(invert(VDUD) %*% Vsand %*% invert(VDUD)) %*% S / nclusters)
    } else {
        if (sas)
          as.numeric(t(S) %*% iVDUD[,testidx] %*% invert( (iVDUD %*%  Vsand  %*% iVDUD)[testidx,testidx] ) %*% iVDUD[testidx,] %*% S)
        else {
            myvar <- Vsand[testidx,testidx] - Vsand[testidx, -testidx] %*% invert(Vsand[-testidx,-testidx]) %*% Vsand[-testidx,testidx]
            as.numeric(t(S[testidx]) %*% invert( myvar ) %*% S[testidx])
        }
    }
}
