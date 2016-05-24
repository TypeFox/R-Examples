

setGeneric("getDesign", function(umf, ...) standardGeneric("getDesign"))
setGeneric("handleNA", function(umf, ...) standardGeneric("handleNA"))



# unmarkedFrame

setMethod("getDesign", "unmarkedFrame",
    function(umf, formula, na.rm=TRUE)
{
    detformula <- as.formula(formula[[2]])
    stateformula <- as.formula(paste("~", formula[3], sep=""))
    detVars <- all.vars(detformula)

    M <- numSites(umf)
    R <- obsNum(umf)

    ## Compute state design matrix
    if(is.null(siteCovs(umf))) {
        siteCovs <- data.frame(placeHolder = rep(1, M))
    } else {
        siteCovs <- siteCovs(umf)
    }
    X.mf <- model.frame(stateformula, siteCovs, na.action = NULL)
    X <- model.matrix(stateformula, X.mf)
    X.offset <- as.vector(model.offset(X.mf))
    if (!is.null(X.offset)) {
        X.offset[is.na(X.offset)] <- 0
    }

    ## Compute detection design matrix
    if(is.null(obsCovs(umf))) {
        obsCovs <- data.frame(placeHolder = rep(1, M*R))
    } else {
        obsCovs <- obsCovs(umf)
    }

    ## Record future column names for obsCovs
    colNames <- c(colnames(obsCovs), colnames(siteCovs))

    ## add site Covariates at observation-level
    obsCovs <- cbind(obsCovs, siteCovs[rep(1:M, each = R),])
    colnames(obsCovs) <- colNames

    ## add observation number if not present
    if(!("obsNum" %in% names(obsCovs))) {
        obsCovs <- cbind(obsCovs, obsNum = as.factor(rep(1:R, M)))
    }

    V.mf <- model.frame(detformula, obsCovs, na.action = NULL)
    V <- model.matrix(detformula, V.mf)
    V.offset <- as.vector(model.offset(V.mf))
    if (!is.null(V.offset)) {
        V.offset[is.na(V.offset)] <- 0
    }

    if (na.rm) {
        out <- handleNA(umf, X, X.offset, V, V.offset)
        y <- out$y
        X <- out$X
        X.offset <- out$X.offset
        V <- out$V
        V.offset <- out$V.offset
        removed.sites <- out$removed.sites
    } else {
        y=getY(umf)
        removed.sites=integer(0)
    }

    return(list(y = y, X = X, X.offset = X.offset, V = V,
                V.offset = V.offset, removed.sites = removed.sites))
})


setMethod("handleNA", "unmarkedFrame",
          function(umf, X, X.offset, V, V.offset)
{
    obsToY <- obsToY(umf)
    if(is.null(obsToY)) stop("obsToY cannot be NULL to clean data.")

    J <- numY(umf)
    R <- obsNum(umf)
    M <- numSites(umf)

    X.long <- X[rep(1:M, each = J),]
    X.long.na <- is.na(X.long)

    V.long.na <- apply(V, 2, function(x) {
        x.mat <- matrix(x, M, R, byrow = TRUE)
        x.mat <- is.na(x.mat)
        x.mat <- x.mat %*% obsToY
        x.long <- as.vector(t(x.mat))
        x.long > 0
        })
    V.long.na <- apply(V.long.na, 1, any)

    y.long <- as.vector(t(getY(umf)))
    y.long.na <- is.na(y.long)

    covs.na <- apply(cbind(X.long.na, V.long.na), 1, any)

    ## are any NA in covs not in y already?
    y.new.na <- covs.na & !y.long.na

    if(sum(y.new.na) > 0) {
        y.long[y.new.na] <- NA
        warning("Some observations have been discarded because corresponding covariates were missing.", call. = FALSE)
    }

    y <- matrix(y.long, M, J, byrow = TRUE)
    sites.to.remove <- apply(y, 1, function(x) all(is.na(x)))

    num.to.remove <- sum(sites.to.remove)
    if(num.to.remove > 0) {
        y <- y[!sites.to.remove, ,drop = FALSE]
        X <- X[!sites.to.remove, ,drop = FALSE]
        X.offset <- X.offset[!sites.to.remove]
        V <- V[!sites.to.remove[rep(1:M, each = R)], ,drop = FALSE]
        V.offset <- V.offset[!sites.to.remove[rep(1:M, each = R)], ]
        warning(paste(num.to.remove,"sites have been discarded because of missing data."), call. = FALSE)
        }

    list(y = y, X = X, X.offset = X.offset, V = V, V.offset = V.offset,
        removed.sites = which(sites.to.remove))
    })


# unmarkedFrameOccuFP
# like the occuFP but there are 3 osbservation formula which are stored in V (true positive detections),
# U (false positive detections), and W (b or probability detetion is certain)


setMethod("getDesign", "unmarkedFrameOccuFP",
          function(umf, detformula,FPformula,Bformula = ~.,stateformula, na.rm=TRUE)
          {

            M <- numSites(umf)
            R <- obsNum(umf)

            ## Compute state design matrix
            if(is.null(siteCovs(umf))) {
              siteCovs <- data.frame(placeHolder = rep(1, M))
            } else {
              siteCovs <- siteCovs(umf)
            }
            X.mf <- model.frame(stateformula, siteCovs, na.action = NULL)
            X <- model.matrix(stateformula, X.mf)
            X.offset <- as.vector(model.offset(X.mf))
            if (!is.null(X.offset)) {
              X.offset[is.na(X.offset)] <- 0
            }

            ## Compute detection design matrix
            if(is.null(obsCovs(umf))) {
              obsCovs <- data.frame(placeHolder = rep(1, M*R))
            } else {
              obsCovs <- obsCovs(umf)
            }

            ## Record future column names for obsCovs
            colNames <- c(colnames(obsCovs), colnames(siteCovs))

            ## add site Covariates at observation-level
            obsCovs <- cbind(obsCovs, siteCovs[rep(1:M, each = R),])
            colnames(obsCovs) <- colNames

            ## add observation number if not present
            if(!("obsNum" %in% names(obsCovs))) {
              obsCovs <- cbind(obsCovs, obsNum = as.factor(rep(1:R, M)))
            }



            V.mf <- model.frame(detformula, obsCovs, na.action = NULL)
            V <- model.matrix(detformula, V.mf)
            V.offset <- as.vector(model.offset(V.mf))
            if (!is.null(V.offset)) {
              V.offset[is.na(V.offset)] <- 0
            }


            U.mf <- model.frame(FPformula, obsCovs, na.action = NULL)
            U <- model.matrix(FPformula, U.mf)
            U.offset <- as.vector(model.offset(U.mf))
            if (!is.null(U.offset)) {
              U.offset[is.na(U.offset)] <- 0
            }

            W.mf <- model.frame(Bformula, obsCovs, na.action = NULL)
            W <- model.matrix(Bformula, W.mf)
            W.offset <- as.vector(model.offset(W.mf))
            if (!is.null(W.offset)) {
              W.offset[is.na(W.offset)] <- 0
            }

            if (na.rm) {
              out <- handleNA(umf, X, X.offset, V,V.offset, U, U.offset, W, W.offset)
              y <- out$y
              X <- out$X
              X.offset <- out$X.offset
              V <- out$V
              V.offset <- out$V.offset
              U <- out$U
              U.offset <- out$U.offset
              U <- out$U
              U.offset <- out$U.offset
              removed.sites <- out$removed.sites
            } else {
              y=getY(umf)
              removed.sites=integer(0)
            }



            return(list(y = y, X = X, X.offset = X.offset, V = V,
                        V.offset = V.offset,U = U, U.offset = U.offset,W = W,
                        W.offset = W.offset, removed.sites = removed.sites))
          })


setMethod("handleNA", "unmarkedFrameOccuFP",
          function(umf, X, X.offset, V, V.offset, U, U.offset, W, W.offset)
          {
            obsToY <- obsToY(umf)
            if(is.null(obsToY)) stop("obsToY cannot be NULL to clean data.")

            J <- numY(umf)
            R <- obsNum(umf)
            M <- numSites(umf)

            X.long <- X[rep(1:M, each = J),]
            X.long.na <- is.na(X.long)

            V.long.na <- apply(V, 2, function(x) {
              x.mat <- matrix(x, M, R, byrow = TRUE)
              x.mat <- is.na(x.mat)
              x.mat <- x.mat %*% obsToY
              x.long <- as.vector(t(x.mat))
              x.long > 0
            })
            V.long.na <- apply(V.long.na, 1, any)

            U.long.na <- apply(U, 2, function(x) {
              x.mat <- matrix(x, M, R, byrow = TRUE)
              x.mat <- is.na(x.mat)
              x.mat <- x.mat %*% obsToY
              x.long <- as.vector(t(x.mat))
              x.long > 0
            })
            U.long.na <- apply(U.long.na, 1, any)

            W.long.na <- apply(W, 2, function(x) {
              x.mat <- matrix(x, M, R, byrow = TRUE)
              x.mat <- is.na(x.mat)
              x.mat <- x.mat %*% obsToY
              x.long <- as.vector(t(x.mat))
              x.long > 0
            })
            W.long.na <- apply(W.long.na, 1, any)

            y.long <- as.vector(t(getY(umf)))
            y.long.na <- is.na(y.long)

            covs.na <- apply(cbind(X.long.na, V.long.na, U.long.na, W.long.na), 1, any)

            ## are any NA in covs not in y already?
            y.new.na <- covs.na & !y.long.na

            if(sum(y.new.na) > 0) {
              y.long[y.new.na] <- NA
              warning("Some observations have been discarded because corresponding covariates were missing.", call. = FALSE)
            }

            y <- matrix(y.long, M, J, byrow = TRUE)
            sites.to.remove <- apply(y, 1, function(x) all(is.na(x)))

            num.to.remove <- sum(sites.to.remove)
            if(num.to.remove > 0) {
              y <- y[!sites.to.remove, ,drop = FALSE]
              X <- X[!sites.to.remove, ,drop = FALSE]
              X.offset <- X.offset[!sites.to.remove]
              V <- V[!sites.to.remove[rep(1:M, each = R)], ,drop = FALSE]
              V.offset <- V.offset[!sites.to.remove[rep(1:M, each = R)], ]
              U <- U[!sites.to.remove[rep(1:M, each = R)], ,drop = FALSE]
              U.offset <- U.offset[!sites.to.remove[rep(1:M, each = R)], ]
              W <- W[!sites.to.remove[rep(1:M, each = R)], ,drop = FALSE]
              W.offset <- W.offset[!sites.to.remove[rep(1:M, each = R)], ]
              warning(paste(num.to.remove,"sites have been discarded because of missing data."), call. = FALSE)
            }

            list(y = y, X = X, X.offset = X.offset, V = V, V.offset = V.offset,
                 U = U, U.offset = U.offset, W = W, W.offset = W.offset,
                 removed.sites = which(sites.to.remove))
          })


# UnmarkedMultFrame




setMethod("getDesign", "unmarkedMultFrame",
    function(umf, formula, na.rm = TRUE) {

    aschar1 <- as.character(formula)
    aschar2 <- as.character(formula[[2]])
    aschar3 <- as.character(formula[[2]][[2]])

    detformula <- as.formula(paste(aschar1[1], aschar1[3]))
    epsformula <- as.formula(paste(aschar2[1], aschar2[3]))
    gamformula <- as.formula(paste(aschar3[1], aschar3[3]))
    psiformula <- as.formula(formula[[2]][[2]][[2]])

    detVars <- all.vars(detformula)

    M <- numSites(umf)
    R <- obsNum(umf)
    nY <- umf@numPrimary
    J <- R / nY

    ## Compute phi design matrices
    if(is.null(umf@yearlySiteCovs)) {
        yearlySiteCovs <- data.frame(placeHolder = rep(1, M*nY))
    } else {
        yearlySiteCovs <- umf@yearlySiteCovs
    }
    ## add siteCovs in so they can be used as well
    if(!is.null(umf@siteCovs)) {
        sC <- umf@siteCovs[rep(1:M, each = nY),,drop=FALSE]
        yearlySiteCovs <- cbind(yearlySiteCovs, sC)
        }

    ## Compute site-level design matrix for psi
    if(is.null(siteCovs(umf)))
        siteCovs <- data.frame(placeHolder = rep(1, M))
    else
        siteCovs <- siteCovs(umf)

    W.mf <- model.frame(psiformula, siteCovs, na.action = NULL)
    if(!is.null(model.offset(W.mf)))
        stop("offsets not currently allowed in colext", call.=FALSE)
    W <- model.matrix(psiformula, W.mf)


    ## Compute detection design matrix
    if(is.null(obsCovs(umf)))
        obsCovs <- data.frame(placeHolder = rep(1, M*R))
    else
        obsCovs <- obsCovs(umf)

    ## add site and yearlysite covariates, which contain siteCovs
    cnames <- c(colnames(obsCovs), colnames(yearlySiteCovs))
    obsCovs <- cbind(obsCovs, yearlySiteCovs[rep(1:(M*nY), each = J),])
    colnames(obsCovs) <- cnames

    ## add observation number if not present
    if(!("obsNum" %in% names(obsCovs)))
        obsCovs <- cbind(obsCovs, obsNum = as.factor(rep(1:R, M)))

    V.mf <- model.frame(detformula, obsCovs, na.action = NULL)
    if(!is.null(model.offset(V.mf)))
        stop("offsets not currently allowed in colext", call.=FALSE)
    V <- model.matrix(detformula, V.mf)

    ## in order to drop factor levels that only appear in last year,
    ## replace last year with NAs and use drop=TRUE
    yearlySiteCovs[seq(nY,M*nY,by=nY),] <- NA
    yearlySiteCovs <- as.data.frame(lapply(yearlySiteCovs, function(x) {
        x[,drop = TRUE]
        }))

    X.mf.gam <- model.frame(gamformula, yearlySiteCovs, na.action = NULL)
    if(!is.null(model.offset(X.mf.gam)))
        stop("offsets not currently allowed in colext", call.=FALSE)
    X.gam <- model.matrix(gamformula, X.mf.gam)
    X.mf.eps <- model.frame(epsformula, yearlySiteCovs, na.action = NULL)
    if(!is.null(model.offset(X.mf.eps)))
        stop("offsets not currently allowed in colext", call.=FALSE)
    X.eps <- model.matrix(epsformula, X.mf.eps)

    if(na.rm)
        out <- handleNA(umf, X.gam, X.eps, W, V)
    else
        out <- list(y=getY(umf), X.gam=X.gam, X.eps=X.eps, W=W, V=V,
            removed.sites=integer(0))

    return(list(y = out$y, X.eps = out$X.eps, X.gam = out$X.gam, W = out$W,
        V = out$V, removed.sites = out$removed.sites))
})






setMethod("handleNA", "unmarkedMultFrame",
    function(umf, X.gam, X.eps, W, V)
{
    obsToY <- obsToY(umf)
    if(is.null(obsToY)) stop("obsToY cannot be NULL to clean data.")

    R <- obsNum(umf)
    M <- numSites(umf)
    nY <- umf@numPrimary
    J <- numY(umf) / nY

    ## treat both X's #######no: and W together
#    X <- cbind(X.gam, X.eps, W[rep(1:M, each = nY), ])
    X <- cbind(X.gam, X.eps)

    X.na <- is.na(X)
    X.na[seq(nY,M*nY,by=nY),] <- FALSE  ## final years are unimportant.
                                        ## not true for W covs!!!
    W.expand <- W[rep(1:M, each=nY),,drop=FALSE]
    W.na <- is.na(W.expand)
    X.na <- cbind(X.na, W.na) # NAs in siteCovs results in removed site

    X.long.na <- X.na[rep(1:(M*nY), each = J),]

    V.long.na <- apply(V, 2, function(x) {
        x.mat <- matrix(x, M, R, byrow = TRUE)
        x.mat <- is.na(x.mat)
        x.mat <- x.mat %*% obsToY
        x.long <- as.vector(t(x.mat))
        x.long > 0
    })
    V.long.na <- apply(V.long.na, 1, any)

    y.long <- as.vector(t(getY(umf)))
    y.long.na <- is.na(y.long)

    # It doesn't make sense to combine X.gam/eps with W here b/c
    # a X.eps does not map correctly to y
    covs.na <- apply(cbind(X.long.na, V.long.na), 1, any)

    ## are any NA in covs not in y already?
    y.new.na <- covs.na & !y.long.na

    if(sum(y.new.na) > 0) {
        y.long[y.new.na] <- NA
        warning("Some observations have been discarded because correspoding covariates were missing.", call. = FALSE)
        }

    y <- matrix(y.long, M, numY(umf), byrow = TRUE)
    sites.to.remove <- apply(y, 1, function(x) all(is.na(x)))
#    Perhaps we need to remove sites that have no data in T=1
#    noDataT1 <- apply(is.na(y[,1:J]), 1, all) #
#    sites.to.remove <- sites.to.remove | noDataT1

    num.to.remove <- sum(sites.to.remove)
    if(num.to.remove > 0) {
        y <- y[!sites.to.remove, ,drop = FALSE]
        X.gam <- X.gam[!sites.to.remove[rep(1:M, each = nY)],,drop = FALSE]
        X.eps <- X.eps[!sites.to.remove[rep(1:M, each = nY)],,drop = FALSE]
        W <- W[!sites.to.remove,, drop = FALSE] # !!! Recent bug fix
        V <- V[!sites.to.remove[rep(1:M, each = R)], ,drop = FALSE]
        warning(paste(num.to.remove,"sites have been discarded because of missing data."), call. = FALSE)
    }
    list(y = y, X.gam = X.gam, X.eps = X.eps, W = W, V = V,
        removed.sites = which(sites.to.remove))
})






# pcountOpen
setMethod("getDesign", "unmarkedFramePCO",
    function(umf, formula, na.rm = TRUE)
{
    aschar1 <- as.character(formula)
    aschar2 <- as.character(formula[[2]])
    aschar3 <- as.character(formula[[2]][[2]])
    aschar4 <- as.character(formula[[2]][[2]][[2]])

    iotaformula <- as.formula(paste(aschar1[1], aschar1[3]))
    pformula <- as.formula(paste(aschar2[1], aschar2[3]))
    omformula <- as.formula(paste(aschar3[1], aschar3[3]))
    gamformula <- as.formula(paste(aschar4[1], aschar4[3]))
    lamformula <- as.formula(formula[[2]][[2]][[2]][[2]])

    y <- getY(umf)
    M <- nrow(y)
    T <- umf@numPrimary
    J <- ncol(y) / T
    delta <- umf@primaryPeriod

    if(is.null(umf@yearlySiteCovs))
        yearlySiteCovs <- data.frame(placeHolder = rep(1, M*T))
    else
        yearlySiteCovs <- umf@yearlySiteCovs

    ## add siteCovs in so they can be used as well
    if(!is.null(umf@siteCovs)) {
        sC <- umf@siteCovs[rep(1:M, each = T),,drop=FALSE]
        yearlySiteCovs <- cbind(yearlySiteCovs, sC)
        }

    if(is.null(siteCovs(umf)))
        siteCovs <- data.frame(placeHolder = rep(1, M))
    else
        siteCovs <- siteCovs(umf)

    Xlam.mf <- model.frame(lamformula, siteCovs, na.action = NULL)
    Xlam <- model.matrix(lamformula, Xlam.mf)
    Xlam.offset <- as.vector(model.offset(Xlam.mf))
    if(!is.null(Xlam.offset))
        Xlam.offset[is.na(Xlam.offset)] <- 0

    if(is.null(obsCovs(umf)))
        obsCovs <- data.frame(placeHolder = rep(1, M*J*T))
    else
        obsCovs <- obsCovs(umf)

    colNames <- c(colnames(obsCovs), colnames(yearlySiteCovs))

    # Add yearlySiteCovs, which contains siteCovs
    obsCovs <- cbind(obsCovs, yearlySiteCovs[rep(1:(M*T), each = J),])
    colnames(obsCovs) <- colNames

    if(!("obsNum" %in% names(obsCovs)))
        obsCovs <- cbind(obsCovs, obsNum = as.factor(rep(1:(J*T), M)))

    # Ignore last year of data
    transCovs <- yearlySiteCovs[-seq(T, M*T, by=T),,drop=FALSE]
    for(i in 1:ncol(transCovs))
        if(is.factor(transCovs[,i]))
            transCovs[,i] <- factor(transCovs[,i]) # drop unused levels

    Xiota.mf <- model.frame(iotaformula, transCovs, na.action = NULL)
    Xiota <- model.matrix(iotaformula, Xiota.mf)
    Xiota.offset <- as.vector(model.offset(Xiota.mf))
    if(!is.null(Xiota.offset))
        Xiota.offset[is.na(Xiota.offset)] <- 0
    Xp.mf <- model.frame(pformula, obsCovs, na.action = NULL)
    Xp <- model.matrix(pformula, Xp.mf)
    Xp.offset <- as.vector(model.offset(Xp.mf))
    if(!is.null(Xp.offset))
        Xp.offset[is.na(Xp.offset)] <- 0
    Xgam.mf <- model.frame(gamformula, transCovs, na.action = NULL)
    Xgam <- model.matrix(gamformula, Xgam.mf)
    Xgam.offset <- as.vector(model.offset(Xgam.mf))
    if(!is.null(Xgam.offset))
        Xgam.offset[is.na(Xgam.offset)] <- 0
    Xom.mf <- model.frame(omformula, transCovs, na.action = NULL)
    Xom <- model.matrix(omformula, Xom.mf)
    Xom.offset <- as.vector(model.offset(Xom.mf))
    if(!is.null(Xom.offset))
        Xom.offset[is.na(Xom.offset)] <- 0

    # determine if gamma, omega, and iota are scalar, vector, or matrix valued
    # Runtime is much faster for scalars and vectors
    Xgo <- cbind(Xgam, Xom, Xiota)
    getGOdims <- function(x) {
        xm <- matrix(x, M, T-1, byrow=TRUE)
#        anyNA <- apply(is.na(xm), 1, any)
#        if(all(anyNA))
#            return("matrix")
#        xm <- xm[!anyNA,] # This is not 100% safe
#        nSites <- nrow(xm)
#        if(all(dim(unique(xm, MARGIN=1)) == c(1, T-1)))
#            return("rowvec")
#        else if(all(dim(unique(xm, MARGIN=2)) == c(nSites, 1)))
#            return("colvec")
#        else return("matrix")
        col.table <- apply(xm, 2, table)
        row.table <- apply(xm, 1, table)
        if(is.vector(col.table) & !is.list(col.table)) {
            return("rowvec")
        } else if(is.vector(row.table) & !is.list(row.table)) {
            return("colvec")
        } else
            return("matrix")
        }
    if(isTRUE(all.equal(gamformula,~1)) & isTRUE(all.equal(omformula, ~1)) &
      isTRUE(all.equal(iotaformula, ~1)))
        go.dims <- "scalar"
    else {
        go.dims.vec <- apply(Xgo, 2, getGOdims)
        if(all(go.dims.vec == "rowvec"))
            go.dims <- "rowvec"
        else if(all(go.dims.vec == "colvec"))
            go.dims <- "matrix" ##"colvec"  ## NOTE: Temporary fix to the problem reported with time-only-varying covariates
        else
            go.dims <- "matrix"
    }

    if(na.rm)
        out <- handleNA(umf, Xlam, Xgam, Xom, Xp, Xiota,
            Xlam.offset, Xgam.offset, Xom.offset, Xp.offset, Xiota.offset,
            delta)
    else {   # delta needs to be formatted first
        ya <- array(y, c(M, J, T))
        yna <- apply(is.na(ya), c(1,3), all)
        delta <- formatDelta(delta, yna)
        out <- list(y=y, Xlam=Xlam, Xgam=Xgam, Xom=Xom, Xp=Xp, Xiota=Xiota,
                    Xlam.offset=Xlam.offset, Xgam.offset=Xgam.offset,
                    Xom.offset=Xom.offset, Xp.offset=Xp.offset,
                    Xiota.offset=Xiota.offset,
                    delta=delta, removed.sites=integer(0))
    }

    return(list(y = out$y, Xlam = out$Xlam, Xgam = out$Xgam,
                Xom = out$Xom, Xp = out$Xp, Xiota = out$Xiota,
                Xlam.offset=Xlam.offset, Xgam.offset=Xgam.offset,
                Xom.offset=Xom.offset, Xp.offset=Xp.offset,
                Xiota.offset=Xiota.offset, delta = out$delta,
                removed.sites = out$removed.sites, go.dims = go.dims))
})








setMethod("handleNA", "unmarkedFramePCO",
    function(umf, Xlam, Xgam, Xom, Xp, Xiota, Xlam.offset, Xgam.offset,
             Xom.offset, Xp.offset, Xiota.offset, delta)
{
	obsToY <- obsToY(umf)
	if(is.null(obsToY)) stop("obsToY cannot be NULL to clean data.")

	M <- numSites(umf)
	T <- umf@numPrimary
	y <- getY(umf)
	J <- ncol(y) / T

	Xlam.long <- Xlam[rep(1:M, each = J*T),]
	Xlam.long.na <- is.na(Xlam.long)

	long.na <- function(x) {
            x.mat <- matrix(x, M, J*T, byrow = TRUE)
            x.mat <- is.na(x.mat)
            x.mat <- x.mat %*% obsToY
            x.long <- as.vector(t(x.mat))
            x.long > 0
        }

        o2y2 <- diag(T)
        o2y2 <- o2y2[-T, -T]

	long.na2 <- function(x) {
            x.mat <- matrix(x, M, T-1, byrow = TRUE)
            x.mat <- is.na(x.mat)
            x.mat <- x.mat %*% o2y2
            x.long <- as.vector(t(x.mat))
            x.long > 0
        }

	Xp.long.na <- apply(Xp, 2, long.na)
	Xp.long.na <- apply(Xp.long.na, 1, any)

        Xgam.long.na <- apply(Xgam, 2, long.na2)
	Xgam.long.na <- apply(Xgam.long.na, 1, any)
	Xom.long.na <- apply(Xom, 2, long.na2)
	Xom.long.na <- apply(Xom.long.na, 1, any)

	y.long <- as.vector(t(y))
	y.long.na <- is.na(y.long)

#  delta.long <- as.vector(t(delta))
#	delta.long.na <- is.na(delta.long)

	covs.na <- apply(cbind(Xlam.long.na, Xp.long.na), 1, any)
	covs.na2 <- apply(cbind(Xgam.long.na, Xom.long.na), 1, any)
	covs.na3 <- rep(covs.na2, each=J)
	# If gamma[1, 1] is NA, remove y[1, 2]
	#common <- 1:(M*J*(T-1))
        ignore <- rep(seq(1, M*J*T, by=J*T), each=J) + 0:(J-1)
        covs.na[-ignore] <- covs.na[-ignore] | covs.na3

	## are any NA in covs not in y already?
	y.new.na <- covs.na & !y.long.na

	if(sum(y.new.na) > 0) {
            y.long[y.new.na] <- NA
            warning("Some observations have been discarded because corresponding covariates were missing.", call. = FALSE)
        }

	y.wide <- matrix(y.long, nrow=M, ncol=J*T, byrow = TRUE)
#	delta <- matrix(delta.long, nrow=M, ncol=T, byrow = TRUE)
	sites.to.remove <- apply(y.wide, 1, function(x) all(is.na(x)))
  # Should also remove sites with no omega and gamma before an observation
  # remove all observations before the one after the first real omega/gamma
#  covs.na2.mat <- matrix(covs.na2, M, T-1, byrow=TRUE)
#  last.y <- apply(y, 1, function(x) max(which(!is.na(y))))
#  last.go <- apply(covs.na2.mat, 1, function(x)
#      if(any(x)) {
#          if(all(x))
#              return(0)
#          else
#              return(max(which(!x)))
#          }
#      else
#          return(T-1))
#  no.go.before.y <- last.y <= last.go

	ya <- array(y.wide, c(M, J, T))
	yna <- apply(is.na(ya), c(1,3), all)
	delta <- formatDelta(delta, yna)

	num.to.remove <- sum(sites.to.remove)
	if(num.to.remove > 0) {
            y.wide <- y.wide[!sites.to.remove, ,drop = FALSE]
            Xlam <- Xlam[!sites.to.remove, ,drop = FALSE]
            Xlam.offset <- Xlam.offset[!sites.to.remove]
            Xgam <- Xgam[!sites.to.remove[rep(1:M, each = T-1)],,
                         drop = FALSE]
            Xgam.offset <- Xgam.offset[!sites.to.remove[rep(1:M, each = T-1)],,
                         drop = FALSE]
            Xom <- Xom[!sites.to.remove[rep(1:M, each = T-1)],,
                       drop = FALSE]
            Xom.offset <- Xom.offset[!sites.to.remove[rep(1:M, each = T-1)],,
                       drop = FALSE]
            Xp <- Xp[!sites.to.remove[rep(1:M, each = J*T)],,
                     drop = FALSE]
            Xp.offset <- Xp.offset[!sites.to.remove[rep(1:M, each = J*T)],,
                     drop = FALSE]
            Xiota <- Xiota[!sites.to.remove[rep(1:M, each = T-1)],,
                       drop = FALSE]
            Xiota.offset <- Xiota.offset[!sites.to.remove[rep(1:M, each = T-1)],,
                       drop = FALSE]
            delta <- delta[!sites.to.remove, ,drop =FALSE]
            warning(paste(num.to.remove, "sites have been discarded because of missing data."), call.=FALSE)
	}

	list(y = y.wide, Xlam = Xlam, Xgam = Xgam, Xom = Xom, Xp = Xp, Xiota = Xiota,
             Xlam.offset=Xlam.offset, Xgam.offset=Xgam.offset,
             Xom.offset=Xom.offset, Xp.offset=Xp.offset,
             Xiota.offset=Xiota.offset, delta = delta,
             removed.sites = which(sites.to.remove))
})










# UnmarkedFrameGMN


setMethod("getDesign", "unmarkedFrameG3",
    function(umf, formula, na.rm = TRUE)
{
    ac1 <- as.character(formula)
    ac2 <- as.character(formula[[2]])

    detformula <- as.formula(paste(ac1[1], ac1[3]))
    phiformula <- as.formula(paste(ac2[1], ac2[3]))
    lamformula <- as.formula(formula[[2]][[2]])

    detVars <- all.vars(detformula)

    M <- numSites(umf)
    T <- umf@numPrimary
    R <- obsNum(umf) # 2*T for double observer sampling
                     # 1*T for distance sampling
                     # nPasses*T for removal sampling

    ## Compute phi design matrices
    if(is.null(umf@yearlySiteCovs)) {
        yearlySiteCovs <- data.frame(placeHolder = rep(1, M*T))
    } else yearlySiteCovs <- umf@yearlySiteCovs

    # add siteCovs in so they can be used as well
    if(!is.null(umf@siteCovs)) {
        sC <- umf@siteCovs[rep(1:M, each = T),,drop=FALSE]
        yearlySiteCovs <- cbind(yearlySiteCovs, sC)
        }

    Xphi.mf <- model.frame(phiformula, yearlySiteCovs, na.action = NULL)
    Xphi <- model.matrix(phiformula, Xphi.mf)
    Xphi.offset <- as.vector(model.offset(Xphi.mf))
    if(!is.null(Xphi.offset)) Xphi.offset[is.na(Xphi.offset)] <- 0

    # Compute site-level design matrix for lambda
    if(is.null(siteCovs(umf))) {
        siteCovs <- data.frame(placeHolder = rep(1, M))
    } else siteCovs <- siteCovs(umf)
    Xlam.mf <- model.frame(lamformula, siteCovs, na.action = NULL)
    Xlam <- model.matrix(lamformula, Xlam.mf)
    Xlam.offset <- as.vector(model.offset(Xlam.mf))
    if(!is.null(Xlam.offset)) Xlam.offset[is.na(Xlam.offset)] <- 0

    # Compute detection design matrix
    if(is.null(obsCovs(umf))) {
        obsCovs <- data.frame(placeHolder = rep(1, M*R))
    } else obsCovs <- obsCovs(umf)

    # add site and yearlysite covariates, which contain siteCovs
    cnames <- c(colnames(obsCovs), colnames(yearlySiteCovs))
    obsCovs <- cbind(obsCovs, yearlySiteCovs[rep(1:(M*T), each = R/T),])
    colnames(obsCovs) <- cnames

    # add observation number if not present
    if(!("obsNum" %in% names(obsCovs)))
        obsCovs <- cbind(obsCovs, obsNum = as.factor(rep(1:R, M)))

    Xdet.mf <- model.frame(detformula, obsCovs, na.action = NULL)
    Xdet <- model.matrix(detformula, Xdet.mf)
    Xdet.offset <- as.vector(model.offset(Xdet.mf))
    if(!is.null(Xdet.offset)) Xdet.offset[is.na(Xdet.offset)] <- 0

    if(na.rm)
        out <- handleNA(umf, Xlam, Xlam.offset, Xphi, Xphi.offset, Xdet,
            Xdet.offset)
    else
        out <- list(y=getY(umf), Xlam=Xlam, Xlam.offset = Xlam.offset,
            Xphi=Xphi, Xphi.offset = Xphi.offset, Xdet=Xdet,
				    removed.sites=integer(0))

    return(list(y = out$y, Xlam = out$Xlam, Xphi = out$Xphi,
                Xdet = out$Xdet,
                Xlam.offset = out$Xlam.offset,
                Xphi.offset = out$Xphi.offset,
                Xdet.offset = out$Xdet.offset,
                removed.sites = out$removed.sites))
})





setMethod("handleNA", "unmarkedFrameG3",
    function(umf, Xlam, Xlam.offset, Xphi, Xphi.offset, Xdet, Xdet.offset)
{

#    browser()

    obsToY <- obsToY(umf)
    if(is.null(obsToY)) stop("obsToY cannot be NULL to clean data.")

    M <- numSites(umf)
    T <- umf@numPrimary
    R <- obsNum(umf)
    J <- numY(umf)/T

    # treat Xphi and Xlam together
    X <- cbind(Xphi, Xlam[rep(1:M, each = T), ])

    X.na <- is.na(X)
    X.long.na <- X.na[rep(1:(M*T), each = J),]

    Xdet.long.na <- apply(Xdet, 2, function(x) {
        x.mat <- matrix(x, M, R, byrow = TRUE)
        x.mat <- is.na(x.mat)
        x.mat <- x.mat %*% obsToY
        x.long <- as.vector(t(x.mat))
        x.long > 0
        })

    Xdet.long.na <- apply(Xdet.long.na, 1, any)

    y.long <- as.vector(t(getY(umf)))
    y.long.na <- is.na(y.long)

    covs.na <- apply(cbind(X.long.na, Xdet.long.na), 1, any)

    ## are any NA in covs not in y already?
    y.new.na <- covs.na & !y.long.na

    if(sum(y.new.na) > 0) {
        y.long[y.new.na] <- NA
        warning("Some observations have been discarded because correspoding covariates were missing.", call. = FALSE)
    }

    y <- matrix(y.long, M, numY(umf), byrow = TRUE)
    sites.to.remove <- apply(y, 1, function(x) all(is.na(x)))

    num.to.remove <- sum(sites.to.remove)
    if(num.to.remove > 0) {
        y <- y[!sites.to.remove,, drop = FALSE]
        Xlam <- Xlam[!sites.to.remove,, drop = FALSE]
        Xlam.offset <- Xlam.offset[!sites.to.remove]
        Xphi <- Xphi[!sites.to.remove[rep(1:M, each = T)],, drop = FALSE]
        Xphi.offset <- Xphi.offset[!sites.to.remove[rep(1:M, each = T)]]
        Xdet <- Xdet[!sites.to.remove[rep(1:M, each = R)],,
                     drop=FALSE]
        Xdet.offset <- Xdet.offset[!sites.to.remove[rep(1:M, each=R)]]
        warning(paste(num.to.remove,
                      "sites have been discarded because of missing data."), call.=FALSE)
    }
    list(y = y, Xlam = Xlam, Xlam.offset = Xlam.offset, Xphi = Xphi,
        Xphi.offset = Xphi.offset, Xdet = Xdet, Xdet.offset = Xdet.offset,
        removed.sites = which(sites.to.remove))
})

