CovOgk <- function(x, niter = 2, beta = 0.9, control)
{
    metodo2 <- function(XX) {

        n <- nrow(XX)
        p <- ncol(XX)

        sigma <- apply(XX, 2, mrob)[2,]
        Y <- XX %*% diag(1/sigma)
        U <- matrix(1, p, p)
        for(i in 1:p)
            for(j in i:p) {
                U[j, i] <- U[i, j] <- vrob(Y[,i], Y[,j])
        }

        diag(U) <- 1
        E <- eigen(U)$vectors
        A <- diag(sigma) %*% E
        Z <- Y %*% E

        restau <- apply(Z, 2, mrob)
        sigma <- as.vector(restau[2,])
        cov <- A %*% diag(sigma^2) %*% t(A)
        loc <- A %*% restau[1,]

        list(cov = cov, center = loc, AA = A, ZZ = Z)
    }

    ## Analize and validate the input parameters ...

    ## If a control object was supplied, take the option parameters from it,
    ##  but if single parameters were passed (not defaults) they will override the
    ##  control object.
    ## The functions 'mrob()' and 'vrob()' can be supplied only via the control
    ##  object. If no control object is passed these function will be taken
    ##  from the default one

    defcontrol <- CovControlOgk()           # default control
    mrob <- defcontrol@mrob
    vrob <- defcontrol@vrob
    smrob <- defcontrol@smrob
    svrob <- defcontrol@svrob
    if(!missing(control)){                  # a control object was supplied
        if(niter == defcontrol@niter)       niter <- control@niter
        if(beta == defcontrol@beta)         beta <- control@beta
        mrob <- control@mrob
        vrob <- control@vrob
        smrob <- control@smrob
        svrob <- control@svrob
    }

    if(is.data.frame(x))
        x <- data.matrix(x)
    else if (!is.matrix(x))
        x <- matrix(x, length(x), 1,
            dimnames = list(names(x), deparse(substitute(x))))

    ## drop all rows with missing values (!!) :
    na.x <- !is.finite(x %*% rep(1, ncol(x)))
    ok <- !na.x
    x <- x[ok, , drop = FALSE]
    dx <- dim(x)
    if(!length(dx))
        stop("All observations have missing values!")
        
    dimn <- dimnames(x)
    n <- dx[1]
    p <- dx[2]
    if(p < 2)
        stop("Need at least 2 columns ")

    if(n <= 0)
        stop("All observations have missing values!")

    call <- match.call()

    ## If the user has supplied own mrob and vrob use the pure R version
    ##  with this functions. Otherwise call the C implementation
    if(!is.null(mrob)){

        ##  iterate two times to obtain OGK2
        first <- metodo2(x)
        cov <- first$cov
        center <- as.vector(first$center)
        ZZ <- first$ZZ
        if(niter >= 2){
            second <- metodo2(first$ZZ)
            cov  <- first$AA %*% second$cov %*% t(first$AA)
            center <- as.vector(first$AA %*% as.vector(second$center))
            ZZ <- second$ZZ
        }

        dimnames(cov) <- list(dimn[[2]], dimn[[2]])
        names(center) <- dimn[[2]]

        ##  compute distances and weights
        ##  do not invert cov to compute the distances, use the transformed data
        ##
        ##  dist2 <- mahalanobis(X, center, cov)
        ##

        musigma <- apply(ZZ,2,mrob)
        ZZ <- sweep(ZZ, 2, musigma[1,])
        ZZ <- sweep(ZZ, 2, musigma[2,], '/')
        dist2 <- rowSums(ZZ^2)

        cdelta <- median(dist2)/qchisq(0.5, p)
        cov <- cdelta * cov

        quantiel <- qchisq(beta, p)
        qq <- (quantiel * median(dist2))/qchisq(0.5, p)
        wt <- ifelse(dist2 < qq, 1, 0)
        sum.wt <- sum(wt)
    } else {
        if(!(smrob %in% c("scaleTau2", "s_mad")))
            stop(paste("Scale function not defined: ", smrob))
        if(!(svrob %in% c("gk", "qc")))
            stop(paste("Bivariate covariance function not defined: ", svrob))

        storage.mode(x) <- "double"
        opw <- .Call("covOPW", x, as.integer(niter), smrob, svrob)

        dimnames(opw$cov) <- list(dimn[[2]], dimn[[2]])
        names(opw$center) <- dimn[[2]]

        dist2 <- opw$distances

        cdelta <- median(dist2)/qchisq(0.5, p)
        cov <- opw$cov <- cdelta * opw$cov
        center <- opw$center

        quantiel <- qchisq(beta, p)
        qq <- (quantiel * median(dist2))/qchisq(0.5, p)
        wt <- ifelse(dist2 < qq, 1, 0)
        sum.wt <- sum(wt)
    }

##  compute the reweighted estimates:  OGK2(0.9)
    wcenter <- colSums(x*wt)/sum.wt
    X <- sqrt(wt) * sweep(x, 2, wcenter)
    wcov <- (t(X) %*% X)/sum.wt

##  Compute consistency correction factor for the reweighted  cov
    qdelta.rew <- qchisq(sum(wt)/n, p)
    cdeltainvers.rew <- pgamma(qdelta.rew/2, p/2 + 1)/(sum(wt)/n)
    cnp2 <- 1/cdeltainvers.rew

 ##   wcov <- cnp2 * wcov

    method="Orthogonalized Gnanadesikan-Kettenring Estimator"
    ans <- new("CovOgk",
               call = call,
               iter=niter,
               crit=1,
               cov=wcov,
               center=wcenter,
               n.obs=n,
               raw.cov=cov,
               raw.center=center,
               raw.mah = dist2,
               raw.wt = wt,
               X = x,
               method=method)
    ans
}
