"swd" <- 
function(sbf) 
{
    if (class(sbf) != "sbf")
        stop("Input must be sbf class.")
        
    site <- cbind(sbf$latlon[, 1] * 0.5 * pi/90, sbf$latlon[, 2] * pi/180)

    grid <- NULL
    grid$theta <- sbf$gridlon * pi/180
    grid$phi <- sbf$gridlat * pi/180  
        
    decom <- mrafield.comp(grid, sbf$coeff, site, sbf$netlab, sbf$eta, c(t(sbf$field)), c(t(sbf$density)))     
    #detail <- mrsfield.comp(grid, sbf$coeff, site, sbf$netlab, sbf$eta)
    
    global <- density <- as.list(NULL)
    for (i in 1:(sbf$nlevels - 1)) {
        global  <- c(global,  list(t(matrix(decom$global[, i],  nrow=length(sbf$gridlat)))))
        density <- c(density, list(t(matrix(decom$density[, i], nrow=length(sbf$gridlat)))))    
    }
    detail <- as.list(NULL)
    detail <- c(detail, list(as.matrix(sbf$field - global[[1]])))
    for (i in 1:(sbf$nlevels - 2))
        detail <- c(detail, list(as.matrix(global[[i]] - global[[i+1]])))
  
    out <- list(obs=sbf$obs, latlon=sbf$latlon, 
        netlab=sbf$netlab, eta=sbf$eta, method=sbf$method, approx=sbf$approx, grid.size=sbf$grid.size, lambda=sbf$lambda, p0=sbf$p0,
        gridlat=sbf$gridlat, gridlon=sbf$gridlon, nlevels=sbf$nlevels, coeff=sbf$coef, field=sbf$field, density1=sbf$density,
        latlim=sbf$latlim, lonlim=sbf$lonlim, global=global, density=density, detail=detail, swcoeff=decom$swcoeff, thresh.info="None")
    
    class(out) <- "swd"
    
    out
}
  
"mrafield.comp" <-
function (grid, coeff, site, netlab, eta, field, density) 
{
    out <- NULL
    J <- length(netlab)
    L <- length(eta)
    t2 <- field
    d1 <- t2 * 0
    ks <- 0
    if (length(density) > 0) {
        t20 <- density
        ks <- 1
    }
    coef <- mracoef.comp(coeff, site, netlab, eta)
    for (l in 1:(L - 1)) {
        beta2 <- coef$beta[, l]
        gamma1 <- coef$gamma[, l]
        temp <- msbf.comp(grid, site, beta2, netlab, eta, l)
        t2 <- cbind(t2, temp$field)
        d1 <- cbind(d1, t2[, l] - temp$field)
        if (ks == 1) 
            t20 <- cbind(t20, temp$density)
    }
    t2 <- t2[, -1]
    d1 <- d1[, -1]
    if (ks == 1) 
        t20 <- t20[, -1]
    out$global <- as.matrix(t2)
    out$local <- as.matrix(d1)
    if (ks == 1) 
        out$density <- as.matrix(t20)
        
    sw.gamma <- coef$gamma
    sw.norm <- (coef$norm)^(0.5)
    
    swcoeff <- as.list(NULL)
    for (l in 1:(L - 1))
        swcoeff <- c(swcoeff, list((sw.gamma[, l] * sw.norm[, l])[netlab == l]))
        
    swcoeff <- c(swcoeff, list(coef$beta[, L-1][netlab == L]))         
    out$swcoeff <- swcoeff
    
    out
}

"mracoef.comp" <-
function (coef, site, netlab, eta) 
{
    out <- NULL
    J <- length(netlab)
    L <- length(eta)
    beta <- rep(0, J)
    gamma <- beta
    norm <- beta
    a1 <- mcov.comp(site, netlab, eta)
    beta1 <- coef
    for (l in 1:(L - 1)) {
        temp <- coefmatrix(beta1, a1, netlab, l)
        nnet <- netlab[netlab >= l]
        dnet <- netlab[netlab == l]
        a2 <- as.matrix(a1[netlab > l, ])
        a2 <- as.matrix(a2[, netlab > l])
        b1p <- a1[netlab == l, netlab > l]
        b1 <- a1[netlab > l, netlab == l]
        if (length(dnet) == 1) 
            b1p <- matrix(b1p, nrow = 1)
        c1 <- a1[netlab == l, netlab == l]
        e1p <- t(qr.coef(qr(a2), t(b1p)))
        beta2 <- temp$beta2
        gamma1 <- temp$gamma1
        temp <- rep(NA, J)
        temp[netlab > l] <- beta2
        beta <- cbind(beta, temp)
        temp <- rep(NA, J)
        temp[netlab == l] <- gamma1
        gamma <- cbind(gamma, temp)
        temp <- rep(NA, J)
        temp[netlab == l] <- diag(c1 - e1p %*% b1)
        norm <- cbind(norm, temp)
        beta1 <- beta2
    }
    out$beta <- as.matrix(beta[, -1])
    out$gamma <- as.matrix(gamma[, -1])
    out$norm <- as.matrix(norm[, -1])
    out
}

"mcov.comp" <-
function (site, netlab, eta) 
{
    site1 <- site[, 1]
    site2 <- site[, 2]
    cov <- cov.comp(site1, site2, eta, netlab)$aa
    cov
}

"coefmatrix" <-
function (beta1, fullcov, netlab, l) 
{
    # computes coef for 
    #       (1) [-e1p,id]:   coefs of SBF in spherical wavelets;
    #       (2) beta2:       coefs of SBF after removing subnet l;
    #       (3) gamma1:       coefs of SW for subnet l;
    # beta1:    coef of SBF from previous SBF representation
    # fullcov:  JxJ covariance matrix of all sites
    # netlab:   J labels of L subnetworks for all sites
    out <- NULL
    nnet <- netlab[netlab >= l]
    dnet <- netlab[netlab == l]
    a2 <- as.matrix(fullcov[netlab > l, ])
    a2 <- as.matrix(a2[, netlab > l])
    b1p <- fullcov[netlab == l, netlab > l]
    b1 <- fullcov[netlab > l, netlab == l]
    if (length(dnet) == 1) 
        b1p <- matrix(b1p, nrow = 1)
    c1 <- fullcov[netlab == l, netlab == l]
    e1p <- t(qr.coef(qr(a2), t(b1p)))
    ii <- diag(rep(1, length(dnet)))
    temp <- cbind(-e1p, ii)
    temp0 <- 0 * temp
    temp0[, nnet == l] <- ii
    temp0[, nnet > l] <- -e1p
    out$wcoef <- temp0
    delta1 <- beta1[nnet > l]
    gamma1 <- beta1[nnet == l]
    beta2 <- delta1 + t(e1p) %*% gamma1
    out$beta2 <- beta2
    out$gamma1 <- gamma1
    out$alpha1 <- t(out$wcoef) %*% gamma1
    out$norm <- c1 - e1p %*% b1
    out
}

"cov.comp" <-
function (site1, site2, eta, netlab) 
{
    JJ <- length(site1)
    L <- length(eta)
    aa <- matrix(0, JJ, JJ)
    temp <- matrix(0, JJ, JJ)
    x <- matrix(0, JJ, JJ)
    pp0 <- 0
    pp <- matrix(0, JJ, JJ)
    eta1 <- 0
    eta2 <- 0
    zz <- .Fortran("cov", PACKAGE = "SpherWave", aa = as.double(aa), temp = as.double(temp), 
        x = as.double(x), pp0 = as.double(pp0), pp = as.double(pp), 
        site1 = as.double(site1), site2 = as.double(site2), eta = as.double(eta), 
        eta1 = as.double(eta1), eta2 = as.double(eta2), netlab = as.integer(netlab), 
        JJ = as.integer(JJ), L = as.integer(L))
        
    list(aa = matrix(zz$aa, JJ, JJ))
}
