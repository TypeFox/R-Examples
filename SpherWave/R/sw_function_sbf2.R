"sbf" <- 
function (obs, latlon, netlab, eta, method, approx=FALSE, grid.size=c(50, 100), lambda=NULL, p0=0, 
    latlim=NULL, lonlim=NULL) 
{
    if (length(obs) != nrow(latlon))
        stop("The number of data must be the same to the number of locations.")
    if (method != "ls" && method != "pls")
        stop("Only least squares method or penalized least squares method is proveided.")
    if (method == "pls" && is.null(lambda))
        stop("Provide smoothing parameter lambda for penalized least squares method.")

    if (is.null(latlim)) latlim <- c(-90, 90)
    if (is.null(lonlim)) lonlim <- c(-180, 180)
 
    site <- cbind(latlon[, 1] * 0.5 * pi/90, latlon[, 2] * pi/180)
    if (approx == FALSE) {
        ssite <- site
        snetlab <- netlab
        seta <- eta
    } else { 
        up <- max(netlab)
        ssite <- NULL
        snetlab <- NULL
        seta <- eta[2:up] 
        for(i in 2:up) {
            ssite <- rbind(ssite, site[netlab == i, ]) 
            snetlab <- c(snetlab, rep(i, nrow(site[netlab == i, ])))
        }
        snetlab <- snetlab - 1
    }

    grid <- mesh(grid.size[2] - 1, grid.size[1] - 1)

    if (method == "ls")
        coef <- as.vector(ls.comp(obs, site, ssite, snetlab, seta)$coef) 
    
    if (method == "pls")
        coef <- as.vector(ridge.comp(obs, site, ssite, snetlab, seta, lambda)$coef)

    if (approx == FALSE)
        tmp <- msbf.comp(grid, site, coef, netlab, eta, p0)
    else 
        tmp <- msbf.comp(grid, ssite, coef, snetlab, seta, p0)
    
    gridlat <- grid$phi * 180/pi
    gridlon <- grid$theta * 180/pi
    
    out <- list(obs=obs, latlon=latlon, 
        netlab=netlab, eta=eta, method=method, approx=approx, grid.size=grid.size, lambda=lambda, p0=p0,
        gridlon=gridlon, gridlat=gridlat, nlevels=length(unique(netlab)), coeff=coef, 
        field=t(matrix(tmp$field, nrow=length(gridlat))), 
        density=t(matrix(tmp$density, nrow=length(gridlat))), latlim=latlim, lonlim=lonlim)
        
    class(out) <- "sbf"
    out
}
 
"mesh" <-
function (M, N) 
{
    theta <- 2 * pi * (1:M - ceiling(M/2))/M
    phi <- pi * (1:N - ceiling(N/2))/N
    
    list(theta = theta, phi = phi)
}


"ls.comp" <-
function (obs, site, ssite, snet, seta) 
{
    site1 <- site[, 1]
    site2 <- site[, 2]
    ssite1 <- ssite[, 1]
    ssite2 <- ssite[, 2]
    gg <- lscoef.comp(site1, site2, ssite1, ssite2, snet, seta)$gg
    temp <- lsfit(gg, obs, intercept = FALSE)
}

"lscoef.comp" <-
function (site1, site2, ssite1, ssite2, snet, seta) 
{
    KK <- length(site1)
    JJ <- length(ssite1)
    p <- length(seta)
    gg <- matrix(0, KK, JJ)
    temp <- matrix(0, JJ, KK)
    x <- matrix(0, JJ, KK)
    zz <- .Fortran("ls", PACKAGE = "SpherWave", gg = as.double(gg), temp = as.double(temp), 
        ssite1 = as.double(ssite1), ssite2 = as.double(ssite2), 
        site1 = as.double(site1), site2 = as.double(site2), x = as.double(x), 
        snet = as.integer(snet), seta = as.double(seta), JJ = as.integer(JJ), 
        KK = as.integer(KK), p = as.integer(p))
        
    list(gg = matrix(zz$gg, KK, JJ))
}

"ridge.comp" <-
function (obs, site, ssite, snet, seta, lam) 
{
    J <- length(snet)
    site1 <- site[, 1]
    site2 <- site[, 2]
    ssite1 <- ssite[, 1]
    ssite2 <- ssite[, 2]
    gg <- lscoef.comp(site1, site2, ssite1, ssite2, snet, seta)$gg
    ggg <- gg.comp(site1, site2, ssite1, ssite2, snet, seta, lam)$gg
    gg0 <- t(gg) %*% gg
    for (j in 1:J) 
        gg0[j, j] <- gg0[j, j] + lam^2
    gg0 <- solve(gg0)
    sumsq <- 0
    for (j in 1:J) 
        sumsq <- sumsq + gg0[j, j]
    temp <- c(obs, rep(0, J))
    temp <- lsfit(ggg, temp, intercept = FALSE)
    out <- temp
    out$sumsq <- sumsq
    out
}

"gg.comp" <-
function (site1, site2, ssite1, ssite2, snet, seta, lam) 
{
    KK <- length(site1)
    JJ <- length(ssite1)
    p <- length(seta)
    gg <- matrix(0, KK + JJ, JJ)
    temp <- matrix(0, JJ, KK)
    x <- matrix(0, JJ, KK)
    zz <- .Fortran("ridge", PACKAGE = "SpherWave", gg = as.double(gg), temp = as.double(temp), 
        ssite1 = as.double(ssite1), ssite2 = as.double(ssite2), 
        site1 = as.double(site1), site2 = as.double(site2), x = as.double(x), 
        snet = as.integer(snet), seta = as.double(seta), JJ = as.integer(JJ), 
        KK = as.integer(KK), p = as.integer(p), lam = as.double(lam))
        
    list(gg = matrix(zz$gg, KK + JJ, JJ))
}


"msbf.comp" <-
function (grid, site, coef, netlab, eta, p0) 
{
    out <- NULL
    point1 <- grid$phi
    point2 <- grid$theta
    site1 <- site[, 1]
    site2 <- site[, 2]
    tmp <- sbf.comp(point1, point2, site1, site2, coef, netlab, eta, p0)
    aa <- tmp$aa
    bb <- tmp$bb
    m <- length(point2)
    mm <- NULL
    for (j in 1:m)
        mm <- c(mm, aa[j, ])

    out$field <- mm
    
    mm <- NULL
    for (j in 1:m)
        mm <- c(mm, bb[j, ])

    out$density <- mm
    out
}

"sbf.comp" <-
function (point1, point2, site1, site2, coef, netlab, eta, p0) 
{
    m <- length(point2)
    n <- length(point1)
    JJ <- length(site1)
    p <- length(eta)
    para <- rep(1, length(coef))
    for (i in 1:p) 
        if (is.na(max(coef[netlab == i]))) {
            coef[netlab == i] <- 0
            para[netlab == i] <- 0
        }

    aa <- matrix(0, m, n)
    bb <- matrix(0, m, n)
    temp <- rep(0, JJ)
    pp <- rep(0, JJ)
    ppp <- rep(0, JJ)
    x <- rep(0, JJ)
    stemp <- 0
    spp <- 0
    beta <- 0
    zz <- .Fortran("sbf", PACKAGE = "SpherWave", aa = as.double(aa), bb = as.double(bb), 
        stemp = as.double(stemp), temp = as.double(temp), spp = as.double(spp), 
        pp = as.double(pp), ppp = as.double(ppp), point1 = as.double(point1), 
        point2 = as.double(point2), site1 = as.double(site1), 
        site2 = as.double(site2), x = as.double(x), coef = as.double(coef), 
        beta = as.double(beta), netlab = as.integer(netlab), 
        para = as.integer(para), eta = as.double(eta), m = as.integer(m), 
        n = as.integer(n), JJ = as.integer(JJ), p = as.integer(p), 
        p0 = as.integer(p0))
        
    list(aa = matrix(zz$aa, m, n), bb = matrix(zz$bb, m, n))
}

"eta.comp" <-
function (netlab) 
{
    n <- length(netlab[netlab == max(netlab)])
    nlevel <- max(netlab)
    first.eta <- bandwidth(n)$h
    first.rho <- -log(first.eta)
    rho <- NULL
    for (l in 1:(nlevel - 1)) 
        rho <- c(rho, first.rho/2^l)
    
    eta <- c(first.eta, exp(-rho))
    
    list(eta = eta[nlevel:1])
}

"bandwidth" <-
function (n) 
{
    nvar <- 1 - (1 - (2/n))^2
    nsig <- sqrt(nvar)
    h2 <- (1 - nsig)/(1 + nsig)
    h <- sqrt(h2)

    list(h=h)
}


"gcv.lambda" <- 
function (obs, latlon, netlab, eta, approx=FALSE, lambda) 
{       
    site <- cbind(latlon[, 1] * 0.5 * pi/90, latlon[, 2] * pi/180)
    if (approx == FALSE) {
        ssite <- site
        snetlab <- netlab
        seta <- eta
    } else { 
        up <- max(netlab)
        ssite <- NULL
        snetlab <- NULL
        seta <- eta[2:up] 
        for(i in 2:up) {
            ssite <- rbind(ssite, site[netlab == i, ]) 
            snetlab <- c(snetlab, rep(i, nrow(site[netlab == i, ])))
        }
        snetlab <- snetlab - 1
    }   
    
    out <- ridge.diacomp(ridge.comp(obs, site, ssite, snetlab, seta, lambda), obs, lambda)
    list(gcv = out$gcv)
}


"ridge.diacomp" <-
function (out.ls, obs, lam) 
{
    resid <- out.ls$residuals
    sumsq <- out.ls$sumsq
    nd <- length(obs)
    rsq <- 1 - sum(resid[1:nd]^2)/sum((obs - mean(obs))^2)
    temp <- ls.diag(out.ls)
    df <- sum(temp$hat) - sumsq * lam^2
    gcv <- mean(resid[1:nd]^2/(1 - df/nd)^2)
    
    list(rsq = rsq, gcv = gcv, df = df)
}

"pk" <-
function (theta, eta) 
{
# normalized poisson kernel as a function of angle
    (1 - eta^2) * (1 - eta)^2/(1 + eta)/(1 - 2 * eta * cos(pi * theta) + eta^2)^(3/2) 
}
 
"sw.plot" <- 
function (sw=NULL, z=NULL, latlon=NULL, latlim=NULL, lonlim=NULL, type="field", nlevel=256, pch=NULL, cex=NULL, ...) 
{
    #oldpar <- par(no.readonly = TRUE)

    if (type != "obs" & type != "network" & type != "field" & type != "swcoeff" & type != "decom" & type != "recon")
        stop("Possible types are obs, network, field, decom and recon.")

    if (!is.null(sw)) latlim <- sw$latlim
    else if (is.null(latlim)) latlim <- c(-90, 90)
 
    if (!is.null(sw)) lonlim <- sw$lonlim
    else if (is.null(lonlim)) lonlim <- c(-180, 180)

    if ((is.null(sw) | class(sw) == "sbf" | class(sw) == "swd") & type == "obs") { 
        if (class(sw) == "sbf" | class(sw) == "swd") {
 
            z <- sw$obs
            latlon <- sw$latlon
        }
      
        old.par <- par(no.readonly = TRUE)

        temp <- imageplot.setup()
        smallplot <- temp$smallplot
        bigplot <- temp$bigplot
                
        timcolors <- tim.colors(nlevel)
 
        rangez <- range(z)    
        boundarys <- seq(rangez[1], rangez[2], length = nlevel + 1)
        
        ix <- 1
        iy <- (boundarys + (boundarys[2] - boundarys[1])/2)[1:nlevel]
        iz <- matrix(iy, nrow = 1, ncol = length(iy)) 
           
        indexcol <- NULL
        for (i in 1:length(z)) {
            tmp <- (1:nlevel)[z[i] >= boundarys[-(nlevel + 1)] & z[i] < boundarys[-1]]
            if(length(tmp) == 0)
                tmp <- nlevel
            indexcol <- c(indexcol, tmp)
        }
        
        if (is.null(pch)) pch <- 20
        if (is.null(cex)) cex <- 1
               
        par(plt = bigplot)
        plot(latlon[1, 2], latlon[1, 1], col = timcolors[indexcol[1]], ylim = latlim, xlim = lonlim, 
            pch = pch, cex = cex, ...)    
        for (i in 2:length(z))
            points(latlon[i, 2], latlon[i, 1], col = timcolors[indexcol[i]], pch = pch, cex = cex, ...) 
        world(add = TRUE, ylim = latlim, xlim = lonlim)     

        big.par <- par(no.readonly = TRUE)
        if ((smallplot[2] < smallplot[1]) | (smallplot[4] < smallplot[3])) {
            par(old.par)
            stop("plot region too small to add legend\n")
        }
                
        par(new = TRUE, pty = "m", plt = smallplot, err = -1)
        image(ix, iy, iz, xaxt = "n", yaxt = "n", xlab = "", ylab = "", col = timcolors)
        axis(4, mgp = c(3, 1, 0), las = 2)
        
        mfg.save <- par()$mfg
        par(big.par)
        par(plt = big.par$plt, xpd = FALSE)
        par(mfg = mfg.save, new = FALSE)
        invisible()
    }   
    if ((is.null(sw) | class(sw) == "sbf" | class(sw) == "swd") & type == "network") {
        if (class(sw) == "sbf" | class(sw) == "swd") {
            z <- sw$netlab
            latlon <- sw$latlon
        }
        
        #network.plot(latlon, z, latlim, lonlim, ...)
        #world(add = TRUE, ylim = latlim, xlim = lonlim)
        
        old.par <- par(no.readonly = TRUE)
        
        temp <- imageplot.setup(legend.mar = 1, legend.width = 1.2 * 5)
        smallplot <- temp$smallplot
        bigplot <- temp$bigplot
                
        nlab <- length(unique(z))
        timcolors <- rev(tim.colors(nlab))

        if (is.null(pch)) pch <- 19
        if (is.null(cex)) cex <- 0.6
                
        par(plt = bigplot)
        tmp <- latlon[z == 1, ]
        plot(tmp[, 2], tmp[, 1], ylim = latlim, xlim = lonlim, pch = pch, cex = cex, col = timcolors[1], ...)        
 
        for (i in 2:nlab) {
            tmp <- latlon[z == i, ]
            points(tmp[, 2], tmp[, 1], pch = pch - 1 + i, cex = cex - 0.1 + i/10, col = timcolors[i], ...)
        }
        world(add = TRUE, ylim = latlim, xlim = lonlim)
        
        big.par <- par(no.readonly = TRUE)
        if ((smallplot[2] < smallplot[1]) | (smallplot[4] < smallplot[3])) {
            par(old.par)
            stop("plot region too small to add legend\n")
        }
                  
        par(new = TRUE, pty = "m", plt = smallplot, err = -1)
        plot(c(1, 2), c(0, 2 * nlab), type = "n", axes = FALSE, xlab = "", ylab = "")
        legend(x = 1, y = (2 * nlab - nlab/2), legend = paste("Level", 1:nlab, sep=" "),
            pch = (pch):(pch - 1 + nlab), cex = 0.8, col = timcolors)
            
        mfg.save <- par()$mfg
        par(big.par)
        par(plt = big.par$plt, xpd = FALSE)
        par(mfg = mfg.save, new = FALSE)
        invisible()
    } 
    if ((class(sw) == "sbf" | class(sw) == "swd") & type == "field") {
        image.plot(sw$gridlon, sw$gridlat, sw$field, nlevel = nlevel,
            zlim=range(c(sw$obs, sw$field)), ylim = latlim, xlim = lonlim, ...)
        world(add = TRUE, ylim = latlim, xlim = lonlim)
    }     
    if (class(sw) == "swd" && type == "swcoeff") {
        par(mfrow=c(sw$nlevels %/% 2 + sw$nlevels %% 2, 2), mar = c(0.3, 1.5, 1.45, 1.5) + 0.1)
        
        timcolors <- tim.colors(nlevel)
        if (is.null(pch)) pch <- 20
        if (is.null(cex)) cex <- 1.2
                
        for (i in 1:sw$nlevels) {
            rangez <- range(sw$swcoeff[[i]])
            boundarys <- seq(rangez[1], rangez[2], length = nlevel + 1) 
            
            ix <- 1
            iy <- (boundarys + (boundarys[2] - boundarys[1])/2)[1:nlevel]
            iz <- matrix(iy, nrow = 1, ncol = length(iy))
        
            indexcol <- NULL
            for (j in 1:length(sw$swcoeff[[i]])) {
                tmp <- (1:nlevel)[sw$swcoeff[[i]][j] >= boundarys[-(nlevel + 1)] & sw$swcoeff[[i]][j] < boundarys[-1]]
                if(length(tmp) == 0)
                    tmp <- nlevel
                indexcol <- c(indexcol, tmp)
            }
        
            old.par <- par(no.readonly = TRUE)

            temp <- imageplot.setup()
            smallplot <- temp$smallplot
            bigplot <- temp$bigplot
            
            #plot(1, type = "n", axes = FALSE, xlab = "", ylab = "", ylim = latlim, xlim = lonlim, ...)
            #par(new = TRUE, pty = "m")
            par(plt = bigplot)
            if (i != sw$nlevels) {
                plot(sw$latlon[sw$netlab == i, 2], sw$latlon[sw$netlab == i, 1], 
                    col = timcolors[indexcol], cex = cex^i, pch = pch, main = paste("Detailed coefficients of level", i),
                    xlab="", ylab="", axes = FALSE, ylim = latlim, xlim = lonlim, ...)
                box()
                #mtext(paste("Level", i), side = 3, line = 0.1, cex = 0.9) 
                world(add = TRUE, ylim = latlim, xlim = lonlim)    
            } else {
                plot(sw$latlon[sw$netlab == i, 2], sw$latlon[sw$netlab == i, 1], 
                    col = timcolors[indexcol], cex = cex^i, pch = pch, main = paste("Global coefficients of level", i),
                    xlab="", ylab="", axes = FALSE, ylim = latlim, xlim = lonlim, ...)
                box()
                #mtext(paste("Level", i), side = 3, line = 0.1, cex = 0.9) 
                world(add = TRUE, ylim = latlim, xlim = lonlim)             
            }
            
            big.par <- par(no.readonly = TRUE)
            if ((smallplot[2] < smallplot[1]) | (smallplot[4] < smallplot[3])) {
                par(old.par)
                stop("plot region too small to add legend\n")
            }           
            
            par(new = TRUE, pty = "m", plt = smallplot, err = -1)
            image(ix, iy, iz, xaxt = "n", yaxt = "n", xlab = "", ylab = "", col = timcolors)
            axis(4, mgp = c(3, 1, 0), las = 2)
        
            mfg.save <- par()$mfg
            par(big.par)
            par(plt = big.par$plt, xpd = FALSE)
            par(mfg = mfg.save, new = FALSE)
            invisible()
        }
                   
        # Alternative 
        #timcolors <- tim.colors(nlevel)
        
        #z <- unlist(sw$swcoeff)
        #rangez <- range(z)
        #boundarys <- seq(rangez[1], rangez[2], length = nlevel + 1)

        #ix <- 1
        #iy <- (boundarys + (boundarys[2] - boundarys[1])/2)[1:nlevel]
        #iz <- matrix(iy, nrow = 1, ncol = length(iy))
        
        #indexcol <- as.list(NULL)
        #for (i in 1:sw$nlevels) {
        #    tmpi <- NULL
        #    for (j in 1:length(sw$swcoeff[[i]])) {
        #        tmp <- (1:nlevel)[sw$swcoeff[[i]][j] >= boundarys[-(nlevel + 1)] & sw$swcoeff[[i]][j] < boundarys[-1]]
        #        if(length(tmp) == 0)
        #            tmp <- nlevel
        #        tmpi <- c(tmpi, tmp)
        #    }
        #    indexcol <- c(indexcol, list(tmpi))
        #}
        
        #smallplot <- c(0.28675617, 0.71324383, 0.04, 0.066) # 0.08098695 0.10711177
        #bigplot12 <- c(0.02612482, 0.50000000 - 0.02612482/2, 0.50000000 + 0.02612482/2, 0.97387518)
        #bigplot34 <- seq(1, 0.09, length = 4) 

        #if (is.null(pch)) pch <- 20
        #if (is.null(cex)) cex <- 1.1
        
        #plot(1, type = "n", axes = FALSE, xlab = "", ylab = "", ylim = latlim, xlim = lonlim, ...)
        #for (i in 1:sw$nlevels) {
        #    if (i %% 2 == 1)
        #        par(new = TRUE, pty = "m", plt = c(bigplot12[1:2], 
        #            bigplot34[c(i %/% 2 + i %% 2 + 1, i %/% 2 + i %% 2)] - c(0, 0.034))) #0.02612482
        #    else
        #        par(new = TRUE, pty = "m", plt = c(bigplot12[3:4], 
        #            bigplot34[c(i %/% 2 + i %% 2 + 1, i %/% 2 + i %% 2)] - c(0, 0.034))) #0.02612482
        #
        #    plot(sw$latlon[sw$netlab == i, 2], sw$latlon[sw$netlab == i, 1], 
        #        col = timcolors[indexcol[[i]]], cex = cex^i, pch = pch,
        #        xlab="", ylab="", axes = FALSE, ylim = latlim, xlim = lonlim, ...)
        #    box()
        #    mtext(paste("Level", i), side = 3, line = 0.1, cex = 0.9) 
        #    world(add = TRUE, ylim = latlim, xlim = lonlim)     
        #}
        
        #par(new = TRUE, pty = "m", plt = smallplot, err = -1)
        #image(iy, ix, t(iz), xaxt = "n", yaxt = "n", xlab = "", ylab = "", col = tim.colors(nlevel))
        #axis(1, mgp = c(3, 1, 0)/2)        
    }        
    if (class(sw) == "swd" & type == "decom") { 
        par(mfrow = c(sw$nlevels - 1, 2), mar = c(0.3, 1.5, 1.4, 1.5) + 0.1)
        rangeglobal <- range(unlist(sw$global))

        for (i in 1:(sw$nlevels - 1)) { 
            image.plot(sw$gridlon, sw$gridlat, sw$global[[i]], axes = FALSE, xlab = "", ylab = paste("Level", i+1), 
                mgp = c(0.5,0,0), nlevel = nlevel, zlim = rangeglobal, ylim = latlim, xlim = lonlim, ...)
            if (i == 1) title("Global Component", line = 0.3)
            world(add = TRUE, ylim = latlim, xlim = lonlim) 
            if (is.list(sw$thresh.info))
                image.plot(sw$gridlon, sw$gridlat, sw$detail[[i]], axes = FALSE, xlab = "", ylab = paste("Level", i), 
                    mgp = c(0.5,0,0), nlevel = nlevel, ylim = latlim, xlim = lonlim, zlim = sw$thresh.info$rangedetail[i, ], ...)
            else
                image.plot(sw$gridlon, sw$gridlat, sw$detail[[i]], axes = FALSE, xlab = "", ylab = paste("Level", i), 
                    mgp = c(0.5,0,0), nlevel = nlevel, ylim = latlim, xlim = lonlim, ...)            
            if (i == 1) title("Local Component", line = 0.3)
            world(add = TRUE, ylim = latlim, xlim = lonlim)
        }

        # Alternative         
        #rangez <- range(c(unlist(sw$global), unlist(sw$detail)))
        #boundarys <- seq(rangez[1], rangez[2], length=nlevel+1)
 
        #ix <- 1
        #iy <- (boundarys + (boundarys[2] - boundarys[1])/2)[1:nlevel]
        #iz <- matrix(iy, nrow = 1, ncol = length(iy))

        #temp <- imageplot.setup(horizon = TRUE)
        #smallplot <- temp$smallplot
        #bigplot <- temp$bigplot
        #bigplot[3] <- bigplot[3] - (bigplot[3] - smallplot[4]) * 0.8
        #smallplot[1:2] <- seq(smallplot[1], smallplot[2], length=5)[c(2,4)]
        
        #bigplot34 <- seq(bigplot[4], bigplot[3], length=sw$nlevels)
        #bigplot12 <- seq(bigplot[1], bigplot[2], length=3)
    
        #plot(1, type="n", axes=FALSE, xlab="", ylab="", ylim=latlim, xlim=lonlim)
        #for (i in 1:(sw$nlevels - 1)) { 
        #    par(new = TRUE, pty="m", plt=c(bigplot12[1:2], bigplot34[c(i+1, i)]))
        #    image(sw$gridlon, sw$gridlat, sw$global[[i]], xlab="", ylab="", axes = FALSE, col=tim.colors(nlevel), 
        #       zlim=rangez, ylim=latlim, xlim=lonlim)
        #    mtext(paste("Level", i), side = 2, line=1)               
        #    world(add = TRUE, ylim=latlim, xlim=lonlim)
        #    par(new = TRUE, plt=c(bigplot12[2:3], bigplot34[c(i+1, i)]))
        #    image(sw$gridlon, sw$gridlat, sw$detail[[i]], xlab="", ylab="", axes = FALSE, col=tim.colors(nlevel),
        #        zlim=rangez, ylim=latlim, xlim=lonlim)
        #    world(add = TRUE, ylim=latlim, xlim=lonlim)
        #}    
        
        #par(new = TRUE, pty = "m", plt = smallplot, err = -1)
        #image(iy, ix, t(iz), xaxt = "n", yaxt = "n", xlab = "", ylab = "", col = tim.colors(nlevel))
        #axis(1, mgp = c(3, 1, 0)/2)
    }
    if (is.null(sw) & type == "recon") {
        image.plot(seq(lonlim[1], lonlim[2], length = nrow(z)), seq(latlim[1], latlim[2], length = ncol(z)), z, 
            nlevel = nlevel, ylim = latlim, xlim = lonlim, ...)
        world(add = TRUE, ylim = latlim, xlim = lonlim)
    }
    #par(oldpar)
}  
