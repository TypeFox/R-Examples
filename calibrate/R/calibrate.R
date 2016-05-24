calibrate <-
function (g, y, tm, Fr, tmlab = tm, tl = 0.05, dt = TRUE, dp = FALSE, 
    lm = TRUE, verb = TRUE, axislab = "", reverse = FALSE, 
    alpha = NULL, labpos = 1, weights = diag(rep(1, length(y))), 
    axiscol = "blue", cex.axislab = 0.75, graphics = TRUE, where = 3, 
    laboffset = c(0, 0), m = matrix(c(0, 0), nrow = 1), markerpos = 3, 
    showlabel = TRUE, lwd = 1, shiftvec = c(0,0), shiftdir = "none", shiftfactor = 1.05 ) 
{
    g <- as.vector(g)
    if (is.matrix(weights)) 
        Dw <- weights
    else if (is.vector(weights)) 
        Dw <- diag(weights)
    else stop("calibrate: weights is not a vector or matrix")
    if (verb > 1) 
        print(Dw)
    if(all(!as.logical(shiftvec))) { # c(0,0) vector
       shiftvec <- switch(shiftdir,
                    right = shiftvector(g,Fr)$dr,
                    left = shiftvector(g,Fr)$dl,
                    none = c(0,0),
                    stop("invalid value for parameter shiftdir"))
     }
    shiftvec <- shiftfactor*shiftvec
    optalpha <- t(g) %*% t(Fr) %*% Dw %*% Fr %*% g/((t(y) %*% 
        Dw %*% Fr %*% g) * (t(g) %*% g))
    optalpha <- optalpha[1, 1]
    useralpha <- NULL
    if (is.null(alpha)) {
        alpha <- optalpha
    }
    else {
        useralpha <- alpha
    }
    Mmean <- matrix(rep(1, length(tm)), ncol = 1) %*% m
    M <- alpha * (tm) %*% t(g) + Mmean
    nrM <- nrow(M)
    M <- M + matrix(rep(1, nrM), ncol = 1) %*% shiftvec
    di <- 1/(alpha * t(g) %*% g)
    yt <- di[1, 1] * Fr %*% g
    e <- y - yt
    Q <- t(e) %*% Dw %*% e
    gos <- 1 - Q/(t(y) %*% Dw %*% y)
    odi <- 1/(optalpha * t(g) %*% g)
    oyt <- odi[1, 1] * Fr %*% g
    oe <- y - oyt
    oQ <- t(oe) %*% Dw %*% oe
    gof <- 1 - oQ/(t(y) %*% Dw %*% y)
    ang <- atan(g[2]/g[1]) * 180/pi
    lengthoneunit <- alpha * sqrt(t(g) %*% g)
    if (verb) {
        cat("---------- Calibration Results for ", axislab, " ")
        for (i in 1:(60 - (38 + nchar(axislab)))) cat("-")
        cat("\n")
        cat("Length of 1 unit of the original variable = ", round(lengthoneunit, 
            digits = 4), " \n")
        cat("Angle                                     = ", round(ang, 
            digits = 2), "degrees\n")
        cat("Optimal calibration factor                = ", round(optalpha, 
            digits = 4), " \n")
        cat("Used calibration factor                   = ", round(alpha, 
            digits = 4), " \n")
        cat("Goodness-of-fit                           = ", round(gof, 
            digits = 4), " \n")
        cat("Goodness-of-scale                         = ", round(gos, 
            digits = 4), " \n")
        cat("------------------------------------------------------------\n")
    }
    Fr2 <- Fr[, 1:2]
    nn <- t(g) %*% g
    scal <- (Fr2 %*% g)/nn[1, 1]
    Dscal <- diag(as.vector(scal))
    Fpr <- Dscal %*% matrix(rep(1, nrow(Fr)), ncol = 1) %*% t(g)
    deltax <- tl * sin(ang * pi/180)
    deltay <- tl * cos(ang * pi/180)
    if (reverse == TRUE) 
        Mn <- cbind(M[, 1] - deltax, M[, 2] + deltay)
    else Mn <- cbind(M[, 1] + deltax, M[, 2] - deltay)
    if (graphics) {
        lines(rbind(M[1, ], M[nrM, ]), col = axiscol, lwd = lwd)
        if (lm) {
            if (reverse == TRUE) 
                text(Mn[, 1], Mn[, 2], tmlab, pos = markerpos, 
                  offset = 0.2, cex = cex.axislab, srt = ang)
            else if (markerpos > 2) 
                text(Mn[, 1], Mn[, 2], tmlab, pos = markerpos - 
                  2, offset = 0.2, cex = cex.axislab, srt = ang)
            else text(Mn[, 1], Mn[, 2], tmlab, pos = markerpos + 
                2, offset = 0.2, cex = cex.axislab, srt = ang)
        }
        nm <- nrow(M)
        if (dt == TRUE) {
            for (i in 1:nm) lines(rbind(M[i, 1:2], Mn[i, 1:2]), 
                col = axiscol, lwd = lwd)
        }
        if (dp) {
            nrFpr <- nrow(Fpr)
            dlines(Fr2 + matrix(rep(1, nrFpr), ncol = 1) %*% 
                m, Fpr + matrix(rep(1, nrFpr), ncol = 1) %*% 
                m + matrix(rep(1, nrFpr), ncol = 1) %*% shiftvec)
        }
        if (showlabel) {
            switch(where, text(M[1, 1] + laboffset[1], M[1, 2] + 
                laboffset[2], axislab, cex = cex.axislab, srt = ang, pos = labpos, adj = c(0.5,0.5)), text(M[round(nrow(M)/2), 
                1] + laboffset[1], M[round(nrow(M)/2), 2] + laboffset[2], 
                axislab, cex = cex.axislab, srt = ang, pos = labpos, 
                                adj = c(0.5,0.5) ), text(M[nrow(M), 1] + laboffset[1], 
                M[nrow(M), 2] + laboffset[2], axislab, cex = cex.axislab, 
                srt = ang, pos = labpos, adj = c(0.5,0.5) ))
        }
    }
    return(list(useralpha = useralpha, gos = gos, optalpha = optalpha, 
        gof = gof, lengthoneunit = lengthoneunit, M = M, Q = Q, 
        ang = ang, shiftvec = shiftvec, yt = yt, e = e, Fpr = Fpr, Mn = Mn))
}

