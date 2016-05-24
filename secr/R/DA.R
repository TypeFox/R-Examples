############################################################################################
## package 'secr'
## DA.R
## last changed 2011 02 15
############################################################################################

write.DA <- function (capthist, buffer, nzeros = 200, units = 1) {

    if (ms(capthist)) {
        lapply(capthist, write.DA, buffer = buffer, nzeros = nzeros, units = units)
    }
    else {
        if (!(detector(traps(capthist)) %in% c('polygon','polygonX')))
            stop ("data should be capthist with polygon detector type")
        if (!inherits(capthist, 'capthist'))
            stop ("requires 'capthist' object")
        if (!((detector(traps(capthist)) == 'polygonX') | (max(abs(capthist))==1)))
            stop ("requires binary data")
        warning ("assuming polygon is a rectangle aligned to x-y axes")
        trp <- traps(capthist)
        if (length(levels(polyID(trp))) > 1)
            stop("cannot export multiple polygons in this format")
        cpt <- matrix(nrow = nrow(capthist), ncol = ncol(capthist))
        cpt[,] <- capthist
        cpt[abs(cpt)>0] <- 1
        obs <- matrix(nrow=nrow(cpt), ncol=ncol(cpt))
        no.obs <- matrix(nrow=nzeros, ncol=ncol(cpt))
        out <- vector(mode='list')
        out$Xl <- (min(trp$x) - buffer) / units
        out$Xu <- (max(trp$x) + buffer) / units
        out$Yl <- (min(trp$y) - buffer) / units
        out$Yu <- (max(trp$y) + buffer) / units
        out$delta <- buffer / units
        out$nind <- nrow(capthist)
        out$nzeros <- nzeros
        out$T <- ncol(cpt)
        out$Y <- rbind(cpt, matrix(0, nrow=nzeros, ncol=ncol(cpt)))
        out$U1 <- obs
        out$U1[cpt > 0] <- xy(capthist)[,1]/units
        out$U1 <- rbind(out$U1, no.obs)
        out$U2 <- obs
        out$U2[cpt > 0] <- xy(capthist)[,2]/units
        out$U2 <- rbind(out$U2, no.obs)
        out
    }
}


read.DA <- function (DAlist, detector = 'polygonX', units = 1, session = 1,
    Y = 'Y', xcoord = 'U1', ycoord = 'U2', xmin = 'Xl', xmax = 'Xu',
    ymin = 'Yl', ymax = 'Yu', buffer = 'delta', verify = TRUE) {

    if (!(detector %in% c('polygon','polygonX')))
        stop ("detector type must be 'polygon' or 'polygonX'")
    ## captures
    Y <- DAlist[[Y]]    ## binary CH matrix
    session <- rep(session, prod(dim(Y)))
    ID <- as.numeric(row(Y))
    occ <- as.numeric(col(Y))
    x <- as.numeric(DAlist[[xcoord]]) * units
    y <- as.numeric(DAlist[[ycoord]]) * units
    capt <- data.frame(session, ID, occ, x, y)
    capt <- capt[!is.na(x),]
    capt <- capt[order(capt$ID, capt$occ),]
    ## detector polygon
    if (is.character(buffer))
        buffer <-  DAlist[[buffer]] * units
    xl <- DAlist[[xmin]] + buffer
    xu <- DAlist[[xmax]] - buffer
    yl <- DAlist[[ymin]] + buffer
    yu <- DAlist[[ymax]] - buffer
    trapdata <- data.frame(
        polyID = c(1,1,1,1,1),
        x = c(xl,xl,xu,xu,xl) * units,
        y = c(yl,yu,yu,yl,yl) * units
    )
    trap <- read.traps(data = trapdata, detector = detector)
    temp <- make.capthist(capt, trap, fmt='XY')
    if (verify) verify(temp)
    temp
}
