areaahulleval <-
function (x) 
{
    if (class(x) != "ahull") {
        cat("Argument is not of class ahull.\n")
        return(invisible())
    }
    ah <- x
    a <- ah$arcs[, 7:8]
    arcsah <- ah$arcs
    aa <- as.numeric(t(a))
    nuevoa <- aa
    miropos <- 2
    newpos <- 3
    fila <- 1
    while (miropos < length(aa)) {
        posigual <- which(aa[miropos] == aa)
        posigual <- posigual[posigual != miropos]
        if (posigual%/%2 == posigual/2) {
            nuevoa[newpos:(newpos + 1)] <- aa[posigual:(posigual - 
                1)]
            nuevoa[posigual:(posigual - 1)] <- aa[newpos:(newpos + 
                1)]
            filan <- posigual/2
            filold <- newpos%/%2 + 1
            filan2 <- arcsah[filan, ]
            filold2 <- arcsah[filold, ]
            arcsah[filan, ] <- filold2
            arcsah[filold, ] <- filan2
        }
        else {
            nuevoa[newpos:(newpos + 1)] <- aa[posigual:(posigual + 
                1)]
            nuevoa[posigual:(posigual + 1)] <- aa[newpos:(newpos + 
                1)]
            filan <- posigual%/%2 + 1
            filold <- newpos%/%2 + 1
            filan2 <- arcsah[filan, ]
            filold2 <- arcsah[filold, ]
            arcsah[filan, ] <- filold2
            arcsah[filold, ] <- filan2
        }
        aa <- nuevoa
        miropos <- miropos + 2
        newpos <- newpos + 2
    }
    arcsah <- cbind(arcsah[, 1:6], matrix(nuevoa, ncol = 2, byrow = T))
    x$arcs <- arcsah
    witharea <- x$arcs[, 3] > 0
    if (sum(witharea) > 0) {
        ind <- x$arcs[witharea, 7:8]
        arcs <- x$arcs[witharea, 1:6]
        row <- 1
        ncomp <- 0
        npoly <- numeric()
        while (row <= dim(ind)[1]) {
            check1 <- ind[row, 1]
            rownew <- match(check1, ind[, 2])
            ncomp <- ncomp + 1
            npoly[row:rownew] <- ncomp
            row <- rownew + 1
        }
        hole <- rep(-1, ncomp)
        area <- rep(0, ncomp)
        for (i in 1:ncomp) {
            compind <- (npoly == i)
            polyind <- unique(as.numeric(ind[compind, ]))
            polypoints <- x$xahull[polyind, ]
            check <- which(arcs[compind, 6] < pi * 0.5)[1]
            if (in.polygon(arcs[compind, 1][check], arcs[compind, 
                2][check], polypoints[, 1], polypoints[, 2])) {
                hole[i] <- 1
            }
            areapoly <- areapl(polypoints)
            areacirc <- arcs[compind, 3]^2 * 0.5 * (2 * arcs[compind, 
                6] - sin(2 * arcs[compind, 6]))
            area[i] <- areapoly + hole[i] * sum(areacirc)
        }
        totalarea <- sum(area[hole == -1]) - sum(area[hole == 
            1])
    }
    else {
        totalarea <- 0
    }
    return(area = totalarea)
}
