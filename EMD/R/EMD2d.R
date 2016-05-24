
extrema2dC <- function(z, nnrow=nrow(z), nncol=ncol(z)) {
    nnrow <- as.integer(nnrow)
    nncol <- as.integer(nncol)
    z <- as.double(z)
    maxindex <- integer(nnrow*nncol)
    nmax <- as.integer(0)
    totalmax <- as.integer(0)    
    minindex <- integer(nnrow*nncol)   
    nmin <- as.integer(0)   
    totalmin <- as.integer(0)      
    thetime <- unix.time(
    out <- .C("extrema2dC",
        z=z,
        nnrow=nnrow,
        nncol=nncol,
        maxindex=maxindex,
        nmax=nmax,
        totalmax=totalmax,
        minindex=minindex,
        nmin=nmin,
        totalmin=totalmin,
        #vertex=integer(nnrow*nncol),
        PACKAGE="EMD"))
    out$time <- thetime
    out$maxindex <- out$maxindex + 1
    out$minindex <- out$minindex + 1
    out$maxindex <- out$maxindex[1:out$totalmax]
    out$minindex <- out$minindex[1:out$totalmin]
    #out$vertex <- out$vertex + 1
        
    out$z <- NULL
    out$nnrow <- NULL
    out$nncol <- NULL
    out$totalmax <- NULL
    out$totalmin <- NULL
       
    return(out)
}


extrema2dVC <- function(z, nnrow=nrow(z), nncol=ncol(z)) {
    nnrow <- as.integer(nnrow)
    nncol <- as.integer(nncol)
    z <- as.double(z)
    maxindex <- integer(nnrow*nncol)
    nmax <- as.integer(0)
    totalmax <- as.integer(0)    
    minindex <- integer(nnrow*nncol)   
    nmin <- as.integer(0)   
    totalmin <- as.integer(0)      
    thetime <- unix.time(
    out <- .C("extrema2dVC",
        z=z,
        nnrow=nnrow,
        nncol=nncol,
        maxindex=maxindex,
        nmax=nmax,
        totalmax=totalmax,
        minindex=minindex,
        nmin=nmin,
        totalmin=totalmin,
        vertex=integer(nnrow*nncol),
        PACKAGE="EMD"))
    out$time <- thetime
    out$maxindex <- out$maxindex + 1
    out$minindex <- out$minindex + 1
    out$maxindex <- out$maxindex[1:out$totalmax]
    out$minindex <- out$minindex[1:out$totalmin]
    out$vertex <- out$vertex + 1
        
    out$z <- NULL
    out$nnrow <- NULL
    out$nncol <- NULL
    out$totalmax <- NULL
    out$totalmin <- NULL
       
    return(out)
}

extractimf2d <-
function (residue, x = NULL, y = NULL, nnrow = nrow(residue), 
    nncol = ncol(residue), tol = sd(c(residue)) * 0.1^2, max.sift = 20, 
    boundary = "reflexive", boundperc = 0.3, sm="none", spar = NULL, weight = NULL, check = FALSE)
{
    if (sm != "none" & (is.null(spar) || spar == 0)) stop("Provide the smoothing parameter.\n")
    if (sm == "none") spar <- 0
    
    if (check) 
        par(mfrow = c(2, 2), oma = c(0, 0, 0, 0), mar = 0.1 + c(1, 0.5, 1.5, 0.5))
    tolN <- 9
    if (is.null(x)) 
        x <- 1:nnrow
    if (is.null(y)) 
        y <- 1:nncol
    input <- residue

    #if (sm == "Tps") {
    #    xx <- NULL
    #    for (i in 1:nncol) xx <- rbind(xx, cbind(x, y[i]))
    #    dimnames(xx) <- list(NULL, NULL)
    #    tmp <- Tps(xx, c(input))
    #    input <- input - matrix(predict(tmp), ncol=nncol)
    #}
    
    if (boundary == "periodic" || boundary == "reflexive") {
        rowext <- round(nnrow * boundperc)
        colext <- round(nncol * boundperc)
        x <- c(sort(seq(x[1], length = rowext + 1, by = x[1] - 
            x[2])[-1]), x, seq(x[nnrow], length = rowext + 1, 
            by = x[2] - x[1])[-1])
        y <- c(sort(seq(y[1], length = colext + 1, by = y[1] - 
            y[2])[-1]), y, seq(y[nncol], length = colext + 1, 
            by = y[2] - y[1])[-1])
        nnrowext <- 2 * rowext + nnrow
        nncolext <- 2 * colext + nncol
        if (boundary == "periodic") {
            input <- cbind(input[, (nncol - colext + 1):nncol], 
                input, input[, 1:colext])
            input <- rbind(input[(nnrow - rowext + 1):nnrow, 
                ], input, input[1:rowext, ])
        }
        if (boundary == "reflexive") {
            input <- cbind(input[, colext:1], input, input[, 
                nncol:(nncol - colext + 1)])
            input <- rbind(input[rowext:1, ], input, input[nnrow:(nnrow - 
                rowext + 1), ])
        }
    }
    else if (boundary == "none") {
        nnrowext <- nnrow
        nncolext <- nncol
    }
    
    xx <- NULL
    for (i in 1:nncolext) xx <- rbind(xx, cbind(x, y[i]))
    dimnames(xx) <- list(NULL, NULL)
    imf <- maxindex <- minindex <- NULL
    
    j <- 1
    repeat {
        tmp <- extrema2dC(input, nnrowext, nncolext)
        if (tmp$nmin <= tolN || tmp$nmax <= tolN) 
            break
        minindex <- unlist(tmp$minindex)
        maxindex <- unlist(tmp$maxindex)
        
        if (sm == "Tps") {
        fmin <- Tps(xx[minindex, ], input[minindex])
        fmax <- Tps(xx[maxindex, ], input[maxindex])
        emin <- predictSurface(fmin, grid.list = list(x, y), extrap = TRUE, lambda = spar)$z
        emax <- predictSurface(fmax, grid.list = list(x, y), extrap = TRUE, lambda = spar)$z
        } else if (sm == "mKrig" || sm == "none") {
        fmin <- mKrig(xx[minindex, ], input[minindex], lambda=spar)
        fmax <- mKrig(xx[maxindex, ], input[maxindex], lambda=spar)        
        emin <- predictSurface(fmin, grid.list = list(x, y), extrap = TRUE)$z
        emax <- predictSurface(fmax, grid.list = list(x, y), extrap = TRUE)$z
        } else if (sm == "locfit") {
        fmin <- locfit(input[minindex] ~ xx[minindex, ], deg=3, kern="gauss", alpha=spar)
        fmax <- locfit(input[maxindex] ~ xx[maxindex, ], deg=3, kern="gauss", alpha=spar)      
        emin <- predict(fmin, xx)
        emax <- predict(fmax, xx)
        } else if (sm == "loess") {
        fmin <- loess(input[minindex] ~ xx[minindex, ], span=spar, degree=2, control = loess.control(surface = "direct"))
        fmax <- loess(input[maxindex] ~ xx[maxindex, ], span=spar, degree=2, control = loess.control(surface = "direct"))
        emin <- predict(fmin, xx)
        emax <- predict(fmax, xx)
        }
        
        #if (sm == "none") {
        #    emin <- predictSurface(fmin, grid.list = list(x, y), extrap = TRUE, lambda = spar)$z
        #    emax <- predictSurface(fmax, grid.list = list(x, y), extrap = TRUE, lambda = spar)$z
        #} else if (sm == "Tps") {
        #    emin <- predictSurface(fmin, grid.list = list(x, y), extrap = TRUE, lambda = spar)$z
        #    emax <- predictSurface(fmax, grid.list = list(x, y), extrap = TRUE, lambda = spar)$z
        #}
        em <- (emin + emax)/2
        if (check) {
            rangez <- range(c(input, emin, emax, em))
            image(x, y, input, main = paste(j, "-th input", sep = ""), 
                xlab = "", ylab = "", zlim = rangez, col = gray(100:0/100), 
                axes = FALSE)
            box()
            image(x, y, emin, main = "lower", xlab = "", ylab = "", 
                zlim = rangez, col = gray(100:0/100), axes = FALSE)
            box()
            image(x, y, emax, main = "upper", xlab = "", ylab = "", 
                zlim = rangez, col = gray(100:0/100), axes = FALSE)
            box()
            image(x, y, em, main = "mean", xlab = "", ylab = "", 
                zlim = rangez, col = gray(100:0/100), axes = FALSE)
            box()
            cat("\a")
            locator(1)
        }
        if (all(abs(em) < tol) || j >= max.sift) {
            imf <- input - em
            if (boundary == "periodic" || boundary == "reflexive") {
                imf <- imf[(rowext + 1):(rowext + nnrow), (colext + 
                  1):(colext + nncol)]
                tmpcol <- ceiling(maxindex/nnrowext)
                tmprow <- maxindex - nnrowext * (tmpcol - 1)
                tmp2 <- (tmprow > rowext & tmprow <= rowext + 
                  nnrow & tmpcol > colext & tmpcol <= colext + 
                  nncol)
                tmprow <- tmprow[tmp2] - rowext
                tmpcol <- tmpcol[tmp2] - colext
                maxindex <- tmprow + nnrow * (tmpcol - 1)
                tmpcol <- ceiling(minindex/nnrowext)
                tmprow <- minindex - nnrowext * (tmpcol - 1)
                tmp2 <- (tmprow > rowext & tmprow <= rowext + 
                  nnrow & tmpcol > colext & tmpcol <= colext + 
                  nncol)
                tmprow <- tmprow[tmp2] - rowext
                tmpcol <- tmpcol[tmp2] - colext
                minindex <- tmprow + nnrow * (tmpcol - 1)
            }
            residue <- residue - imf
            break
        }
        input <- input - em
        j <- j + 1
    }
    list(imf = imf, residue = residue, maxindex = maxindex, minindex = minindex, niter = j)
}

emd2d <-
function (z, x = NULL, y = NULL, tol = sd(c(z)) * 0.1^2, max.sift = 20, 
    boundary = "reflexive", boundperc = 0.3, max.imf = 5, sm = "none", smlevels = 1, spar=NULL,
    weight = NULL, plot.imf = FALSE) 
{
    if (sm != "none") {
        if (is.null(spar)) stop("Provide the smoothing parameter.\n")
        else if (length(spar) == 1)
            spar <- rep(spar, length(smlevels))  
        else if (length(smlevels) != length(spar))
            stop("Provide the smoothing parameter.")
    }
  
    nnrow <- nrow(z)
    nncol <- ncol(z)
    if (is.null(x)) 
        x <- 1:nnrow
    if (is.null(y)) 
        y <- 1:nncol
    if (plot.imf) 
        par(mfrow = c(3, 1), oma = c(0, 0, 0, 0), mar = 0.1 + c(1, 0.5, 1.5, 0.5))
    residue <- z
    imf <- maxindex <- minindex <- as.list(NULL)
    j <- 1
    repeat {
        if (j > max.imf) 
            break
        if ((any(j == smlevels)) & sm != "none")
            tmp <- extractimf2d(residue, x, y, nnrow, nncol, tol, max.sift, boundary, boundperc, sm, spar, NULL, check = FALSE)
        else tmp <- extractimf2d(residue, x, y, nnrow, nncol, tol, max.sift, boundary, boundperc, "none", 0, NULL, check = FALSE)
        if (is.null(tmp$imf)) 
            break
        if (plot.imf) {
            rangez <- range(c(residue, tmp$imf, tmp$residue))
            image(x, y, residue, main = paste(j - 1, "-th residue", 
                sep = ""), xlab = "", ylab = "", zlim = rangez, 
                col = gray(100:0/100), axes = FALSE)
            box()
            image(x, y, tmp$imf, main = paste(j, "-th imf", sep = ""), 
                xlab = "", ylab = "", zlim = rangez, col = gray(100:0/100), 
                axes = FALSE)
            box()
            image(x, y, tmp$residue, main = paste(j, "-th residue", 
                sep = ""), xlab = "", ylab = "", zlim = rangez, 
                col = gray(100:0/100), axes = FALSE)
            box()
            cat("\a")
            locator(1)
        }
        imf <- c(imf, list(tmp$imf))
        maxindex <- c(maxindex, list(tmp$maxindex))
        minindex <- c(minindex, list(tmp$minindex))
        residue <- tmp$residue
        j <- j + 1
    }
    list(imf = imf, residue = residue, maxindex = maxindex, minindex = minindex, 
        nimf = j - 1)
}


imageEMD <- function(z=z, emdz, extrema=FALSE, ...) {

    if(extrema)
        par(mfrow=c(emdz$nimf+1, 2), oma=c(0,0,0,0), mar=0.1+c(1,0.5,1.5,0.5))
    if(!extrema) 
        par(mfrow=c(emdz$nimf+2, 1), oma=c(0,0,0,0), mar=0.1+c(1,0.5,1.5,0.5)) 

    nnrow <- nrow(z); nncol <- ncol(z); ndata <- nncol * nnrow        
    rangez <- range(c(z, emdz$imf, emdz$residue))
    image(z, main="Image", xlab="", ylab="", zlim=rangez, axes=FALSE, ...)
    if(extrema)
        image(emdz$residue, main="residue", xlab="", ylab="", zlim=rangez, axes=FALSE, ...)
    for(i in 1:emdz$nimf) {
        image(emdz$imf[[i]], main=paste("imf ", i, sep=""), xlab="", ylab="", zlim=rangez,
                    col=gray(100:0/100), axes=FALSE)
        if(extrema) {
            tmp <- matrix(rep(50, ndata), ncol=nncol)
            tmp[emdz$maxindex[[i]]] <- 100
            tmp[emdz$minindex[[i]]] <- 0
            image(tmp, main=paste("extrema ", i, sep=""), xlab="", ylab="", col=gray(0:100/100), axes=FALSE)
        }
    }
    if(!extrema)
        image(emdz$residue, main="residue", xlab="", ylab="", zlim=rangez, col=gray(100:0/100), axes=FALSE)
}
