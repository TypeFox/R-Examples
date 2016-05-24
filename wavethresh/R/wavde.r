"Chires5" <-
function(x, tau=1, J, filter.number=10, family="DaubLeAsymm", nT=20)

# data x
# fine tuning parameter tau
# resolution level J
# family and filter.number specify the scaling function to be used
# nT is the number of iterations performed in the Daubechies-Lagarias algorithm

{

# calculate support of father wavelet

    sup <- support(filter.number, family)
    sup <- c(sup$phi.lh, sup$phi.rh)

# extract filter coefficients

    filcf <- filter.select(filter.number, family)$H

# calculate primary resolution

    p <- tau * 2^J

# calculate bounds on translation

    kmin <- ceiling(p*min(x)-sup[2])
    kmax <- floor(p*max(x)-sup[1])

# create vector to put estimated coefficients in

    chat <- rep(0, kmax-kmin+1)

# call C code!

    error <- 0

    ans <- .C("SFDE5",
        x = as.double(x),
        nx = as.integer(length(x)),
        p = as.double(p),
        filter = as.double(filcf),
        nf = as.integer(2*filter.number - 1),
        prec = as.integer(nT),
        chat = as.double(chat),
        kmin = as.integer(kmin),
        kmax = as.integer(kmax),
        philh = as.double(sup[1]),
        phirh = as.double(sup[2]),
	error = as.integer(error), PACKAGE = "wavethresh")

    if (ans$error != 0)
	stop(paste("PLDF2 function returned error code:", ans$error))


    filter <- list(filter.number=filter.number, family=family)
    res <- list(p=p, tau=tau, J=J)
    list(coef=ans$chat, klim=c(ans$kmin, ans$kmax), p=ans$p, filter=filter,
        n=length(x), res=res)
}

"Chires6" <-
function(x, tau=1, J, filter.number=10, family="DaubLeAsymm", nT=20)

# data x
# fine tuning parameter tau
# resolution level J
# family and filter.number specify the scaling function to be used
# nT is the number of iterations performed in the Daubechies-Lagarias algorithm

{

# calculate support of father wavelet

    sup <- support(filter.number, family)
    sup <- c(sup$phi.lh, sup$phi.rh)

# extract filter coefficients

    filcf <- filter.select(filter.number, family)$H

# calculate primary resolution

    p <- tau * 2^J

# calculate bounds on translation

    kmin <- ceiling(p*min(x)-sup[2])
    kmax <- floor(p*max(x)-sup[1])

# create output vector/matrix

    ncoef <- kmax-kmin+1
    chat <- rep(0, ncoef)
    covar <- matrix(0, nrow=ncoef, ncol=(2*filter.number-1))

# call C code!

    error <- 0

    ans <- .C("SFDE6",
        x = as.double(x),
        nx = as.integer(length(x)),
        p = as.double(p),
        filter = as.double(filcf),
        nf = as.integer(2*filter.number - 1),
        prec = as.integer(nT),
        chat = as.double(chat),
        covar = as.double(covar),
        kmin = as.integer(kmin),
        kmax = as.integer(kmax),
        philh = as.double(sup[1]),
        phirh = as.double(sup[2]),
	error = as.integer(error), PACKAGE = "wavethresh")

    if (ans$error != 0)
	stop(paste("PLDF2 function returned error code:", ans$error))


    filter <- list(filter.number=filter.number, family=family)
    res <- list(p=p, tau=tau, J=J)
    list(coef=ans$chat, covar=matrix(ans$covar, nrow=ncoef),
        klim=c(ans$kmin, ans$kmax), p=ans$p, filter=filter,
        n=length(x), res=res)
}

"dclaw" <-
function(x)

{
    den <- dnorm(x)/2
    for(i in 0:4){
        den <- den + dnorm(x, mean=(i/2-1), sd=1/10)/10
    }

    den
}

"dencvwd" <-
function(hrproj, filter.number=hrproj$filter$filter.number,
    family=hrproj$filter$family, type="wavelet", bc="zero", firstk=hrproj$klim,
    RetFather=TRUE, verbose=FALSE)

{
    image <- hrproj$covar

# Select wavelet filter

    filter <- filter.select(filter.number = filter.number, family = family)
    Csize <- nrow(image)

# Set-up first/last database

    if(is.null(firstk))
        firstk <- c(0, Csize-1)
    if(verbose == TRUE)
        cat("...selected\nFirst/last database...")
    fl.dbase <- first.last.dh(LengthH = length(filter$H), DataLength = Csize, 
        bc = bc, type = type, firstk = firstk)
    first.last.c <- fl.dbase$first.last.c
    first.last.d <- fl.dbase$first.last.d

    nlev <- nrow(first.last.d)

# Set up answer list

    image.decomp <- list(nlevels = nlev, fl.dbase = fl.dbase, filter = filter,
        type = type, bc = bc, date = date())

    if(verbose == TRUE) cat("...built\n")

# Ok, go round loop doing decompositions

    nbc <- switch(bc,
        periodic = 1,
        symmetric = 2,
        zero = 3)
    if(is.null(nbc))
        stop("Unknown boundary handling")
    if(type == "station" && bc == "symmetric")
        stop("Cannot do stationary transform with symmetric boundary conditions"
            )
    ntype <- switch(type,
        wavelet = 1,
        station = 2)
    if(is.null(ntype)) stop("Unknown type of transform")

# Load up original image

    smoothed <- as.vector(image)
    if(verbose == TRUE) {
        cat(bc, " boundary handling\n")
        cat("Decomposing...")
    }
    for(level in seq(nrow(first.last.d), 1, -1)) {
        if(verbose == TRUE)
            cat(level - 1, "")
        LengthCin <- first.last.c[level+1, 2] - first.last.c[level+1, 1] + 1
        LengthCout <- first.last.c[level, 2] - first.last.c[level, 1] + 1
        LengthDout <- first.last.d[level, 2] - first.last.d[level, 1] + 1
        ImCC <- rep(0, (LengthCout * (2*filter.number-1)))
        ImDD <- rep(0, (LengthDout * (2*filter.number-1)))
        error <- 0
        z <- .C("StoDCDS",
            C = as.double(smoothed),
            Csize = as.integer(LengthCin),
            firstCin = as.integer(first.last.c[level + 1, 1]),
            H = as.double(filter$H),
            LengthH = as.integer(length(filter$H)),
            LengthCout = as.integer(LengthCout),
            firstCout = as.integer(first.last.c[level, 1]),
            lastCout = as.integer(first.last.c[level, 2]),
            LengthDout = as.integer(LengthDout),
            firstDout = as.integer(first.last.d[level, 1]),
            lastDout = as.integer(first.last.d[level, 2]),
            ImCC = as.double(ImCC),
            ImDD = as.double(ImDD),
            nbc = as.integer(nbc),
            ntype = as.integer(ntype),
            error = as.integer(error), PACKAGE = "wavethresh")
        error <- z$error
        if(error != 0) {
            cat("Error was ", error, "\n")
            stop("Error reported")
        }
        smoothed <- z$ImCC
        if(RetFather == TRUE) {
            nm <- lt.to.name(level - 1, "CC")
            image.decomp[[nm]] <- matrix(z$ImCC, nrow=LengthCout)
        }
        nm <- lt.to.name(level - 1, "DD")
        image.decomp[[nm]] <- matrix(z$ImDD, nrow=LengthDout)
    }
    if(verbose == TRUE)
        cat("\nReturning answer...\n")
    image.decomp$w0Lconstant <- smoothed
    image.decomp$bc <- bc
    image.decomp$date <- date()
    class(image.decomp) <- "imwd"

    l <- list(C=NULL, D=rep(0, fl.dbase$ntotal.d),
                nlevels=nrow(fl.dbase$first.last.d), fl.dbase=fl.dbase,
                filter=filter, type=type, bc=bc, date=date())
    class(l) <- "wd"
    for(level in 1:nlevelsWT(l)) {
        covar <- image.decomp[[lt.to.name(level - 1, "DD")]]
        l <- putD.wd(l, level-1, covar[,1], boundary=TRUE)
    }
 
    l

}

"denplot" <-
function(wr, coef, nT=20, lims, n=50)

# smoothed high level coefficients wr
# coef is the output from denproj for this analysis
# nT is the number of iterations performed in the Daubechies-Lagarias algorithm
# estimate is evaluated at n points between lims

{
    p <- coef$res$p
    filter <- coef$filter

# calculate support of father wavelet

    sup <- support(filter$filter.number, filter$family)
    sup <- c(sup$phi.lh, sup$phi.rh)

# extract filter coefficients

    filcf <- filter.select(filter$filter.number, filter$family)$H

# create grid for drawing density estimate and vector to put values in

    gx <- seq(lims[1], lims[2], length=n)
    gy <- c(rep(0, length(gx)))

# find range of high resolution coefficients

    kmin <- coef$klim[1]
    kmax <- coef$klim[2]

# call C code!

    error <- 0
    ans <- .C("PLDE2",
        C = as.double(wr),
        p = as.double(p),
        filter = as.double(filcf),
        nf = as.integer(2*filter$filter.number - 1),
        prec = as.integer(nT),
        kmin = as.integer(kmin),
        kmax = as.integer(kmax),
        gx = as.double(gx),
        gy = as.double(gy),
        ng = as.integer(n),
        philh = as.double(sup[1]),
        phirh = as.double(sup[2]),
	error = as.integer(error), PACKAGE = "wavethresh")

    if (ans$error != 0)
	stop(paste("PLDF2 function returned error code:", ans$error))


    list(x=ans$gx, y=ans$gy)
}

"denproj" <-
function(x, tau=1, J, filter.number=10, family="DaubLeAsymm",
    covar=FALSE, nT=20)

# data x
# fine tuning parameter tau
# resolution level J
# family and filter.number specify the scaling function to be used
# covar - logical variable indicating whether covariances should be calculated
# nT is the number of iterations performed in the Daubechies-Lagarias algorithm

{
    if(covar)
        ans <- Chires6(x, tau, J, filter.number, family, nT)
    else
        ans <- Chires5(x, tau, J, filter.number, family, nT)

    ans
}

"denwd" <-
function(coef)

{
    wd.dh(coef$coef, filter.number=coef$filter$filter.number,
        family=coef$filter$family, bc="zero", firstk=coef$klim)

}

"denwr" <-
function(wd, start.level=0, verbose=FALSE, bc=wd$bc, return.object=FALSE,
    filter.number=wd$filter$filter.number, family=wd$filter$family)
{
    if(IsEarly(wd)) {
        ConvertMessage()
        stop()
    }
    if(verbose == TRUE)
        cat("Argument checking...")

# Check class of wd

    if(verbose == TRUE)
        cat("Argument checking\n")
    ctmp <- class(wd)
    if(is.null(ctmp))
        stop("wd has no class")
    else if(ctmp != "wd")
        stop("wd is not of class wd")
    if(start.level < 0)
        stop("start.level must be nonnegative")
    if(start.level >= nlevelsWT(wd))
        stop("start.level must be less than the number of levels")
    if(is.null(wd$filter$filter.number))
        stop("NULL filter.number for wd")
    if(bc != wd$bc)
        warning("Boundary handling is different to original")
    if(wd$type == "station")
        stop("Use convert to generate wst object and then AvBasis or InvBasis"
            )
    type <- wd$type
    filter <- filter.select(filter.number = filter.number, family = family)
    LengthH <- length(filter$H) 

# Build the reconstruction first/last database

    if(verbose == TRUE)
        cat("...done\nFirst/last database...")
    r.first.last.c <- wd$fl.dbase$first.last.c[(start.level+1):(nlevelsWT(wd)+1), ]
    ntotal <- r.first.last.c[1,3] + r.first.last.c[1,2] -
                                               r.first.last.c[1,1] + 1
    names(ntotal) <- NULL
    C <- accessC(wd, level = start.level, boundary = TRUE)
    C <- c(rep(0, length = (ntotal - length(C))), C)
    nlevels <- nlevelsWT(wd) - start.level
    error <- 0

# Load object code

    if(verbose == TRUE)
        cat("...built\n")
    if(verbose == TRUE) {
        cat("Reconstruction...")
        error <- 1
    }
    ntype <- switch(type,
        wavelet = 1,
        station = 2)
    if(is.null(ntype))
        stop("Unknown type of decomposition")
    nbc <- switch(bc,
        periodic = 1,
        symmetric = 2,
        zero = 3)
    if(is.null(nbc))
        stop("Unknown boundary handling")
    if(!is.complex(wd$D)) {
        wavelet.reconstruction <- .C("waverecons_dh",
            C = as.double(C),
            D = as.double(wd$D),
            H = as.double(filter$H),
            LengthH = as.integer(LengthH),
            nlevels = as.integer(nlevels),
            firstC = as.integer(r.first.last.c[, 1]),
            lastC = as.integer(r.first.last.c[, 2]),
            offsetC = as.integer(r.first.last.c[, 3]),
            firstD = as.integer(wd$fl.dbase$first.last.d[, 1]),
            lastD = as.integer(wd$fl.dbase$first.last.d[, 2]),
            offsetD = as.integer(wd$fl.dbase$first.last.d[, 3]),
            ntype = as.integer(ntype),
            nbc = as.integer(nbc),
            error = as.integer(error), PACKAGE = "wavethresh")
    }
    if(verbose == TRUE)
        cat("done\n")
    error <- wavelet.reconstruction$error
    if(error != 0) {
        cat("Error code returned from waverecons: ", error, "\n")
        stop("waverecons returned error")
    }
    fl.dbase <- list(first.last.c=r.first.last.c,
        ntotal=wavelet.reconstruction$LengthC,
        first.last.d=wd$fl.dbase$ first.last.d, ntotal.d=wd$fl.dbase$ntotal.d)
    if(!is.complex(wd$D)) {
        l <- list(C=wavelet.reconstruction$C, D=wavelet.reconstruction$D,
                    fl.dbase=fl.dbase, nlevels=nlevelsWT(wavelet.reconstruction),
                    filter=filter, type=type, bc=bc, date=date())
    }
    class(l) <- "wd"
    if(return.object == TRUE)
        return(l)
    else {
        if(bc == "zero")
            return(accessC(l, boundary = TRUE))
        else return(accessC(l))
    }
    stop("Shouldn't get here\n")
}

"first.last.dh" <-
function(LengthH, DataLength, type = "wavelet",
    bc = "periodic", firstk=c(0, DataLength-1))
{
    if(type == "station" && bc != "periodic")
        stop("Can only do periodic boundary conditions with station")
    if(type != "station" && type != "wavelet")
        stop("Type can only be wavelet or station")

    if(bc=="periodic" || bc=="symmetric") {
        levels <- log(DataLength)/log(2)
        first.last.c <- matrix(0, nrow = levels + 1, ncol = 3,
            dimnames = list(NULL, c("First", "Last", "Offset")))
        first.last.d <- matrix(0, nrow = levels, ncol = 3,
            dimnames = list(NULL, c("First", "Last", "Offset")))
    }

    if(bc == "periodic") {
# Periodic boundary correction
        if(type == "wavelet") {
            first.last.c[, 1] <- rep(0, levels + 1)
            first.last.c[, 2] <- 2^(0:levels) - 1
            first.last.c[, 3] <- rev(c(0, cumsum(rev(1 +
                first.last.c[, 2]))[1:levels]))
            first.last.d[, 1] <- rep(0, levels)
            first.last.d[, 2] <- 2^(0:(levels - 1)) - 1
            first.last.d[, 3] <- rev(c(0, cumsum(rev(1 +
                first.last.d[, 2]))[1:(levels - 1)]))
            ntotal <- 2 * DataLength - 1
            ntotal.d <- DataLength - 1
        }
        else if(type == "station") {
            first.last.c[, 1] <- rep(0, levels + 1)
            first.last.c[, 2] <- 2^levels - 1
            first.last.c[, 3] <- rev(c(0, cumsum(rev(1 +
                first.last.c[, 2]))[1:levels]))
            first.last.d[, 1] <- rep(0, levels)
            first.last.d[, 2] <- 2^levels - 1
            first.last.d[, 3] <- rev(c(0, cumsum(rev(1 +
                first.last.d[, 2]))[1:(levels - 1)]))
            ntotal <- (levels + 1) * 2^levels
            ntotal.d <- levels * 2^levels
        }
    }

    else if(bc == "symmetric") {
# Symmetric boundary reflection
        first.last.c[levels + 1, 1] <- 0
        first.last.c[levels + 1, 2] <- DataLength - 1
        first.last.c[levels + 1, 3] <- 0
        ntotal <- first.last.c[levels + 1, 2] - first.last.c[levels + 1,1] + 1
        ntotal.d <- 0
        for(i in levels:1) {
            first.last.c[i, 1] <- trunc(0.5 * (1 - LengthH +
                first.last.c[i + 1, 1]))
            first.last.c[i, 2] <- trunc(0.5 * first.last.c[i + 1, 2])
            first.last.c[i, 3] <- first.last.c[i + 1, 3] +
                first.last.c[i + 1, 2] - first.last.c[i + 1, 1] + 1
            first.last.d[i, 1] <- trunc(0.5 * (first.last.c[i + 1, 1] - 1))
            first.last.d[i, 2] <- trunc(0.5 * (first.last.c[i + 1,
                2] + LengthH - 2))
            if(i != levels) {
                first.last.d[i, 3] <- first.last.d[i + 1, 3] +
                    first.last.d[i + 1, 2] - first.last.d[i + 1, 1] + 1
            }
            ntotal <- ntotal + first.last.c[i, 2] - first.last.c[i,1] + 1
            ntotal.d <- ntotal.d + first.last.d[i, 2] -  first.last.d[i, 1] + 1
        }
    }

    else if(bc=="zero") {
        first.c <- firstk[1]
        last.c <- firstk[2]
        offset.c <- 0
        first.d <- NULL
        last.d <- NULL
        offset.d <- 0
        ntotal <- last.c - first.c + 1
        ntotal.d <- 0

        while( (first.c[1] > 2 - LengthH || first.c[1] < 1 - LengthH) ||
                (last.c[1] > 0 || last.c[1] < -1) ) {
            first.c <- c(ceiling(0.5*(first.c[1] - LengthH + 1)), first.c)
            last.c <- c(floor(0.5*last.c[1]), last.c)
            offset.c <- c(offset.c[1] + last.c[2] - first.c[2] +1, offset.c)
            ntotal <- ntotal + last.c[1] - first.c[1] + 1

            first.d <- c(ceiling(0.5*(first.c[2]-1)), first.d)
            last.d <- c(floor(0.5*(last.c[2] + LengthH - 2)), last.d)
            if(length(first.d) > 1)
                offset.d <- c(offset.d[1] + last.d[2] - 
                    first.d[2] + 1, offset.d)
            ntotal.d <- ntotal.d + last.d[1] - first.d[1] +1
            
        }
        
        first.last.c <- matrix(c(first.c, last.c, offset.c), ncol=3, 
            dimnames=list(NULL, c("First", "Last", "Offset")))
        first.last.d <- matrix(c(first.d, last.d, offset.d), ncol=3,
            dimnames=list(NULL, c("First", "Last", "Offset")))
    }
        
    else {
        stop("Unknown boundary correction method")
    }
    names(ntotal) <- NULL
    names(ntotal.d) <- NULL
    list(first.last.c = first.last.c, ntotal = ntotal, first.last.d =
        first.last.d, ntotal.d = ntotal.d)
}

"pclaw" <-
function(q)

{
    prob <- pnorm(q)/2
    for(i in 0:4){
        prob <- prob + pnorm(q, mean=(i/2-1), sd=1/10)/10
    }

    prob
}

"plotdenwd" <-
function(wd, xlabvals, xlabchars, ylabchars, first.level=0,
    top.level=nlevelsWT(wd) - 1,
    main="Wavelet Decomposition Coefficients", scaling="global", rhlab=FALSE,
    sub, NotPlotVal=0.005, xlab="Translate",
    ylab="Resolution Level", aspect="Identity", ...)

{
    ctmp <- class(wd)
    if(is.null(ctmp))
        stop("wd has no class")
    else if(ctmp != "wd")
        stop("wd is not of class wd")

    levels <- nlevelsWT(wd)
    nlevels <- levels - first.level
    cfac <- top.level - (levels-1)

    sfac <- rep(2, nlevels) ^ c((nlevels-1):0)
    first <- wd$fl.dbase$first.last.d[(first.level+1):levels,1]
    first <- first * sfac + (sfac-1)/2
    last <- wd$fl.dbase$first.last.d[(first.level+1):levels,2]
    last <- last * sfac + (sfac-1)/2
    xrange <- c(floor(min(first)), ceiling(max(last)))
    type <- wd$type
    if(type == "wavelet")
         n <- 2^(levels-2)

    if(missing(sub))
        sub <- paste(switch(type,
            wavelet = "Standard transform",
            station = "Nondecimated transform"), wd$filter$name)
    if(aspect != "Identity")
        sub <- paste(sub, "(", aspect, ")")
    plot(c(xrange[1], xrange[1], xrange[2], xrange[2]),
        c(0, nlevels+1, nlevels+1, 0), type="n", xlab=xlab,
        ylab=ylab, main=main, yaxt="n", xaxt="n", sub=sub, ...)

    yll <- top.level:(first.level+cfac)
    if(missing(ylabchars))
        axis(2, at = 1:(nlevels), labels = yll)
    else if(length(ylabchars) != nlevels)
        stop(paste("Should have ", nlevels, " entries in ylabchars"))
    else axis(2, at = 1:(nlevels), labels = ylabchars)

    if(missing(xlabchars)) {
        if(missing(xlabvals)) {
            if(type == "wavelet") {
                if(wd$bc != "zero") {
                    axx <- c(0, 2^(levels - 3), 2^(levels - 2),
                        2^(levels - 2) + 2^(levels - 3), 2^(levels - 1))
                }
                else {
                    jrange <- floor(logb(abs(xrange), 2))
                    xlabr <- sign(xrange) * 2^jrange
                    xsp <- diff(xlabr)
                    axx <- xlabr[1] + c(0, xsp/4, xsp/2, 3*xsp/4, xsp)
                    if((xlabr[2]+xsp/4) <= xrange[2])
                        axx <- c(axx, xlabr[2]+xsp/4)
                    if((xlabr[1]-xsp/4) >= xrange[1])
                        axx <- c(xlabr[1]-xsp/4, axx)
                }
            }
            else axx <- c(0, 2^(levels - 2), 2^(levels - 1),
                2^(levels - 1) + 2^(levels - 2), 2^levels)
            if(is.null(tsp(wd)))
                axis(1, at = axx)
            else {
                v <- seq(from=tsp(wd)["start"], by=tsp(wd)["deltat"],
                    length=n)
                if(type == "wavelet")
                    atl <- 2 * v
                else atl <- v
                atl <- pretty(atl, n = 4)
                ats <- (n * atl)/(max(atl) - min(atl))
                axis(1, at = ats, labels = atl)
            }
        }
        else {
            lx <- pretty(xlabvals, n = 4)
            cat("lx is ", lx, "\n")
            if(lx[1] < min(xlabvals))
                lx[1] <- min(xlabvals)
            if(lx[length(lx)] > max(xlabvals))
                lx[length(lx)] <- max(xlabvals)
            cat("lx is ", lx, "\n")
            xix <- NULL
            for(i in 1:length(lx)) {
                u <- (xlabvals - lx[i])^2
                xix <- c(xix, (1:length(u))[u == min(u)])
            }
            axx <- xix
            if(type == "wavelet")
                axx <- xix/2
            axl <- signif(lx, digits = 2)
            axis(1, at = axx, labels = axl)
        }
    }
    else axis(1, at = xlabvals, labels = xlabchars)



    x <- 1:n
    height <- 1
    first.last.d <- wd$fl.dbase$first.last.d
    axr <- NULL
    if(scaling == "global") {
        my <- 0
        for(i in ((levels - 1):first.level)) {
            y <- accessD(wd, i, boundary=TRUE, aspect = aspect)
            my <- max(c(my, abs(y)))
        }
    }
    if(scaling == "compensated") {
        my <- 0
        for(i in ((levels - 1):first.level)) {
            y <- accessD(wd, i, boundary=TRUE, aspect = aspect) * 2^(i/2)
            my <- max(c(my, abs(y)))
        }
    }
    if(scaling == "super") {
        my <- 0
        for(i in ((levels - 1):first.level)) {
            y <- accessD(wd, i, boundary=TRUE, aspect = aspect) * 2^i
            my <- max(c(my, abs(y)))
        }
    }
    shift <- 1
    for(i in ((levels - 1):first.level)) {
        y <- accessD(wd, i, boundary=TRUE, aspect = aspect)
        if(type == "wavelet")
            n <- first.last.d[i+1,2]-first.last.d[i+1,1]+1
        else {
            y <- y[c((n - shift + 1):n, 1:(n - shift))]
            shift <- shift * 2
        }
        xplot <- seq(from=first[i-first.level+1], to=last[i-first.level+1], by=2^(nlevels-(i-first.level)-1))
        ly <- length(y)
        if(scaling == "by.level")
            my <- max(abs(y))
        if(scaling == "compensated")
            y <- y * 2^(i/2)
        if(scaling == "super")
            y <- y * 2^i
        if(my == 0) {
            y <- rep(0, length(y))
        }
        else y <- (0.5 * y)/my
        axr <- c(axr, my)
        if(max(abs(y)) > NotPlotVal)
            segments(xplot, height, xplot, height + y)
        if(i != first.level) {
            if(type == "wavelet") {
#                x1 <- x[seq(1, n - 1, 2)]
#                x2 <- x[seq(2, n, 2)]
#                x <- (x1 + x2)/2
#                x <- 1:n
            }
            height <- height + 1
        }
    }
    if(rhlab == TRUE)
        axis(4, at = 1:length(axr), labels = signif(axr, 3))
    axr
}

"rclaw" <-
function(n)

{
    nx <- rnorm(n)
    p <- runif(n)
    oldx <- nx
    nx[p<=0.5] <- nx[p<=0.5]/10 + (trunc(p[p<=0.5] * 10)/2 -1)
    nx
}

"wd.dh" <-
function(data, filter.number = 10, family = "DaubLeAsymm",
    type = "wavelet", bc = "periodic", firstk=NULL, verbose = FALSE)
{
    if(verbose == TRUE)
        cat("wd: Argument checking...")
    if(!is.atomic(data))
        stop("Data is not atomic")
    DataLength <- length(data)

# Check that we have a power of 2 data elements if not using zero bcs
    if(bc=="periodic" || bc=="symmetric") {
        nlevels <- nlevelsWT(data)
        if(is.na(nlevels)) stop("Data length is not power of two")
    }

# Check for correct type
    if(type != "wavelet" && type != "station")
        stop("Unknown type of wavelet decomposition")
    if(type == "station" && bc != "periodic")
        stop("Can only do periodic boundary conditions with station")

# Select the appropriate filter
    if(verbose == TRUE)
        cat("...done\nFilter...")
    filter <- filter.select(filter.number = filter.number, family = family)

# Build the first/last database
    if(verbose == TRUE)
        cat("...selected\nFirst/last database...")
    fl.dbase <- first.last.dh(LengthH = length(filter$H), DataLength =
        DataLength, type = type, bc = bc, firstk = firstk)

# Find number of levels in zero bc case
    if(bc=="zero")
        nlevels <- nrow(fl.dbase$first.last.d)

# Save time series attribute if there is one
    dtsp <- tsp(data)

# Put in the data
    C <- rep(0, fl.dbase$ntotal)
    C[1:DataLength] <- data
    if(verbose == TRUE)
        error <- 1
    else error <- 0
    if(verbose == TRUE) cat("built\n")

# Compute the decomposition
    if(verbose == TRUE)
        cat("Decomposing...\n")
    nbc <- switch(bc,
        periodic = 1,
        symmetric = 2,
        zero = 3)
    if(is.null(nbc))
        stop("Unknown boundary condition")
    ntype <- switch(type,
        wavelet = 1,
        station = 2)
    if(is.null(filter$G)) {
        wavelet.decomposition <- .C("wavedecomp_dh",
            C = as.double(C),
            D = as.double(rep(0, fl.dbase$ntotal.d)),
            H = as.double(filter$H),
            LengthH = as.integer(length(filter$H)),
            nlevels = as.integer(nlevels),
            firstC = as.integer(fl.dbase$first.last.c[, 1]),
            lastC = as.integer(fl.dbase$first.last.c[, 2]),
            offsetC = as.integer(fl.dbase$first.last.c[, 3]),
            firstD = as.integer(fl.dbase$first.last.d[, 1]),
            lastD = as.integer(fl.dbase$first.last.d[, 2]),
            offsetD = as.integer(fl.dbase$first.last.d[, 3]),
            ntype = as.integer(ntype),
            nbc = as.integer(nbc),
            error = as.integer(error), PACKAGE = "wavethresh")
    }
    if(verbose == TRUE)
        cat("done\n")
    error <- wavelet.decomposition$error
    if(error != 0) {
        cat("Error ", error, " occured in wavedecomp\n")
        stop("Error")
    }
    if(is.null(filter$G)) {
        l <- list(C = wavelet.decomposition$C, D =
            wavelet.decomposition$D, nlevels =
            nlevelsWT(wavelet.decomposition), fl.dbase = fl.dbase,
            filter = filter, type = type, bc = bc, date = date())
    }
    class(l) <- "wd"
    if(!is.null(dtsp))
        tsp(l) <- dtsp
    l
}
