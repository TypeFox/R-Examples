histniche <- function(x, pr, type = c("h", "l"), adjust = 1,
                      Acol, Ucol, Aborder, Uborder, Alwd = 1, Ulwd = 1,
                      ylim, ncla=15, ...)
{
    ## Verifications
    type <- match.arg(type)
    if (nrow(x)!=length(pr))
        stop("x should have the same length as pr")
    glou <- missing(ylim)

    ## Graphical settings
    if (missing(Acol)) {
        Acol <- NULL
        Acolf <- "white"
        Acold <- "black"
    }
    else {
        Acold <- Acol
        Acolf <- Acol
    }
    if (missing(Aborder))
        Aborder <- "black"
    if (missing(Ucol)) {
        Ucol <- gray(0.8)
        Ucold <- gray(0.8)
    }
    else
        Ucold <- Ucol
    if (missing(Uborder))
        Uborder <- gray(0.8)

    ## The type of variable for which a histogram is wanted
    clas <- rep("", ncol(x))
    for (j in 1:ncol(x)) {
        w1 <- "q"
        if (is.factor(x[, j]))
            w1 <- "f"
        clas[j] <- w1
    }
    if (any(clas == "f") & type == "l")
        warning("Type = 'l' is not possible for factors, type = 'h' used instead.\n")

    ## Again graphical settings
    if (ncol(x)>1) {
        old.par <- par(no.readonly = TRUE)
        on.exit(par(old.par))
        par(mar = c(0.5, 0.5, 2, 0.5))
        par(mfrow = rev(n2mfrow(ncol(x))))
    }

    ## The function used for each histogram
    f1 <- function(j) {

        ## Use and availability
        tmpU <- rep(x[, j], pr)
        tmpA <- x[, j]
        name <- names(x)[j]


        ## For factor maps: a barplot
        if (clas[j] == "f") {
            par(mar = c(3, 0.5, 2, 0.5))
            mat <- t(cbind(table(tmpA), table(tmpU)))
            mat <- lapply(1:2, function(i) mat[i, ]/sum(mat[i,]))
            mat <- rbind(mat[[1]], mat[[2]])
            max <- max(mat)
            max <- max + max/20
            if (glou)
                ylim <- c(0, max)
            barplot(mat, col = c(Acolf, Ucol), border = c(Aborder, Uborder),
                    ylim = ylim, main = name, ylab = NULL, axes = FALSE,
                    beside = TRUE, ...)
            par(mar = c(0.5, 0.5, 2, 0.5))
        }
        else {

            ## for continuous maps: either a histogram...
            if (type == "h") {
                xrange <- range(tmpA)
                H <- hist(tmpU, plot = FALSE, br = seq(range(tmpA)[1],
                                              range(tmpA)[2], length = ncla))
                G <- hist(tmpA, plot = FALSE, br = seq(range(tmpA)[1],
                                              range(tmpA)[2], length = ncla))
                if (glou)
                    ylim <- c(0, max(H$density, G$density))
                plot(H, freq = FALSE, col = Ucol, border = Uborder,
                     xlim = xrange, ylim = ylim, main = name, xlab = NULL,
                     ylab = "Density", axes = FALSE, ...)
                plot(G, freq = FALSE, col = Acol, border = Aborder, add = TRUE)
            }

            ## ... or a smoothing of the density
            if (type == "l") {
                densA <- density(tmpA, adjust = adjust)
                densU <- density(tmpU, adjust = adjust, from = min(densA$x),
                                 to = max(densA$x))
                max <- max(densU$y, densA$y)
                max <- max + max/20
                if (glou)
                    ylim <- c(0, max)
                plot(densU, col = Ucol, ylim = ylim, type = "l",
                     lwd = Ulwd, main = name, xlab = NULL, ylab = "Density",
                     axes = FALSE, ...)
                lines(rep(mean(tmpU), 2), c(0, densU$y[512 - sum(densU$x >
                                                                 mean(tmpU))]),
                      col = Ucol, lty = 2, lwd = Ulwd)
                lines(densA, col = Acold, lwd = Alwd)
                lines(rep(mean(tmpA), 2), c(0, densA$y[512 - sum(densA$x >
                                                                 mean(tmpA))]),
                      col = Acold, lty = 2, lwd = Alwd)
            }
        }
        box()
    }


    f2 <- function(j) {

        ## Use and availability
        tmpU <- rep(x[,1], pr)
        tmpA <- x[,1]


        ## For factor maps: a barplot
        if (clas[j] == "f") {
            mat <- t(cbind(table(tmpA), table(tmpU)))
            mat <- lapply(1:2, function(i) mat[i, ]/sum(mat[i,]))
            mat <- rbind(mat[[1]], mat[[2]])
            max <- max(mat)
            max <- max + max/20
            if (glou)
                ylim <- c(0, max)
            barplot(mat, col = c(Acolf, Ucol), border = c(Aborder, Uborder),
                    ylim = ylim, beside=TRUE,...)
        }
        else {
            ## for continuous maps: either a histogram...
            if (type == "h") {
                xrange <- range(tmpA)
                H <- hist(tmpU, plot = FALSE, breaks = seq(range(tmpA)[1],
                                              range(tmpA)[2], length = ncla))
                G <- hist(tmpA, plot = FALSE, breaks = seq(range(tmpA)[1],
                                              range(tmpA)[2], length = ncla))
                if (glou)
                    ylim <- c(0, max(H$density, G$density))
                plot(H, freq = FALSE, col = Ucol, border = Uborder,
                     xlim = xrange, ylim = ylim, ...)
                plot(G, freq = FALSE, col = Acol, border = Aborder, add = TRUE)
            }

            ## ... or a smoothing of the density
            if (type == "l") {
                densA <- density(tmpA, adjust = adjust)
                densU <- density(tmpU, adjust = adjust, from = min(densA$x),
                                 to = max(densA$x))
                max <- max(densU$y, densA$y)
                max <- max + max/20
                if (glou)
                    ylim <- c(0, max)
                plot(densU, col = Ucol, ylim = ylim, type = "l",
                     lwd = Ulwd, ...)
                lines(rep(mean(tmpU), 2), c(0, densU$y[512 - sum(densU$x >
                                                                 mean(tmpU))]),
                      col = Ucol, lty = 2, lwd = Ulwd)
                lines(densA, col = Acold, lwd = Alwd)
                lines(rep(mean(tmpA), 2), c(0, densA$y[512 - sum(densA$x >
                                                                 mean(tmpA))]),
                      col = Acold, lty = 2, lwd = Alwd)
            }
        }
    }



    ## And we apply this function to each variable
    if (ncol(x)>1) {
        lapply(1:ncol(x), f1)
    } else {
        f2(1)
    }

    return(invisible(NULL))
}

