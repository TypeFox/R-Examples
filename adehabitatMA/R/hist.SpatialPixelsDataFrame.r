hist.SpatialPixelsDataFrame <- function (x, type = c("h", "l", "b"),
                                         adjust = 1, col, border,
                                         lwd = 1, ...)
{
    type <- match.arg(type)
    if (is(x, "SpatialGrid"))
        fullgrid(x) = FALSE
    if (!inherits(x, "SpatialPixelsDataFrame"))
        stop("should be an object of class \"SpatialPixelsDataFramex\"")
    gr <- gridparameters(x)
    if (nrow(gr) > 2)
        stop("x should be defined in two dimensions")
    if ((gr[1, 2] - gr[2, 2])> get(".adeoptions",
                                   envir=.adehabitatMAEnv)$epsilon)
        stop("the cellsize should be the same in x and y directions")

    if (missing(col)) {
        col <- NULL
        cold <- "black"
    }
    else cold <- col
    if (missing(border))
        border <- "black"
    tab <- as.data.frame(x)
    tab <- tab[,-c((ncol(tab)-1):ncol(tab))]
    clas <- rep("", ncol(tab))
    for (j in 1:ncol(tab)) {
        w1 <- "q"
        if (is.factor(tab[, j]))
            w1 <- "f"
        clas[j] <- w1
    }
    if (any(clas == "f") & type == "l")
        warning("Type = 'l' is not possible for factors, type = 'h' used instead.\n")
    if (any(clas == "f") & type == "b")
        warning("Type = 'b' is not possible for factors, type = 'h' used instead.\n")
    old.par <- par(no.readonly = TRUE)
    on.exit(par(old.par))
    par(mar = c(3, 0.5, 2, 0.5))
    par(mfrow = rev(n2mfrow(ncol(tab))))
    f1 <- function(j) {
        tmpA <- tab[, j]
        name <- names(tab)[j]
        if (clas[j] == "f") {
            max <- max(table(tmpA))
            max <- max + max/20
            ylim <- c(0, max)
            barplot(unclass(summary(tmpA[!is.na(tmpA)])), ylim = ylim,
                border = border, col = col, main = name, ylab = NULL,
                axes = FALSE, ...)
        }
        else {
            if (type == "h") {
                xrange <- range(tmpA)
                G <- hist(tmpA, plot = FALSE)
                plot(G, freq = FALSE, border = border, col = col,
                  main = name, xlab = NULL, ylab = NULL, axes = FALSE,
                  ...)
                axis(side = 1)
            }
            if (type == "l") {
                dens <- density(tmpA, adjust = adjust, na.rm = TRUE)
                plot(dens, col = cold, type = "l", lwd = lwd,
                  main = name, xlab = NULL, ylab = "Density",
                  axes = FALSE, ...)
                mean <- mean(tmpA, na.rm = TRUE)
                lines(rep(mean, 2), c(0, dens$y[512 - sum(dens$x >
                  mean)]), col = cold, lty = 2, lwd = lwd)
                axis(side = 1)
            }
            if (type == "b") {
                xrange <- range(tmpA)
                G <- hist(tmpA, plot = FALSE)
                plot(G, freq = FALSE, border = border, col = col,
                  main = name, xlab = NULL, ylab = NULL, axes = FALSE,
                  ...)
                dens <- density(tmpA, adjust = adjust, na.rm = TRUE)
                lines(dens, col = cold, type = "l", lwd = lwd)
                mean <- mean(tmpA, na.rm = TRUE)
                lines(rep(mean, 2), c(0, dens$y[512 - sum(dens$x >
                  mean)]), col = cold, lty = 2, lwd = lwd)
                axis(side = 1)
            }
        }
    box()
    }
    lapply(1:ncol(tab), f1)
    return(invisible(NULL))
}
