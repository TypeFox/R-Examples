plotlocus <-
function (reg, loc, winscale = 0.5, effectscol = c("red", "blue"), 
    referencecol = "green", xlab = "Phenotype", ylab = "Genotypic values", 
    main = c(paste("Linear regression : Locus", loc)), ylim = NA, 
    ...) 
{
    "freqloc" <- function(reg, loc) {
        ss <- Sloc(reference = reg$reference, i = loc, genZ = reg$genZ)
        return(solve(ss)[1, ])
    }
    "meangenot" <- function(reg, loc) {
        ss <- Sloc(reference = reg$reference, i = loc, genZ = reg$genZ)
        return(sum(solve(ss)[1, ] * rep(c(0, 1, 2))))
    }
    Gval <- recomputeG(reg, loc)[, 1]
    Effect <- isolateLoc(reg, loc, remove0 = T)
    freq <- freqloc(reg, loc)
    xR <- meangenot(reg, loc)
    yR <- Effect[1]
    a <- Effect[2]
    y0 <- yR - a * xR
    miny <- min(Gval, y0, (a * 2 + y0))
    maxy <- max(Gval, y0, (a * 2 + y0))
    marge <- abs(maxy - miny) * 0.1
    if (is.na(ylim[1])) {
        ylim <- c(miny - marge, maxy + marge)
    }
    plot(c(0, 1, 2), Gval, pch = 20, xlim = c(0 - winscale, 2 + 
        winscale), ylim = ylim, cex.main = 1, cex.axis = 0.7, 
        cex = (freq + 0.15) * 7, xaxt = "n", xlab = xlab, ylab = ylab, 
        main = main, ...)
    axis(1, at = seq(0, 2, 1), labels = c(paste(tolower(LETTERS[loc]), 
        tolower(LETTERS[loc])), paste(tolower(LETTERS[loc]), 
        LETTERS[loc]), paste(LETTERS[loc], LETTERS[loc])), cex.axis = 0.8, 
        tick = F, las = 1, padj = -2)
    abline(y0, a)
    segments(1, Gval[2], 1, 0.5 * Gval[1] + 0.5 * Gval[3], col = effectscol[1])
    segments(0, Gval[1], 2, Gval[3], lty = "dashed")
    text(1 + 0.15, (((a * 1 + y0) + Gval[2])/2), "d", col = effectscol[1])
    segments(2, a * 2 + y0, 2.2, a * 2 + y0, lty = "dashed")
    segments(0, y0, 2.2, y0, lty = "dashed")
    arrows(2.2, y0, 2.2, a * 2 + y0, col = effectscol[2], code = 3, 
        length = 0.1)
    text(2.35, ((y0 + (a * 2 + y0))/2), "2a", col = effectscol[2])
    points(xR, yR, pch = 18, cex = 2, col = referencecol)
}
