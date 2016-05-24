plottwoloci <-
function (x, loc1, loc2, shift = 0.02, stderr = TRUE, GPcol = c("red", 
    "green", "blue"), arrowcol = aperm(matrix(rep(GPcol, 3), 
    3, 3)), ylab = "Genotypic values", xlab = "Phenotype", ylim = NA, 
    main = paste("GPmap:", "Loc", loc1, "/", "Loc", loc2, "(", 
        err.name, ")"), ...) 
{
    if (stderr == TRUE) {
        err <- 1
        err.name = "stderr"
    }
    else {
        err <- 1.96
        err.name <- "IC95%"
    }
    Gv <- recomputeG(x, c(loc1, loc2))
    segp1 <- matrix((Gv[, 1] - err * (Gv[, 2]/2)), 3, 3)
    segp2 <- matrix((Gv[, 1] + err * (Gv[, 2]/2)), 3, 3)
    matGval <- matrix(Gv[, 1], 3, 3)
    lab1 <- c(paste(tolower(LETTERS[loc1]), tolower(LETTERS[loc1])), 
        paste(tolower(LETTERS[loc1]), LETTERS[loc1]), paste(LETTERS[loc1], 
            LETTERS[loc1]))
    lab2 <- c(paste(tolower(LETTERS[loc2]), tolower(LETTERS[loc2])), 
        paste(tolower(LETTERS[loc2]), LETTERS[loc2]), paste(LETTERS[loc2], 
            LETTERS[loc2]))
    z <- matrix(c(0, 1, 2), 3, ((3^2)/3))
    c <- c(0, shift, -shift)
    z <- aperm(apply(z, 1, function(x) {
        x + c
    }))
    if (is.na(ylim[1])) {
        ylim <- c(floor(min(segp1)), ceiling(max(segp2)) * 1.12)
    }
    matplot(z, matGval, xaxt = "n", xlab = xlab, col = GPcol, 
        type = "b", cex.main = 1, cex.axis = 0.7, pch = 3, ylim = ylim, 
        main = main, ylab = ylab, ...)
    axis(1, at = seq(0, 2, 1), labels = lab1, cex.axis = 0.8, 
        tick = FALSE, las = 1, padj = -2)
    arrows(z, segp1, z, segp2, code = 3, length = 0.02, angle = 90, 
        col = arrowcol, cex.lab = 0.3)
    legend("top", legend = c(if (x$nloc > 1) {
        lab2
    }), pch = c("_"), col = GPcol, bty = "n", cex = 0.8, pt.cex = 0.8)
}
