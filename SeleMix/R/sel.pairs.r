#   A matrix of scatterplots with outliers and influentials units highlighted is produced
sel.pairs <- function (x, outl = rep(0, nrow(x)), sel = rep(0, nrow(x)), labs = NULL,
    log = TRUE, legend=TRUE )
{
    if (!inherits(x, c("matrix", "data.frame")))
        stop("data must be supplied as matrix or data frame ")
    testo_legenda <- c("Only Outlier           ", 
                       "Only Influential       ", 
                       "Outlier and Influential")
    titolo = "Selective Editing - outliers and influential errors"
    nvar1 <- ncol(x)
    if (is.null(labs)) {
        if (length(colnames(x) > 0))
            labs <- colnames(x)
        else labs <- paste("x", 1:nvar1, sep = "")
    }
    if (log == TRUE) {
        x[x == 0] <- 1e-07
        x <- log(x)
    }

    par(mfrow = c(nvar1, nvar1), mar = c(3, 2, 2, 3), oma = c(0, 0, 3, 0))    
        
    for (j in 1:nvar1) {
        for (i in 1:nvar1) {
            if (i == j) {
                boxplot(x[, i], main = NULL, col = "peachpuff",
                  horizontal = TRUE
#                , names = labs[i], show.names=TRUE
                )
                aa <- x[outl == 1 & sel == 0, i]
                points(aa, rep(1, length(aa)), pch = 21, cex = 1.2,
                  col = "blue4", bg = "blue")
                
                aa <- x[outl == 1 & sel == 1, i]
                points(aa, rep(1, length(aa)), pch = 22, cex = 1.2,
                  col = "cyan4", bg = "cyan")
                aa <- x[sel == 1 & outl == 0, i]
                points(aa, rep(1, length(aa)), pch = 24, cex = 1.2,
                  col = "red4", bg = "red") 
                par(mgp=c(2,1,0))
                title(xlab = labs[i], col.lab = "red3", cex=1.2, font.lab=2)
                par(mgp=c(3,1,0))

                if (i == 1 & j == 1 & legend == TRUE) {
                  legend("topleft", legend = testo_legenda,
                    pch = c(21, 24, 22), col = c("blue4","red4", "cyan4"),
                    pt.bg = c("blue", "red","cyan"), xjust = 0,
                    yjust = 1,  
                    cex = 2/3)
               }

            }
            else {
                plot(x[, i], x[, j], xlab = labs[i], ylab = labs[j],
                  pch = 21, col = "lightgrey")
                appo <- cbind(x[, i], x[, j])
                points(appo[outl == 1 & sel == 0, , drop = FALSE], pch = 21,
                  cex = 1.2, col = "blue4", bg = "blue")
                points(appo[sel == 1 & outl == 1, , drop = FALSE], pch = 22,
                  cex = 1.2, col = "cyan4", bg = "cyan")
                points(appo[outl == 0 & sel == 1, , drop = FALSE], pch = 24,
                  cex = 1.2, col = "red4",  bg = "red")

                }

        }
    }
    mtext(titolo, outer = TRUE)    
}  
