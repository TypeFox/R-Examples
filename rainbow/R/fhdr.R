fhdr = function (data, alpha = c(0.01, 0.5), label = TRUE, xlab, ylab, 
    plotlegend, legendpos, ncol, projmethod, ...) 
{
    y = t(data$y)
    x = data$x
	if(projmethod == "PCAproj")
	{
		sco = PCAproj(y, k = 2, center = median)$scores
	}
	if(projmethod == "rapca")
	{
		sco = fdpca(x, data$y)$coeff[,2:3]
		rownames(sco) = 1:ncol(data$y)
	}
    ylim = range(y, na.rm = TRUE)
    band = Hscv.diag(sco, binned = TRUE)
    if (any(diag(band) < 10^(-30))) {
        stop("Computationally singular due to at least one of the diagonal elements of bandwidth matrix is very close to 0.")
    }
    else {
        den <- kde(x = sco, H = 0.8 * band)
        hdr1 <- hdrcde::hdr.2d(sco[, 1], sco[, 2], prob = alpha,
            list(x = den$eval.points[[1]], y = den$eval.points[[2]],
                z = den$estimate))
        index <- (hdr1$fxy < min(hdr1$falpha))
        outlier <- which(as.vector(index))
        index <- (hdr1$fxy > hdr1$falpha[1])
        out <- which(as.vector(index))
        outcurve <- y[out, ]
        maximum2 <- apply(outcurve, 2, max)
        minimum2 <- apply(outcurve, 2, min)
        index <- (hdr1$fxy > hdr1$falpha[2])
        inside <- as.matrix(which(as.vector(index)))
        insideone <- which(pam(inside, 2)$clustering == 2)
        insidetwo <- which(pam(inside, 2)$clustering == 1)
        insideone <- y[inside[insideone, ], ]
        insidetwo <- y[inside[insidetwo, ], ]
        maximum3 <- apply(insideone, 2, max)
        minimum3 <- apply(insideone, 2, min)
        maximum4 <- apply(insidetwo, 2, max)
        minimum4 <- apply(insidetwo, 2, min)
        dist <- (sco[, 1] - hdr1$mode[1])^2 + (sco[, 2] - hdr1$mode[2])^2
        center <- order(dist)[1]
        centercurve <- y[center, ]
        n <- length(outlier)
        x <- data$x
        plot(c(x, rev(x)), c(maximum2, rev(minimum2)), type = "n", 
            main = "", ylim = ylim, xlab = xlab, ylab = ylab, 
            ...)
        polygon(c(x, rev(x)), c(maximum2, rev(minimum2)), border = FALSE, 
            col = "light gray", ylim = ylim, ...)
        polygon(c(x, rev(x)), c(maximum3, rev(minimum3)), border = FALSE, 
            col = "dark gray", ...)
        polygon(c(x, rev(x)), c(maximum4, rev(minimum4)), border = FALSE, 
            col = "dark gray", ...)
        lines(fts(x, centercurve), col = "black", ...)
        if (n > 0) {
            outliercurve <- y[outlier, ]
            if (n == 1) {
                lines(fts(x, as.matrix(outliercurve)), col = rainbow(n), 
                  ...)
            }
            if (n > 1) {
                lines(fts(x, t(outliercurve)), col = rainbow(n), 
                  ...)
            }
            if (plotlegend == TRUE) {
                legend(legendpos, c(colnames(data$y)[outlier]), 
                  col = rainbow(n), lty = 1, ncol = ncol, ...)
            }
        }
        return(outlier)
    }
}

