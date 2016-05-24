`plot.cca.pertables` <-
function (x, pch = 18, ...) 
{
    require(vegan)
    sites.raw <- scores(x$raw$cca.raw, display = "sites")
    biplot.raw <- scores(x$raw$cca.raw, display = "bp")
    sign.axis1 <- ifelse(sites.raw[, 1][abs(sites.raw[, 1]) == 
        max(abs(sites.raw[, 1]))] > 0, 1, 0)
    sign.axis2 <- ifelse(sites.raw[, 2][abs(sites.raw[, 2]) == 
        max(abs(sites.raw[, 2]))] > 0, 1, 0)
    for (i in 1:length(x$simulation$sites)) {
        x$simulation$sites[[i]][, 1] <- x$simulation$sites[[i]][, 
            1] * ifelse(ifelse(x$simulation$sites[[i]][rownames(x$simulation$sites[[i]]) == 
            attr(sign.axis1, "names"), 1] > 0, 1, 0) == sign.axis1, 
            1, -1)
        x$simulation$sites[[i]][, 2] <- x$simulation$sites[[i]][, 
            2] * ifelse(ifelse(x$simulation$sites[[i]][rownames(x$simulation$sites[[i]]) == 
            attr(sign.axis2, "names"), 2] > 0, 1, 0) == sign.axis2, 
            1, -1)
    }
    allsites <- do.call("rbind", x$simulation$sites)
    factor <- data.frame(CCA1 = paste(rownames(allsites), "CCA1"), 
        CCA2 = paste(rownames(allsites), "CCA2"))
    sign.axis1 <- ifelse(biplot.raw[, 1][abs(biplot.raw[, 1]) == 
        max(abs(biplot.raw[, 1]))] > 0, 1, 0)
    sign.axis2 <- ifelse(biplot.raw[, 2][abs(biplot.raw[, 2]) == 
        max(abs(biplot.raw[, 2]))] > 0, 1, 0)
    for (i in 1:length(x$simulation$biplot)) {
        x$simulation$biplot[[i]][, 1] <- x$simulation$biplot[[i]][, 
            1] * ifelse(ifelse(x$simulation$biplot[[i]][rownames(x$simulation$biplot[[i]]) == 
            attr(sign.axis1, "names"), 1] > 0, 1, 0) == sign.axis1, 
            1, -1)
        x$simulation$biplot[[i]][, 2] <- x$simulation$biplot[[i]][, 
            2] * ifelse(ifelse(x$simulation$biplot[[i]][rownames(x$simulation$biplot[[i]]) == 
            attr(sign.axis2, "names"), 2] > 0, 1, 0) == sign.axis2, 
            1, -1)
    }
    plot(x$raw$cca.raw, type = "n", ...)
    sapply(x$simulation$biplot, function(arg1) apply(arg1, 1, 
        function(arg2) arrows(0, 0, arg2[[1]], arg2[[2]], col = "grey", 
            length = 0.05, lty = 3)))
    sapply(x$simulation$sites, function(x) points(x, pch = pch, 
        col = "lightblue"))
    ordiellipse(allsites, factor[, 1], conf = 0.99, lty = 3, 
        col = "blue")
    points(sites.raw, pch = 3, cex = 2, col = "blue")
    text(x$raw$cca.raw, display = "bp", col = "red")
    text(x$raw$cca.raw, display = "sites", col = "blue", cex = 0.8, 
        pos = 2)
}

