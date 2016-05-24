show.simdata <-
function (dataset, plot.coord = FALSE, file.plot.coord = NA, 
    plot.tess = FALSE, file.plot.tess = NA, plot.freq.grid = FALSE, 
    file.plot.freq.grid = NA, loc.grid = 1, plot.freq.indiv = FALSE, 
    file.plot.freq.indiv = NA, loc.indiv = 1, zlim.freq = c(0, 
        1), plot.gen = FALSE, file.plot.gen = NA) 
{
    nindiv <- nrow(dataset$genotypes)
    nloc <- length(dataset$allele.numbers)
    if (plot.coord == TRUE) {
        if (is.na(file.plot.coord)) {
            dev.new()
        }
        else pdf(file.plot.coord)
        plot(dataset$coord.indiv[, 1], dataset$coord.indiv[, 
            2], xlab = "x dataset$coordinates", ylab = "y dataset$coordinates", 
            asp = )
        points(dataset$coord.nuclei[, 1], dataset$coord.nuclei[, 
            2], col = 2)
        text(dataset$coord.nuclei[, 1], dataset$coord.nuclei[, 
            2], 1:dataset$number.nuclei, pos = 1, col = 2, pch = 2, 
            cex = 1.2)
        text(dataset$coord.indiv[, 1], dataset$coord.indiv[, 
            2], dataset$color.nuclei[dataset$nearest.nucleus.indiv], 
            pos = 2)
        title(main = "Location of individuals")
        if (!is.na(file.plot.coord)) 
            dev.off()
    }
    if (plot.tess == TRUE) {
        if (is.na(file.plot.tess)) {
            dev.new()
        }
        else pdf(file.plot.tess)
        image(seq(from = dataset$coord.lim[1], to = dataset$coord.lim[2], 
            length = dataset$npix[1]), seq(from = dataset$coord.lim[3], 
            to = dataset$coord.lim[4], length = dataset$npix[2]), 
            matrix(nrow = dataset$npix[1], ncol = dataset$npix[2], 
                dataset$color.nuclei[dataset$nearest.nucleus.grid], 
                byrow = FALSE), xlab = "", ylab = "", col = terrain.colors(dataset$npop), 
            asp = 1)
        points(dataset$coord.nuclei[1, ], dataset$coord.nuclei[2, 
            ], col = 1, pch = ".", cex = 1.5)
        if (!is.na(file.plot.tess)) 
            dev.off()
    }
    if (plot.freq.grid == TRUE) {
        for (iloc in loc.grid) {
            for (iall in 1:(dataset$allele.numbers[iloc])) {
                if (is.na(file.plot.freq.grid)) {
                  dev.new()
                }
                else {
                  pdf(paste(substr(file.plot.freq.grid, start = 1, 
                    stop = nchar(file.plot.freq.grid) - 3), iall, 
                    ".ps", sep = ""))
                }
                FF <- rep(-999, prod(dataset$npix))
                for (ipop in 1:dataset$npop) {
                  ff <- dataset$freq.grid[ipop, , iloc, iall]
                  cc <- dataset$color.nuclei[dataset$nearest.nucleus.grid]
                  FF[cc == ipop] <- ff[cc == ipop]
                }
                image(seq(from = dataset$coord.lim[1], to = dataset$coord.lim[2], 
                  length = dataset$npix[1]), seq(from = dataset$coord.lim[3], 
                  to = dataset$coord.lim[4], length = dataset$npix[2]), 
                  matrix(nrow = dataset$npix[1], ncol = dataset$npix[2], 
                    FF, byrow = FALSE), col = heat.colors(500), 
                  xlab = "", ylab = "", zlim = zlim.freq, asp = 1)
                title(paste("Frequencies of allele #", iall, 
                  "at locus #", iloc))
                if (!is.na(file.plot.freq.grid)) 
                  dev.off()
            }
        }
    }
    if (plot.freq.indiv == TRUE) {
        for (iloc in 1:loc.grid) {
            for (iall in 1:(dataset$allele.numbers[iloc])) {
                FF <- rep(NA, nindiv)
                for (ipop in 1:dataset$npop) {
                  ff <- dataset$freq.indiv[ipop, , iloc, iall]
                  cc <- dataset$color.nuclei[dataset$nearest.nucleus.indiv]
                  FF[cc == ipop] <- ff[cc == ipop]
                }
                if (is.na(file.plot.freq.indiv)) {
                  dev.new()
                }
                else pdf(paste(substr(file.plot.freq.indiv, start = 1, 
                  stop = nchar(file.plot.freq.indiv) - 3), iall, 
                  ".ps", sep = ""))
                look <- as.image(x = dataset$coord.indiv, Z = FF)
                image.plot(look, main = paste("Field of frequencies for locus #", 
                  iloc, "allele #", iall), asp = 1)
                points(dataset$coord.nuclei[1, ], dataset$coord.nuclei[2, 
                  ], col = 2, cex = 2, lwd = 3)
                if (!is.na(file.plot.freq.indiv)) 
                  dev.off()
            }
        }
    }
    if (plot.gen == TRUE) {
        nindiv <- ncol(dataset$coord.indiv)
        for (iloc in 1:loc.grid) {
            if (is.na(file.plot.gen)) {
                dev.new()
            }
            else pdf(paste(substr(file.plot.gen, start = 1, stop = nchar(file.plot.gen) - 
                3), iloc, ".ps", sep = ""))
            plot(dataset$coord.indiv[, 1], dataset$coord.indiv[, 
                2], type = "n", cex = 1, lwd = 1, xlab = "", 
                ylab = "", asp = 1)
            title(paste("Genotypes at locus #", iloc))
            text(dataset$coord.indiv[, 1], dataset$coord.indiv[, 
                2], dataset$genotypes[, 2 * (iloc - 1) + 1], 
                cex = 1, lwd = 2, pos = 2, adj = 0)
            text(dataset$coord.indiv[, 1], dataset$coord.indiv[, 
                2], dataset$genotypes[, 2 * (iloc - 1) + 2], 
                cex = 1, lwd = 2, pos = 4, col = 2, adj = 0.5, 
                offset = -0.5)
            if (!is.na(file.plot.gen)) 
                dev.off()
        }
    }
}
