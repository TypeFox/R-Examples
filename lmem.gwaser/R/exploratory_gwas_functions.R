#########################################
#' Principal Component Analysis.
#'
#'Performs Principal Component Analysis of marker data from
#'an object of cross class created by the gwas.cross function.
#'
#' @usage pca.analysis(crossobj, p.val)
#'
#' @param crossobj An object of class = cross obtained from the gwas.cross
#' function from this package, or the read.cross
#' function from r/qtl package (Broman and Sen, 2009).
#' This file contains phenotypic means, genotypic marker score, and genetic map.
#'
#' @param p.val Alpha level (a number) to identify the number of significant axis
#'
#' @return A PCA plot with two principal components and a scree
#' plot for all significant axes indicating the proportion of
#' the variance explained by each marker.
#'
#' @references Comadran J, Thomas W, van Eeuwijk F, Ceccarelli S, Grando S,
#'             Stanca A, Pecchioni N, Akar T, Al-Yassin A,
#'             Benbelkacem A, Ouabbou H, Bort J, Romagosa I,
#'             Hackett C, Russell J (2009) Patterns of genetic diversity
#'             and linkage disequilibrium in a highly structured
#'             Hordeum vulgare association-mapping population for
#'             the Mediterranean basin. Theor Appl Genet 119:175-187
#'
#'             Becker, R. A., Chambers, J. M. and Wilks, A. R. (1988) The New S
#'             Language. Wadsworth & Brooks/Cole.
#'
#'             Mardia, K. V., J. T. Kent, and J. M. Bibby (1979) Multivariate
#'             Analysis, London: Academic Press.
#'
#'             Venables, W. N. and B. D. Ripley (2002) Modern Applied Statistics
#'             with S, Springer-Verlag.
#'
#' @author Lucia Gutierrez
#'
#' @details Performs two plots.
#'
#' @note In gwas.memq function, the pca.anlysis function is already included.
#'
#' @seealso gwas.analysis
#'
#' @import qtl
#' @import pastecs
#' @import LDheatmap
#' @import genetics
#'
#' @export
#'
#' @examples
#' \dontrun{
#' data (QA_geno)
#' data (QA_map)
#' data (QA_pheno)
#'
#' P.data <- QA_pheno
#' G.data <- QA_geno
#' map.data <- QA_map
#'
#' cross.data <- gwas.cross (P.data, G.data, map.data,
#' cross='gwas', heterozygotes=FALSE)
#' summary (cross.data)
#'
#' pca <- pca.analysis(crossobj=cross.data, p.val=0.05)
#'
#'}
pca.analysis <- function(crossobj, p.val) {
    scores <- NULL
    for (i in 1:nchr(crossobj)) {
        a <- paste("crossobj$geno$'", i, "'$data", sep = "")
        s <- eval(parse(text = a))
        scores <- cbind(scores, s)
    }
    scores[scores == 1] <- 0
    # Impute missing values with means
    average <- NULL
    for (i in 1:dim(scores)[2]) {
        a <- mean(scores[, i], na.rm = TRUE)
        average <- cbind(average, a)
    }
    for (i in 1:dim(scores)[1]) {
        for (j in 1:dim(scores)[2]) {
            if (is.na(scores[i, j]) == TRUE)
                scores[i, j] <- average[1, j]
        }
    }
    pca.analysis1 <- prcomp(scores, scale = TRUE)
    meff1 <- nind(crossobj) - 1
    ngeno <- nind(crossobj)
    lambda <- summary(pca.analysis1)$importance[2, ]
    lambda <- lambda[1:meff1]
    # loop to test whether each axis is significant or not
    TM <- NULL
    x <- 10
    while (x > 0.9792895) {
        # This is using the nominal p-value of 0.05 from TW statistic
        meff <- length(lambda)
        sum.lambda <- sum(lambda, na.rm = TRUE)
        sum.lambda2 <- sum(lambda^2, na.rm = TRUE)
        neff <- ((ngeno + 1) * (sum.lambda^2)) /
          (((ngeno - 1) * sum.lambda2) - (sum.lambda^2))
        mu <- ((sqrt(neff - 1) + sqrt(meff))^2)/(neff)
        sigma <- ((sqrt(neff - 1) + sqrt(meff)) / neff) * (((1 / (sqrt(neff - 1)))
          + (1/sqrt(meff))) ^ (1 / 3))
        L1 <- (meff * lambda[1]) / sum.lambda
        x <- (L1 - mu) / sigma
        TM <- c(TM, x)
        lambda <- lambda[2:length(lambda)]
    }
    n.signif <- length(TM) - 1

    plot(summary(pca.analysis1)$importance[2, ],
      ylim = c(0, 1), xlim = c(0, n.signif), col = "red",
      pch = 4, cex = 0.5, ylab = "Proportion of variance",
      main = "Proportion of variance explained by each significant axis")
     lines(summary(pca.analysis1)$importance[2, ], col = "red")
     points(summary(pca.analysis1)$importance[3, ], col = "blue", pch = 4, cex = 0.5)
     lines(summary(pca.analysis1)$importance[3, ], col = "blue")
    legend("topleft", col = c("blue", "red"), lty = c(1, 1),
      legend = c("Cummulative variance explained",
      "Proportion of total variance explained"))

    plot(pca.analysis1$x[, 1], pca.analysis1$x[, 2],
      xlab = paste("PC1 (", (summary(
        pca.analysis1)$importance[2, 1] * 100), " %)",
        sep = ""), ylab = paste("PC2 (",
          (summary(pca.analysis1)$importance[2, 2] * 100), " % )",
          sep = ""), main = "Scores of PC1 and PC2")

    out <- NULL
    out$number.pca <- n.signif
    out$prop.variance <- summary(pca.analysis1)$importance[, 1:n.signif]
    out$scores <- pca.analysis1$x[, 1:n.signif]
    out
}

#########################################
#' Linkage Disequilibrium heatmap plot
#'
#'Performs a Linkage Disequilibrium heatmap plot for the GWAS analysis.
#'Non-random association of markers (linkage disequilibrium) are estimated
#'as Lewont in \eqn{'/c}s \eqn{D'/c} (Lewontin\eqn{'/c}s, 1964) with
#'the LD function of the genetics package (Warnes and Leisch, 2005) and
#'isualized with the LD.heatmap package (Shin et al., 2015). \eqn{D'/c}
#'is estimated as: \eqn{D'/c = D/DMax} where D= pAB - pApB,
#'DMax=Min (pApb, papB), and pA is the probability of the A allele for marker 1,
#'pa=1-pA, pB is the probability of the B allele for marker 2, pb=1-pB,
#'and pAB is the probability of AB alleles.
#'
#' @usage linkdis.plots (crossobj, heterozygotes, chr)
#'
#' @param crossobj An object of class = cross obtained from the gwas.cross
#' function from this package, or the read.cross
#' function from r/qtl package (Broman and Sen, 2009).
#' This file contains phenotypic means, genotypic marker score, and genetic map.
#'
#' @param heterozygotes Logical value indicating whether heterozygotes are present.
#'
#' @param chr A vector containing chromosome number to use.
#'
#' @return Return a Linkage Disequilibrium heatmap plot.
#'
#' @references Comadran J, Thomas W, van Eeuwijk F, Ceccarelli S, Grando S,
#'             Stanca A, Pecchioni N, Akar T, Al-Yassin A,
#'             Benbelkacem A, Ouabbou H, Bort J, Romagosa I,
#'             Hackett C, Russell J (2009) Patterns of genetic diversity
#'             and linkage disequilibrium in a highly structured
#'             Hordeum vulgare association-mapping population for
#'             the Mediterranean basin. Theor Appl Genet 119:175-187
#'
#'             Warnes, G; Leisch, F. 2005. Genetics: Population genetics R package
#'             1.2.0. Lewontin, R. 1964. The interaction of selection and linkage.
#'              I. General Considerations: Heterotic models. Genetics 49: 49-67.
#'
#' @author Lucia Gutierrez
#'
#' @details The function returns the LD.heatmap for the chromosomes selected.
#'
#' @note When large data sets are being used, linkdis.plots is encourage to be
#' performed for each chromosome separately.
#'
#' @seealso gwas.analysis
#'
#' @import qtl
#' @import pastecs
#' @import genetics
#' @import LDheatmap
#'
#' @export
#'
#' @examples
#' \dontrun{
#' data (QA_geno)
#' data (QA_map)
#' data (QA_pheno)
#'
#' P.data <- QA_pheno
#' G.data <- QA_geno
#' map.data <- QA_map
#'
#' cross.data <- gwas.cross (P.data, G.data, map.data,
#' cross='gwas', heterozygotes=FALSE)
#' summary (cross.data)
#'
#'LD.plots
#'
#'linkdis.plots(crossobj = cross.data, heterozygotes = FALSE, chr = c('1', '7'))
#'}
#'
linkdis.plots <- function(crossobj, heterozygotes, chr) {

    # construct matrix in correct format
    data <- NULL

    for (i in 1:length(chr)) {
        a <- paste("crossobj$geno$'", chr[i], "'$data", sep = "")
        p1 <- eval(parse(text = a))
        data <- cbind(data, p1)
    }

    if (heterozygotes == "FALSE") {
        data[data == 1] <- "A/A"
        data[data == 2] <- "B/B"
    }

    ###### Check this part!!!!
    if (heterozygotes == "TRUE") {
        data[data == 1] <- "A/A"
        data[data == 2] <- "B/B"
        data[data == 3] <- "A/B"
    }

    genos <- genotype(data[, 1])
    for (i in 2:dim(data)[2]) {
        g <- genotype(data[, i])
        genos <- data.frame(genos, g)
    }

    #### add marker names calculate LD
    ld <- LD(genos)

    # plot LD heatmap
    plot.hm <- LDheatmap(ld$"D'", color = colorRampPalette(c("red", "yellow"))(50))
    plot(plot.hm)
}
